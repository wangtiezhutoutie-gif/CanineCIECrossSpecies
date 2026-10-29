#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
gsea_root <- file.path(project_root, "06_pseudobulk", "gsea")
gene_set_dir <- file.path(gsea_root, "gene_sets")
mapping_dir <- file.path(gsea_root, "mapping")
full7_root <- file.path(gsea_root, "full7")
leading_dir <- file.path(gsea_root, "leading_edges")
figure_dir <- file.path(gsea_root, "figures")
dir.create(leading_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

lineages <- c("Epithelial", "T", "Myeloid")
collection_specs <- list(
  Hallmark = list(
    file = file.path(gene_set_dir, "hallmark_dog_ortholog_projected_membership.tsv.gz"),
    outdir = "hallmark", tag = "Hallmark"
  ),
  Reactome = list(
    file = file.path(gene_set_dir, "reactome_canis_familiaris_gene_sets.tsv.gz"),
    outdir = "reactome", tag = "Reactome"
  ),
  GO_BP = list(
    file = file.path(gene_set_dir, "go_bp_canis_familiaris_gene_sets.tsv.gz"),
    outdir = "go_bp", tag = "GOBP"
  )
)

min_size <- 15L
max_size <- 500L
fgsea_eps <- 1e-10
base_seed <- 20260824L

wrap_text <- function(x, width = 48L) {
  vapply(x, function(z) paste(strwrap(gsub("_", " ", z), width = width), collapse = "\n"), character(1))
}

select_top_both <- function(dt, n_each = 12L) {
  pos <- dt[NES > 0][order(padj, -abs(NES))][seq_len(min(.N, n_each))]
  neg <- dt[NES < 0][order(padj, -abs(NES))][seq_len(min(.N, n_each))]
  rbindlist(list(pos, neg), fill = TRUE)
}

plot_overview <- function(result, lineage, collection, filename) {
  show <- select_top_both(copy(result), if (collection == "Hallmark") 15L else 12L)
  if (nrow(show) == 0L) return(invisible(NULL))
  show[, plot_label := wrap_text(pathway_name)]
  show[, plot_label := factor(plot_label, levels = plot_label[order(NES)])]
  show[, fdr_group := fifelse(padj < 0.05, "FDR < 0.05", fifelse(padj < 0.10, "FDR < 0.10", "FDR >= 0.10"))]
  p <- ggplot(show, aes(x = NES, y = plot_label, color = direction, size = -log10(pmax(padj, 1e-300)), alpha = fdr_group)) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_point() +
    scale_color_manual(values = c(CIE_higher = "#B33A3A", CIE_lower = "#2E63A3")) +
    scale_alpha_manual(values = c("FDR < 0.05" = 1, "FDR < 0.10" = 0.8, "FDR >= 0.10" = 0.45)) +
    labs(title = paste(lineage, collection, "Full7 pre-ranked GSEA"), x = "Normalized enrichment score (NES)", y = NULL, size = "-log10 FDR", color = NULL, alpha = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), panel.grid.major.y = element_blank(), legend.position = "bottom")
  ggsave(filename, p, width = 9.2, height = if (collection == "Hallmark") 10 else 8.5, dpi = 180)
}

all_results <- list()
all_leading <- list()

for (lineage_index in seq_along(lineages)) {
  lineage <- lineages[[lineage_index]]
  rank_path <- file.path(mapping_dir, paste0(lineage, "_Full7_collapsed_rank.csv"))
  rank_dt <- fread(rank_path)
  stopifnot(!anyDuplicated(rank_dt$target_geneset_id))
  setorder(rank_dt, -rank_metric, feature_id)
  stats <- rank_dt$rank_metric
  names(stats) <- as.character(rank_dt$target_geneset_id)
  if (any(!is.finite(stats))) stop("Non-finite rank metric in ", lineage)

  for (collection_index in seq_along(collection_specs)) {
    collection <- names(collection_specs)[[collection_index]]
    spec <- collection_specs[[collection]]
    members <- unique(fread(spec$file), by = c("pathway_id", "target_geneset_id"))
    members[, target_geneset_id := as.character(target_geneset_id)]
    pathway_info <- unique(members[, .(
      pathway_id, pathway_name, mapping_type, source_version, source_species, target_species
    )], by = "pathway_id")
    pathways <- split(members$target_geneset_id, members$pathway_id)
    set.seed(base_seed + lineage_index * 100L + collection_index)
    message("Running ", lineage, " / ", collection, " with ", length(stats), " ranked genes and ", length(pathways), " source pathways...")
    fg <- as.data.table(fgseaMultilevel(
      pathways = pathways,
      stats = stats,
      minSize = min_size,
      maxSize = max_size,
      eps = fgsea_eps,
      scoreType = "std",
      nproc = 1L
    ))
    setnames(fg, "pathway", "pathway_id")
    fg[, padj := p.adjust(pval, method = "BH")]
    result <- merge(fg, pathway_info, by = "pathway_id", all.x = TRUE, sort = FALSE)
    result[, `:=`(
      lineage = lineage,
      collection = collection,
      leading_edge_n = lengths(leadingEdge),
      direction = fifelse(NES > 0, "CIE_higher", "CIE_lower")
    )]

    leading_rows <- list()
    for (i in seq_len(nrow(result))) {
      ids <- as.character(result$leadingEdge[[i]])
      if (length(ids) == 0L) next
      detail <- rank_dt[match(ids, as.character(rank_dt$target_geneset_id))]
      detail[, `:=`(
        lineage = lineage,
        collection = collection,
        pathway_id = result$pathway_id[[i]],
        pathway_name = result$pathway_name[[i]],
        pathway_NES = result$NES[[i]],
        pathway_padj = result$padj[[i]],
        pathway_direction = result$direction[[i]],
        mapping_type = result$mapping_type[[i]],
        source_version = result$source_version[[i]]
      )]
      leading_rows[[length(leading_rows) + 1L]] <- detail[, .(
        lineage, collection, pathway_id, pathway_name, pathway_NES, pathway_padj,
        pathway_direction, target_geneset_id, target_symbol, feature_id,
        gene_symbol, unique_var_name, logFC, statistic_F, PValue, FDR,
        rank_metric, mapping_type, source_version
      )]
    }
    leading <- rbindlist(leading_rows, fill = TRUE)
    leading_key <- leading[, .(
      leading_edge_genes = paste(unique(na.omit(target_symbol[nzchar(target_symbol)])), collapse = ";"),
      leading_edge_feature_ids = paste(unique(feature_id), collapse = ";")
    ), by = pathway_id]
    result <- merge(result, leading_key, by = "pathway_id", all.x = TRUE, sort = FALSE)
    result[is.na(leading_edge_genes), leading_edge_genes := ""]
    result[is.na(leading_edge_feature_ids), leading_edge_feature_ids := ""]
    result[, leadingEdge := NULL]
    setcolorder(result, c(
      "lineage", "collection", "pathway_id", "pathway_name", "ES", "NES", "pval", "padj",
      "size", "leading_edge_n", "leading_edge_genes", "leading_edge_feature_ids", "direction",
      "mapping_type", "source_version", "source_species", "target_species", "log2err"
    ))
    result[, abs_NES_sort := abs(NES)]
    setorder(result, padj, -abs_NES_sort, pathway_id)
    result[, abs_NES_sort := NULL]

    outdir <- file.path(full7_root, spec$outdir)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    result_path <- file.path(outdir, paste0(lineage, "_", spec$tag, "_GSEA.csv"))
    leading_path <- file.path(leading_dir, paste0(lineage, "_", spec$tag, "_leading_edges.csv"))
    fwrite(result, result_path)
    fwrite(leading, leading_path)
    plot_overview(result, lineage, collection, file.path(figure_dir, paste0(lineage, "_", spec$tag, "_NES_overview.png")))
    message("Completed ", lineage, " / ", collection, ": ", nrow(result), " pathways tested; FDR<0.05=", result[padj < 0.05, .N], ".")
    all_results[[paste(lineage, collection, sep = "|")]] <- result
    all_leading[[paste(lineage, collection, sep = "|")]] <- leading
  }

  combined_result <- rbindlist(all_results[grepl(paste0("^", lineage, "\\|"), names(all_results))], fill = TRUE)
  combined_leading <- rbindlist(all_leading[grepl(paste0("^", lineage, "\\|"), names(all_leading))], fill = TRUE)
  focus <- combined_result[order(padj, -abs(NES)), head(.SD, 2L), by = collection]
  effect <- merge(
    combined_leading,
    focus[, .(collection, pathway_id)],
    by = c("collection", "pathway_id")
  )
  effect <- effect[order(pathway_padj, -abs(rank_metric)), head(.SD, 12L), by = .(collection, pathway_id)]
  if (nrow(effect) > 0L) {
    effect[, gene_label := fifelse(!is.na(target_symbol) & nzchar(target_symbol), target_symbol, feature_id)]
    effect[, facet_label := paste(collection, wrap_text(pathway_name, 38), sep = ": ")]
    p <- ggplot(effect, aes(x = rank_metric, y = reorder(gene_label, rank_metric), color = pathway_direction)) +
      geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
      geom_point(size = 2) +
      facet_wrap(~facet_label, scales = "free_y", ncol = 2) +
      scale_color_manual(values = c(CIE_higher = "#B33A3A", CIE_lower = "#2E63A3")) +
      labs(title = paste(lineage, "leading-edge rank effects"), x = "Phase3 rank metric: sign(logFC) x sqrt(F)", y = NULL, color = NULL) +
      theme_bw(base_size = 9) +
      theme(plot.title = element_text(face = "bold"), panel.grid.major.y = element_blank(), legend.position = "bottom")
    ggsave(file.path(figure_dir, paste0(lineage, "_leading_edge_effect.png")), p, width = 11, height = 10, dpi = 180)
  }
}

software <- data.table(
  software = c("R", "fgsea", "msigdbr", "babelgene", "data.table", "dplyr", "ggplot2", "Reactome data", "MSigDB data", "GO ontology/annotation"),
  version = c(
    R.version.string,
    as.character(packageVersion("fgsea")),
    as.character(packageVersion("msigdbr")),
    as.character(packageVersion("babelgene")),
    as.character(packageVersion("data.table")),
    as.character(packageVersion("dplyr")),
    as.character(packageVersion("ggplot2")),
    "v97", "2026.1.Hs", "2026-08-05 / CANLF GAF 2026-07-28"
  )
)
fwrite(software, file.path(project_root, "00_metadata", "PHASE4_SOFTWARE_VERSIONS.csv"))
message("All Full7 Phase4 GSEA analyses completed.")
