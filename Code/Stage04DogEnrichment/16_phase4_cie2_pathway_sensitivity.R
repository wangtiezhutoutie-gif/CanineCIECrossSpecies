#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
gsea_root <- file.path(project_root, "06_pseudobulk", "gsea")
ranked_dir <- file.path(project_root, "06_pseudobulk", "ranked_lists")
gene_sens_dir <- file.path(project_root, "06_pseudobulk", "sensitivity")
mapping_dir <- file.path(gsea_root, "mapping")
gene_set_dir <- file.path(gsea_root, "gene_sets")
sens_dir <- file.path(gsea_root, "sensitivity")
common_out <- file.path(sens_dir, "common_universe_gsea")
figure_dir <- file.path(gsea_root, "figures")
dir.create(sens_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(common_out, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

lineages <- c("Epithelial", "T", "Myeloid")
collection_specs <- list(
  Hallmark = list(file = file.path(gene_set_dir, "hallmark_dog_ortholog_projected_membership.tsv.gz"), tag = "Hallmark"),
  Reactome = list(file = file.path(gene_set_dir, "reactome_canis_familiaris_gene_sets.tsv.gz"), tag = "Reactome"),
  GO_BP = list(file = file.path(gene_set_dir, "go_bp_canis_familiaris_gene_sets.tsv.gz"), tag = "GOBP")
)
min_size <- 15L
max_size <- 500L
fgsea_eps <- 1e-10
base_seed <- 20260824L

mapping <- fread(file.path(mapping_dir, "phase4_feature_to_geneset_id.csv"))
mapping <- mapping[!is.na(target_geneset_id) & nzchar(as.character(target_geneset_id)), .(
  feature_id, target_geneset_id = as.character(target_geneset_id), target_symbol
)]

flatten_result <- function(fg, info, lineage, collection, model) {
  result <- merge(fg, info, by = "pathway_id", all.x = TRUE, sort = FALSE)
  result[, `:=`(
    lineage = lineage,
    collection = collection,
    model = model,
    padj = p.adjust(pval, method = "BH"),
    direction = fifelse(NES > 0, "CIE_higher", "CIE_lower"),
    leading_edge_n = lengths(leadingEdge),
    leading_edge_target_ids = vapply(leadingEdge, function(x) paste(as.character(x), collapse = ";"), character(1))
  )]
  result[, leadingEdge := NULL]
  setorder(result, padj, pathway_id)
  result
}

jaccard <- function(a, b) {
  aa <- if (is.na(a) || !nzchar(a)) character() else unique(strsplit(a, ";", fixed = TRUE)[[1L]])
  bb <- if (is.na(b) || !nzchar(b)) character() else unique(strsplit(b, ";", fixed = TRUE)[[1L]])
  union_n <- length(union(aa, bb))
  if (union_n == 0L) return(NA_real_)
  length(intersect(aa, bb)) / union_n
}

top20_jaccard <- function(dt) {
  full <- dt[!is.na(NES_full)][order(padj_full, -abs(NES_full), pathway_id), head(pathway_id, 20L)]
  no <- dt[!is.na(NES_no_cie2)][order(padj_no_cie2, -abs(NES_no_cie2), pathway_id), head(pathway_id, 20L)]
  if (length(union(full, no)) == 0L) return(NA_real_)
  length(intersect(full, no)) / length(union(full, no))
}

all_comparisons <- list()
all_summaries <- list()

for (lineage_index in seq_along(lineages)) {
  lineage <- lineages[[lineage_index]]
  full <- fread(file.path(ranked_dir, paste0(lineage, "_ranked.csv")))
  no <- fread(file.path(gene_sens_dir, paste0(lineage, "_no_cie2_edger.csv")))
  no[, rank_metric_no := sign(logFC) * sqrt(F)]
  setnames(
    full,
    c("gene_symbol", "unique_var_name", "logFC", "statistic_F", "PValue", "FDR", "rank_metric"),
    c("gene_symbol_full", "unique_var_name_full", "logFC_full", "F_full", "PValue_full", "FDR_full", "rank_metric_full")
  )
  setnames(
    no,
    c("gene_symbol", "unique_var_name", "logFC", "F", "PValue", "FDR"),
    c("gene_symbol_no", "unique_var_name_no", "logFC_no", "F_no", "PValue_no", "FDR_no")
  )
  common <- merge(full, no, by = "feature_id", all = FALSE)
  common_features_n <- nrow(common)
  common <- merge(common, mapping, by = "feature_id", all = FALSE)
  common[, abs_rank_full := abs(rank_metric_full)]
  setorder(common, target_geneset_id, -abs_rank_full, PValue_full, feature_id)
  common[, collision_rank := seq_len(.N), by = target_geneset_id]
  common[, kept_after_collision := collision_rank == 1L]
  collapse_audit <- common[,.N, by = target_geneset_id][N > 1L]
  common_kept <- common[kept_after_collision == TRUE]
  setorder(common_kept, -rank_metric_full, feature_id)
  fwrite(common, file.path(sens_dir, paste0(lineage, "_common_universe_rank_mapping.csv")))

  full_stats <- common_kept$rank_metric_full
  no_stats <- common_kept$rank_metric_no
  names(full_stats) <- common_kept$target_geneset_id
  names(no_stats) <- common_kept$target_geneset_id
  stopifnot(identical(names(full_stats), names(no_stats)), !anyDuplicated(names(full_stats)))
  message(
    lineage, ": Full7 features=", nrow(full), ", No-cie2 features=", nrow(no),
    ", common features=", common_features_n, ", common mapped target genes=", length(full_stats),
    ", collision groups=", nrow(collapse_audit), "."
  )

  lineage_comparisons <- list()
  for (collection_index in seq_along(collection_specs)) {
    collection <- names(collection_specs)[[collection_index]]
    spec <- collection_specs[[collection]]
    members <- unique(fread(spec$file), by = c("pathway_id", "target_geneset_id"))
    members[, target_geneset_id := as.character(target_geneset_id)]
    info <- unique(members[, .(pathway_id, pathway_name, mapping_type, source_version)], by = "pathway_id")
    pathways <- split(members$target_geneset_id, members$pathway_id)
    seed <- base_seed + lineage_index * 100L + collection_index
    message("Running common-universe sensitivity: ", lineage, " / ", collection, " / Full7...")
    set.seed(seed)
    fg_full <- as.data.table(fgseaMultilevel(
      pathways = pathways, stats = full_stats, minSize = min_size, maxSize = max_size,
      eps = fgsea_eps, scoreType = "std", nproc = 1L
    ))
    setnames(fg_full, "pathway", "pathway_id")
    full_result <- flatten_result(fg_full, info, lineage, collection, "Full7_common")
    message("Running common-universe sensitivity: ", lineage, " / ", collection, " / No-cie2...")
    set.seed(seed)
    fg_no <- as.data.table(fgseaMultilevel(
      pathways = pathways, stats = no_stats, minSize = min_size, maxSize = max_size,
      eps = fgsea_eps, scoreType = "std", nproc = 1L
    ))
    setnames(fg_no, "pathway", "pathway_id")
    no_result <- flatten_result(fg_no, info, lineage, collection, "No_cie2_common")

    collection_dir <- file.path(common_out, tolower(collection))
    dir.create(collection_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(full_result, file.path(collection_dir, paste0(lineage, "_", spec$tag, "_Full7_common_GSEA.csv")))
    fwrite(no_result, file.path(collection_dir, paste0(lineage, "_", spec$tag, "_No_cie2_common_GSEA.csv")))

    full_cmp <- full_result[, .(
      pathway_id, pathway_name, mapping_type, source_version,
      ES_full = ES, NES_full = NES, pval_full = pval, padj_full = padj,
      size_full = size, leading_edge_n_full = leading_edge_n,
      leading_edge_full = leading_edge_target_ids, direction_full = direction
    )]
    no_cmp <- no_result[, .(
      pathway_id, pathway_name_no = pathway_name,
      ES_no_cie2 = ES, NES_no_cie2 = NES, pval_no_cie2 = pval, padj_no_cie2 = padj,
      size_no_cie2 = size, leading_edge_n_no_cie2 = leading_edge_n,
      leading_edge_no_cie2 = leading_edge_target_ids, direction_no_cie2 = direction
    )]
    cmp <- merge(full_cmp, no_cmp, by = "pathway_id", all = TRUE)
    cmp[is.na(pathway_name), pathway_name := pathway_name_no]
    cmp[, pathway_name_no := NULL]
    cmp[, `:=`(
      lineage = lineage,
      collection = collection,
      delta_NES = NES_no_cie2 - NES_full,
      same_direction = !is.na(NES_full) & !is.na(NES_no_cie2) & sign(NES_full) == sign(NES_no_cie2),
      leading_edge_jaccard = mapply(jaccard, leading_edge_full, leading_edge_no_cie2),
      common_feature_universe_n = common_features_n,
      common_target_gene_universe_n = length(full_stats),
      pathway_size_change = size_no_cie2 - size_full
    )]
    cmp[, sensitivity_flag := fifelse(
      is.na(NES_full) | is.na(NES_no_cie2), "Not_assessable",
      fifelse(
        !same_direction, "Direction_flip",
        fifelse(abs(delta_NES) >= 0.5 | (!is.na(leading_edge_jaccard) & leading_edge_jaccard < 0.25),
                "Sensitive_to_cie2", "Stable_to_cie2")
      )
    )]
    setcolorder(cmp, c(
      "lineage", "collection", "pathway_id", "pathway_name", "NES_full", "padj_full",
      "NES_no_cie2", "padj_no_cie2", "delta_NES", "same_direction",
      "leading_edge_jaccard", "sensitivity_flag"
    ))
    cmp[, abs_delta_NES_sort := abs(delta_NES)]
    setorder(cmp, sensitivity_flag, -abs_delta_NES_sort, pathway_id)
    cmp[, abs_delta_NES_sort := NULL]
    all_comparisons[[paste(lineage, collection, sep = "|")]] <- cmp
    lineage_comparisons[[collection]] <- cmp

    comparable <- cmp[!is.na(NES_full) & !is.na(NES_no_cie2)]
    sig05 <- comparable[padj_full < 0.05]
    sig10 <- comparable[padj_full < 0.10]
    all_summaries[[paste(lineage, collection, sep = "|")]] <- data.table(
      lineage = lineage,
      collection = collection,
      full7_tested_features = nrow(full),
      no_cie2_tested_features = nrow(no),
      common_feature_universe_n = common_features_n,
      common_mapped_target_genes_n = length(full_stats),
      common_collapsed_features_n = nrow(common) - nrow(common_kept),
      pathways_full = nrow(full_result),
      pathways_no_cie2 = nrow(no_result),
      pathways_comparable = nrow(comparable),
      pathways_not_assessable = cmp[sensitivity_flag == "Not_assessable", .N],
      nes_spearman = if (nrow(comparable) > 2L) cor(comparable$NES_full, comparable$NES_no_cie2, method = "spearman") else NA_real_,
      nes_pearson = if (nrow(comparable) > 2L) cor(comparable$NES_full, comparable$NES_no_cie2, method = "pearson") else NA_real_,
      direction_concordance = mean(comparable$same_direction),
      full_sig05_n = nrow(sig05),
      full_sig05_direction_retained_fraction = if (nrow(sig05)) mean(sig05$same_direction) else NA_real_,
      full_sig05_fdr05_and_direction_retained_fraction = if (nrow(sig05)) mean(sig05$same_direction & sig05$padj_no_cie2 < 0.05) else NA_real_,
      full_sig10_n = nrow(sig10),
      full_sig10_direction_retained_fraction = if (nrow(sig10)) mean(sig10$same_direction) else NA_real_,
      full_sig10_fdr10_and_direction_retained_fraction = if (nrow(sig10)) mean(sig10$same_direction & sig10$padj_no_cie2 < 0.10) else NA_real_,
      top20_pathway_jaccard = top20_jaccard(comparable),
      median_abs_delta_NES = median(abs(comparable$delta_NES)),
      max_abs_delta_NES = max(abs(comparable$delta_NES)),
      median_leading_edge_jaccard = median(comparable$leading_edge_jaccard, na.rm = TRUE),
      stable_n = comparable[sensitivity_flag == "Stable_to_cie2", .N],
      sensitive_n = comparable[sensitivity_flag == "Sensitive_to_cie2", .N],
      direction_flip_n = comparable[sensitivity_flag == "Direction_flip", .N],
      pathway_size_change_max_abs = max(abs(comparable$pathway_size_change), na.rm = TRUE),
      note = if (lineage == "Myeloid") "cie_2 Myeloid has 27 cells; interpret sensitivity under low-depth limitation." else "Common-universe sensitivity only; Full7 remains primary."
    )
    message("Completed sensitivity ", lineage, " / ", collection, ": comparable=", nrow(comparable), ", flips=", comparable[sensitivity_flag == "Direction_flip", .N], ".")
  }

  plot_dt <- rbindlist(lineage_comparisons, fill = TRUE)
  plot_dt <- plot_dt[!is.na(NES_full) & !is.na(NES_no_cie2)]
  p <- ggplot(plot_dt, aes(x = NES_full, y = NES_no_cie2, color = sensitivity_flag)) +
    geom_abline(slope = 1, intercept = 0, color = "grey65", linewidth = 0.4) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_point(alpha = 0.7, size = 1.3) +
    facet_wrap(~collection) +
    scale_color_manual(values = c(
      Stable_to_cie2 = "#2E7D62", Sensitive_to_cie2 = "#D28E2D",
      Direction_flip = "#B33A3A", Not_assessable = "#777777"
    )) +
    coord_equal() +
    labs(
      title = paste(lineage, "pathway sensitivity to cie_2"),
      subtitle = "Both models use the identical mapped-gene universe",
      x = "Full7 common-universe NES", y = "No-cie2 common-universe NES", color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
  ggsave(file.path(figure_dir, paste0(lineage, "_Full7_vs_No_cie2_NES_scatter.png")), p, width = 10, height = 4.5, dpi = 180)
}

comparisons <- rbindlist(all_comparisons, fill = TRUE)
summaries <- rbindlist(all_summaries, fill = TRUE)
fwrite(comparisons, file.path(sens_dir, "cie2_pathway_sensitivity.csv"))
fwrite(summaries, file.path(sens_dir, "cie2_pathway_sensitivity_summary.csv"))

parameter_path <- file.path(project_root, "00_metadata", "PHASE4_GSEA_PARAMETERS.csv")
parameters <- fread(parameter_path)
extra <- data.table(
  parameter = c("delta_NES_definition", "top_pathway_overlap_definition", "leading_edge_overlap_definition"),
  value = c(
    "NES_no_cie2 - NES_full on the identical common mapped-gene universe",
    "Jaccard of top 20 pathways ordered by padj then abs(NES)",
    "Jaccard of NCBI Gene IDs in paired leading edges"
  )
)
parameters <- unique(rbindlist(list(parameters, extra), fill = TRUE), by = "parameter")
fwrite(parameters, parameter_path)
message("All cie_2 common-universe pathway sensitivity analyses completed.")
