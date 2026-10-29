#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
base <- file.path(project_root, "07_human_crohn", "GSE266616")
pb_root <- file.path(base, "pseudobulk")
gsea_root <- file.path(base, "gsea")
gene_set_dir <- file.path(gsea_root, "gene_sets")
primary_dir <- file.path(gsea_root, "full")
leading_dir <- file.path(gsea_root, "leading_edges")
figure_dir <- file.path(base, "figures", "gsea")
sensitivity_dir <- file.path(gsea_root, "treatment_sensitivity")
secondary_dir <- file.path(gsea_root, "secondary_paired")
metadata_dir <- file.path(gsea_root, "metadata")
log_path <- file.path(project_root, "logs", "phase7_human_gsea.log")

for (path in c(primary_dir, leading_dir, figure_dir, sensitivity_dir, secondary_dir, metadata_dir, dirname(log_path))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

lineages <- c("Epithelial", "T", "Myeloid")
collection_specs <- list(
  Hallmark = file.path(gene_set_dir, "HUMAN_HALLMARK_GENE_SETS.tsv.gz"),
  Reactome = file.path(gene_set_dir, "HUMAN_REACTOME_GENE_SETS.tsv.gz"),
  GO_BP = file.path(gene_set_dir, "HUMAN_GO_BP_GENE_SETS.tsv.gz")
)
min_size <- 15L
max_size <- 500L
fgsea_eps <- 1e-10
base_seed <- 20260825L

atomic_write_csv <- function(frame, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  fwrite(frame, tmp)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Atomic output rename failed: ", path)
}

read_rank <- function(path) {
  rank <- fread(path)
  required <- c("unique_gene_id", "gene_id", "gene_symbol", "logFC", "statistic_F", "PValue", "FDR", "rank_metric")
  if (!all(required %in% names(rank))) stop("Rank file schema failure: ", path)
  rank <- rank[!is.na(gene_symbol) & nzchar(gene_symbol) & is.finite(rank_metric)]
  rank[, abs_rank := abs(rank_metric)]
  setorder(rank, -abs_rank, PValue, unique_gene_id)
  rank <- unique(rank, by = "gene_symbol")
  rank[, abs_rank := NULL]
  setorder(rank, -rank_metric, gene_symbol)
  if (anyDuplicated(rank$gene_symbol)) stop("Symbol collapse failed: ", path)
  stats <- rank$rank_metric
  names(stats) <- rank$gene_symbol
  list(rank = rank, stats = stats)
}

load_memberships <- function(collection) {
  path <- collection_specs[[collection]]
  if (!file.exists(path)) stop("Missing Human-native gene sets: ", path)
  members <- unique(fread(path), by = c("pathway_id", "gene_symbol"))
  required <- c("collection", "pathway_id", "pathway_name", "gene_symbol", "source_version", "source_species")
  if (!all(required %in% names(members))) stop("Gene-set membership schema failure: ", path)
  if (any(members$source_species != "Homo sapiens")) stop("Non-Human source entered Phase7 GSEA: ", collection)
  members
}

select_top_both <- function(dt, n_each = 12L) {
  pos <- dt[NES > 0][order(padj, -abs(NES))][seq_len(min(.N, n_each))]
  neg <- dt[NES < 0][order(padj, -abs(NES))][seq_len(min(.N, n_each))]
  rbindlist(list(pos, neg), fill = TRUE)
}

wrap_label <- function(x, width = 48L) {
  vapply(x, function(z) paste(strwrap(gsub("_", " ", z), width = width), collapse = "\n"), character(1))
}

plot_overview <- function(result, lineage, collection, analysis_id, output) {
  show <- select_top_both(copy(result), if (collection == "Hallmark") 15L else 12L)
  if (!nrow(show)) return(invisible(NULL))
  show[, plot_label := factor(wrap_label(pathway_name), levels = wrap_label(pathway_name)[order(NES)])]
  show[, fdr_group := fifelse(padj < 0.05, "FDR < 0.05", fifelse(padj < 0.10, "FDR < 0.10", "FDR >= 0.10"))]
  p <- ggplot(show, aes(x = NES, y = plot_label, color = direction, size = -log10(pmax(padj, 1e-300)), alpha = fdr_group)) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_point() +
    scale_color_manual(values = c(ActiveCrohn_higher = "#B33A3A", ActiveCrohn_lower = "#2E63A3", Active_higher = "#B33A3A", Active_lower = "#2E63A3")) +
    scale_alpha_manual(values = c("FDR < 0.05" = 1, "FDR < 0.10" = 0.8, "FDR >= 0.10" = 0.45)) +
    labs(title = paste(lineage, collection, analysis_id), x = "Normalized enrichment score (NES)", y = NULL, size = "-log10 FDR", color = NULL, alpha = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), panel.grid.major.y = element_blank(), legend.position = "bottom")
  ggsave(output, p, width = 9.2, height = if (collection == "Hallmark") 10 else 8.5, dpi = 180)
}

run_one <- function(rank_path, lineage, collection, analysis_id, primary = FALSE, make_plot = FALSE) {
  rank_info <- read_rank(rank_path)
  members <- load_memberships(collection)
  pathway_info <- unique(members[, .(pathway_id, pathway_name, source_version, source_species)], by = "pathway_id")
  pathways <- split(members$gene_symbol, members$pathway_id)
  set.seed(base_seed + match(lineage, lineages) * 100L + match(collection, names(collection_specs)))
  fg <- as.data.table(fgseaMultilevel(
    pathways = pathways,
    stats = rank_info$stats,
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
    analysis_id = analysis_id,
    leading_edge_n = lengths(leadingEdge),
    leading_edge = vapply(leadingEdge, function(x) paste(as.character(x), collapse = ";"), character(1)),
    direction = ifelse(NES > 0, if (primary) "ActiveCrohn_higher" else "Active_higher", if (primary) "ActiveCrohn_lower" else "Active_lower")
  )]
  leading <- rbindlist(lapply(seq_len(nrow(result)), function(i) {
    genes <- as.character(result$leadingEdge[[i]])
    if (!length(genes)) return(NULL)
    detail <- rank_info$rank[match(genes, gene_symbol)]
    detail[, `:=`(
      lineage = lineage, collection = collection, analysis_id = analysis_id,
      pathway_id = result$pathway_id[[i]], pathway_name = result$pathway_name[[i]],
      pathway_NES = result$NES[[i]], pathway_padj = result$padj[[i]],
      pathway_direction = result$direction[[i]]
    )]
    detail[, .(lineage, collection, analysis_id, pathway_id, pathway_name, pathway_NES, pathway_padj,
               pathway_direction, gene_symbol, unique_gene_id, gene_id, logFC, statistic_F, PValue, FDR, rank_metric)]
  }), fill = TRUE)
  result[, leadingEdge := NULL]
  setcolorder(result, c("lineage", "collection", "analysis_id", "pathway_id", "pathway_name", "ES", "NES", "pval", "padj", "size", "leading_edge_n", "leading_edge", "direction", "source_version", "source_species", "log2err"))
  result[, abs_NES_sort := abs(NES)]
  setorder(result, padj, -abs_NES_sort, pathway_id)
  result[, abs_NES_sort := NULL]
  if (make_plot) {
    safe_analysis <- gsub("[^A-Za-z0-9]+", "_", analysis_id)
    plot_overview(result, lineage, collection, analysis_id, file.path(figure_dir, paste0("Human_", lineage, "_", collection, "_", safe_analysis, "_NES.png")))
  }
  list(result = result, leading = leading)
}

run_primary <- function() {
  results <- list()
  for (lineage in lineages) {
    rank_path <- file.path(pb_root, "ranked_lists", paste0("Human_", lineage, "_ActiveCrohn_vs_Healthy_ranked.csv"))
    if (!file.exists(rank_path)) stop("Missing primary Human rank file: ", rank_path)
    for (collection in names(collection_specs)) {
      analysis_id <- "Primary_ActiveCrohn_vs_Healthy"
      one <- run_one(rank_path, lineage, collection, analysis_id, primary = TRUE, make_plot = TRUE)
      atomic_write_csv(one$result, file.path(primary_dir, paste0("Human_", lineage, "_", collection, "_GSEA.csv")))
      atomic_write_csv(one$leading, file.path(leading_dir, paste0("Human_", lineage, "_", collection, "_leading_edges.csv")))
      results[[paste(lineage, collection, sep = "::")]] <- one$result
      message("Primary ", lineage, " / ", collection, ": tested=", nrow(one$result), " FDR05=", sum(one$result$padj < 0.05))
    }
  }
  results
}

run_generic_ranks <- function(input_dir, output_dir, pattern, analysis_prefix) {
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  results <- list()
  if (!length(files)) return(results)
  for (path in files) {
    name <- basename(path)
    match <- regexec("^Human_(Epithelial|T|Myeloid)_(.+)_ranked\\.csv$", name)
    parts <- regmatches(name, match)[[1]]
    if (length(parts) != 3L) next
    lineage <- parts[[2]]
    comparison <- parts[[3]]
    analysis_id <- paste0(analysis_prefix, comparison)
    for (collection in names(collection_specs)) {
      one <- run_one(path, lineage, collection, analysis_id, primary = FALSE, make_plot = FALSE)
      out <- file.path(output_dir, paste0("Human_", lineage, "_", comparison, "_", collection, "_GSEA.csv"))
      atomic_write_csv(one$result, out)
      results[[paste(lineage, comparison, collection, sep = "::")]] <- one$result
    }
  }
  results
}

treatment_label <- function(primary_row, sensitivity_results) {
  matches <- rbindlist(lapply(sensitivity_results, function(x) x[
    lineage == primary_row$lineage & collection == primary_row$collection & pathway_id == primary_row$pathway_id
  ]), fill = TRUE)
  if (!nrow(matches)) return("NOT_EXECUTED")
  expected_sign <- sign(primary_row$NES)
  same <- sign(matches$NES) == expected_sign
  sig_same <- same & matches$padj < 0.05
  if (all(same) && all(sig_same)) return(paste0("ALL_", nrow(matches), "_FDR05_SAME_DIRECTION"))
  if (all(same)) return(paste0("ALL_", nrow(matches), "_SAME_DIRECTION_DIRECTIONALLY_STABLE"))
  paste0("DIRECTION_CHANGE_IN_", sum(!same), "_OF_", nrow(matches), "_SENSITIVITIES")
}

freeze_programs <- function(primary_results, sensitivity_results) {
  combined <- rbindlist(primary_results, fill = TRUE)
  programs <- combined[padj < 0.05]
  programs[, treatment_sensitivity := vapply(seq_len(.N), function(i) treatment_label(programs[i], sensitivity_results), character(1))]
  programs[, `:=`(
    program = pathway_name,
    pathway = pathway_id,
    FDR = padj,
    primary_or_secondary = "PRIMARY",
    notes = "Unbiased Human-native GSEA; frozen before PHASE6_HUMAN_FORWARD_LIST was read"
  )]
  frozen <- programs[, .(lineage, program, collection, pathway, NES, FDR, direction, leading_edge, treatment_sensitivity, primary_or_secondary, notes)]
  frozen[, abs_NES_sort := abs(NES)]
  setorder(frozen, lineage, FDR, -abs_NES_sort, collection, pathway)
  frozen[, abs_NES_sort := NULL]
  atomic_write_csv(frozen, file.path(gsea_root, "HUMAN_DISEASE_PROGRAMS.csv"))
  frozen
}

write_gate <- function(primary_results, programs, sensitivity_results, secondary_results) {
  combined <- rbindlist(primary_results, fill = TRUE)
  criteria <- data.table(
    criterion = c(
      "Nine unbiased primary lineage-by-collection GSEA results exist",
      "Human-native gene space only",
      "Both NES directions retained in full results",
      "At least one FDR<0.05 program per primary lineage",
      "Human program table frozen before Dog-forward observation",
      "Treatment sensitivity GSEA executed",
      "Secondary paired GSEA handled according to eligible rank availability"
    ),
    pass = c(
      length(primary_results) == 9L,
      all(combined$source_species == "Homo sapiens"),
      all(combined[, any(NES > 0) && any(NES < 0), by = lineage]$V1),
      all(programs[, .N, by = lineage]$N > 0) && uniqueN(programs$lineage) == 3L,
      file.exists(file.path(gsea_root, "HUMAN_DISEASE_PROGRAMS.csv")),
      length(sensitivity_results) > 0L,
      TRUE
    ),
    detail = c(
      paste0("results=", length(primary_results)),
      paste(unique(combined$source_species), collapse = ";"),
      paste(combined[, paste0(lineage, "=+", sum(NES > 0), "/-", sum(NES < 0)), by = lineage]$V1, collapse = ";"),
      paste(programs[, paste0(lineage, "=", .N), by = lineage]$V1, collapse = ";"),
      "HUMAN_DISEASE_PROGRAMS.csv",
      paste0("sensitivity lineage-collection results=", length(sensitivity_results)),
      paste0("secondary lineage-collection results=", length(secondary_results))
    )
  )
  criteria[, status := ifelse(pass, "PASS", "FAIL")]
  atomic_write_csv(criteria, file.path(metadata_dir, "PHASE7C_GSEA_GATE.csv"))
  if (!all(criteria$pass)) stop("Phase7C GSEA gate failed")
}

main <- function() {
  sink(log_path, split = TRUE)
  on.exit(sink(), add = TRUE)
  cat("Phase7 Human-native GSEA\nStart:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
  primary_results <- run_primary()
  sensitivity_results <- run_generic_ranks(
    file.path(pb_root, "sensitivity", "ranked_lists"), sensitivity_dir,
    "^Human_(Epithelial|T|Myeloid)_.+_ranked\\.csv$", "TreatmentSensitivity_"
  )
  secondary_results <- run_generic_ranks(
    file.path(pb_root, "paired_ranked_lists"), secondary_dir,
    "^Human_(Epithelial|T|Myeloid)_.+_ranked\\.csv$", "SecondaryPaired_"
  )
  programs <- freeze_programs(primary_results, sensitivity_results)
  write_gate(primary_results, programs, sensitivity_results, secondary_results)
  software <- data.table(
    software = c("R", "fgsea", "data.table", "ggplot2"),
    version = c(as.character(getRversion()), as.character(packageVersion("fgsea")), as.character(packageVersion("data.table")), as.character(packageVersion("ggplot2"))),
    role = c("runtime", "multilevel pre-ranked GSEA", "tabular I/O", "figures")
  )
  atomic_write_csv(software, file.path(project_root, "00_metadata", "PHASE7_GSEA_SOFTWARE_VERSIONS.csv"))
  cat("End:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
}

status <- 0L
tryCatch(main(), error = function(error) {
  message("FATAL: ", conditionMessage(error))
  status <<- 1L
})
quit(save = "no", status = status)
