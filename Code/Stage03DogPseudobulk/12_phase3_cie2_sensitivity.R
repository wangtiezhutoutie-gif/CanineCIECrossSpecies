#!/usr/bin/env Rscript

# Pre-specified targeted sensitivity analysis excluding cie_2 only.

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
pb_root <- file.path(project_root, "06_pseudobulk")
matrix_dir <- file.path(pb_root, "matrices")
metadata_dir <- file.path(pb_root, "metadata")
de_dir <- file.path(pb_root, "de_results")
sensitivity_dir <- file.path(pb_root, "sensitivity")
figure_dir <- file.path(pb_root, "figures")
log_path <- file.path(project_root, "logs", "phase3_cie2_sensitivity.log")
dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)

log_connection <- file(log_path, open = "wt", encoding = "UTF-8")
sink(log_connection, type = "output", split = TRUE)
sink(log_connection, type = "message")

cleanup_log <- function() {
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(log_connection), silent = TRUE)
}

atomic_write_csv <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  write.csv(x, temporary, row.names = FALSE, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) stop("Could not atomically write: ", path)
}

open_png <- function(path, width = 1500, height = 1100) png(path, width = width, height = height, res = 150)

main <- function() {
  set.seed(20240824)
  start_time <- Sys.time()
  suppressPackageStartupMessages(library(edgeR))
  cat("Phase3 targeted cie_2 sensitivity\n")
  cat("Start:", format(start_time, tz = "Asia/Shanghai"), "\n")
  cat("edgeR:", as.character(packageVersion("edgeR")), "\n")
  cat("Random seed: 20240824 (no stochastic inferential step)\n")
  cat("Primary remains Full7; this model temporarily excludes only GSM8031701_cie_2\n")
  cat("Same workflow and design: ~condition, Healthy reference, conditionCIE coefficient\n")

  sample_metadata <- read.csv(file.path(metadata_dir, "pseudobulk_sample_metadata.csv"), check.names = FALSE)
  annotation <- read.csv(file.path(metadata_dir, "pseudobulk_gene_annotation.csv"), check.names = FALSE)
  rownames(annotation) <- annotation$feature_id
  lineages <- c("Epithelial", "T", "Myeloid")
  gene_level_rows <- list()
  summary_rows <- list()
  threshold <- 0.25
  top_n <- 100L
  excluded_donor <- "GSM8031701_cie_2"

  for (lineage in lineages) {
    cat("\n===", lineage, "===\n")
    count_frame <- read.csv(gzfile(file.path(matrix_dir, paste0("pseudobulk_", lineage, "_counts.csv.gz"))), check.names = FALSE)
    feature_ids <- count_frame$feature_id
    counts <- as.matrix(count_frame[, -1, drop = FALSE])
    storage.mode(counts) <- "integer"
    rownames(counts) <- feature_ids
    counts <- counts[, colnames(counts) != excluded_donor, drop = FALSE]
    eligible_flag <- tolower(as.character(sample_metadata$eligible)) == "true"
    meta <- sample_metadata[sample_metadata$lineage == lineage & eligible_flag, , drop = FALSE]
    meta <- meta[meta$donor_id != excluded_donor, , drop = FALSE]
    meta <- meta[match(colnames(counts), meta$donor_id), , drop = FALSE]
    condition <- factor(meta$condition, levels = c("Healthy", "CIE"))
    if (nrow(meta) != 6 || sum(condition == "Healthy") != 3 || sum(condition == "CIE") != 3) stop("No-cie2 sample mismatch: ", lineage)
    design <- model.matrix(~ condition)
    rownames(design) <- meta$donor_id
    cat("Design matrix:\n")
    print(design)

    y0 <- DGEList(counts = counts, samples = data.frame(donor_id = meta$donor_id, condition = condition, row.names = meta$donor_id))
    keep <- filterByExpr(y0, design = design)
    y <- y0[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y, method = "TMM")
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE)
    qlf <- glmQLFTest(fit, coef = "conditionCIE")
    tab <- topTags(qlf, n = Inf, sort.by = "none")$table
    no_cie2 <- data.frame(
      feature_id = rownames(tab),
      gene_symbol = annotation[rownames(tab), "original_gene_symbol"],
      unique_var_name = annotation[rownames(tab), "unique_var_name"],
      logFC = tab$logFC,
      logCPM = tab$logCPM,
      F = tab$F,
      PValue = tab$PValue,
      FDR = tab$FDR,
      direction = ifelse(tab$logFC > 0, "CIE_higher", ifelse(tab$logFC < 0, "CIE_lower", "No_change")),
      lineage = lineage,
      n_healthy = 3L,
      n_cie = 3L,
      stringsAsFactors = FALSE
    )
    no_cie2 <- no_cie2[order(no_cie2$PValue, -abs(no_cie2$logFC)), ]
    atomic_write_csv(no_cie2, file.path(sensitivity_dir, paste0(lineage, "_no_cie2_edger.csv")))

    full <- read.csv(file.path(de_dir, paste0(lineage, "_CIE_vs_Healthy_edger.csv")), check.names = FALSE)
    common <- merge(
      full[, c("feature_id", "gene_symbol", "unique_var_name", "logFC", "F", "PValue", "FDR")],
      no_cie2[, c("feature_id", "logFC", "F", "PValue", "FDR")],
      by = "feature_id", suffixes = c("_full7", "_no_cie2"), sort = FALSE
    )
    common$lineage <- lineage
    common$direction_concordant <- sign(common$logFC_full7) == sign(common$logFC_no_cie2)
    common$meaningful_full7 <- abs(common$logFC_full7) >= threshold
    common$meaningful_no_cie2 <- abs(common$logFC_no_cie2) >= threshold
    common$meaningful_in_both <- common$meaningful_full7 & common$meaningful_no_cie2
    common$meaningful_direction_concordant <- common$direction_concordant & common$meaningful_in_both
    common$large_effect_sign_flip <- common$meaningful_in_both & !common$direction_concordant
    common$delta_logFC_no_cie2_minus_full7 <- common$logFC_no_cie2 - common$logFC_full7
    gene_level_rows[[lineage]] <- common

    top_full <- head(common$feature_id[order(abs(common$logFC_full7), decreasing = TRUE)], top_n)
    top_no <- head(common$feature_id[order(abs(common$logFC_no_cie2), decreasing = TRUE)], top_n)
    top_overlap <- length(intersect(top_full, top_no))
    both_meaningful <- common$meaningful_in_both
    meaningful_concordance <- if (sum(both_meaningful) > 0) mean(common$direction_concordant[both_meaningful]) else NA_real_
    summary_rows[[lineage]] <- data.frame(
      lineage = lineage,
      full7_tested_genes = nrow(full),
      no_cie2_tested_genes = nrow(no_cie2),
      common_tested_genes = nrow(common),
      logFC_spearman = cor(common$logFC_full7, common$logFC_no_cie2, method = "spearman"),
      logFC_pearson = cor(common$logFC_full7, common$logFC_no_cie2, method = "pearson"),
      direction_concordance = mean(common$direction_concordant),
      meaningful_effect_threshold = threshold,
      genes_meaningful_in_both = sum(both_meaningful),
      meaningful_direction_concordance = meaningful_concordance,
      top_effect_n = top_n,
      top_effect_overlap_n = top_overlap,
      top_effect_overlap_fraction = top_overlap / top_n,
      large_effect_sign_flips = sum(common$large_effect_sign_flip),
      median_abs_delta_logFC = median(abs(common$delta_logFC_no_cie2_minus_full7)),
      max_abs_delta_logFC = max(abs(common$delta_logFC_no_cie2_minus_full7)),
      stringsAsFactors = FALSE
    )

    point_colours <- ifelse(common$large_effect_sign_flip, "#B3261E", ifelse(common$meaningful_in_both, "#5C78A8", "#BBBBBB"))
    open_png(file.path(figure_dir, paste0(lineage, "_cie2_sensitivity_logFC.png")))
    plot(common$logFC_full7, common$logFC_no_cie2, pch = 19, cex = 0.42, col = point_colours,
         xlab = "Full7 logFC (CIE - Healthy)", ylab = "No-cie2 logFC (CIE - Healthy)",
         main = paste(lineage, "targeted cie_2 sensitivity"))
    abline(a = 0, b = 1, lty = 2, col = "black")
    abline(h = c(-threshold, threshold), v = c(-threshold, threshold), lty = 3, col = "grey55")
    legend("topleft", legend = c("Both |logFC| >= 0.25", "Large-effect sign flip", "Other"),
           col = c("#5C78A8", "#B3261E", "#BBBBBB"), pch = 19, bty = "n", cex = 0.8)
    dev.off()

    cat("Genes full/no/common:", nrow(full), "/", nrow(no_cie2), "/", nrow(common), "\n")
    print(summary_rows[[lineage]])
  }

  gene_level <- do.call(rbind, gene_level_rows)
  rownames(gene_level) <- NULL
  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL
  atomic_write_csv(gene_level, file.path(sensitivity_dir, "cie2_sensitivity_gene_level.csv"))
  atomic_write_csv(summary, file.path(sensitivity_dir, "cie2_sensitivity_summary.csv"))
  cat("\nThreshold 0.25 is a pre-specified internal descriptive effect threshold, not a significance cutoff.\n")
  cat("No other donor was excluded; complete LODO was not run.\n")
  cat("End:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
  cat("Elapsed seconds:", as.numeric(difftime(Sys.time(), start_time, units = "secs")), "\n")
}

status <- 0L
tryCatch(
  main(),
  error = function(e) {
    cat("FATAL ERROR:", conditionMessage(e), "\n")
    traceback(2)
    status <<- 1L
  }
)
cleanup_log()
quit(save = "no", status = status)
