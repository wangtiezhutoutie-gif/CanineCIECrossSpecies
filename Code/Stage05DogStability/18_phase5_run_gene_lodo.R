#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
pb_root <- file.path(project_root, "06_pseudobulk")
matrix_dir <- file.path(pb_root, "matrices")
pb_meta_dir <- file.path(pb_root, "metadata")
lodo_root <- file.path(pb_root, "lodo")
gene_root <- file.path(lodo_root, "gene_level")
meta_root <- file.path(lodo_root, "metadata")
qc_figure_root <- file.path(lodo_root, "figures", "qc")
log_path <- file.path(project_root, "logs", "phase5_gene_lodo.log")
for (path in c(gene_root, meta_root, qc_figure_root, dirname(log_path))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

log_connection <- file(log_path, open = "wt", encoding = "UTF-8")
sink(log_connection, type = "output", split = TRUE)
sink(log_connection, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(log_connection), silent = TRUE)
}, add = TRUE)

atomic_write_csv <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  write.csv(x, temporary, row.names = FALSE, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) stop("Could not atomically write: ", path)
}

lineages <- c("Epithelial", "T", "Myeloid")
donor_table <- data.frame(
  donor_id = c(
    "GSM8031698_h_1", "GSM8031699_h_2", "GSM8031700_h_3",
    "GSM8031701_cie_2", "GSM8031702_cie_3", "GSM8031703_cie_4", "GSM8031704_cie_6"
  ),
  lodo_run = c("LODO_H1", "LODO_H2", "LODO_H3", "LODO_CIE2", "LODO_CIE3", "LODO_CIE4", "LODO_CIE6"),
  short_label = c("drop_h1", "drop_h2", "drop_h3", "drop_cie2", "drop_cie3", "drop_cie4", "drop_cie6"),
  condition = c("Healthy", "Healthy", "Healthy", "CIE", "CIE", "CIE", "CIE"),
  stringsAsFactors = FALSE
)

sample_metadata <- read.csv(file.path(pb_meta_dir, "pseudobulk_sample_metadata.csv"), check.names = FALSE)
annotation <- read.csv(file.path(pb_meta_dir, "pseudobulk_gene_annotation.csv"), check.names = FALSE)
if (anyDuplicated(annotation$feature_id)) stop("Duplicated feature IDs in frozen pseudobulk annotation")
rownames(annotation) <- annotation$feature_id

status_rows <- list()
normalization_rows <- list()
qc_rows <- list()

cat("Phase5 donor-level gene LODO\n")
cat("Start:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
cat("R:", R.version.string, "\n")
cat("edgeR:", as.character(packageVersion("edgeR")), "\n")
cat("limma:", as.character(packageVersion("limma")), "\n")
cat("Model: ~condition; Healthy reference; robust QL; every run refits from raw counts\n")

for (lineage in lineages) {
  matrix_path <- file.path(matrix_dir, paste0("pseudobulk_", lineage, "_counts.csv.gz"))
  count_frame <- read.csv(gzfile(matrix_path), check.names = FALSE)
  feature_ids <- count_frame$feature_id
  counts_all <- as.matrix(count_frame[, -1, drop = FALSE])
  storage.mode(counts_all) <- "integer"
  rownames(counts_all) <- feature_ids
  eligible <- tolower(as.character(sample_metadata$eligible)) == "true"
  lineage_meta <- sample_metadata[sample_metadata$lineage == lineage & eligible, , drop = FALSE]
  if (nrow(lineage_meta) != 7L || !setequal(colnames(counts_all), lineage_meta$donor_id)) {
    stop("Frozen pseudobulk donor mismatch: ", lineage)
  }
  lineage_dir <- file.path(gene_root, lineage)
  ranked_dir <- file.path(lineage_dir, "ranked")
  dir.create(lineage_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ranked_dir, recursive = TRUE, showWarnings = FALSE)

  for (drop_index in seq_len(nrow(donor_table))) {
    drop <- donor_table[drop_index, ]
    run_start <- Sys.time()
    status_key <- paste(lineage, drop$lodo_run, sep = "|")
    cat("\n===", lineage, drop$lodo_run, "drop", drop$donor_id, "===\n")
    run_status <- tryCatch({
      keep_donors <- setdiff(colnames(counts_all), drop$donor_id)
      counts <- counts_all[, keep_donors, drop = FALSE]
      meta <- lineage_meta[match(keep_donors, lineage_meta$donor_id), , drop = FALSE]
      if (any(is.na(meta$donor_id))) stop("Metadata alignment failed")
      condition <- factor(meta$condition, levels = c("Healthy", "CIE"))
      n_healthy <- sum(condition == "Healthy")
      n_cie <- sum(condition == "CIE")
      if (!((n_healthy == 2L && n_cie == 4L) || (n_healthy == 3L && n_cie == 3L))) {
        stop("Unexpected remaining group sizes")
      }
      design <- model.matrix(~condition)
      rownames(design) <- meta$donor_id

      y0 <- DGEList(
        counts = counts,
        samples = data.frame(donor_id = meta$donor_id, condition = condition, row.names = meta$donor_id)
      )
      keep_gene <- filterByExpr(y0, design = design)
      if (sum(keep_gene) < 1000L) stop("Too few genes after filterByExpr: ", sum(keep_gene))
      y <- y0[keep_gene, , keep.lib.sizes = FALSE]
      y <- calcNormFactors(y, method = "TMM")
      y <- estimateDisp(y, design, robust = TRUE)
      fit <- glmQLFit(y, design, robust = TRUE)
      if (!all(is.finite(fit$coefficients))) stop("Non-finite edgeR coefficients")
      qlf <- glmQLFTest(fit, coef = "conditionCIE")
      raw <- topTags(qlf, n = Inf, sort.by = "none")$table
      if (nrow(raw) != sum(keep_gene)) stop("edgeR result row mismatch")

      result <- data.frame(
        feature_id = rownames(raw),
        gene_symbol = annotation[rownames(raw), "original_gene_symbol"],
        unique_var_name = annotation[rownames(raw), "unique_var_name"],
        logFC = raw$logFC,
        logCPM = raw$logCPM,
        F = raw$F,
        PValue = raw$PValue,
        FDR = raw$FDR,
        direction = ifelse(raw$logFC > 0, "CIE_higher", ifelse(raw$logFC < 0, "CIE_lower", "No_change")),
        rank_metric = sign(raw$logFC) * sqrt(pmax(raw$F, 0)),
        lineage = lineage,
        lodo_run = drop$lodo_run,
        short_label = drop$short_label,
        dropped_donor = drop$donor_id,
        dropped_condition = drop$condition,
        n_healthy = n_healthy,
        n_cie = n_cie,
        stringsAsFactors = FALSE
      )
      result <- result[order(result$PValue, -abs(result$logFC), result$feature_id), ]
      out_path <- file.path(lineage_dir, paste0(lineage, "_drop_", drop$donor_id, ".csv"))
      atomic_write_csv(result, out_path)
      ranked <- result[order(result$rank_metric, decreasing = TRUE, na.last = TRUE), c(
        "feature_id", "gene_symbol", "unique_var_name", "logFC", "F", "PValue", "FDR",
        "rank_metric", "lineage", "lodo_run", "short_label", "dropped_donor"
      )]
      names(ranked)[names(ranked) == "F"] <- "statistic_F"
      rank_path <- file.path(ranked_dir, paste0(lineage, "_drop_", drop$donor_id, "_ranked.csv"))
      atomic_write_csv(ranked, rank_path)

      norm <- data.frame(
        lineage = lineage, lodo_run = drop$lodo_run, short_label = drop$short_label,
        dropped_donor = drop$donor_id, donor_id = colnames(y), condition = as.character(condition),
        raw_library_size = y$samples$lib.size, normalization_factor = y$samples$norm.factors,
        effective_library_size = y$samples$lib.size * y$samples$norm.factors,
        stringsAsFactors = FALSE
      )
      normalization_rows[[paste(status_key, "norm")]] <- norm

      log_cpm <- cpm(y, log = TRUE, prior.count = 2)
      variable_order <- order(apply(log_cpm, 1, var), decreasing = TRUE)
      pca_features <- variable_order[seq_len(min(500L, length(variable_order)))]
      pca <- prcomp(t(log_cpm[pca_features, , drop = FALSE]), center = TRUE, scale. = FALSE)
      pca_var <- pca$sdev^2 / sum(pca$sdev^2)
      mds <- plotMDS(y, plot = FALSE)
      if (is.null(mds$x) || is.null(mds$y)) {
        classical <- cmdscale(dist(t(log_cpm)), k = 2)
        mds_x <- classical[, 1]
        mds_y <- classical[, 2]
      } else {
        mds_x <- mds$x
        mds_y <- mds$y
      }
      colours <- ifelse(condition == "Healthy", "#2F7EBB", "#D45A42")
      png(file.path(qc_figure_root, paste0(lineage, "_", drop$short_label, "_PCA_MDS.png")), width = 1800, height = 850, res = 150)
      par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
      plot(pca$x[, 1], pca$x[, 2], col = colours, pch = 19, cex = 1.4,
           xlab = sprintf("PC1 (%.1f%%)", 100 * pca_var[1]), ylab = sprintf("PC2 (%.1f%%)", 100 * pca_var[2]),
           main = paste(lineage, drop$short_label, "PCA"))
      text(pca$x[, 1], pca$x[, 2], labels = meta$donor_id, pos = 3, cex = 0.65)
      plot(mds_x, mds_y, col = colours, pch = 19, cex = 1.4, xlab = "MDS1", ylab = "MDS2",
           main = paste(lineage, drop$short_label, "MDS"))
      text(mds_x, mds_y, labels = meta$donor_id, pos = 3, cex = 0.65)
      dev.off()

      qc_rows[[status_key]] <- data.frame(
        lineage = lineage, lodo_run = drop$lodo_run, short_label = drop$short_label,
        dropped_donor = drop$donor_id, dropped_condition = drop$condition,
        n_healthy = n_healthy, n_cie = n_cie, genes_before = nrow(counts), genes_after = sum(keep_gene),
        common_dispersion = y$common.dispersion,
        biological_coefficient_variation = sqrt(y$common.dispersion),
        median_tagwise_dispersion = median(y$tagwise.dispersion, na.rm = TRUE),
        min_norm_factor = min(y$samples$norm.factors), max_norm_factor = max(y$samples$norm.factors),
        min_library_size = min(y$samples$lib.size), max_library_size = max(y$samples$lib.size),
        finite_fit_coefficients = all(is.finite(fit$coefficients)),
        tested_genes = nrow(result), model_status = "SUCCESS",
        stringsAsFactors = FALSE
      )
      cat("SUCCESS; groups", n_healthy, "Healthy vs", n_cie, "CIE; genes", nrow(counts), "->", sum(keep_gene),
          "; BCV", sqrt(y$common.dispersion), "; norm range", paste(range(y$samples$norm.factors), collapse = "-"), "\n")
      list(status = "SUCCESS", error = "", n_healthy = n_healthy, n_cie = n_cie,
           genes_before = nrow(counts), genes_after = sum(keep_gene), output = out_path, ranked = rank_path)
    }, error = function(e) {
      cat("LODO_RUN_FAILED:", conditionMessage(e), "\n")
      list(status = "LODO_RUN_FAILED", error = conditionMessage(e), n_healthy = NA_integer_, n_cie = NA_integer_,
           genes_before = nrow(counts_all), genes_after = NA_integer_, output = "", ranked = "")
    })
    status_rows[[status_key]] <- data.frame(
      lineage = lineage, lodo_run = drop$lodo_run, short_label = drop$short_label,
      dropped_donor = drop$donor_id, dropped_condition = drop$condition,
      status = run_status$status, error = run_status$error,
      n_healthy = run_status$n_healthy, n_cie = run_status$n_cie,
      genes_before = run_status$genes_before, genes_after = run_status$genes_after,
      output_file = run_status$output, ranked_file = run_status$ranked,
      elapsed_seconds = as.numeric(difftime(Sys.time(), run_start, units = "secs")),
      stringsAsFactors = FALSE
    )
  }
}

status <- do.call(rbind, status_rows)
normalization <- do.call(rbind, normalization_rows)
qc <- do.call(rbind, qc_rows)
atomic_write_csv(status, file.path(meta_root, "LODO_RUN_STATUS.csv"))
atomic_write_csv(normalization, file.path(meta_root, "LODO_NORMALIZATION_FACTORS.csv"))
atomic_write_csv(qc, file.path(meta_root, "LODO_MODEL_QC.csv"))

software <- data.frame(
  software = c("R", "edgeR", "limma"),
  version = c(R.version.string, as.character(packageVersion("edgeR")), as.character(packageVersion("limma"))),
  role = c("Phase5 gene LODO runtime", "LODO donor-level QL DE", "edgeR dependency/QC"),
  stringsAsFactors = FALSE
)
atomic_write_csv(software, file.path(project_root, "00_metadata", "PHASE5_SOFTWARE_VERSIONS_R.csv"))
cat("\nCompleted gene LODO:", sum(status$status == "SUCCESS"), "/", nrow(status), "successful runs.\n")
if (any(status$status != "SUCCESS")) quit(status = 2L)
