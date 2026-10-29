#!/usr/bin/env Rscript

# Phase 3 primary donor-level pseudobulk DE with edgeR quasi-likelihood.

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)

pb_root <- file.path(project_root, "06_pseudobulk")
matrix_dir <- file.path(pb_root, "matrices")
metadata_dir <- file.path(pb_root, "metadata")
de_dir <- file.path(pb_root, "de_results")
figure_dir <- file.path(pb_root, "figures")
ranked_dir <- file.path(pb_root, "ranked_lists")
log_path <- file.path(project_root, "logs", "phase3_edger.log")
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ranked_dir, recursive = TRUE, showWarnings = FALSE)

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

open_png <- function(path, width = 1500, height = 1100) {
  png(path, width = width, height = height, res = 150)
}

robust_flag <- function(values, side = c("two", "low", "high"), threshold = 3) {
  side <- match.arg(side)
  center <- median(values, na.rm = TRUE)
  scale <- mad(values, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale) || scale == 0) return(rep(FALSE, length(values)))
  score <- (values - center) / scale
  if (side == "two") abs(score) > threshold else if (side == "low") score < -threshold else score > threshold
}

lineage_slug <- function(lineage) gsub("[^A-Za-z0-9]+", "_", lineage)

main <- function() {
  set.seed(20240824)
  start_time <- Sys.time()
  cat("Phase3 edgeR primary analysis\n")
  cat("Start:", format(start_time, tz = "Asia/Shanghai"), "\n")
  cat("Project root:", project_root, "\n")
  cat("Random seed: 20240824 (no stochastic inferential step)\n")

  suppressPackageStartupMessages(library(edgeR))
  suppressPackageStartupMessages(library(limma))
  cat("R version:", R.version.string, "\n")
  cat("edgeR version:", as.character(packageVersion("edgeR")), "\n")
  cat("limma version:", as.character(packageVersion("limma")), "\n")
  cat("Model: ~ condition; Healthy reference; coefficient conditionCIE\n")
  cat("Workflow: DGEList -> filterByExpr -> TMM -> estimateDisp -> glmQLFit(robust=TRUE) -> glmQLFTest\n")

  sample_metadata_path <- file.path(metadata_dir, "pseudobulk_sample_metadata.csv")
  annotation_path <- file.path(metadata_dir, "pseudobulk_gene_annotation.csv")
  if (!file.exists(sample_metadata_path) || !file.exists(annotation_path)) stop("Missing Phase3 pseudobulk metadata")
  sample_metadata <- read.csv(sample_metadata_path, check.names = FALSE)
  annotation <- read.csv(annotation_path, check.names = FALSE)
  if (anyDuplicated(annotation$feature_id)) stop("Feature IDs are duplicated in gene annotation")
  rownames(annotation) <- annotation$feature_id

  lineages <- c("Epithelial", "T", "Myeloid")
  filter_rows <- list()
  normalization_rows <- list()
  qc_rows <- list()

  for (lineage in lineages) {
    cat("\n===", lineage, "===\n")
    slug <- lineage_slug(lineage)
    matrix_path <- file.path(matrix_dir, paste0("pseudobulk_", slug, "_counts.csv.gz"))
    if (!file.exists(matrix_path)) stop("Missing matrix: ", matrix_path)
    count_frame <- read.csv(gzfile(matrix_path), check.names = FALSE)
    if (ncol(count_frame) != 8 || nrow(count_frame) != 20257) stop("Unexpected matrix dimensions for ", lineage)
    if (anyDuplicated(count_frame$feature_id)) stop("Duplicated matrix feature_id for ", lineage)
    feature_ids <- count_frame$feature_id
    counts <- as.matrix(count_frame[, -1, drop = FALSE])
    storage.mode(counts) <- "integer"
    rownames(counts) <- feature_ids
    if (any(counts < 0L) || any(is.na(counts))) stop("Invalid counts for ", lineage)

    eligible_flag <- tolower(as.character(sample_metadata$eligible)) == "true"
    meta <- sample_metadata[sample_metadata$lineage == lineage & eligible_flag, , drop = FALSE]
    if (nrow(meta) != 7 || length(unique(meta$donor_id)) != 7) stop("Primary lineage does not have 7 eligible donors: ", lineage)
    if (!setequal(colnames(counts), meta$donor_id)) stop("Count columns and metadata donors differ: ", lineage)
    meta <- meta[match(colnames(counts), meta$donor_id), , drop = FALSE]
    condition <- factor(meta$condition, levels = c("Healthy", "CIE"))
    if (sum(condition == "Healthy") != 3 || sum(condition == "CIE") != 4) stop("Unexpected condition counts: ", lineage)
    design <- model.matrix(~ condition)
    rownames(design) <- meta$donor_id
    cat("Design matrix:\n")
    print(design)

    y_unfiltered <- DGEList(counts = counts, samples = data.frame(
      donor_id = meta$donor_id,
      condition = condition,
      row.names = meta$donor_id
    ))
    keep <- filterByExpr(y_unfiltered, design = design)
    y <- y_unfiltered[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y, method = "TMM")
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE)
    qlf <- glmQLFTest(fit, coef = "conditionCIE")
    raw_table <- topTags(qlf, n = Inf, sort.by = "none")$table
    if (nrow(raw_table) != sum(keep)) stop("Test table row count mismatch: ", lineage)

    result <- data.frame(
      feature_id = rownames(raw_table),
      gene_symbol = annotation[rownames(raw_table), "original_gene_symbol"],
      unique_var_name = annotation[rownames(raw_table), "unique_var_name"],
      logFC = raw_table$logFC,
      logCPM = raw_table$logCPM,
      F = raw_table$F,
      PValue = raw_table$PValue,
      FDR = raw_table$FDR,
      direction = ifelse(raw_table$logFC > 0, "CIE_higher", ifelse(raw_table$logFC < 0, "CIE_lower", "No_change")),
      lineage = lineage,
      n_healthy = 3L,
      n_cie = 4L,
      rank_metric = sign(raw_table$logFC) * sqrt(pmax(raw_table$F, 0)),
      stringsAsFactors = FALSE
    )
    result <- result[order(result$PValue, -abs(result$logFC), na.last = TRUE), ]
    out_path <- file.path(de_dir, paste0(slug, "_CIE_vs_Healthy_edger.csv"))
    atomic_write_csv(result, out_path)
    atomic_write_csv(result[result$FDR < 0.05, ], file.path(de_dir, paste0(slug, "_CIE_vs_Healthy_edger_FDR05.csv")))
    atomic_write_csv(result[result$FDR < 0.10, ], file.path(de_dir, paste0(slug, "_CIE_vs_Healthy_edger_FDR10.csv")))

    ranked <- result[order(result$rank_metric, decreasing = TRUE, na.last = TRUE), c(
      "feature_id", "gene_symbol", "unique_var_name", "logFC", "F", "PValue", "FDR", "rank_metric"
    )]
    names(ranked)[names(ranked) == "F"] <- "statistic_F"
    atomic_write_csv(ranked, file.path(ranked_dir, paste0(slug, "_ranked.csv")))

    filter_rows[[lineage]] <- data.frame(
      lineage = lineage,
      genes_before = nrow(counts),
      genes_after = sum(keep),
      genes_removed = sum(!keep),
      filter_method = "edgeR::filterByExpr",
      parameters = "default arguments; lineage-specific design=~condition",
      stringsAsFactors = FALSE
    )
    norm_frame <- data.frame(
      lineage = lineage,
      donor_id = colnames(y),
      condition = as.character(condition),
      raw_library_size = y$samples$lib.size,
      normalization_factor = y$samples$norm.factors,
      effective_library_size = y$samples$lib.size * y$samples$norm.factors,
      stringsAsFactors = FALSE
    )
    normalization_rows[[lineage]] <- norm_frame

    log_cpm <- cpm(y, log = TRUE, prior.count = 2)
    variable_order <- order(apply(log_cpm, 1, var), decreasing = TRUE)
    pca_features <- variable_order[seq_len(min(500, length(variable_order)))]
    pca <- prcomp(t(log_cpm[pca_features, , drop = FALSE]), center = TRUE, scale. = FALSE)
    pca_variance <- pca$sdev^2 / sum(pca$sdev^2)
    correlation <- cor(log_cpm, method = "pearson")
    correlation_without_self <- correlation
    diag(correlation_without_self) <- NA_real_
    min_correlation <- apply(correlation_without_self, 2, min, na.rm = TRUE)
    pca_distance <- sqrt(pca$x[, 1]^2 + pca$x[, 2]^2)
    log_library <- log2(norm_frame$raw_library_size)
    flag_low_library <- robust_flag(log_library, side = "low")
    flag_low_library_ratio <- norm_frame$raw_library_size < 0.25 * median(norm_frame$raw_library_size)
    flag_low_correlation <- robust_flag(min_correlation, side = "low")
    flag_pca_distance <- robust_flag(pca_distance, side = "high")

    top_de_ids <- head(result$feature_id[order(abs(result$rank_metric), decreasing = TRUE)], min(100L, nrow(result)))
    top_de_expression <- log_cpm[top_de_ids, , drop = FALSE]
    top_de_contribution <- vapply(seq_len(ncol(top_de_expression)), function(j) {
      same_group <- condition == condition[j]
      group_mean <- rowMeans(top_de_expression[, same_group, drop = FALSE])
      median(abs(top_de_expression[, j] - group_mean))
    }, numeric(1))
    flag_top_de_contribution <- robust_flag(top_de_contribution, side = "high")

    mds <- plotMDS(y, plot = FALSE)
    mds_x <- mds$x
    mds_y <- mds$y
    if (is.null(mds_x) || is.null(mds_y)) {
      classical_mds <- cmdscale(dist(t(log_cpm)), k = 2)
      mds_x <- classical_mds[, 1]
      mds_y <- classical_mds[, 2]
    }
    mds_distance <- sqrt(mds_x^2 + mds_y^2)
    flag_mds_distance <- robust_flag(mds_distance, side = "high")
    qc_rows[[lineage]] <- data.frame(
      lineage = lineage,
      donor_id = meta$donor_id,
      condition = as.character(condition),
      n_cells = meta$n_cells,
      raw_library_size = norm_frame$raw_library_size,
      n_detected_features = meta$n_detected_features,
      PCA1 = pca$x[, 1],
      PCA2 = pca$x[, 2],
      MDS1 = mds_x,
      MDS2 = mds_y,
      min_sample_correlation = min_correlation,
      pca_distance = pca_distance,
      mds_distance = mds_distance,
      top100_de_within_group_deviation = top_de_contribution,
      low_library_flag = flag_low_library,
      low_library_ratio_flag = flag_low_library_ratio,
      low_correlation_flag = flag_low_correlation,
      pca_distance_flag = flag_pca_distance,
      mds_distance_flag = flag_mds_distance,
      top_de_contribution_flag = flag_top_de_contribution,
      descriptive_outlier_flag = flag_low_library | flag_low_library_ratio | flag_low_correlation |
        flag_pca_distance | flag_mds_distance | flag_top_de_contribution,
      note = ifelse(meta$donor_id == "GSM8031701_cie_2", "Pre-specified cie_2 focus; retained in primary model", "Retained in primary model"),
      stringsAsFactors = FALSE
    )

    condition_colours <- ifelse(condition == "Healthy", "#2F7EBB", "#D45A42")
    point_shapes <- ifelse(meta$donor_id == "GSM8031701_cie_2", 17, 19)

    open_png(file.path(figure_dir, paste0(slug, "_library_and_detected_genes.png")), width = 1700, height = 900)
    par(mfrow = c(1, 2), mar = c(10, 4, 3, 1))
    barplot(meta$library_size, names.arg = meta$donor_id, las = 2, col = condition_colours,
            ylab = "Raw pseudobulk library size", main = paste(lineage, "library size"))
    barplot(meta$n_detected_features, names.arg = meta$donor_id, las = 2, col = condition_colours,
            ylab = "Detected features", main = paste(lineage, "detected features"))
    dev.off()

    open_png(file.path(figure_dir, paste0(slug, "_PCA.png")))
    plot(pca$x[, 1], pca$x[, 2], col = condition_colours, pch = point_shapes, cex = 1.5,
         xlab = sprintf("PC1 (%.1f%%)", 100 * pca_variance[1]),
         ylab = sprintf("PC2 (%.1f%%)", 100 * pca_variance[2]), main = paste(lineage, "pseudobulk PCA"))
    text(pca$x[, 1], pca$x[, 2], labels = meta$donor_id, pos = 3, cex = 0.68)
    legend("topright", legend = c("Healthy", "CIE", "cie_2"), col = c("#2F7EBB", "#D45A42", "#D45A42"), pch = c(19, 19, 17), bty = "n")
    dev.off()

    open_png(file.path(figure_dir, paste0(slug, "_MDS.png")))
    plotMDS(y, col = condition_colours, pch = point_shapes, labels = meta$donor_id,
            main = paste(lineage, "edgeR MDS"))
    legend("topright", legend = c("Healthy", "CIE", "cie_2"), col = c("#2F7EBB", "#D45A42", "#D45A42"), pch = c(19, 19, 17), bty = "n")
    dev.off()

    open_png(file.path(figure_dir, paste0(slug, "_sample_correlation.png")), width = 1500, height = 1300)
    heatmap(correlation, symm = TRUE, margins = c(12, 12), col = colorRampPalette(c("#2F4B7C", "white", "#B3261E"))(101),
            main = paste(lineage, "logCPM sample correlation"))
    dev.off()

    open_png(file.path(figure_dir, paste0(slug, "_BCV.png")))
    plotBCV(y, main = paste(lineage, "biological coefficient of variation"))
    dev.off()

    open_png(file.path(figure_dir, paste0(slug, "_QL_dispersion.png")))
    plotQLDisp(fit, main = paste(lineage, "QL dispersion"))
    dev.off()

    open_png(file.path(figure_dir, paste0(slug, "_MA.png")))
    md_status <- ifelse(raw_table$FDR < 0.05 & raw_table$logFC > 0, 1L,
                        ifelse(raw_table$FDR < 0.05 & raw_table$logFC < 0, -1L, 0L))
    names(md_status) <- rownames(raw_table)
    plotMD(qlf, status = md_status,
           main = paste(lineage, "CIE vs Healthy MA"))
    abline(h = 0, lty = 2, col = "grey40")
    dev.off()

    volcano_colours <- ifelse(result$FDR < 0.05 & result$logFC > 0, "#D45A42",
                              ifelse(result$FDR < 0.05 & result$logFC < 0, "#2F7EBB", "#B8B8B8"))
    open_png(file.path(figure_dir, paste0(slug, "_volcano.png")))
    plot(result$logFC, -log10(pmax(result$PValue, .Machine$double.xmin)), pch = 19, cex = 0.45,
         col = volcano_colours, xlab = "log2 fold-change (CIE - Healthy)", ylab = "-log10(PValue)",
         main = paste(lineage, "CIE-associated DE"))
    abline(v = 0, lty = 2, col = "grey40")
    dev.off()

    cat("Cells:", sum(meta$n_cells), "\n")
    cat("Genes before/after filter:", nrow(counts), "/", sum(keep), "\n")
    cat("TMM normalization factors:\n")
    print(norm_frame)
    cat("FDR <0.05:", sum(result$FDR < 0.05), "; FDR <0.10:", sum(result$FDR < 0.10), "\n")
    cat("Descriptive sample outlier flags:", paste(qc_rows[[lineage]]$donor_id[qc_rows[[lineage]]$descriptive_outlier_flag], collapse = "; "), "\n")
  }

  atomic_write_csv(do.call(rbind, filter_rows), file.path(metadata_dir, "gene_filter_summary.csv"))
  atomic_write_csv(do.call(rbind, normalization_rows), file.path(metadata_dir, "normalization_factors.csv"))
  atomic_write_csv(do.call(rbind, qc_rows), file.path(metadata_dir, "pseudobulk_sample_qc.csv"))

  environment_lines <- c(
    paste("Phase3 environment captured", format(Sys.time(), tz = "Asia/Shanghai")),
    paste("R", R.version.string),
    paste("Platform", R.version$platform),
    paste("edgeR", as.character(packageVersion("edgeR"))),
    paste("limma", as.character(packageVersion("limma"))),
    paste("locfit", as.character(packageVersion("locfit"))),
    paste("statmod", as.character(packageVersion("statmod"))),
    paste("R executable", file.path(R.home("bin"), "R.exe")),
    paste("R library paths", paste(.libPaths(), collapse = "; ")),
    "Primary model: ~ condition (Healthy reference)",
    "Primary framework: edgeR quasi-likelihood, robust=TRUE",
    "Random seed: 20240824"
  )
  writeLines(environment_lines, file.path(project_root, "00_metadata", "phase3_environment.txt"), useBytes = TRUE)
  software <- data.frame(
    software = c("R", "edgeR", "limma", "locfit", "statmod"),
    version = c(as.character(getRversion()), as.character(packageVersion("edgeR")), as.character(packageVersion("limma")),
                as.character(packageVersion("locfit")), as.character(packageVersion("statmod"))),
    role = c("statistical runtime", "pseudobulk quasi-likelihood DE", "edgeR dependency/plotting", "robust dispersion dependency", "robust QL dependency"),
    stringsAsFactors = FALSE
  )
  atomic_write_csv(software, file.path(project_root, "00_metadata", "PHASE3_SOFTWARE_VERSIONS.csv"))
  cat("\nEnd:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
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
