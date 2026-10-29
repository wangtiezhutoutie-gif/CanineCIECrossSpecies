#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
base <- file.path(project_root, "07_human_crohn", "GSE266616")
pb_root <- file.path(base, "pseudobulk")
matrix_dir <- file.path(pb_root, "matrices")
metadata_dir <- file.path(pb_root, "metadata")
de_dir <- file.path(pb_root, "de_results")
ranked_dir <- file.path(pb_root, "ranked_lists")
sensitivity_dir <- file.path(pb_root, "sensitivity")
sensitivity_ranked_dir <- file.path(sensitivity_dir, "ranked_lists")
paired_dir <- file.path(pb_root, "paired_de_results")
paired_ranked_dir <- file.path(pb_root, "paired_ranked_lists")
figure_dir <- file.path(base, "figures", "pseudobulk")
log_path <- file.path(project_root, "logs", "phase7_edger.log")

for (path in c(de_dir, ranked_dir, sensitivity_dir, sensitivity_ranked_dir, paired_dir, paired_ranked_dir, figure_dir, dirname(log_path))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

atomic_write_csv <- function(frame, path) {
  temp <- paste0(path, ".tmp")
  write.csv(frame, temp, row.names = FALSE, quote = TRUE, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temp, path)) stop("Could not atomically publish ", path)
}

open_png <- function(path, width = 1300, height = 1000) {
  png(path, width = width, height = height, res = 150)
}

read_count_matrix <- function(path) {
  if (!file.exists(path)) stop("Missing count matrix: ", path)
  frame <- read.csv(gzfile(path), check.names = FALSE)
  # Pandas emits an UTF-8 BOM in this controlled CSV; the portable R locale
  # cannot decode that prefix reliably, so freeze the audited first column.
  names(frame)[1] <- "unique_gene_id"
  required <- c("unique_gene_id", "gene_symbol")
  if (!all(required %in% names(frame))) stop("Missing gene identity columns in ", path)
  if (anyDuplicated(frame$unique_gene_id)) stop("Duplicated unique_gene_id in ", path)
  counts <- as.matrix(frame[, setdiff(names(frame), required), drop = FALSE])
  storage.mode(counts) <- "integer"
  rownames(counts) <- frame$unique_gene_id
  if (anyNA(counts) || any(counts < 0L)) stop("Invalid counts in ", path)
  list(frame = frame, counts = counts)
}

fit_edger <- function(counts, meta, design, coef_name, lineage, analysis_name) {
  if (!setequal(colnames(counts), meta$sample_column)) stop("Count columns and metadata sample identity differ: ", analysis_name)
  meta <- meta[match(colnames(counts), meta$sample_column), , drop = FALSE]
  rownames(design) <- meta$sample_column
  y0 <- DGEList(counts = counts)
  keep <- filterByExpr(y0, design = design)
  if (sum(keep) < 1000L) stop("Too few genes after filterByExpr for ", analysis_name)
  y <- y0[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = coef_name)
  raw <- topTags(qlf, n = Inf, sort.by = "none")$table
  if (nrow(raw) != sum(keep)) stop("edgeR result row mismatch for ", analysis_name)
  list(meta = meta, keep = keep, y = y, fit = fit, qlf = qlf, raw = raw)
}

result_frame <- function(fit_result, gene_annotation, lineage, analysis_name, positive_label, negative_label, n_reference, n_test) {
  raw <- fit_result$raw
  annotation <- gene_annotation[match(rownames(raw), gene_annotation$unique_gene_id), , drop = FALSE]
  if (anyNA(annotation$unique_gene_id)) stop("Gene annotation lookup failed for ", analysis_name)
  result <- data.frame(
    unique_gene_id = rownames(raw),
    gene_id = annotation$gene_id,
    gene_symbol = annotation$gene_symbol,
    logFC = raw$logFC,
    logCPM = raw$logCPM,
    F = raw$F,
    PValue = raw$PValue,
    FDR = raw$FDR,
    rank_metric = sign(raw$logFC) * sqrt(pmax(raw$F, 0)),
    direction = ifelse(raw$logFC > 0, positive_label, ifelse(raw$logFC < 0, negative_label, "No_change")),
    lineage = lineage,
    analysis = analysis_name,
    n_reference = n_reference,
    n_test = n_test,
    stringsAsFactors = FALSE
  )
  result[order(result$PValue, -abs(result$logFC), na.last = TRUE), ]
}

ranked_frame <- function(result) {
  ranked <- result[order(result$rank_metric, decreasing = TRUE, na.last = TRUE), c(
    "unique_gene_id", "gene_id", "gene_symbol", "logFC", "F", "PValue", "FDR", "rank_metric"
  )]
  names(ranked)[names(ranked) == "F"] <- "statistic_F"
  ranked
}

safe_slug <- function(value) gsub("[^A-Za-z0-9]+", "_", value)

primary_plots <- function(lineage, fit_result, result, meta) {
  y <- fit_result$y
  fit <- fit_result$fit
  qlf <- fit_result$qlf
  log_cpm <- cpm(y, log = TRUE, prior.count = 2)
  variable_order <- order(apply(log_cpm, 1, var), decreasing = TRUE)
  pca_features <- variable_order[seq_len(min(1000L, length(variable_order)))]
  pca <- prcomp(t(log_cpm[pca_features, , drop = FALSE]), center = TRUE, scale. = FALSE)
  condition <- factor(meta$condition, levels = c("Healthy", "ActiveCrohn"))
  colours <- ifelse(condition == "Healthy", "#2878B5", "#C74440")

  qc <- data.frame(
    lineage = lineage,
    patient_id = meta$patient_id,
    sample_id = meta$sample_id,
    condition = meta$condition,
    n_cells = meta$n_cells,
    raw_library_size = y$samples$lib.size,
    normalization_factor = y$samples$norm.factors,
    effective_library_size = y$samples$lib.size * y$samples$norm.factors,
    detected_genes = meta$detected_genes,
    PCA1 = pca$x[, 1],
    PCA2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )

  open_png(file.path(figure_dir, paste0("Human_", lineage, "_library_size.png")), 1800, 900)
  par(mar = c(10, 5, 3, 1))
  barplot(qc$raw_library_size, names.arg = qc$patient_id, las = 2, col = colours,
          ylab = "Raw pseudobulk library size", main = paste("Human", lineage, "library size"))
  dev.off()

  open_png(file.path(figure_dir, paste0("Human_", lineage, "_PCA.png")))
  plot(qc$PCA1, qc$PCA2, col = colours, pch = 19, cex = 1.2,
       xlab = "PCA1", ylab = "PCA2", main = paste("Human", lineage, "pseudobulk PCA"))
  text(qc$PCA1, qc$PCA2, labels = qc$patient_id, pos = 3, cex = 0.55)
  legend("topright", legend = c("Healthy", "Active Crohn"), col = c("#2878B5", "#C74440"), pch = 19, bty = "n")
  dev.off()

  open_png(file.path(figure_dir, paste0("Human_", lineage, "_MDS.png")))
  plotMDS(y, col = colours, pch = 19, labels = meta$patient_id, main = paste("Human", lineage, "edgeR MDS"))
  legend("topright", legend = c("Healthy", "Active Crohn"), col = c("#2878B5", "#C74440"), pch = 19, bty = "n")
  dev.off()

  correlation <- cor(log_cpm, method = "pearson")
  open_png(file.path(figure_dir, paste0("Human_", lineage, "_sample_correlation.png")), 1800, 1600)
  heatmap(correlation, symm = TRUE, margins = c(10, 10), col = colorRampPalette(c("#26456E", "white", "#A3333D"))(101),
          main = paste("Human", lineage, "logCPM sample correlation"))
  dev.off()

  open_png(file.path(figure_dir, paste0("Human_", lineage, "_BCV.png")))
  plotBCV(y, main = paste("Human", lineage, "biological coefficient of variation"))
  dev.off()

  open_png(file.path(figure_dir, paste0("Human_", lineage, "_QL_dispersion.png")))
  plotQLDisp(fit, main = paste("Human", lineage, "QL dispersion"))
  dev.off()

  open_png(file.path(figure_dir, paste0("Human_", lineage, "_MA.png")))
  status <- ifelse(result$FDR < 0.05 & result$logFC > 0, 1L,
                   ifelse(result$FDR < 0.05 & result$logFC < 0, -1L, 0L))
  names(status) <- result$unique_gene_id
  plotMD(qlf, status = status, main = paste("Human", lineage, "Active Crohn vs Healthy"))
  abline(h = 0, lty = 2, col = "grey40")
  dev.off()

  open_png(file.path(figure_dir, paste0("Human_", lineage, "_volcano.png")))
  colours_de <- ifelse(result$FDR < 0.05 & result$logFC > 0, "#C74440",
                       ifelse(result$FDR < 0.05 & result$logFC < 0, "#2878B5", "#B8B8B8"))
  plot(result$logFC, -log10(pmax(result$PValue, .Machine$double.xmin)), pch = 19, cex = 0.4,
       col = colours_de, xlab = "log2 fold-change (Active Crohn - Healthy)", ylab = "-log10(PValue)",
       main = paste("Human", lineage, "differential expression"))
  abline(v = 0, lty = 2, col = "grey40")
  dev.off()
  qc
}

subset_for_sensitivity <- function(meta, sensitivity) {
  healthy <- meta$condition == "Healthy"
  active <- meta$condition == "ActiveCrohn"
  if (grepl("_active_Crohn_TI_vs_Healthy_TI$", sensitivity)) {
    class_name <- sub("_active_Crohn_TI_vs_Healthy_TI$", "", sensitivity)
    return(meta[healthy | (active & meta$biologic_class == class_name), , drop = FALSE])
  }
  if (grepl("^exclude_", sensitivity)) {
    excluded <- sub("^exclude_", "", sensitivity)
    return(meta[healthy | (active & meta$biologic_class != excluded), , drop = FALSE])
  }
  stop("Unknown treatment sensitivity: ", sensitivity)
}

run_primary_and_treatment <- function(gene_annotation) {
  primary_meta <- read.csv(file.path(metadata_dir, "HUMAN_PSEUDOBULK_SAMPLE_METADATA.csv"), check.names = FALSE)
  names(primary_meta)[1] <- "analysis"
  plan <- read.csv(file.path(project_root, "00_metadata", "HUMAN_TREATMENT_SENSITIVITY_PLAN.csv"), check.names = FALSE)
  names(plan)[1] <- "sensitivity"
  lineages <- c("Epithelial", "T", "Myeloid")
  filter_rows <- list()
  qc_rows <- list()
  treatment_summary <- list()
  primary_results <- list()

  for (lineage in lineages) {
    matrix_info <- read_count_matrix(file.path(matrix_dir, paste0("human_", lineage, "_counts.csv.gz")))
    counts <- matrix_info$counts
    meta <- primary_meta[primary_meta$lineage == lineage & primary_meta$eligible == "YES", , drop = FALSE]
    meta$sample_column <- meta$patient_id
    if (nrow(meta) != 27L || length(unique(meta$patient_id)) != 27L) stop("Expected 27 eligible Primary donors for ", lineage)
    meta <- meta[match(colnames(counts), meta$sample_column), , drop = FALSE]
    condition <- factor(meta$condition, levels = c("Healthy", "ActiveCrohn"))
    if (sum(condition == "Healthy") != 16L || sum(condition == "ActiveCrohn") != 11L) stop("Primary condition counts changed for ", lineage)
    design <- model.matrix(~ condition)
    fit_result <- fit_edger(counts, meta, design, "conditionActiveCrohn", lineage, "Primary")
    result <- result_frame(fit_result, gene_annotation, lineage, "Primary_ActiveCrohn_vs_Healthy",
                           "ActiveCrohn_higher", "ActiveCrohn_lower", 16L, 11L)
    primary_results[[lineage]] <- result
    atomic_write_csv(result, file.path(de_dir, paste0("Human_", lineage, "_ActiveCrohn_vs_Healthy.csv")))
    atomic_write_csv(result[result$FDR < 0.05, ], file.path(de_dir, paste0("Human_", lineage, "_ActiveCrohn_vs_Healthy_FDR05.csv")))
    atomic_write_csv(ranked_frame(result), file.path(ranked_dir, paste0("Human_", lineage, "_ActiveCrohn_vs_Healthy_ranked.csv")))
    qc_rows[[lineage]] <- primary_plots(lineage, fit_result, result, fit_result$meta)
    filter_rows[[lineage]] <- data.frame(
      lineage = lineage,
      genes_before = nrow(counts),
      genes_after = sum(fit_result$keep),
      genes_removed = sum(!fit_result$keep),
      filter_method = "edgeR::filterByExpr",
      design = "~condition",
      stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(plan))) {
      if (plan$executable[i] != "YES") next
      sensitivity <- plan$sensitivity[i]
      subset <- subset_for_sensitivity(meta, sensitivity)
      active_n <- sum(subset$condition == "ActiveCrohn")
      if (active_n != plan$active_crohn_donor_n[i]) stop("Treatment plan donor count mismatch: ", sensitivity)
      sub_counts <- counts[, subset$sample_column, drop = FALSE]
      subset <- subset[match(colnames(sub_counts), subset$sample_column), , drop = FALSE]
      sub_condition <- factor(subset$condition, levels = c("Healthy", "ActiveCrohn"))
      sub_design <- model.matrix(~ sub_condition)
      fit_sens <- fit_edger(sub_counts, subset, sub_design, "sub_conditionActiveCrohn", lineage, sensitivity)
      sens_result <- result_frame(fit_sens, gene_annotation, lineage, sensitivity,
                                  "ActiveCrohn_higher", "ActiveCrohn_lower", 16L, active_n)
      slug <- safe_slug(sensitivity)
      atomic_write_csv(sens_result, file.path(sensitivity_dir, paste0("Human_", lineage, "_", slug, ".csv")))
      atomic_write_csv(ranked_frame(sens_result), file.path(sensitivity_ranked_dir, paste0("Human_", lineage, "_", slug, "_ranked.csv")))
      common <- merge(
        result[, c("unique_gene_id", "logFC", "FDR")],
        sens_result[, c("unique_gene_id", "logFC", "FDR")],
        by = "unique_gene_id",
        suffixes = c("_primary", "_sensitivity")
      )
      treatment_summary[[paste(lineage, sensitivity, sep = "::")]] <- data.frame(
        lineage = lineage,
        sensitivity = sensitivity,
        active_crohn_donor_n = active_n,
        tested_genes = nrow(sens_result),
        spearman_logFC = cor(common$logFC_primary, common$logFC_sensitivity, method = "spearman"),
        sign_concordance = mean(sign(common$logFC_primary) == sign(common$logFC_sensitivity)),
        primary_FDR05 = sum(result$FDR < 0.05),
        sensitivity_FDR05 = sum(sens_result$FDR < 0.05),
        interpretation_limit = "Treatment-confounding sensitivity; not a causal treatment-effect analysis",
        stringsAsFactors = FALSE
      )
    }
  }
  atomic_write_csv(do.call(rbind, filter_rows), file.path(metadata_dir, "HUMAN_EDGER_FILTER_SUMMARY.csv"))
  atomic_write_csv(do.call(rbind, qc_rows), file.path(metadata_dir, "HUMAN_PSEUDOBULK_SAMPLE_QC.csv"))
  atomic_write_csv(do.call(rbind, treatment_summary), file.path(sensitivity_dir, "HUMAN_TREATMENT_GENE_SENSITIVITY_SUMMARY.csv"))
  primary_results
}

run_paired <- function(gene_annotation) {
  plan <- read.csv(file.path(metadata_dir, "HUMAN_PAIRED_SECONDARY_PLAN.csv"), check.names = FALSE)
  names(plan)[1] <- "comparison"
  meta_all <- read.csv(file.path(metadata_dir, "HUMAN_PAIRED_PSEUDOBULK_SAMPLE_METADATA.csv"), check.names = FALSE)
  names(meta_all)[1] <- "comparison"
  summary_rows <- list()
  lineages <- c("Epithelial", "T", "Myeloid")
  for (i in seq_len(nrow(plan))) {
    if (plan$executable[i] != "YES") next
    comparison <- plan$comparison[i]
    target_state <- plan$target_state[i]
    for (lineage in lineages) {
      matrix_info <- read_count_matrix(file.path(matrix_dir, paste0("human_", lineage, "_", comparison, "_counts.csv.gz")))
      counts <- matrix_info$counts
      meta <- meta_all[meta_all$comparison == comparison & meta_all$lineage == lineage, , drop = FALSE]
      meta <- meta[meta$eligible_ge20_cells == "YES", , drop = FALSE]
      complete_donors <- names(which(table(meta$patient_id) == 2L))
      meta <- meta[meta$patient_id %in% complete_donors, , drop = FALSE]
      if (length(complete_donors) < 3L) {
        summary_rows[[paste(comparison, lineage, sep = "::")]] <- data.frame(
          comparison = comparison, lineage = lineage, paired_donors = length(complete_donors), status = "NO_DESCRIPTIVE_ONLY_AFTER_CELL_COVERAGE",
          tested_genes = 0L, FDR05 = 0L, stringsAsFactors = FALSE
        )
        next
      }
      counts <- counts[, meta$sample_id, drop = FALSE]
      meta$sample_column <- meta$sample_id
      meta <- meta[match(colnames(counts), meta$sample_column), , drop = FALSE]
      meta$patient_id <- factor(meta$patient_id)
      meta$state <- factor(meta$state, levels = c(target_state, "Active"))
      design <- model.matrix(~ patient_id + state, data = meta)
      fit_result <- fit_edger(counts, meta, design, "stateActive", lineage, comparison)
      result <- result_frame(fit_result, gene_annotation, lineage, comparison,
                             "Active_higher", "Active_lower", length(complete_donors), length(complete_donors))
      atomic_write_csv(result, file.path(paired_dir, paste0("Human_", lineage, "_", comparison, ".csv")))
      atomic_write_csv(ranked_frame(result), file.path(paired_ranked_dir, paste0("Human_", lineage, "_", comparison, "_ranked.csv")))
      summary_rows[[paste(comparison, lineage, sep = "::")]] <- data.frame(
        comparison = comparison, lineage = lineage, paired_donors = length(complete_donors), status = "EXECUTED_DONOR_BLOCKED",
        tested_genes = nrow(result), FDR05 = sum(result$FDR < 0.05), stringsAsFactors = FALSE
      )
    }
  }
  if (length(summary_rows)) {
    atomic_write_csv(do.call(rbind, summary_rows), file.path(paired_dir, "HUMAN_PAIRED_DE_SUMMARY.csv"))
  } else {
    atomic_write_csv(data.frame(comparison = character(), lineage = character(), paired_donors = integer(), status = character(), tested_genes = integer(), FDR05 = integer()),
                     file.path(paired_dir, "HUMAN_PAIRED_DE_SUMMARY.csv"))
  }
}

write_platform_decision <- function() {
  audit <- read.csv(file.path(project_root, "00_metadata", "GSE266616_PRIMARY_PLATFORM_AUDIT.csv"), check.names = FALSE)
  names(audit)[1] <- "platform"
  nextseq <- audit[audit$platform == "NextSeq 500", , drop = FALSE]
  decision <- data.frame(
    sensitivity = "~platform+condition",
    executable = "NO_PRIMARY_NEXTSEQ_N1_PER_CONDITION",
    design_estimable_algebraically = "YES",
    sparsity = paste(paste(nextseq$primary_condition, nextseq$donor_n, sep = "="), collapse = ";"),
    interpretation = "Primary remains ~condition; no claim of platform adjustment because the secondary platform has one donor per condition",
    stringsAsFactors = FALSE
  )
  atomic_write_csv(decision, file.path(metadata_dir, "HUMAN_PLATFORM_SENSITIVITY_DECISION.csv"))
}

main <- function() {
  sink(log_path, split = TRUE)
  on.exit(sink(), add = TRUE)
  cat("Phase7 Human edgeR analysis\n")
  cat("Start:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
  set.seed(20240825)
  suppressPackageStartupMessages(library(edgeR))
  suppressPackageStartupMessages(library(limma))
  cat("R:", R.version.string, "\n")
  cat("edgeR:", as.character(packageVersion("edgeR")), " limma:", as.character(packageVersion("limma")), "\n")
  gene_annotation <- read.csv(file.path(metadata_dir, "HUMAN_PSEUDOBULK_GENE_ANNOTATION.csv"), check.names = FALSE)
  names(gene_annotation)[1] <- "unique_gene_id"
  if (anyDuplicated(gene_annotation$unique_gene_id)) stop("Duplicated Human gene annotation IDs")
  run_primary_and_treatment(gene_annotation)
  run_paired(gene_annotation)
  write_platform_decision()
  software <- data.frame(
    software = c("R", "edgeR", "limma", "locfit", "statmod"),
    version = c(as.character(getRversion()), as.character(packageVersion("edgeR")), as.character(packageVersion("limma")),
                as.character(packageVersion("locfit")), as.character(packageVersion("statmod"))),
    role = c("statistical runtime", "donor-level quasi-likelihood DE", "edgeR dependency/plotting", "robust dispersion dependency", "robust QL dependency"),
    stringsAsFactors = FALSE
  )
  atomic_write_csv(software, file.path(project_root, "00_metadata", "PHASE7_EDGER_SOFTWARE_VERSIONS.csv"))
  cat("End:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
}

status <- 0L
tryCatch(main(), error = function(error) {
  message("FATAL: ", conditionMessage(error))
  message("CALL: ", paste(deparse(conditionCall(error)), collapse = " "))
  calls <- sys.calls()
  for (call in tail(calls, 12L)) message("STACK: ", paste(deparse(call), collapse = " "))
  status <<- 1L
})
quit(save = "no", status = status)
