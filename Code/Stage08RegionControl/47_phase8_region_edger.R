#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
phase8_root <- file.path(project_root, "08_region_control")
base <- file.path(phase8_root, "human_AC_TI")
matrix_dir <- file.path(base, "pseudobulk", "matrices")
metadata_path <- file.path(phase8_root, "metadata", "HUMAN_REGION_PSEUDOBULK_METADATA.csv")
gene_path <- file.path(project_root, "00_metadata", "GSE266616_HUMAN_GENE_MASTER.csv")
de_dir <- file.path(base, "de_results")
ranked_dir <- file.path(base, "ranked_lists")
figure_dir <- file.path(phase8_root, "figures", "edgeR")
summary_path <- file.path(phase8_root, "summary", "HUMAN_REGION_DE_SUMMARY.csv")
sample_qc_path <- file.path(phase8_root, "metadata", "HUMAN_REGION_EDGER_SAMPLE_QC.csv")
gate_path <- file.path(phase8_root, "metadata", "PHASE8A_REGION_EDGER_GATE.csv")
log_paths <- c(
  file.path(project_root, "logs", "phase8_region_analysis.log"),
  file.path(phase8_root, "logs", "phase8_region_analysis.log")
)

for (path in c(de_dir, ranked_dir, figure_dir, dirname(summary_path), dirname(gate_path), dirname(log_paths))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

suppressPackageStartupMessages({
  library(edgeR)
})

log_line <- function(message) {
  line <- sprintf("[%s] %s", format(Sys.time(), tz = "UTC", usetz = TRUE), message)
  message(line)
  for (path in log_paths) cat(line, "\n", file = path, append = TRUE)
}

atomic_write_csv <- function(frame, path) {
  temp <- paste0(path, ".tmp")
  write.csv(frame, temp, row.names = FALSE, quote = TRUE, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temp, path)) stop("Could not atomically publish ", path)
}

read_count_matrix <- function(path) {
  if (!file.exists(path)) stop("Missing matrix: ", path)
  frame <- read.csv(gzfile(path), check.names = FALSE)
  names(frame)[1:2] <- c("unique_gene_id", "gene_symbol")
  if (anyDuplicated(frame$unique_gene_id)) stop("Duplicated unique gene IDs: ", path)
  counts <- as.matrix(frame[, -(1:2), drop = FALSE])
  storage.mode(counts) <- "integer"
  rownames(counts) <- frame$unique_gene_id
  if (anyNA(counts) || any(counts < 0L)) stop("Invalid counts: ", path)
  list(frame = frame, counts = counts)
}

fit_one <- function(analysis, lineage, metadata, gene_annotation) {
  matrix_path <- file.path(matrix_dir, paste0(analysis, "_", lineage, "_counts.csv.gz"))
  matrix_info <- read_count_matrix(matrix_path)
  counts <- matrix_info$counts
  meta <- metadata[metadata$analysis == analysis & metadata$lineage == lineage & metadata$eligible == "YES", , drop = FALSE]
  if (!setequal(colnames(counts), meta$sample_column)) stop("Matrix/metadata mismatch: ", analysis, " ", lineage)
  meta <- meta[match(colnames(counts), meta$sample_column), , drop = FALSE]

  if (analysis == "Healthy_AC_vs_TI") {
    tab <- table(meta$patient_id, meta$region)
    paired_ids <- rownames(tab)[rowSums(tab == 1L) == 2L]
    meta <- meta[meta$patient_id %in% paired_ids, , drop = FALSE]
    counts <- counts[, meta$sample_column, drop = FALSE]
    meta <- meta[match(colnames(counts), meta$sample_column), , drop = FALSE]
    patient <- factor(meta$patient_id)
    region <- factor(meta$region, levels = c("TI", "AC"))
    design <- model.matrix(~ patient + region)
    coefficient <- "regionAC"
    positive_label <- "AC_higher"
    negative_label <- "AC_lower"
    group_evidence <- paste0("paired_donors=", length(unique(meta$patient_id)))
    if (length(unique(meta$patient_id)) != 15L) stop("Expected 15 Healthy pairs for ", lineage)
  } else if (analysis == "ActiveCrohn_AC_vs_TI") {
    region <- factor(meta$region, levels = c("TI", "AC"))
    design <- model.matrix(~ region)
    coefficient <- "regionAC"
    positive_label <- "AC_higher"
    negative_label <- "AC_lower"
    region_tab <- table(region)
    group_evidence <- paste0("TI=", region_tab[["TI"]], ";AC=", region_tab[["AC"]])
    if (min(region_tab) < 3L) stop("Active regional group below n=3: ", lineage)
  } else if (analysis == "AC_ActiveCrohn_vs_Healthy") {
    condition <- factor(meta$condition, levels = c("Healthy", "ActiveCrohn"))
    design <- model.matrix(~ condition)
    coefficient <- "conditionActiveCrohn"
    positive_label <- "ActiveCrohn_higher"
    negative_label <- "ActiveCrohn_lower"
    condition_tab <- table(condition)
    group_evidence <- paste0("Healthy=", condition_tab[["Healthy"]], ";ActiveCrohn=", condition_tab[["ActiveCrohn"]])
    if (min(condition_tab) < 3L) stop("AC disease group below n=3: ", lineage)
  } else {
    stop("Unknown analysis: ", analysis)
  }

  rownames(design) <- meta$sample_column
  if (!(coefficient %in% colnames(design))) stop("Missing coefficient: ", coefficient)
  y0 <- DGEList(counts = counts)
  keep <- filterByExpr(y0, design = design)
  if (sum(keep) < 1000L) stop("Too few tested genes: ", analysis, " ", lineage)
  y <- y0[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = coefficient)
  raw <- topTags(qlf, n = Inf, sort.by = "none")$table
  annotation <- gene_annotation[match(rownames(raw), gene_annotation$unique_gene_id), , drop = FALSE]
  if (anyNA(annotation$unique_gene_id)) stop("Gene annotation lookup failed")
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
    analysis = analysis,
    design = paste(colnames(design), collapse = "+"),
    coefficient = coefficient,
    group_evidence = group_evidence,
    stringsAsFactors = FALSE
  )
  result <- result[order(result$PValue, -abs(result$logFC), na.last = TRUE), ]
  output <- file.path(de_dir, paste0(analysis, "_", lineage, "_edgeR.csv"))
  atomic_write_csv(result, output)
  ranked <- result[order(result$rank_metric, decreasing = TRUE, na.last = TRUE), c(
    "unique_gene_id", "gene_id", "gene_symbol", "logFC", "F", "PValue", "FDR", "rank_metric"
  )]
  names(ranked)[names(ranked) == "F"] <- "statistic_F"
  atomic_write_csv(ranked, file.path(ranked_dir, paste0(analysis, "_", lineage, "_ranked.csv")))

  log_cpm <- cpm(y, log = TRUE, prior.count = 2)
  variable <- order(apply(log_cpm, 1, var), decreasing = TRUE)[seq_len(min(1000L, nrow(log_cpm)))]
  pca <- prcomp(t(log_cpm[variable, , drop = FALSE]), center = TRUE, scale. = FALSE)
  colour <- if (analysis == "AC_ActiveCrohn_vs_Healthy") {
    ifelse(meta$condition == "Healthy", "#2878B5", "#C74440")
  } else {
    ifelse(meta$region == "TI", "#2878B5", "#C74440")
  }
  png(file.path(figure_dir, paste0(analysis, "_", lineage, "_PCA.png")), width = 1400, height = 1100, res = 150)
  plot(pca$x[, 1], pca$x[, 2], pch = 19, col = colour, xlab = "PCA1", ylab = "PCA2", main = paste(analysis, lineage))
  text(pca$x[, 1], pca$x[, 2], labels = meta$patient_id, pos = 3, cex = 0.55)
  dev.off()
  png(file.path(figure_dir, paste0(analysis, "_", lineage, "_volcano.png")), width = 1400, height = 1100, res = 150)
  de_colour <- ifelse(result$FDR < 0.05 & result$logFC > 0, "#C74440", ifelse(result$FDR < 0.05 & result$logFC < 0, "#2878B5", "#B8B8B8"))
  plot(result$logFC, -log10(pmax(result$PValue, .Machine$double.xmin)), pch = 19, cex = 0.4, col = de_colour,
       xlab = "logFC (positive per registered direction)", ylab = "-log10(PValue)", main = paste(analysis, lineage))
  abline(v = 0, lty = 2, col = "grey40")
  dev.off()

  qc <- data.frame(
    analysis = analysis,
    lineage = lineage,
    patient_id = meta$patient_id,
    sample_id = meta$sample_id,
    sample_column = meta$sample_column,
    region = meta$region,
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
  summary <- data.frame(
    analysis = analysis,
    lineage = lineage,
    design = if (analysis == "Healthy_AC_vs_TI") "~patient_id + region" else if (analysis == "ActiveCrohn_AC_vs_TI") "~region" else "~condition",
    coefficient = coefficient,
    group_evidence = group_evidence,
    tested_genes = nrow(result),
    FDR05 = sum(result$FDR < 0.05),
    positive_FDR05 = sum(result$FDR < 0.05 & result$logFC > 0),
    negative_FDR05 = sum(result$FDR < 0.05 & result$logFC < 0),
    positive_direction = positive_label,
    negative_direction = negative_label,
    stringsAsFactors = FALSE
  )
  log_line(sprintf("EDGER %s %s %s genes=%d FDR05=%d", analysis, lineage, group_evidence, nrow(result), summary$FDR05))
  list(summary = summary, qc = qc)
}

metadata <- read.csv(metadata_path, check.names = FALSE)
names(metadata)[1] <- "analysis"
gene_annotation <- read.csv(gene_path, check.names = FALSE)
names(gene_annotation)[1] <- "unique_gene_id"
if (!("gene_id" %in% names(gene_annotation))) gene_annotation$gene_id <- gene_annotation$unique_gene_id
if (anyDuplicated(gene_annotation$unique_gene_id)) stop("Duplicated frozen Human gene IDs")

analyses <- c("Healthy_AC_vs_TI", "ActiveCrohn_AC_vs_TI", "AC_ActiveCrohn_vs_Healthy")
lineages <- c("Epithelial", "T", "Myeloid")
summaries <- list()
qcs <- list()
for (analysis in analyses) {
  for (lineage in lineages) {
    result <- fit_one(analysis, lineage, metadata, gene_annotation)
    key <- paste(analysis, lineage, sep = "::")
    summaries[[key]] <- result$summary
    qcs[[key]] <- result$qc
  }
}
summary_frame <- do.call(rbind, summaries)
qc_frame <- do.call(rbind, qcs)
rownames(summary_frame) <- NULL
rownames(qc_frame) <- NULL
atomic_write_csv(summary_frame, summary_path)
atomic_write_csv(qc_frame, sample_qc_path)

checks <- data.frame(
  criterion = c(
    "Nine preregistered edgeR QL fits completed",
    "Healthy region comparison remained paired",
    "Active Crohn region comparison remained explicitly unpaired",
    "AC disease effect remained secondary and donor-level",
    "Every output retains all tested genes"
  ),
  status = c(
    ifelse(nrow(summary_frame) == 9L, "PASS", "FAIL"),
    ifelse(all(summary_frame$group_evidence[summary_frame$analysis == "Healthy_AC_vs_TI"] == "paired_donors=15"), "PASS", "FAIL"),
    ifelse(all(grepl("TI=11;AC=", summary_frame$group_evidence[summary_frame$analysis == "ActiveCrohn_AC_vs_TI"])), "PASS", "FAIL"),
    ifelse(all(grepl("Healthy=16;ActiveCrohn=", summary_frame$group_evidence[summary_frame$analysis == "AC_ActiveCrohn_vs_Healthy"])), "PASS", "FAIL"),
    ifelse(all(summary_frame$tested_genes >= 1000L), "PASS", "FAIL")
  ),
  evidence = c(
    paste0("fits=", nrow(summary_frame)),
    "15 paired donors per lineage; ~patient_id+region",
    paste(unique(summary_frame$group_evidence[summary_frame$analysis == "ActiveCrohn_AC_vs_TI"]), collapse = ";"),
    paste(unique(summary_frame$group_evidence[summary_frame$analysis == "AC_ActiveCrohn_vs_Healthy"]), collapse = ";"),
    paste0("min_tested=", min(summary_frame$tested_genes))
  ),
  stringsAsFactors = FALSE
)
checks$gate <- if (all(checks$status == "PASS")) "GREEN_REGION_EDGER_COMPLETE" else "YELLOW_REGION_EDGER_LIMITED"
atomic_write_csv(checks, gate_path)
if (!all(checks$status == "PASS")) stop("Phase8 edgeR gate failed")
log_line(paste("EDGER_GATE", checks$gate[1]))
