#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(fgsea)
  library(data.table)
})

options(stringsAsFactors = FALSE)

ROOT <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
OUT <- file.path(ROOT, "WORKING", "12_SELECTION_AWARE_NULL")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SELECTED_PATH <- file.path(ROOT, "10_cross_species_baseline", "summary", "PHASE10_PROGRAM_BASELINE_MATRIX.csv")
META_PATH <- file.path(ROOT, "07_human_crohn", "GSE266616", "pseudobulk", "metadata", "HUMAN_PSEUDOBULK_SAMPLE_METADATA.csv")
ANNOTATION_PATH <- file.path(ROOT, "07_human_crohn", "GSE266616", "pseudobulk", "metadata", "HUMAN_PSEUDOBULK_GENE_ANNOTATION.csv")
FILTER_SUMMARY_PATH <- file.path(ROOT, "07_human_crohn", "GSE266616", "pseudobulk", "metadata", "HUMAN_EDGER_FILTER_SUMMARY.csv")
MATRIX_DIR <- file.path(ROOT, "07_human_crohn", "GSE266616", "pseudobulk", "matrices")
RANK_DIR <- file.path(ROOT, "10_cross_species_baseline", "pathway_effects", "paired_ranks")
GENESET_DIR <- file.path(ROOT, "07_human_crohn", "GSE266616", "gsea", "gene_sets")

LINEAGES <- c("Epithelial", "T", "Myeloid")
COLLECTIONS <- c("Hallmark", "Reactome", "GO_BP")
GENESET_PATHS <- c(
  Hallmark = file.path(GENESET_DIR, "HUMAN_HALLMARK_GENE_SETS.tsv.gz"),
  Reactome = file.path(GENESET_DIR, "HUMAN_REACTOME_GENE_SETS.tsv.gz"),
  GO_BP = file.path(GENESET_DIR, "HUMAN_GO_BP_GENE_SETS.tsv.gz")
)
SEED <- 20261017L
PILOT_N <- 10L
TARGET_N <- 1000L
MAX_PROJECTED_SECONDS <- 3600
MIN_SIZE <- 15L
MAX_SIZE <- 500L

atomic_fwrite <- function(frame, path) {
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  fwrite(frame, tmp, sep = "	", quote = FALSE, na = "NA")
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not publish ", path)
}

read_count_matrix <- function(path) {
  if (!file.exists(path)) stop("Missing count matrix: ", path)
  frame <- read.csv(gzfile(path), check.names = FALSE)
  names(frame)[1] <- "unique_gene_id"
  required <- c("unique_gene_id", "gene_symbol")
  if (!all(required %in% names(frame))) stop("Count schema failure: ", path)
  if (anyDuplicated(frame$unique_gene_id)) stop("Duplicate unique_gene_id: ", path)
  counts <- as.matrix(frame[, setdiff(names(frame), required), drop = FALSE])
  storage.mode(counts) <- "integer"
  rownames(counts) <- frame$unique_gene_id
  if (anyNA(counts) || any(counts < 0L)) stop("Invalid counts: ", path)
  counts
}

load_gene_sets <- function() {
  lapply(GENESET_PATHS, function(path) {
    x <- fread(path)
    required <- c("pathway_id", "pathway_name", "gene_symbol")
    if (!all(required %in% names(x))) stop("Gene-set schema failure: ", path)
    unique(x[, .(pathway_id, pathway_name, gene_symbol)])
  })
}

prepare_lineage <- function(lineage, selected, metadata, annotation, filter_summary, source_sets) {
  lineage_value <- lineage
  message("PREP_START ", lineage)
  counts <- read_count_matrix(file.path(MATRIX_DIR, paste0("human_", lineage, "_counts.csv.gz")))
  meta <- metadata[lineage == lineage_value & eligible == "YES"]
  if (nrow(meta) != 27L || uniqueN(meta$patient_id) != 27L) stop("Expected 27 donors for ", lineage)
  if (!setequal(colnames(counts), meta$patient_id)) stop("Count/metadata identity mismatch for ", lineage)
  setkey(meta, patient_id)
  meta <- meta[colnames(counts)]
  observed_condition <- factor(meta$condition, levels = c("Healthy", "ActiveCrohn"))
  if (sum(observed_condition == "Healthy") != 16L || sum(observed_condition == "ActiveCrohn") != 11L) {
    stop("Observed group-size mismatch for ", lineage)
  }
  observed_design <- model.matrix(~ observed_condition)
  y0 <- DGEList(counts = counts)
  keep <- filterByExpr(y0, design = observed_design)
  expected_after <- filter_summary[lineage == lineage_value, genes_after]
  if (length(expected_after) != 1L || sum(keep) != expected_after) {
    stop("Frozen filter count mismatch for ", lineage, ": observed=", sum(keep), " expected=", paste(expected_after, collapse = ","))
  }
  y <- y0[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  ann <- annotation[match(rownames(y$counts), unique_gene_id)]
  if (anyNA(ann$unique_gene_id)) stop("Gene annotation failure for ", lineage)
  filtered_map <- data.table(unique_gene_id = rownames(y$counts), human_gene_id = ann$gene_id)

  rank_path <- file.path(RANK_DIR, paste0(lineage, "_Human_common_one2one_rank.csv"))
  frozen_rank <- fread(rank_path)
  if (anyDuplicated(frozen_rank$human_gene_id)) stop("Frozen rank identity duplicate: ", lineage)
  gene_map <- merge(
    filtered_map,
    frozen_rank[, .(human_gene_id, human_symbol)],
    by = "human_gene_id",
    all = FALSE
  )
  if (nrow(gene_map) != nrow(frozen_rank)) {
    stop("Frozen common-universe mapping count mismatch for ", lineage, ": ", nrow(gene_map), " versus ", nrow(frozen_rank))
  }
  if (anyDuplicated(gene_map$unique_gene_id) || anyDuplicated(gene_map$human_gene_id)) {
    stop("Common-universe mapping collision for ", lineage)
  }

  selected_lineage <- selected[lineage == lineage_value & testable == TRUE]
  pathway_list <- list()
  pathway_meta <- list()
  for (i in seq_len(nrow(selected_lineage))) {
    row <- selected_lineage[i]
    collection <- row$representative_collection
    pathway_id <- row$representative_pathway
    pathway_id_value <- pathway_id
    source <- source_sets[[collection]]
    symbol_map <- unique(frozen_rank[, .(human_gene_id, gene_symbol = human_symbol)])
    membership <- merge(
      unique(source[pathway_id == pathway_id_value, .(pathway_id, pathway_name, gene_symbol)]),
      symbol_map,
      by = "gene_symbol",
      allow.cartesian = TRUE
    )
    ids <- sort(unique(intersect(membership$human_gene_id, frozen_rank$human_gene_id)))
    if (length(ids) < MIN_SIZE || length(ids) > MAX_SIZE) {
      stop("Selected pathway size outside frozen limits: ", lineage, " / ", pathway_id, " n=", length(ids))
    }
    key <- paste(lineage, collection, pathway_id, sep = "::")
    message("PREP_PATH ", key)
    pathway_list[[key]] <- ids
    pathway_meta[[key]] <- data.table(
      unit_key = key,
      lineage = lineage,
      program = row$program,
      collection = collection,
      pathway_id = pathway_id,
      dog_NES = row$dog_common1to1_NES,
      frozen_human_NES = row$human_common1to1_NES,
      pathway_size = length(ids)
    )
  }

  list(
    lineage = lineage,
    counts = counts,
    meta = meta,
    observed_condition = observed_condition,
    y = y,
    gene_map = gene_map,
    frozen_rank = frozen_rank,
    pathways = pathway_list,
    pathway_meta = rbindlist(pathway_meta),
    genes_before = nrow(counts),
    genes_after = sum(keep)
  )
}

fit_rank <- function(prepared, labels) {
  condition <- factor(labels, levels = c("Healthy", "ActiveCrohn"))
  design <- model.matrix(~ condition)
  y <- estimateDisp(prepared$y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "conditionActiveCrohn")
  raw <- topTags(qlf, n = Inf, sort.by = "none")$table
  result <- data.table(
    unique_gene_id = rownames(raw),
    logFC = raw$logFC,
    F = raw$F,
    rank_metric = sign(raw$logFC) * sqrt(pmax(raw$F, 0))
  )
  rank <- merge(prepared$gene_map, result, by = "unique_gene_id", all = FALSE)
  if (nrow(rank) != nrow(prepared$frozen_rank)) {
    stop("Permuted common-rank size mismatch for ", prepared$lineage)
  }
  if (anyDuplicated(rank$human_gene_id)) stop("Permuted common-rank duplicate for ", prepared$lineage)
  rank <- rank[is.finite(rank_metric)]
  setorder(rank, -rank_metric, human_gene_id)
  stats <- rank$rank_metric
  names(stats) <- rank$human_gene_id
  stats
}

enrichment_signs <- function(prepared, stats, permutation_id) {
  output <- copy(prepared[["pathway_meta"]])
  pathway_map <- prepared[["pathways"]]
  pathway_keys <- as.character(output[["unit_key"]])
  if (permutation_id == 0L) message("DEBUG pathway_map class=", paste(class(pathway_map), collapse=","), "; typeof=", typeof(pathway_map), "; keys=", paste(names(pathway_map), collapse="|"))
  human_es <- vapply(pathway_keys, function(k) {
    ids <- pathway_map[[k]]
    if (is.null(ids)) stop("Missing prepared pathway key: ", k)
    selected_stats <- sort(match(ids, names(stats)))
    selected_stats <- selected_stats[!is.na(selected_stats)]
    if (length(selected_stats) < MIN_SIZE || length(selected_stats) > MAX_SIZE) {
      stop("Permuted pathway size changed: ", k)
    }
    as.numeric(calcGseaStat(
      stats = stats,
      selectedStats = selected_stats,
      gseaParam = 1,
      returnAllExtremes = FALSE
    ))
  }, numeric(1))
  output[, human_ES := human_es]
  output[, human_sign := sign(human_ES)]
  output[, dog_sign := sign(dog_NES)]
  output[, concordant := human_sign == dog_sign]
  output[, permutation_id := permutation_id]
  output
}
generate_assignments <- function(participants, n, seed) {
  set.seed(seed)
  nova <- participants[platform == "Illumina NovaSeq 6000", patient_id]
  nextseq <- participants[platform == "NextSeq 500", patient_id]
  if (length(nova) != 25L || length(nextseq) != 2L) stop("Platform block size mismatch")
  assignments <- list()
  observed_key <- paste(sort(participants[condition == "ActiveCrohn", patient_id]), collapse = "|")
  seen <- new.env(hash = TRUE, parent = emptyenv())
  i <- 1L
  while (length(assignments) < n) {
    active <- sort(c(sample(nova, 10L, replace = FALSE), sample(nextseq, 1L, replace = FALSE)))
    key <- paste(active, collapse = "|")
    if (exists(key, envir = seen, inherits = FALSE)) next
    assign(key, TRUE, envir = seen)
    labels <- ifelse(participants$patient_id %in% active, "ActiveCrohn", "Healthy")
    assignments[[length(assignments) + 1L]] <- list(
      permutation_id = i,
      labels = setNames(labels, participants$patient_id),
      assignment_key = key,
      equals_observed = identical(key, observed_key)
    )
    i <- i + 1L
  }
  assignments
}

run_assignment <- function(assignment, prepared_by_lineage, participant_order) {
  pieces <- lapply(LINEAGES, function(lineage) {
    prepared <- prepared_by_lineage[[lineage]]
    labels <- assignment$labels[prepared$meta$patient_id]
    stats <- fit_rank(prepared, labels)
    enrichment_signs(prepared, stats, assignment$permutation_id)
  })
  program <- rbindlist(pieces)
  stats <- data.table(
    permutation_id = assignment$permutation_id,
    assignment_key = assignment$assignment_key,
    equals_observed = assignment$equals_observed,
    concordant_n = sum(program$concordant),
    denominator = nrow(program),
    concordance = mean(program$concordant),
    epithelial_concordant_n = sum(program[lineage == "Epithelial", concordant]),
    epithelial_n = sum(program$lineage == "Epithelial"),
    t_concordant_n = sum(program[lineage == "T", concordant]),
    t_n = sum(program$lineage == "T"),
    myeloid_concordant_n = sum(program[lineage == "Myeloid", concordant]),
    myeloid_n = sum(program$lineage == "Myeloid")
  )
  list(stat = stats, program = program)
}

main <- function() {
  start_time <- Sys.time()
  selected <- fread(SELECTED_PATH)
  selected[, testable := is.finite(dog_common1to1_NES) & is.finite(human_common1to1_NES)]
  if (nrow(selected) != 25L || sum(selected$testable) != 23L) stop("Selected-set testability mismatch")

  metadata <- fread(META_PATH)
  setnames(metadata, names(metadata)[1], "analysis")
  metadata <- metadata[eligible == "YES"]
  annotation <- fread(ANNOTATION_PATH)
  setnames(annotation, names(annotation)[1], "unique_gene_id")
  if (anyDuplicated(annotation$unique_gene_id) || anyDuplicated(annotation$gene_id)) {
    stop("Human annotation identity duplication")
  }
  filter_summary <- fread(FILTER_SUMMARY_PATH)
  source_sets <- load_gene_sets()

  participants <- unique(metadata[, .(patient_id, condition, platform)])
  if (nrow(participants) != 27L || anyDuplicated(participants$patient_id)) stop("Participant identity failure")
  setorder(participants, patient_id)
  block_table <- participants[, .N, by = .(platform, condition)]
  expected_blocks <- data.table(
    platform = c("Illumina NovaSeq 6000", "Illumina NovaSeq 6000", "NextSeq 500", "NextSeq 500"),
    condition = c("ActiveCrohn", "Healthy", "ActiveCrohn", "Healthy"),
    N = c(10L, 15L, 1L, 1L)
  )
  setorder(block_table, platform, condition)
  setorder(expected_blocks, platform, condition)
  if (!identical(block_table, expected_blocks)) stop("Exchangeability block audit failed")

  prepared_by_lineage <- setNames(
    lapply(LINEAGES, function(lineage) {
      prepare_lineage(lineage, selected, metadata, annotation, filter_summary, source_sets)
    }),
    LINEAGES
  )

  observed_pieces <- lapply(LINEAGES, function(lineage) {
    prepared <- prepared_by_lineage[[lineage]]
    stats <- fit_rank(prepared, prepared$observed_condition)
    enrichment_signs(prepared, stats, 0L)
  })
  observed <- rbindlist(observed_pieces)
  observed[, frozen_human_sign := sign(frozen_human_NES)]
  observed[, human_sign_reproduced := human_sign == frozen_human_sign]
  observed[, observed_concordance_reproduced := human_sign == dog_sign]
  atomic_fwrite(observed, file.path(OUT, "NULL_B_OBSERVED_SIGN_REPRODUCTION.tsv"))
  if (nrow(observed) != 23L || !all(observed$human_sign_reproduced)) {
    stop("Observed-label sign reproduction failed")
  }
  observed_s <- mean(observed$observed_concordance_reproduced)
  frozen_s <- mean(sign(selected[testable == TRUE, dog_common1to1_NES]) == sign(selected[testable == TRUE, human_common1to1_NES]))
  if (!isTRUE(all.equal(observed_s, frozen_s, tolerance = 0))) {
    stop("Observed primary statistic reproduction failed")
  }

  assignments <- generate_assignments(participants, TARGET_N, SEED)
  pilot_start <- Sys.time()
  pilot_results <- lapply(assignments[seq_len(PILOT_N)], run_assignment, prepared_by_lineage, participants$patient_id)
  pilot_seconds <- as.numeric(difftime(Sys.time(), pilot_start, units = "secs"))
  projected_seconds <- pilot_seconds / PILOT_N * TARGET_N
  pilot_stats <- rbindlist(lapply(pilot_results, function(x) x$stat))
  pilot_program <- rbindlist(lapply(pilot_results, function(x) x$program))
  atomic_fwrite(pilot_stats, file.path(OUT, "NULL_B_PILOT_STATISTICS.tsv"))
  atomic_fwrite(pilot_program, file.path(OUT, "NULL_B_PILOT_PROGRAM_SIGNS.tsv"))

  feasible <- is.finite(projected_seconds) && projected_seconds <= MAX_PROJECTED_SECONDS
  if (!feasible) {
    summary <- data.table(
      status = "STOPPED_AFTER_PILOT_RUNTIME_INFEASIBLE",
      observed_concordant_n = sum(observed$observed_concordance_reproduced),
      observed_denominator = nrow(observed),
      observed_concordance = observed_s,
      pilot_n = PILOT_N,
      pilot_seconds = pilot_seconds,
      projected_seconds_for_1000 = projected_seconds,
      valid_permutations = PILOT_N,
      empirical_upper_tail = NA_real_,
      seed = SEED
    )
    atomic_fwrite(summary, file.path(OUT, "NULL_B_SUMMARY.tsv"))
    quit(save = "no", status = 0L)
  }

  remaining_results <- lapply(assignments[(PILOT_N + 1L):TARGET_N], run_assignment, prepared_by_lineage, participants$patient_id)
  all_results <- c(pilot_results, remaining_results)
  permutation_stats <- rbindlist(lapply(all_results, function(x) x$stat))
  permutation_program <- rbindlist(lapply(all_results, function(x) x$program))
  if (nrow(permutation_stats) != TARGET_N || uniqueN(permutation_stats$assignment_key) != TARGET_N) {
    stop("Permutation count or uniqueness failure")
  }
  tail_count <- sum(permutation_stats$concordance >= observed_s)
  empirical <- (1 + tail_count) / (1 + nrow(permutation_stats))
  summary <- data.table(
    status = "COMPLETE",
    observed_concordant_n = sum(observed$observed_concordance_reproduced),
    observed_denominator = nrow(observed),
    observed_concordance = observed_s,
    observed_signs_reproduced_n = sum(observed$human_sign_reproduced),
    pilot_n = PILOT_N,
    pilot_seconds = pilot_seconds,
    projected_seconds_for_1000 = projected_seconds,
    valid_permutations = nrow(permutation_stats),
    unique_assignments = uniqueN(permutation_stats$assignment_key),
    null_mean = mean(permutation_stats$concordance),
    null_sd = sd(permutation_stats$concordance),
    null_q025 = quantile(permutation_stats$concordance, 0.025),
    null_q05 = quantile(permutation_stats$concordance, 0.05),
    null_q50 = quantile(permutation_stats$concordance, 0.50),
    null_q95 = quantile(permutation_stats$concordance, 0.95),
    null_q975 = quantile(permutation_stats$concordance, 0.975),
    observed_minus_null_mean = observed_s - mean(permutation_stats$concordance),
    upper_tail_count = tail_count,
    empirical_upper_tail = empirical,
    seed = SEED,
    total_runtime_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
    r_version = as.character(getRversion()),
    edgeR_version = as.character(packageVersion("edgeR")),
    fgsea_version = as.character(packageVersion("fgsea")),
    data_table_version = as.character(packageVersion("data.table"))
  )
  atomic_fwrite(permutation_stats, file.path(OUT, "NULL_B_PERMUTATION_STATISTICS.tsv"))
  atomic_fwrite(permutation_program, file.path(OUT, "NULL_B_PERMUTATION_PROGRAM_SIGNS.tsv"))
  atomic_fwrite(summary, file.path(OUT, "NULL_B_SUMMARY.tsv"))
  print(summary)
}

status <- 0L
tryCatch(
  main(),
  error = function(e) {
    message("FATAL: ", conditionMessage(e))
    for (cl in tail(sys.calls(), 20L)) message("CALL: ", paste(deparse(cl), collapse = " "))
    status <<- 1L
  }
)
quit(save = "no", status = status)
