#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
pb_root <- file.path(project_root, "06_pseudobulk")
lodo_root <- file.path(pb_root, "lodo")
gsea_root <- file.path(lodo_root, "gsea")
leading_root <- file.path(lodo_root, "leading_edges")
meta_root <- file.path(lodo_root, "metadata")
phase4_gsea <- file.path(pb_root, "gsea")
phase4_mapping <- file.path(phase4_gsea, "mapping", "phase4_feature_to_geneset_id.csv")
log_path <- file.path(project_root, "logs", "phase5_gsea_lodo.log")
for (path in c(gsea_root, leading_root, meta_root, dirname(log_path))) dir.create(path, recursive = TRUE, showWarnings = FALSE)

log_connection <- file(log_path, open = "wt", encoding = "UTF-8")
sink(log_connection, type = "output", split = TRUE)
sink(log_connection, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(log_connection), silent = TRUE)
}, add = TRUE)

lineages <- c("Epithelial", "T", "Myeloid")
collection_specs <- list(
  Hallmark = list(
    file = file.path(phase4_gsea, "gene_sets", "hallmark_dog_ortholog_projected_membership.tsv.gz"),
    tag = "Hallmark"
  ),
  Reactome = list(
    file = file.path(phase4_gsea, "gene_sets", "reactome_canis_familiaris_gene_sets.tsv.gz"),
    tag = "Reactome"
  ),
  GO_BP = list(
    file = file.path(phase4_gsea, "gene_sets", "go_bp_canis_familiaris_gene_sets.tsv.gz"),
    tag = "GOBP"
  )
)
min_size <- 15L
max_size <- 500L
fgsea_eps <- 1e-10
base_seed <- 20260824L

atomic_fwrite <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  fwrite(x, temporary)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) stop("Could not atomically write: ", path)
}

jaccard_strings <- function(a, b) {
  aa <- if (is.na(a) || !nzchar(a)) character() else unique(strsplit(a, ";", fixed = TRUE)[[1L]])
  bb <- if (is.na(b) || !nzchar(b)) character() else unique(strsplit(b, ";", fixed = TRUE)[[1L]])
  union_n <- length(union(aa, bb))
  if (union_n == 0L) return(NA_real_)
  length(intersect(aa, bb)) / union_n
}

top20_jaccard <- function(full, lodo) {
  f <- full[order(padj, -abs(NES), pathway_id), head(pathway_id, 20L)]
  l <- lodo[order(padj, -abs(NES), pathway_id), head(pathway_id, 20L)]
  if (length(union(f, l)) == 0L) return(NA_real_)
  length(intersect(f, l)) / length(union(f, l))
}

flatten_fgsea <- function(fg, info, lineage, collection, model, lodo_run, short_label, dropped_donor) {
  setnames(fg, "pathway", "pathway_id")
  fg[, padj := p.adjust(pval, method = "BH")]
  result <- merge(fg, info, by = "pathway_id", all.x = TRUE, sort = FALSE)
  result[, `:=`(
    lineage = lineage, collection = collection, model = model,
    lodo_run = lodo_run, short_label = short_label, dropped_donor = dropped_donor,
    direction = fifelse(NES > 0, "CIE_higher", "CIE_lower"),
    leading_edge_n = lengths(leadingEdge),
    leading_edge_target_ids = vapply(leadingEdge, function(x) paste(as.character(x), collapse = ";"), character(1))
  )]
  result[, leadingEdge := NULL]
  result[, abs_NES_sort := abs(NES)]
  setorder(result, padj, -abs_NES_sort, pathway_id)
  result[, abs_NES_sort := NULL]
  result
}

mapping <- fread(phase4_mapping)
mapping <- mapping[!is.na(target_geneset_id) & nzchar(as.character(target_geneset_id)), .(
  feature_id, target_geneset_id = as.character(target_geneset_id), target_symbol
)]
run_status <- fread(file.path(meta_root, "LODO_RUN_STATUS.csv"))
if (nrow(run_status) != 21L || any(run_status$status != "SUCCESS")) stop("Gene LODO status is not 21/21 SUCCESS")

collection_cache <- list()
for (collection in names(collection_specs)) {
  spec <- collection_specs[[collection]]
  members <- unique(fread(spec$file), by = c("pathway_id", "target_geneset_id"))
  members[, target_geneset_id := as.character(target_geneset_id)]
  info <- unique(members[, .(pathway_id, pathway_name, mapping_type, source_version)], by = "pathway_id")
  collection_cache[[collection]] <- list(info = info, pathways = split(members$target_geneset_id, members$pathway_id))
  cat("Frozen", collection, ":", length(collection_cache[[collection]]$pathways), "pathways\n")
}

status_rows <- list()
qc_rows <- list()
cat("Phase5 common-universe LODO GSEA\n")
cat("Start:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
cat("fgsea:", as.character(packageVersion("fgsea")), "\n")
cat("Parameters: fgseaMultilevel; minSize=15; maxSize=500; eps=1e-10; scoreType=std; nproc=1\n")

for (lineage_index in seq_along(lineages)) {
  lineage <- lineages[[lineage_index]]
  lineage_value <- lineage
  full <- fread(file.path(pb_root, "ranked_lists", paste0(lineage, "_ranked.csv")))
  setnames(full,
           c("gene_symbol", "unique_var_name", "logFC", "statistic_F", "PValue", "FDR", "rank_metric"),
           c("gene_symbol_full", "unique_var_name_full", "logFC_full", "F_full", "PValue_full", "FDR_full", "rank_metric_full"))
  lineage_runs <- run_status[which(run_status[["lineage"]] == lineage_value), ]
  if (nrow(lineage_runs) != 7L || any(lineage_runs[["lineage"]] != lineage)) {
    stop("LODO run lineage subset failed for ", lineage)
  }
  for (run_index in seq_len(nrow(lineage_runs))) {
    run <- lineage_runs[run_index]
    lodo <- fread(run$ranked_file)
    if (!all(lodo$lineage == lineage)) stop("Ranked file lineage mismatch: ", run$ranked_file)
    setnames(lodo,
             c("gene_symbol", "unique_var_name", "logFC", "statistic_F", "PValue", "FDR", "rank_metric"),
             c("gene_symbol_lodo", "unique_var_name_lodo", "logFC_lodo", "F_lodo", "PValue_lodo", "FDR_lodo", "rank_metric_lodo"))
    common <- merge(full, lodo, by = "feature_id", all = FALSE)
    common_feature_n <- nrow(common)
    common <- merge(common, mapping, by = "feature_id", all = FALSE)
    common[, abs_rank_full := abs(rank_metric_full)]
    setorder(common, target_geneset_id, -abs_rank_full, PValue_full, feature_id)
    common[, collision_rank := seq_len(.N), by = target_geneset_id]
    common[, kept_after_collision := collision_rank == 1L]
    common_kept <- common[kept_after_collision == TRUE]
    setorder(common_kept, -rank_metric_full, feature_id)
    full_stats <- common_kept$rank_metric_full
    lodo_stats <- common_kept$rank_metric_lodo
    names(full_stats) <- common_kept$target_geneset_id
    names(lodo_stats) <- common_kept$target_geneset_id
    stopifnot(identical(names(full_stats), names(lodo_stats)), !anyDuplicated(names(full_stats)))

    pair_dir <- file.path(gsea_root, lineage, run$short_label)
    rank_audit_dir <- file.path(meta_root, "common_rank_mappings", lineage)
    dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(rank_audit_dir, recursive = TRUE, showWarnings = FALSE)
    atomic_fwrite(common, file.path(rank_audit_dir, paste0(lineage, "_", run$short_label, "_common_rank_mapping.csv")))
    cat("\n===", lineage, run$lodo_run, "drop", run$dropped_donor,
        "; common features", common_feature_n, "; target genes", length(full_stats), "===\n")

    for (collection_index in seq_along(collection_specs)) {
      collection <- names(collection_specs)[[collection_index]]
      spec <- collection_specs[[collection]]
      cached <- collection_cache[[collection]]
      key <- paste(lineage, run$lodo_run, collection, sep = "|")
      run_start <- Sys.time()
      pair_status <- tryCatch({
        seed <- base_seed + lineage_index * 10000L + run_index * 100L + collection_index
        cat("Running", collection, "Full7-common...")
        set.seed(seed)
        fg_full <- as.data.table(fgseaMultilevel(
          pathways = cached$pathways, stats = full_stats, minSize = min_size, maxSize = max_size,
          eps = fgsea_eps, scoreType = "std", nproc = 1L
        ))
        full_result <- flatten_fgsea(
          fg_full, cached$info, lineage, collection, "Full7_common",
          run$lodo_run, run$short_label, run$dropped_donor
        )
        cat("Running", collection, "LODO-common...")
        set.seed(seed)
        fg_lodo <- as.data.table(fgseaMultilevel(
          pathways = cached$pathways, stats = lodo_stats, minSize = min_size, maxSize = max_size,
          eps = fgsea_eps, scoreType = "std", nproc = 1L
        ))
        lodo_result <- flatten_fgsea(
          fg_lodo, cached$info, lineage, collection, "LODO_common",
          run$lodo_run, run$short_label, run$dropped_donor
        )

        prefix <- paste0(lineage, "_", run$short_label, "_", spec$tag)
        full_path <- file.path(pair_dir, paste0(prefix, "_Full7_common_GSEA.csv"))
        lodo_path <- file.path(pair_dir, paste0(prefix, "_LODO_common_GSEA.csv"))
        atomic_fwrite(full_result, full_path)
        atomic_fwrite(lodo_result, lodo_path)

        full_cmp <- full_result[, .(
          pathway_id, pathway_name, mapping_type, source_version,
          pair_full_ES = ES, pair_full_NES = NES, pair_full_pval = pval, pair_full_padj = padj,
          pair_full_size = size, pair_full_direction = direction,
          pair_full_leading_edge_n = leading_edge_n, pair_full_leading_edge_target_ids = leading_edge_target_ids
        )]
        lodo_cmp <- lodo_result[, .(
          pathway_id, pathway_name_lodo = pathway_name,
          lodo_ES = ES, lodo_NES = NES, lodo_pval = pval, lodo_padj = padj,
          lodo_size = size, lodo_direction = direction,
          lodo_leading_edge_n = leading_edge_n, lodo_leading_edge_target_ids = leading_edge_target_ids
        )]
        comparison <- merge(full_cmp, lodo_cmp, by = "pathway_id", all = TRUE)
        comparison[is.na(pathway_name), pathway_name := pathway_name_lodo]
        comparison[, pathway_name_lodo := NULL]
        comparison[, `:=`(
          lineage = lineage, collection = collection, lodo_run = run$lodo_run,
          short_label = run$short_label, dropped_donor = run$dropped_donor,
          common_feature_universe_n = common_feature_n,
          common_target_gene_universe_n = length(full_stats),
          delta_NES = lodo_NES - pair_full_NES,
          same_direction = !is.na(lodo_NES) & !is.na(pair_full_NES) & sign(lodo_NES) == sign(pair_full_NES),
          sign_flip = !is.na(lodo_NES) & !is.na(pair_full_NES) & sign(lodo_NES) != sign(pair_full_NES),
          strong_effect_flip = !is.na(lodo_NES) & !is.na(pair_full_NES) & sign(lodo_NES) != sign(pair_full_NES) &
            abs(lodo_NES) >= 1.5 & abs(pair_full_NES) >= 1.5,
          leading_edge_jaccard = mapply(jaccard_strings, pair_full_leading_edge_target_ids, lodo_leading_edge_target_ids),
          pathway_size_change = lodo_size - pair_full_size
        )]
        comparison_path <- file.path(pair_dir, paste0(prefix, "_PAIR_COMPARISON.csv"))
        atomic_fwrite(comparison, comparison_path)

        leading_dir <- file.path(leading_root, lineage)
        dir.create(leading_dir, recursive = TRUE, showWarnings = FALSE)
        leading_output <- comparison[, .(
          lineage, collection, lodo_run, short_label, dropped_donor, pathway_id, pathway_name,
          pair_full_NES, lodo_NES, pair_full_leading_edge_n, lodo_leading_edge_n,
          pair_full_leading_edge_target_ids, lodo_leading_edge_target_ids, leading_edge_jaccard
        )]
        atomic_fwrite(leading_output, file.path(leading_dir, paste0(prefix, "_LEADING_EDGE_COMPARISON.csv")))

        comparable <- comparison[!is.na(pair_full_NES) & !is.na(lodo_NES)]
        qc_rows[[key]] <- data.frame(
          lineage = lineage, collection = collection, lodo_run = run$lodo_run,
          short_label = run$short_label, dropped_donor = run$dropped_donor,
          full7_tested_features = nrow(full), lodo_tested_features = nrow(lodo),
          common_feature_universe_n = common_feature_n, common_mapping_candidates_n = nrow(common),
          common_target_gene_universe_n = length(full_stats), collapsed_features_n = nrow(common) - nrow(common_kept),
          pathways_full = nrow(full_result), pathways_lodo = nrow(lodo_result),
          pathways_comparable = nrow(comparable), pathways_not_comparable = nrow(comparison) - nrow(comparable),
          top20_jaccard = top20_jaccard(full_result, lodo_result),
          nes_spearman = cor(comparable$pair_full_NES, comparable$lodo_NES, method = "spearman"),
          nes_pearson = cor(comparable$pair_full_NES, comparable$lodo_NES, method = "pearson"),
          direction_concordance = mean(comparable$same_direction),
          max_abs_pathway_size_change = max(abs(comparable$pathway_size_change), na.rm = TRUE),
          status = "SUCCESS", stringsAsFactors = FALSE
        )
        cat("SUCCESS", nrow(comparable), "pathways; Spearman",
            qc_rows[[key]]$nes_spearman, "; direction", qc_rows[[key]]$direction_concordance, "\n")
        list(status = "SUCCESS", error = "", full_path = full_path, lodo_path = lodo_path, comparison_path = comparison_path)
      }, error = function(e) {
        cat("LODO_GSEA_RUN_FAILED:", conditionMessage(e), "\n")
        list(status = "LODO_GSEA_RUN_FAILED", error = conditionMessage(e), full_path = "", lodo_path = "", comparison_path = "")
      })
      status_rows[[key]] <- data.frame(
        lineage = lineage, collection = collection, lodo_run = run$lodo_run,
        short_label = run$short_label, dropped_donor = run$dropped_donor,
        status = pair_status$status, error = pair_status$error,
        full7_common_file = pair_status$full_path, lodo_common_file = pair_status$lodo_path,
        comparison_file = pair_status$comparison_path,
        elapsed_seconds = as.numeric(difftime(Sys.time(), run_start, units = "secs")),
        stringsAsFactors = FALSE
      )
    }
  }
}

status <- rbindlist(status_rows, fill = TRUE)
qc <- rbindlist(qc_rows, fill = TRUE)
atomic_fwrite(status, file.path(meta_root, "LODO_GSEA_RUN_STATUS.csv"))
atomic_fwrite(qc, file.path(meta_root, "LODO_GSEA_PAIR_QC.csv"))
software <- data.table(
  software = c("R", "fgsea", "data.table"),
  version = c(R.version.string, as.character(packageVersion("fgsea")), as.character(packageVersion("data.table"))),
  role = c("Phase5 GSEA LODO runtime", "Common-universe preranked GSEA", "I/O and aggregation")
)
atomic_fwrite(software, file.path(project_root, "00_metadata", "PHASE5_SOFTWARE_VERSIONS_GSEA_R.csv"))
cat("\nCompleted paired GSEA:", sum(status$status == "SUCCESS"), "/", nrow(status), "successful lineage/donor/collection pairs.\n")
if (any(status$status != "SUCCESS")) quit(status = 2L)
