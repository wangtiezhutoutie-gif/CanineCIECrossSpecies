#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
phase8_root <- file.path(project_root, "08_region_control")
base <- file.path(phase8_root, "human_AC_TI")
ranked_dir <- file.path(base, "ranked_lists")
gsea_dir <- file.path(base, "gsea")
metadata_dir <- file.path(phase8_root, "metadata")
summary_dir <- file.path(phase8_root, "summary")
figure_dir <- file.path(phase8_root, "figures")
core_registry_path <- file.path(metadata_dir, "REGION_CORE_PROGRAM_REGISTRY.csv")
dog_registry_path <- file.path(metadata_dir, "DOG_FORWARD_REGION_REPRESENTATIVE_FREEZE.csv")
core_output <- file.path(summary_dir, "REGION_PROGRAM_STABILITY.csv")
matrix_output <- file.path(summary_dir, "REGION_INTERPRETATION_MATRIX.csv")
gsea_summary_output <- file.path(summary_dir, "HUMAN_REGION_REGISTERED_GSEA_SUMMARY.csv")
gate_path <- file.path(metadata_dir, "PHASE8A_REGION_PROGRAM_GATE.csv")
log_paths <- c(
  file.path(project_root, "logs", "phase8_region_analysis.log"),
  file.path(phase8_root, "logs", "phase8_region_analysis.log")
)

for (path in c(gsea_dir, metadata_dir, summary_dir, figure_dir, dirname(log_paths))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

gene_set_dir <- file.path(project_root, "07_human_crohn", "GSE266616", "gsea", "gene_sets")
collection_paths <- list(
  Hallmark = file.path(gene_set_dir, "HUMAN_HALLMARK_GENE_SETS.tsv.gz"),
  Reactome = file.path(gene_set_dir, "HUMAN_REACTOME_GENE_SETS.tsv.gz"),
  GO_BP = file.path(gene_set_dir, "HUMAN_GO_BP_GENE_SETS.tsv.gz")
)
analyses <- c("Healthy_AC_vs_TI", "ActiveCrohn_AC_vs_TI", "AC_ActiveCrohn_vs_Healthy")
lineages <- c("Epithelial", "T", "Myeloid")
base_seed <- 20260826L

log_line <- function(message) {
  line <- sprintf("[%s] %s", format(Sys.time(), tz = "UTC", usetz = TRUE), message)
  message(line)
  for (path in log_paths) cat(line, "\n", file = path, append = TRUE)
}

atomic_write_csv <- function(frame, path) {
  temp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  fwrite(frame, temp)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temp, path)) stop("Could not atomically publish ", path)
}

read_rank <- function(path) {
  rank <- fread(path)
  required <- c("unique_gene_id", "gene_id", "gene_symbol", "logFC", "statistic_F", "PValue", "FDR", "rank_metric")
  if (!all(required %in% names(rank))) stop("Rank schema failure: ", path)
  rank <- rank[!is.na(gene_symbol) & nzchar(gene_symbol) & is.finite(rank_metric)]
  rank[, abs_rank := abs(rank_metric)]
  setorder(rank, -abs_rank, PValue, unique_gene_id)
  rank <- unique(rank, by = "gene_symbol")
  rank[, abs_rank := NULL]
  stats <- rank$rank_metric
  names(stats) <- rank$gene_symbol
  stats
}

load_collection <- function(collection) {
  path <- collection_paths[[collection]]
  if (is.null(path) || !file.exists(path)) stop("Missing registered collection: ", collection)
  members <- unique(fread(path), by = c("pathway_id", "gene_symbol"))
  if (any(members$source_species != "Homo sapiens")) stop("Non-Human gene set entered Phase8")
  members
}

core_registry <- fread(core_registry_path)
dog_registry <- fread(dog_registry_path)
if (nrow(core_registry) != 36L || nrow(dog_registry) != 25L) stop("Frozen program registry cardinality changed")
members_cache <- lapply(names(collection_paths), load_collection)
names(members_cache) <- names(collection_paths)

# data.table evaluates column names before function arguments. Use this wrapper to avoid name collision.
registered_for_lineage <- function(lineage_value) {
  core <- core_registry[lineage == lineage_value, .(collection, pathway_id)]
  dog <- dog_registry[lineage == lineage_value, .(
    collection = representative_collection,
    pathway_id = representative_pathway_id
  )]
  unique(rbindlist(list(core, dog)))
}

run_registered <- function(analysis, lineage) {
  rank_path <- file.path(ranked_dir, paste0(analysis, "_", lineage, "_ranked.csv"))
  stats <- read_rank(rank_path)
  requested <- registered_for_lineage(lineage)
  outputs <- list()
  for (collection in unique(requested$collection)) {
    # Avoid data.table name collision using a scalar alias.
    collection_value <- collection
    ids <- requested[collection == collection_value, pathway_id]
    members <- members_cache[[collection_value]][pathway_id %in% ids]
    info <- unique(members[, .(pathway_id, pathway_name, source_version, source_species)], by = "pathway_id")
    pathways <- split(members$gene_symbol, members$pathway_id)
    missing <- setdiff(ids, names(pathways))
    if (length(missing)) stop("Registered pathways absent from frozen memberships: ", paste(missing, collapse = ";"))
    set.seed(base_seed + match(analysis, analyses) * 100L + match(lineage, lineages) * 10L + match(collection_value, names(collection_paths)))
    fg <- as.data.table(fgseaMultilevel(
      pathways = pathways,
      stats = stats,
      minSize = 15L,
      maxSize = 500L,
      eps = 1e-10,
      scoreType = "std",
      nproc = 1L
    ))
    setnames(fg, "pathway", "pathway_id")
    if (!setequal(fg$pathway_id, ids)) stop("One or more registered pathways were not assessable: ", analysis, " ", lineage, " ", collection_value)
    result <- merge(fg, info, by = "pathway_id", all.x = TRUE, sort = FALSE)
    result[, `:=`(
      analysis = analysis,
      lineage = lineage,
      collection = collection_value,
      leading_edge_n = lengths(leadingEdge),
      leading_edge = vapply(leadingEdge, function(value) paste(as.character(value), collapse = ";"), character(1))
    )]
    result[, leadingEdge := NULL]
    outputs[[collection_value]] <- result
  }
  combined <- rbindlist(outputs, fill = TRUE)
  combined[, registered_FDR := p.adjust(pval, method = "BH")]
  combined[, direction := ifelse(
    NES > 0,
    ifelse(analysis == "AC_ActiveCrohn_vs_Healthy", "ActiveCrohn_higher", "AC_higher"),
    ifelse(analysis == "AC_ActiveCrohn_vs_Healthy", "ActiveCrohn_lower", "AC_lower")
  )]
  setcolorder(combined, c(
    "analysis", "lineage", "collection", "pathway_id", "pathway_name", "ES", "NES", "pval", "registered_FDR",
    "size", "leading_edge_n", "leading_edge", "direction", "source_version", "source_species", "log2err"
  ))
  combined[, abs_NES_sort := abs(NES)]
  setorder(combined, registered_FDR, -abs_NES_sort, collection, pathway_id)
  combined[, abs_NES_sort := NULL]
  output <- file.path(gsea_dir, paste0(analysis, "_", lineage, "_REGISTERED_GSEA.csv"))
  atomic_write_csv(combined, output)
  log_line(sprintf("REGISTERED_GSEA %s %s pathways=%d FDR05=%d", analysis, lineage, nrow(combined), sum(combined$registered_FDR < 0.05)))
  combined
}

all_results <- list()
for (analysis in analyses) {
  for (lineage in lineages) {
    key <- paste(analysis, lineage, sep = "::")
    all_results[[key]] <- run_registered(analysis, lineage)
  }
}
combined_results <- rbindlist(all_results, fill = TRUE)
gsea_summary <- combined_results[, .(
  registered_pathways = .N,
  registered_FDR05 = sum(registered_FDR < 0.05),
  positive_NES = sum(NES > 0),
  negative_NES = sum(NES < 0)
), by = .(analysis, lineage)]
atomic_write_csv(gsea_summary, gsea_summary_output)

lookup_result <- function(analysis_value, lineage_value, collection_value, pathway_value) {
  match <- combined_results[
    analysis == analysis_value & lineage == lineage_value & collection == collection_value & pathway_id == pathway_value
  ]
  if (nrow(match) != 1L) stop("Registered GSEA lookup failed: ", analysis_value, " ", lineage_value, " ", pathway_value)
  match[1]
}

core_rows <- list()
for (i in seq_len(nrow(core_registry))) {
  frozen <- core_registry[i]
  ac_disease <- lookup_result("AC_ActiveCrohn_vs_Healthy", frozen$lineage, frozen$collection, frozen$pathway_id)
  healthy_region <- lookup_result("Healthy_AC_vs_TI", frozen$lineage, frozen$collection, frozen$pathway_id)
  active_region <- lookup_result("ActiveCrohn_AC_vs_TI", frozen$lineage, frozen$collection, frozen$pathway_id)
  same <- sign(frozen$phase7_TI_NES) == sign(ac_disease$NES)
  status <- if (same && ac_disease$registered_FDR < 0.05) {
    "ROBUST_FDR"
  } else if (same) {
    "ROBUST_DIRECTIONAL_LIMITED"
  } else if (ac_disease$registered_FDR < 0.05) {
    "REGION_SENSITIVE_FDR"
  } else {
    "REGION_UNKNOWN_OPPOSITE_WEAK"
  }
  core_rows[[i]] <- data.table(
    lineage = frozen$lineage,
    program = frozen$program,
    collection = frozen$collection,
    pathway = frozen$pathway_id,
    TI_NES = frozen$phase7_TI_NES,
    TI_FDR = frozen$phase7_TI_FDR,
    TI_direction = frozen$phase7_TI_direction,
    AC_NES = ac_disease$NES,
    AC_registered_FDR = ac_disease$registered_FDR,
    AC_direction = ac_disease$direction,
    same_direction = ifelse(same, "YES", "NO"),
    opposite_direction = ifelse(same, "NO", "YES"),
    healthy_AC_vs_TI_NES = healthy_region$NES,
    healthy_AC_vs_TI_registered_FDR = healthy_region$registered_FDR,
    activeCrohn_AC_vs_TI_NES = active_region$NES,
    activeCrohn_AC_vs_TI_registered_FDR = active_region$registered_FDR,
    region_status = status,
    notes = "AC FDR is BH-adjusted within the preregistered pathway family; AC disease effect is secondary and cannot replace Phase7 TI Primary"
  )
}
core_stability <- rbindlist(core_rows)
atomic_write_csv(core_stability, core_output)

matrix_rows <- list()
for (i in seq_len(nrow(dog_registry))) {
  frozen <- dog_registry[i]
  ac_disease <- lookup_result("AC_ActiveCrohn_vs_Healthy", frozen$lineage, frozen$representative_collection, frozen$representative_pathway_id)
  healthy_region <- lookup_result("Healthy_AC_vs_TI", frozen$lineage, frozen$representative_collection, frozen$representative_pathway_id)
  active_region <- lookup_result("ActiveCrohn_AC_vs_TI", frozen$lineage, frozen$representative_collection, frozen$representative_pathway_id)
  same <- sign(frozen$phase7_TI_NES) == sign(ac_disease$NES)
  risk <- if (same && ac_disease$registered_FDR < 0.05) {
    "LOW_REGION_RISK"
  } else if (same) {
    "MODERATE_DIRECTION_ONLY"
  } else if (ac_disease$registered_FDR < 0.05) {
    "HIGH_REGION_RISK"
  } else {
    "REGION_UNKNOWN"
  }
  priority <- switch(risk,
    LOW_REGION_RISK = "PROCEED_REGION_ROBUST",
    MODERATE_DIRECTION_ONLY = "PROCEED_WITH_REGION_CAUTION",
    HIGH_REGION_RISK = "HIGH_REGION_CAUTION",
    REGION_UNKNOWN = "REQUIRES_INDEPENDENT_REGION_CONTROL"
  )
  matrix_rows[[i]] <- data.table(
    lineage = frozen$lineage,
    program = frozen$program,
    Dog_CIE_direction = frozen$Dog_CIE_direction,
    FRE_direction = frozen$FRE_direction,
    representative_collection = frozen$representative_collection,
    representative_pathway = frozen$representative_pathway_id,
    Human_TI_direction = frozen$phase7_TI_direction,
    Human_TI_NES = frozen$phase7_TI_NES,
    Human_TI_FDR = frozen$phase7_TI_FDR,
    Human_AC_direction = ac_disease$direction,
    Human_AC_NES = ac_disease$NES,
    Human_AC_registered_FDR = ac_disease$registered_FDR,
    TI_AC_consistency = ifelse(same, "SAME_DIRECTION", "OPPOSITE_DIRECTION"),
    healthy_AC_vs_TI_NES = healthy_region$NES,
    healthy_AC_vs_TI_registered_FDR = healthy_region$registered_FDR,
    activeCrohn_AC_vs_TI_NES = active_region$NES,
    activeCrohn_AC_vs_TI_registered_FDR = active_region$registered_FDR,
    region_risk = risk,
    Phase10_priority = priority,
    scope_limit = "Region robustness only; not a cross-species conservation claim"
  )
}
interpretation <- rbindlist(matrix_rows)
atomic_write_csv(interpretation, matrix_output)

plot_data <- melt(
  core_stability,
  id.vars = c("lineage", "program"),
  measure.vars = c("TI_NES", "AC_NES", "healthy_AC_vs_TI_NES", "activeCrohn_AC_vs_TI_NES"),
  variable.name = "contrast",
  value.name = "NES"
)
plot_data[, contrast := factor(contrast, levels = c("TI_NES", "AC_NES", "healthy_AC_vs_TI_NES", "activeCrohn_AC_vs_TI_NES"), labels = c("TI disease", "AC disease", "Healthy AC-TI", "Active Crohn AC-TI"))]
plot_data[, label := paste(lineage, program, sep = " | ")]
plot_data[, label := factor(label, levels = rev(unique(label)))]
p <- ggplot(plot_data, aes(x = contrast, y = label, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2878B5", mid = "white", high = "#C74440", midpoint = 0) +
  labs(title = "Phase8A preregistered Human region-program NES", x = NULL, y = NULL, fill = "NES") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 7), panel.grid = element_blank())
ggsave(file.path(figure_dir, "PHASE8_REGION_CORE_PROGRAM_HEATMAP.png"), p, width = 11, height = 14, dpi = 180)

checks <- data.table(
  criterion = c(
    "Nine registered GSEA analyses completed",
    "Only preregistered program families were tested",
    "Thirty-six core lineage-program rows classified",
    "Twenty-five Dog-forward representatives classified",
    "No Dog-Human conservation claim encoded"
  ),
  status = c(
    ifelse(nrow(gsea_summary) == 9L, "PASS", "FAIL"),
    ifelse(all(combined_results$pathway_id %in% unique(c(core_registry$pathway_id, dog_registry$representative_pathway_id))), "PASS", "FAIL"),
    ifelse(nrow(core_stability) == 36L && !anyNA(core_stability$region_status), "PASS", "FAIL"),
    ifelse(nrow(interpretation) == 25L && !anyNA(interpretation$region_risk), "PASS", "FAIL"),
    ifelse(all(interpretation$scope_limit == "Region robustness only; not a cross-species conservation claim"), "PASS", "FAIL")
  ),
  evidence = c(
    paste0("analyses=", nrow(gsea_summary)),
    paste0("unique_pathways=", uniqueN(combined_results$pathway_id)),
    paste0("core_rows=", nrow(core_stability)),
    paste0("Dog_forward_rows=", nrow(interpretation)),
    "Phase10 claims explicitly prohibited"
  )
)
checks[, gate := ifelse(all(status == "PASS"), "GREEN_REGION_PROGRAM_ASSESSMENT_COMPLETE", "YELLOW_REGION_PROGRAM_LIMITED")]
atomic_write_csv(checks, gate_path)
if (!all(checks$status == "PASS")) stop("Phase8 region program gate failed")
log_line(paste("REGION_PROGRAM_GATE", checks$gate[1], "core=36 Dog_forward=25"))
