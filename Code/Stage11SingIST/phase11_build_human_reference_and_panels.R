project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
lib <- file.path(project_root, "11_singIST", "environment", "r_library")
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(scuttle)
  library(SummarizedExperiment)
})

set.seed(20260826)
human_dir <- file.path(project_root, "11_singIST", "human_reference")
panel_dir <- file.path(project_root, "11_singIST", "superpathways")
metadata_dir <- file.path(project_root, "00_metadata")
dir.create(human_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(project_root, "logs", "phase11_human_reference.log")
zz <- file(log_file, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({
  sink(type = "message"); sink(type = "output"); close(zz)
}, add = TRUE)

cat("PHASE11 HUMAN REFERENCE AND PANEL FREEZE\n")
cat("started=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")

lineages <- c("Epithelial", "T", "Myeloid")
pb_root <- file.path(project_root, "07_human_crohn", "GSE266616", "pseudobulk")
phase7_meta <- fread(file.path(pb_root, "metadata", "HUMAN_PSEUDOBULK_SAMPLE_METADATA.csv"))
phase7_meta <- phase7_meta[
  analysis == "Primary_ActiveCrohn_vs_Healthy" & lineage %in% lineages
]
stopifnot(nrow(phase7_meta) == 81L)
stopifnot(all(table(phase7_meta$lineage) == 27L))
stopifnot(all(table(phase7_meta$lineage, phase7_meta$condition)[, "Healthy"] == 16L))
stopifnot(all(table(phase7_meta$lineage, phase7_meta$condition)[, "ActiveCrohn"] == 11L))

fread_gz_frozen <- function(path, sep = "auto") {
  tmp <- tempfile(fileext = if (sep == "\t") ".tsv" else ".csv")
  input <- gzfile(path, open = "rb")
  output <- file(tmp, open = "wb")
  repeat {
    buf <- readBin(input, what = "raw", n = 1024L * 1024L)
    if (!length(buf)) break
    writeBin(buf, output)
  }
  close(input); close(output)
  on.exit(unlink(tmp), add = TRUE)
  if (sep == "auto") fread(tmp) else fread(tmp, sep = sep)
}

read_symbol_counts <- function(lineage) {
  path <- file.path(pb_root, "matrices", paste0("human_", lineage, "_counts.csv.gz"))
  tab <- fread_gz_frozen(path)
  required <- c("unique_gene_id", "gene_symbol")
  stopifnot(all(required %in% names(tab)))
  sample_cols <- setdiff(names(tab), required)
  mat <- as.matrix(tab[, ..sample_cols])
  storage.mode(mat) <- "double"
  symbols <- trimws(as.character(tab$gene_symbol))
  keep <- !is.na(symbols) & nzchar(symbols) & symbols != "NA"
  invalid_library <- colSums(mat[!keep, , drop = FALSE])
  if (any(invalid_library != 0)) {
    stop("Non-zero counts have missing Human gene symbols in ", lineage)
  }
  # HGNC symbols must be unique columns in the singIST training matrix. Summing
  # duplicate-symbol features preserves every raw count and the library size.
  collapsed <- rowsum(mat[keep, , drop = FALSE], symbols[keep], reorder = FALSE)
  list(
    counts = collapsed,
    donor = sample_cols,
    original_feature_n = nrow(tab),
    unique_symbol_n = nrow(collapsed),
    duplicate_feature_n = sum(duplicated(symbols[keep])),
    source = path
  )
}

blocks <- setNames(lapply(lineages, read_symbol_counts), lineages)
donors <- blocks[[1]]$donor
stopifnot(all(vapply(blocks, function(x) identical(x$donor, donors), logical(1))))
stopifnot(length(donors) == 27L)
genes <- unique(unlist(lapply(blocks, function(x) rownames(x$counts)), use.names = FALSE))

align_block <- function(x, genes) {
  out <- matrix(0, nrow = length(genes), ncol = ncol(x$counts),
                dimnames = list(genes, x$donor))
  out[rownames(x$counts), ] <- x$counts
  out
}
aligned <- lapply(blocks, align_block, genes = genes)
combined_counts <- do.call(cbind, Map(function(mat, lineage) {
  colnames(mat) <- paste(lineage, colnames(mat), sep = "_")
  mat
}, aligned, names(aligned)))

profile_meta <- rbindlist(lapply(lineages, function(lineage) {
  lineage_value <- lineage
  m <- copy(phase7_meta[lineage == lineage_value])
  m <- m[match(donors, patient_id)]
  stopifnot(identical(m$patient_id, donors))
  data.table(
    profile_id = paste(lineage, m$patient_id, sep = "_"),
    donor = m$patient_id,
    sample_id = m$sample_id,
    condition = m$condition,
    celltype = lineage,
    cl_id = c(Epithelial = "CL:0000066", T = "CL:0000084",
              Myeloid = "CL:0000763")[[lineage]],
    n_cells = as.integer(m$n_cells),
    phase7_raw_library_size = as.numeric(m$library_size)
  )
}))
stopifnot(identical(colnames(combined_counts), profile_meta$profile_id))
raw_lib <- colSums(combined_counts)
if (!isTRUE(all.equal(as.numeric(raw_lib), profile_meta$phase7_raw_library_size,
                      tolerance = 0))) {
  cmp <- data.table(profile_id = profile_meta$profile_id,
                    computed = as.numeric(raw_lib),
                    phase7 = profile_meta$phase7_raw_library_size,
                    delta = as.numeric(raw_lib) - profile_meta$phase7_raw_library_size)
  print(cmp[order(-abs(delta))][1:10])
  stop("Phase7 raw-count library size mismatch after symbol collapse")
}

sce <- SingleCellExperiment(
  assays = list(counts = combined_counts),
  colData = S4Vectors::DataFrame(profile_meta)
)
sce <- scuttle::computeLibraryFactors(sce)
sce <- scuttle::logNormCounts(
  sce, assay.type = "counts", size.factors = sizeFactors(sce),
  pseudo_count = 1, center_size_factors = FALSE
)
lognorm <- t(as.matrix(assay(sce, "logcounts")))
rownames(lognorm) <- profile_meta$profile_id
stopifnot(all(is.finite(lognorm)))
stopifnot(identical(rownames(lognorm), profile_meta$profile_id))

profile_meta[, raw_library_size := as.numeric(raw_lib)]
profile_meta[, size_factor := as.numeric(sizeFactors(sce))]
profile_meta[, normalization := paste0(
  "singIST1.0.2 SCE route: raw donor-celltype sum; ",
  "scuttle::computeLibraryFactors; scuttle::logNormCounts(",
  "pseudo_count=1,center_size_factors=FALSE)"
)]
profile_meta[, gene_count := ncol(lognorm)]
profile_meta[, detected_gene_count := colSums(combined_counts > 0)]
profile_meta[, base_class := "Healthy"]
profile_meta[, target_class := "ActiveCrohn"]
profile_meta[, tissue := "terminal_ileum"]

pb_path <- file.path(human_dir, "HUMAN_SINGIST_PSEUDOBULK_LOGNORM.rds")
saveRDS(lognorm, pb_path, compress = "gzip")
fwrite(profile_meta,
       file.path(human_dir, "HUMAN_SINGIST_PSEUDOBULK_METADATA.csv"))
roundtrip <- readRDS(pb_path)
if (!isTRUE(all.equal(lognorm, roundtrip, tolerance = 0))) {
  stop("Human pseudobulk RDS round-trip failed")
}

# Freeze the exact Human-driven Hallmark panel without reading Dog singIST data.
gsea_dir <- file.path(project_root, "07_human_crohn", "GSE266616", "gsea")
hallmark_results <- rbindlist(lapply(lineages, function(lineage) {
  fread(file.path(gsea_dir, "full", paste0("Human_", lineage,
                                           "_Hallmark_GSEA.csv")))
}), fill = TRUE)
hallmark_sig <- hallmark_results[is.finite(padj) & padj < 0.05]
panel_a <- hallmark_sig[, .(
  pathway_name = pathway_name[1],
  significant_lineages = paste(lineages[lineages %in% unique(lineage)],
                               collapse = ";"),
  significant_lineage_n = uniqueN(lineage),
  minimum_human_FDR = min(padj),
  max_abs_human_NES = max(abs(NES)),
  source_version = source_version[1],
  source_species = source_species[1]
), by = .(collection, pathway_id)]
setorder(panel_a, pathway_id)
panel_a[, panel_A := TRUE]

# Freeze all 25 registered Phase10 programs and their original representatives.
phase10 <- fread(file.path(project_root, "10_cross_species_baseline", "summary",
                           "PHASE10_PROGRAM_BASELINE_MATRIX.csv"))
stopifnot(nrow(phase10) == 25L)
panel_b <- phase10[, .(
  lineage, program,
  collection = representative_collection,
  pathway_id = representative_pathway,
  phase10_baseline_status = baseline_status,
  region_risk, phase10_priority,
  dog_lodo_stability,
  human_treatment_status = human_treatment_stability,
  fre_status = fre_summary_status
)]
panel_b[, panel_B := TRUE]

membership_files <- c(
  Hallmark = file.path(gsea_dir, "gene_sets", "HUMAN_HALLMARK_GENE_SETS.tsv.gz"),
  Reactome = file.path(gsea_dir, "gene_sets", "HUMAN_REACTOME_GENE_SETS.tsv.gz"),
  GO_BP = file.path(gsea_dir, "gene_sets", "HUMAN_GO_BP_GENE_SETS.tsv.gz")
)
members <- rbindlist(lapply(membership_files, fread_gz_frozen, sep = "\t"),
                     fill = TRUE)
members <- unique(members, by = c("collection", "pathway_id", "gene_symbol"))

panel_b <- merge(
  panel_b,
  members[, .(pathway_name = pathway_name[1], candidate_gene_n = uniqueN(gene_symbol),
              source_version = source_version[1], source_species = source_species[1]),
          by = .(collection, pathway_id)],
  by = c("collection", "pathway_id"), all.x = TRUE, sort = FALSE
)
panel_b[, testability := fifelse(
  is.na(candidate_gene_n), "NOT_TESTABLE_FROZEN_GENESET_NOT_FOUND", "PENDING_OBSERVED_GENE_AUDIT"
)]
panel_b[, lineage_order__ := match(lineage, lineages)]
setorder(panel_b, lineage_order__, program)
panel_b[, lineage_order__ := NULL]

# One model per unique pathway. Duplicates across Panel A/B or registered
# lineages are fitted only once, while panel membership is retained.
all_keys <- unique(rbindlist(list(
  panel_a[, .(collection, pathway_id)],
  panel_b[, .(collection, pathway_id)]
)))
unified <- merge(
  all_keys,
  members[, .(pathway_name = pathway_name[1], source_version = source_version[1],
              source_species = source_species[1], candidate_gene_n = uniqueN(gene_symbol)),
          by = .(collection, pathway_id)],
  by = c("collection", "pathway_id"), all.x = TRUE
)
unified[, panel_A := pathway_id %in% panel_a$pathway_id]
unified[, panel_B := pathway_id %in% panel_b$pathway_id]
unified[, panel_A_significant_lineages := panel_a$significant_lineages[
  match(pathway_id, panel_a$pathway_id)]]
unified[, registered_programs := vapply(pathway_id, function(p) {
  z <- panel_b[pathway_id == p]
  if (!nrow(z)) return(NA_character_)
  paste(paste(z$lineage, z$program, sep = ": "), collapse = "; ")
}, character(1))]

# Candidate-to-observed feature audit is Human-only and fixed before any Dog
# recapitulation. A feature is observed in a block if raw pseudobulk count is >0
# in at least one of its 27 donors and log-normalized variance is positive.
definition <- list()
obs_rows <- list()
for (i in seq_len(nrow(unified))) {
  key <- unified[i]
  candidate <- members[
    collection == key$collection & pathway_id == key$pathway_id,
    unique(gene_symbol)
  ]
  gs_blocks <- setNames(vector("list", length(lineages)), lineages)
  for (celltype in lineages) {
    idx <- which(as.character(profile_meta$celltype) == celltype)
    observed <- intersect(candidate, colnames(lognorm))
    observed <- observed[colSums(combined_counts[observed,
                                                  profile_meta$celltype == celltype,
                                                  drop = FALSE]) > 0]
    if (length(observed)) {
      gene_idx <- match(observed, colnames(lognorm))
      keep_idx <- !is.na(gene_idx)
      observed <- observed[keep_idx]
      gene_idx <- gene_idx[keep_idx]
      variances <- apply(lognorm[idx, gene_idx, drop = FALSE], 2, stats::var)
      observed <- observed[is.finite(variances) & variances > 0]
    }
    gs_blocks[[celltype]] <- observed
    obs_rows[[length(obs_rows) + 1L]] <- data.table(
      collection = key$collection, pathway_id = key$pathway_id,
      pathway_name = key$pathway_name, celltype = celltype,
      candidate_gene_n = length(candidate), human_observed_n = length(observed)
    )
  }
  definition[[key$pathway_id]] <- list(
    collection = key$collection,
    pathway_id = key$pathway_id,
    pathway_name = key$pathway_name,
    source_version = key$source_version,
    candidate_genes = candidate,
    gene_sets_celltype = gs_blocks
  )
}
obs_audit <- rbindlist(obs_rows)
testability <- obs_audit[, .(
  minimum_human_observed_n = min(human_observed_n),
  observed_all_three_blocks = all(human_observed_n >= 2L)
), by = .(collection, pathway_id)]
unified <- merge(unified, testability, by = c("collection", "pathway_id"), all.x = TRUE)
unified[, testability := fifelse(
  is.na(candidate_gene_n), "NOT_TESTABLE_FROZEN_GENESET_NOT_FOUND",
  fifelse(observed_all_three_blocks, "TESTABLE", "NOT_TESTABLE_LT2_OBSERVED_GENES")
)]
panel_b[, `:=`(
  minimum_human_observed_n = unified$minimum_human_observed_n[
    match(pathway_id, unified$pathway_id)],
  testability = unified$testability[match(pathway_id, unified$pathway_id)]
)]

fwrite(panel_a, file.path(panel_dir, "HUMAN_DRIVEN_HALLMARK_PANEL.csv"))
fwrite(panel_b, file.path(panel_dir, "REGISTERED_PROGRAM_PANEL.csv"))
fwrite(unified, file.path(panel_dir, "PHASE11_UNIFIED_SUPERPATHWAY_PANEL.csv"))
fwrite(obs_audit, file.path(panel_dir, "HUMAN_OBSERVED_FEATURE_AUDIT.csv"))
saveRDS(definition, file.path(panel_dir, "PHASE11_SUPERPATHWAY_DEFINITIONS.rds"),
        compress = "gzip")

# Hyperparameters are frozen before the first formal project fit.
quantile_values <- seq(0.1, 0.95, by = 0.3)
hyper <- data.table(
  parameter = c(
    "quantile_values_per_block", "quantile_combination_count", "number_PLS",
    "folds_CV", "repetition_CV", "measure", "Method",
    "expected_measure_increase", "maxiter", "npermut", "nsubsampling",
    "validation_type", "global_significance_full", "CIP.GIP_significance_full",
    "random_seed", "parallel", "parallel_backend", "workers"
  ),
  value = c(
    paste(quantile_values, collapse = ";"), length(quantile_values)^3,
    3, 1, 1, "B_accuracy", "NULL_package_default", 0.005, 100, 1000, 100,
    "jackknife", FALSE, FALSE, 20260826, TRUE,
    "BiocParallel::SnowParam", 3
  ),
  source = c(
    "singIST1.0.2 vignette LOOCV example lines 255-270",
    "Cartesian product across Epithelial/T/Myeloid",
    "Phase11 protocol and singIST1.0.2 vignette LOOCV example",
    "Phase11 protocol: LOOCV", "singIST1.0.2 LOOCV requirement",
    "singIST1.0.2 fitOptimal default and Phase11 protocol",
    "singIST1.0.2 fitOptimal default", "singIST1.0.2 fitOptimal default",
    "singIST1.0.2 fitOptimal default", "Phase11 protocol", "Phase11 protocol",
    "singIST1.0.2 default; LOOCV small-sample validation",
    "Phase11 protocol", "Phase11 protocol", "Phase11 protocol",
    "singIST1.0.2 vignette", "singIST1.0.2 vignette cross-platform backend",
    "Preregistered local worker count"
  )
)
fwrite(hyper, file.path(metadata_dir, "PHASE11_SINGIST_HYPERPARAMETERS.csv"))

cat("donors=", length(donors), "\n", sep = "")
cat("healthy=", uniqueN(profile_meta[condition == "Healthy", donor]), "\n", sep = "")
cat("active_crohn=", uniqueN(profile_meta[condition == "ActiveCrohn", donor]), "\n", sep = "")
cat("profiles=", nrow(lognorm), "\n", sep = "")
cat("genes=", ncol(lognorm), "\n", sep = "")
for (lineage in lineages) {
  cat(lineage, "_original_features=", blocks[[lineage]]$original_feature_n,
      " unique_symbols=", blocks[[lineage]]$unique_symbol_n,
      " duplicate_features=", blocks[[lineage]]$duplicate_feature_n, "\n", sep = "")
}
cat("human_driven_hallmark_n=", nrow(panel_a), "\n", sep = "")
cat("registered_program_n=", nrow(panel_b), "\n", sep = "")
cat("unique_superpathway_n=", nrow(unified), "\n", sep = "")
cat("testable_superpathway_n=", sum(unified$testability == "TESTABLE"), "\n", sep = "")
cat("completed=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("PHASE11_HUMAN_REFERENCE_COMPLETE\n")
