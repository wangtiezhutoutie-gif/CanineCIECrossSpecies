project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
lib <- file.path(project_root, "11_singIST", "environment", "r_library")
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(singIST)
  library(data.table)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

set.seed(20260826)
singist_root <- file.path(project_root, "11_singIST")
dog_dir <- file.path(singist_root, "dog_model")
human_dir <- file.path(singist_root, "human_reference")
model_dir <- file.path(human_dir, "models")
primary_dir <- file.path(singist_root, "primary", "method_faithful")
sens_dir <- file.path(singist_root, "donor_aware_sensitivity")
checkpoint_dir <- file.path(singist_root, "primary", "checkpoints")
dir.create(primary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sens_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
safe_pathway_stem <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
log_file <- file.path(project_root, "logs", "phase11_primary.log")
zz <- file(log_file, open = "at")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({sink(type = "message"); sink(type = "output"); close(zz)}, add = TRUE)

cat("\nPHASE11 PRIMARY TRACKS\n")
cat("started=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")

fit_status <- fread(file.path(human_dir, "HUMAN_REFERENCE_MODEL_FIT_STATUS.csv"))
if (nrow(fit_status) != 46L || any(fit_status$status != "SUCCESS")) {
  stop("All 46 frozen Human reference models must succeed before final primary tracks")
}
setorder(fit_status, pathway_id)
panel <- fread(file.path(singist_root, "superpathways",
                         "PHASE11_UNIFIED_SUPERPATHWAY_PANEL.csv"))
defs <- readRDS(file.path(singist_root, "superpathways",
                          "PHASE11_SUPERPATHWAY_DEFINITIONS.rds"))
mapping <- readRDS(file.path(dog_dir, "DOG_SINGIST_MAPPING_ORGANISM.rds"))
method_fc <- readRDS(file.path(dog_dir, "DOG_METHOD_FAITHFUL_FC_LIST.rds"))
donor_fc <- readRDS(file.path(dog_dir, "DOG_DONOR_AWARE_FC_LIST.rds"))
lineages <- c("Epithelial", "T", "Myeloid")

# Official Ensembl 104 Phase10 mapping. The stable-ID table has 14,945 Primary
# one-to-one rows. Thirty-five rows lack a usable Human symbol and cannot be
# represented in the symbol-based singIST model; no remapping or update occurs.
master <- fread(file.path(project_root, "10_cross_species_baseline", "mapping",
                          "DOG_HUMAN_ORTHOLOG_MASTER.csv"))
primary_master <- master[tolower(as.character(primary_one2one)) == "true"]
stopifnot(nrow(primary_master) == 14945L)
cache <- primary_master[
  !is.na(human_ortholog_symbol) & nzchar(human_ortholog_symbol) &
    human_ortholog_symbol != "NA_NOT_AVAILABLE" &
    !is.na(dog_unique_var_name) & nzchar(dog_unique_var_name),
  .(input_gene = human_ortholog_symbol, output_gene = dog_unique_var_name,
    input_ensembl = human_ortholog_release104_id,
    output_ensembl = dog_source_gene_id,
    source_release = dog_source_ensembl_release,
    homology_type = homology_type_release104)
]
if (anyDuplicated(cache$input_gene) || anyDuplicated(cache$output_gene)) {
  stop("Usable Ensembl104 symbol cache is not one-to-one")
}
stopifnot(nrow(cache) == 14910L)
fwrite(cache, file.path(singist_root, "method_gate",
                        "FROZEN_ENSEMBL104_ORTHOLOGY_CACHE.csv"))
fwrite(data.table(
  mapping_level = c("Phase10 stable-ID Primary1:1", "singIST symbol-usable cache",
                    "excluded missing Human symbol"),
  pair_n = c(nrow(primary_master), nrow(cache), nrow(primary_master) - nrow(cache)),
  release = 104L, assembly = "CanFam3.1",
  status = c("FROZEN_INPUT", "OFFICIAL_CACHE_API_INPUT",
             "EXCLUDED_IDENTIFIER_NOT_REPRESENTABLE"),
  notes = c(
    "No Phase11 remapping", "Unique Human symbols and Dog frozen feature names",
    "Human symbol is Phase10 placeholder NA_NOT_AVAILABLE"
  )
), file.path(singist_root, "method_gate", "FROZEN_ORTHOLOGY_CACHE_AUDIT.csv"))

dog_counts <- assay(mapping$counts, "counts")
dog_detected <- setNames(lapply(lineages, function(ct) {
  idx <- mapping$counts$celltype_cluster == ct
  rownames(mapping$counts)[Matrix::rowSums(dog_counts[, idx, drop = FALSE]) > 0]
}), lineages)

reference_info <- function(model) {
  c_ref <- derive_contributions(model, model$model_fit$predictor_block)
  base <- which(model$superpathway_input$sample_class ==
                  model$superpathway_input$base_class)
  target <- which(model$superpathway_input$sample_class ==
                    model$superpathway_input$target_class)
  omega <- median(c_ref$superpathway_score[target]) -
    median(c_ref$superpathway_score[base])
  gamma <- vapply(seq_len(nrow(c_ref$celltype_contribution)), function(i) {
    median(c_ref$celltype_contribution[i, target]) -
      median(c_ref$celltype_contribution[i, base])
  }, numeric(1))
  names(gamma) <- rownames(c_ref$celltype_contribution)
  list(omega = omega, gamma = gamma)
}

run_one <- function(pathway_id, model, fc, track, ref_info) {
  out <- singISTrecapitulations(
    mapping, model, orthology_cache = cache, FC_list = fc
  )
  out$superpathway$pathway <- pathway_id
  out$celltype$pathway <- pathway_id
  out$gene$pathway <- pathway_id
  out$superpathway$track <- track
  out$celltype$track <- track
  out$gene$track <- track
  out$superpathway$reference_shift <- ref_info$omega
  out$celltype$reference_shift <- ref_info$gamma[
    match(out$celltype$celltype, names(ref_info$gamma))]
  out$superpathway$reference_shift_status <- ifelse(
    !is.finite(ref_info$omega) || abs(ref_info$omega) < .Machine$double.eps,
    "NOT_INTERPRETABLE_ZERO_OR_NONFINITE",
    ifelse(abs(ref_info$omega) < 1e-8,
           "REFERENCE_SHIFT_NUMERICALLY_UNSTABLE", "INTERPRETABLE")
  )
  out$celltype$reference_shift_status <- ifelse(
    !is.finite(out$celltype$reference_shift) |
      abs(out$celltype$reference_shift) < .Machine$double.eps,
    "NOT_INTERPRETABLE_ZERO_OR_NONFINITE",
    ifelse(abs(out$celltype$reference_shift) < 1e-8,
           "REFERENCE_SHIFT_NUMERICALLY_UNSTABLE", "INTERPRETABLE")
  )
  out
}

for (i in seq_len(nrow(fit_status))) {
  fs <- fit_status[i]
  pathway_id <- fs$pathway_id
  cp <- file.path(checkpoint_dir,
                  paste0(safe_pathway_stem(pathway_id), "_PRIMARY_TRACKS.rds"))
  if (file.exists(cp)) {
    cat("SKIP_CHECKPOINT ", i, "/", nrow(fit_status), " ", pathway_id, "\n", sep = "")
    next
  }
  cat("PRIMARY_START ", i, "/", nrow(fit_status), " ", pathway_id, "\n", sep = "")
  model <- readRDS(fs$model_path)
  ref_info <- reference_info(model)
  a <- run_one(pathway_id, model, method_fc,
               "METHOD_FAITHFUL_SINGIST", ref_info)
  b <- run_one(pathway_id, model, donor_fc,
               "DONOR_AWARE_FC_SENSITIVITY", ref_info)
  saveRDS(list(method_faithful = a, donor_aware = b, reference = ref_info),
          cp, compress = "gzip")
  cat("PRIMARY_SUCCESS ", pathway_id,
      " A=", a$superpathway$recapitulation,
      " B=", b$superpathway$recapitulation, "\n", sep = "")
  flush.console()
}

checkpoints <- file.path(checkpoint_dir,
                         paste0(safe_pathway_stem(fit_status$pathway_id),
                                "_PRIMARY_TRACKS.rds"))
stopifnot(all(file.exists(checkpoints)))
res <- lapply(checkpoints, readRDS)

super_a <- rbindlist(lapply(res, function(x) as.data.table(x$method_faithful$superpathway)),
                     fill = TRUE)
cell_a <- rbindlist(lapply(res, function(x) as.data.table(x$method_faithful$celltype)),
                    fill = TRUE)
gene_a <- rbindlist(lapply(res, function(x) as.data.table(x$method_faithful$gene)),
                    fill = TRUE)
super_b <- rbindlist(lapply(res, function(x) as.data.table(x$donor_aware$superpathway)),
                     fill = TRUE)
cell_b <- rbindlist(lapply(res, function(x) as.data.table(x$donor_aware$celltype)),
                    fill = TRUE)
gene_b <- rbindlist(lapply(res, function(x) as.data.table(x$donor_aware$gene)),
                    fill = TRUE)

annot <- panel[, .(pathway = pathway_id, collection, pathway_name,
                   panel_A, panel_B, panel_A_significant_lineages,
                   registered_programs)]
valid <- fit_status[, .(
  pathway = pathway_id,
  human_reference_PValue = pvalue_global_significance,
  human_reference_FDR = FDR_global_significance,
  reference_model_status,
  selected_PLS, selected_quantiles
)]
add_annotations <- function(x) {
  x <- merge(x, annot, by = "pathway", all.x = TRUE, sort = FALSE)
  merge(x, valid, by = "pathway", all.x = TRUE, sort = FALSE)
}
super_a <- add_annotations(super_a)
cell_a <- add_annotations(cell_a)
gene_a <- add_annotations(gene_a)
super_b <- add_annotations(super_b)
cell_b <- add_annotations(cell_b)
gene_b <- add_annotations(gene_b)

# Add Human→Dog frozen gene identity to contribution tables.
gene_a[, dog_ortholog := cache$output_gene[match(gene, cache$input_gene)]]
gene_b[, dog_ortholog := cache$output_gene[match(gene, cache$input_gene)]]
gene_a[, ortholog_mapped := !is.na(dog_ortholog)]
gene_b[, ortholog_mapped := !is.na(dog_ortholog)]

# Official gene-contribution additivity audit uses unrounded cell-type values.
gene_sum <- gene_a[, .(gene_contribution_sum = sum(contribution)),
                   by = .(pathway, celltype)]
add_audit <- merge(
  gene_sum, cell_a[, .(pathway, celltype,
                       celltype_recapitulation = recapitulation)],
  by = c("pathway", "celltype"), all = TRUE
)
add_audit[, absolute_difference := abs(gene_contribution_sum -
                                         celltype_recapitulation)]
add_audit[, tolerance := 1e-8]
add_audit[, status := fifelse(is.finite(absolute_difference) &
                                absolute_difference <= tolerance,
                              "PASS", "ERROR")]
fwrite(add_audit, file.path(primary_dir,
                            "GENE_CONTRIBUTION_ADDITIVITY_AUDIT.csv"))
if (any(add_audit$status != "PASS", na.rm = TRUE))
  stop("Gene contribution additivity failure")

# Explicit per-pathway, per-block orthology and expression coverage.
coverage_rows <- list()
for (pathway_id in fit_status$pathway_id) {
  pid <- pathway_id
  model <- readRDS(fit_status[pathway_id == pid]$model_path)
  d <- defs[[pathway_id]]
  for (b in seq_along(lineages)) {
    ct <- lineages[[b]]
    human_obs <- model$model_fit$observed_gene_sets[[b]]
    hits <- cache[input_gene %in% human_obs]
    present <- unique(hits$output_gene[hits$output_gene %in%
                                         rownames(mapping$counts)])
    detected <- unique(hits$output_gene[hits$output_gene %in% dog_detected[[ct]]])
    detected_human <- unique(hits$input_gene[hits$output_gene %in% detected])
    coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
      pathway = pathway_id, collection = d$collection,
      pathway_name = d$pathway_name, celltype = ct,
      pathway_gene_n = length(d$candidate_genes),
      human_observed_n = length(human_obs),
      dog_observed_n = length(detected),
      one2one_observed_n = length(detected_human),
      orthology_coverage = 100 * length(detected_human) / length(human_obs),
      one2one_object_presence_n = length(present),
      object_presence_coverage = 100 * length(present) / length(human_obs),
      mapping_source = "Phase10 Ensembl104 CanFam3.1 ortholog_one2one",
      mapping_pair_n_stable_id = 14945L,
      mapping_pair_n_symbol_usable = 14910L
    )
  }
}
coverage <- rbindlist(coverage_rows)

# Correct a data.table non-standard-evaluation ambiguity above if encountered by
# rebuilding coverage from the deterministic model-path map.
if (nrow(coverage) != 46L * 3L) stop("Orthology coverage row count mismatch")

fwrite(super_a, file.path(primary_dir, "SUPERPATHWAY_RECAPITULATION.csv"))
fwrite(cell_a, file.path(primary_dir, "CELLTYPE_RECAPITULATION.csv"))
fwrite(gene_a, file.path(primary_dir, "GENE_CONTRIBUTIONS.csv.gz"),
       compress = "gzip")
fwrite(coverage, file.path(primary_dir,
                           "ORTHOLOGY_COVERAGE_BY_SUPERPATHWAY.csv"))

super_b[, recap_level := "SUPERPATHWAY"]
super_b[, celltype := NA_character_]
cell_b[, recap_level := "CELLTYPE"]
donor_recap <- rbindlist(list(super_b, cell_b), fill = TRUE)
fwrite(donor_recap, file.path(sens_dir,
                              "SINGIST_DONOR_AWARE_RECAPITULATION.csv"))
fwrite(gene_b, file.path(sens_dir,
                         "SINGIST_DONOR_AWARE_GENE_CONTRIBUTIONS.csv.gz"),
       compress = "gzip")

sa <- super_a[, .(pathway, collection, pathway_name,
                  method_faithful_recap = recapitulation)]
sb <- super_b[, .(pathway, donor_aware_recap = recapitulation)]
sens_super <- merge(sa, sb, by = "pathway")
sens_super[, `:=`(
  recap_level = "SUPERPATHWAY", celltype = NA_character_,
  same_sign = sign(method_faithful_recap) == sign(donor_aware_recap),
  delta_recap = donor_aware_recap - method_faithful_recap
)]
ca <- cell_a[, .(pathway, collection, pathway_name, celltype,
                 method_faithful_recap = recapitulation,
                 method_faithful_orthology_coverage = orthology)]
cb <- cell_b[, .(pathway, celltype, donor_aware_recap = recapitulation,
                 donor_aware_orthology_coverage = orthology)]
sens_cell <- merge(ca, cb, by = c("pathway", "celltype"))
sens_cell[, `:=`(
  recap_level = "CELLTYPE",
  same_sign = sign(method_faithful_recap) == sign(donor_aware_recap),
  delta_recap = donor_aware_recap - method_faithful_recap
)]
sens <- rbindlist(list(sens_super, sens_cell), fill = TRUE)
sens[, notes := paste0(
  "Track A uses precomputed official diff_expressed; Track B replaces only ",
  "p_val/p_val_adj with Phase3 Full7 edgeR PValue/FDR"
)]
fwrite(sens, file.path(sens_dir, "SINGIST_FC_SOURCE_SENSITIVITY.csv"))

# Update numerical-stability flags in the formal validity table without dropping
# weak models or unstable raw values.
stability <- super_a[, .(
  pathway_id = pathway, reference_shift, reference_shift_status
)]
validity <- merge(fit_status, stability, by = "pathway_id", all.x = TRUE)
validity[, numeric_stability := reference_shift_status]
fwrite(validity, file.path(primary_dir, "HUMAN_REFERENCE_MODEL_VALIDITY.csv"))

cat("superpathway_rows=", nrow(super_a), "\n", sep = "")
cat("celltype_rows=", nrow(cell_a), "\n", sep = "")
cat("gene_rows=", nrow(gene_a), "\n", sep = "")
cat("max_additivity_error=", max(add_audit$absolute_difference), "\n", sep = "")
cat("coverage_range=", paste(range(coverage$orthology_coverage), collapse = ";"), "\n", sep = "")
cat("completed=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("PHASE11_PRIMARY_TRACKS_COMPLETE\n")
