project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
lib <- file.path(project_root, "11_singIST", "environment", "r_library")
Sys.setenv(R_LIBS_USER = lib)
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(singIST)
  library(data.table)
  library(BiocParallel)
})

seed <- 20260826L
set.seed(seed)
human_dir <- file.path(project_root, "11_singIST", "human_reference")
model_dir <- file.path(human_dir, "models")
input_dir <- file.path(human_dir, "model_inputs")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(human_dir, "HUMAN_REFERENCE_MODEL_FIT_STATUS.csv")
validity_path <- file.path(project_root, "11_singIST", "primary",
                           "method_faithful", "HUMAN_REFERENCE_MODEL_VALIDITY.csv")
dir.create(dirname(validity_path), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(project_root, "logs", "phase11_model_fit.log")
zz <- file(log_file, open = "at")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({
  try(BiocParallel::register(BiocParallel::SerialParam(), default = TRUE), silent = TRUE)
  sink(type = "message"); sink(type = "output"); close(zz)
}, add = TRUE)

cat("\nPHASE11 FORMAL HUMAN MODEL FIT RUN\n")
cat("started=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")

pb <- readRDS(file.path(human_dir, "HUMAN_SINGIST_PSEUDOBULK_LOGNORM.rds"))
meta <- fread(file.path(human_dir, "HUMAN_SINGIST_PSEUDOBULK_METADATA.csv"))
panel <- fread(file.path(project_root, "11_singIST", "superpathways",
                         "PHASE11_UNIFIED_SUPERPATHWAY_PANEL.csv"))
defs <- readRDS(file.path(project_root, "11_singIST", "superpathways",
                          "PHASE11_SUPERPATHWAY_DEFINITIONS.rds"))
panel <- panel[testability == "TESTABLE"]
setorder(panel, pathway_id)
stopifnot(nrow(panel) == 46L)

lineages <- c("Epithelial", "T", "Myeloid")
donor_class <- unique(meta[, .(donor, condition)])
stopifnot(nrow(donor_class) == 27L)
setorder(donor_class, donor)
sample_id <- donor_class$donor
sample_class <- donor_class$condition
stopifnot(sum(sample_class == "Healthy") == 16L)
stopifnot(sum(sample_class == "ActiveCrohn") == 11L)
stopifnot(identical(sample_id, sort(sample_id)))

# Vignette LOOCV example grid: seq(0.1, 0.95, by=0.3) independently per block.
q <- seq(0.1, 0.95, by = 0.3)
qgrid <- as.matrix(expand.grid(Epithelial = q, T = q, Myeloid = q))
storage.mode(qgrid) <- "double"
stopifnot(nrow(qgrid) == 27L, ncol(qgrid) == 3L)
hp <- create_hyperparameters(
  quantile_comb_table = qgrid, outcome_type = "binary",
  number_PLS = as.integer(3), folds_CV = as.integer(1),
  repetition_CV = as.integer(1)
)

# Keep one reproducible cross-platform backend. Every pathway gets an explicit
# RNG seed so checkpoint/resume does not change its permutation stream.
backend <- BiocParallel::SnowParam(
  workers = 3, type = "SOCK", progressbar = FALSE,
  RNGseed = seed, exportglobals = FALSE
)
BiocParallel::register(backend, default = TRUE)

if (file.exists(status_path)) {
  status <- fread(status_path)
} else {
  status <- data.table()
}

write_status <- function(row) {
  if (nrow(status)) status <<- status[pathway_id != row$pathway_id]
  status <<- rbindlist(list(status, row), fill = TRUE)
  setorder(status, pathway_id)
  fwrite(status, status_path)
}

safe_pathway_stem <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

for (i in seq_len(nrow(panel))) {
  p <- panel[i]
  safe_stem <- safe_pathway_stem(p$pathway_id)
  model_path <- file.path(model_dir, paste0(safe_stem, ".rds"))
  input_path <- file.path(input_dir, paste0(safe_stem, "_input.rds"))
  existing <- if (nrow(status)) status[pathway_id == p$pathway_id] else data.table()
  if (file.exists(model_path) && nrow(existing) == 1L &&
      existing$status == "SUCCESS") {
    cat("SKIP_SUCCESS ", i, "/", nrow(panel), " ", p$pathway_id, "\n", sep = "")
    next
  }

  start <- Sys.time()
  cat("FIT_START ", i, "/", nrow(panel), " ", p$pathway_id,
      " at ", format(start, "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
  flush.console()
  result <- tryCatch({
    d <- defs[[p$pathway_id]]
    stopifnot(!is.null(d))
    genes <- unique(unlist(d$gene_sets_celltype, use.names = FALSE))
    genes <- intersect(genes, colnames(pb))
    pb_small <- pb[, genes, drop = FALSE]
    # matrixToBlock sorts donor row names. The response vector is therefore
    # explicitly supplied in the same alphabetical donor order.
    required_rows <- unlist(lapply(lineages, function(ct)
      paste(ct, sample_id, sep = "_")), use.names = FALSE)
    pb_small <- pb_small[required_rows, , drop = FALSE]
    # singIST 1.0.2 validates pathway metadata against an older c2/CP-only
    # namespace and cannot represent Hallmark or GO_BP identifiers directly.
    # This transparent alias is metadata-only; frozen identifiers and genes are
    # retained as extra fields and in every external result table.
    pathway_alias <- paste0(
      "REACTOME_PHASE11_ALIAS_",
      gsub("[^A-Za-z0-9]+", "_", d$pathway_id)
    )
    pathway <- create_pathway(
      standard_name = pathway_alias, dbsource = "REACTOME",
      collection = "c2", subcollection = "CP"
    )
    pathway$frozen_pathway_id <- d$pathway_id
    pathway$frozen_collection <- d$collection
    pathway$frozen_pathway_name <- d$pathway_name
    pathway$frozen_source_version <- d$source_version
    superpathway <- create_superpathway(
      pathway_info = pathway, celltypes = lineages,
      gene_sets_celltype = unname(d$gene_sets_celltype[lineages])
    )
    input <- create_superpathway_input(
      superpathway_info = superpathway, hyperparameters_info = hp,
      pseudobulk_lognorm = pb_small, sample_id = sample_id,
      sample_class = sample_class, base_class = "Healthy",
      target_class = "ActiveCrohn"
    )
    pre <- matrixToBlock(input)
    stopifnot(identical(rownames(pre$block_predictor), sample_id))
    expected_y <- as.numeric(sample_class == "ActiveCrohn")
    stopifnot(identical(as.numeric(pre$matrix_response[, 1]), expected_y))
    saveRDS(input, input_path, compress = "gzip")

    set.seed(seed + i)
    fit <- fitOptimal(
      input, parallel = TRUE, measure = "B_accuracy", Method = NULL,
      expected_measure_increase = 0.005, maxiter = 100,
      global_significance_full = FALSE,
      CIP.GIP_significance_full = FALSE,
      npermut = 1000, type = "jackknife", nsubsampling = 100
    )
    if (!identical(rownames(fit$model_fit$predictor_block), sample_id)) {
      stop("Fitted predictor row order differs from frozen donor order")
    }
    if (!identical(as.numeric(fit$model_fit$response_matrix[, 1]), expected_y)) {
      stop("Fitted response direction/order mismatch")
    }
    saveRDS(fit, model_path, compress = "gzip")
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    qfit <- fit$hyperparameters_fit$quantile_comb_table
    row <- data.table(
      pathway_id = p$pathway_id, collection = p$collection,
      pathway_name = p$pathway_name, status = "SUCCESS",
      panel_A = p$panel_A, panel_B = p$panel_B,
      observed_epithelial_n = length(fit$model_fit$observed_gene_sets[[1]]),
      observed_T_n = length(fit$model_fit$observed_gene_sets[[2]]),
      observed_myeloid_n = length(fit$model_fit$observed_gene_sets[[3]]),
      selected_PLS = fit$hyperparameters_fit$number_PLS,
      selected_quantiles = paste(as.numeric(qfit), collapse = ";"),
      pvalue_global_significance = fit$model_validation$pvalue_global_significance,
      base_class = "Healthy", target_class = "ActiveCrohn",
      sample_order_verified = TRUE, npermut = 1000L,
      random_seed = seed + i, elapsed_seconds = elapsed,
      model_path = model_path, error = NA_character_
    )
    cat("FIT_SUCCESS ", p$pathway_id, " p=",
        fit$model_validation$pvalue_global_significance,
        " elapsed_s=", round(elapsed, 1), "\n", sep = "")
    row
  }, error = function(e) {
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    cat("FIT_ERROR ", p$pathway_id, " ", conditionMessage(e), "\n", sep = "")
    data.table(
      pathway_id = p$pathway_id, collection = p$collection,
      pathway_name = p$pathway_name, status = "FAILED",
      panel_A = p$panel_A, panel_B = p$panel_B,
      observed_epithelial_n = NA_integer_, observed_T_n = NA_integer_,
      observed_myeloid_n = NA_integer_, selected_PLS = NA_integer_,
      selected_quantiles = NA_character_, pvalue_global_significance = NA_real_,
      base_class = "Healthy", target_class = "ActiveCrohn",
      sample_order_verified = FALSE, npermut = 1000L,
      random_seed = seed + i, elapsed_seconds = elapsed,
      model_path = model_path, error = conditionMessage(e)
    )
  })
  write_status(result)
  flush.console()
}

status <- fread(status_path)
# Drop first so that a stale logical column from an interrupted/older summary
# cannot force numeric adjusted p values to TRUE/FALSE by reference.
if ("FDR_global_significance" %in% names(status)) {
  status[, FDR_global_significance := NULL]
}
status[, FDR_global_significance := NA_real_]
status[status == "SUCCESS", FDR_global_significance :=
         p.adjust(pvalue_global_significance, method = "BH")]
status[, reference_model_status := fifelse(
  status != "SUCCESS", "MODEL_FAILED",
  fifelse(FDR_global_significance < 0.10,
          "REFERENCE_MODEL_VALIDATED", "REFERENCE_MODEL_WEAK")
)]
status[, numeric_stability := fifelse(
  status != "SUCCESS", "NOT_ASSESSABLE", "PENDING_RECAPITULATION_AUDIT")]
setorder(status, pathway_id)
fwrite(status, status_path)
fwrite(status, validity_path)

cat("success_n=", sum(status$status == "SUCCESS"), "\n", sep = "")
cat("failed_n=", sum(status$status == "FAILED"), "\n", sep = "")
cat("validated_n=", sum(status$reference_model_status == "REFERENCE_MODEL_VALIDATED"), "\n", sep = "")
cat("weak_n=", sum(status$reference_model_status == "REFERENCE_MODEL_WEAK"), "\n", sep = "")
cat("completed=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("PHASE11_HUMAN_MODEL_FIT_COMPLETE\n")

if (any(status$status == "FAILED")) quit(status = 2L)
