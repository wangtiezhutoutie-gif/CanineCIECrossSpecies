project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
lib <- file.path(project_root, "11_singIST", "environment", "r_library")
Sys.setenv(R_LIBS_USER = lib)
.libPaths(c(lib, .libPaths()))
library(data.table)

human_dir <- file.path(project_root, "11_singIST", "human_reference")
model_dir <- file.path(human_dir, "models")
input_dir <- file.path(human_dir, "model_inputs")
status_path <- file.path(human_dir, "HUMAN_REFERENCE_MODEL_FIT_STATUS.csv")
validity_path <- file.path(project_root, "11_singIST", "primary",
                           "method_faithful", "HUMAN_REFERENCE_MODEL_VALIDITY.csv")
checkpoint_dir <- file.path(project_root, "11_singIST", "primary", "checkpoints")

status <- fread(status_path)
safe_stem <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

for (i in seq_len(nrow(status))) {
  id <- status$pathway_id[i]
  old_model <- file.path(model_dir, paste0(id, ".rds"))
  old_input <- file.path(input_dir, paste0(id, "_input.rds"))
  new_model <- file.path(model_dir, paste0(safe_stem(id), ".rds"))
  new_input <- file.path(input_dir, paste0(safe_stem(id), "_input.rds"))
  if (!file.exists(new_model)) {
    if (!file.exists(old_model)) stop("Missing legacy model: ", old_model)
    saveRDS(readRDS(old_model), new_model, compress = "gzip")
  }
  if (!file.exists(new_input)) {
    if (!file.exists(old_input)) stop("Missing legacy input: ", old_input)
    saveRDS(readRDS(old_input), new_input, compress = "gzip")
  }
  old_checkpoint <- file.path(checkpoint_dir,
                              paste0(id, "_PRIMARY_TRACKS.rds"))
  new_checkpoint <- file.path(checkpoint_dir,
                              paste0(safe_stem(id), "_PRIMARY_TRACKS.rds"))
  if (file.exists(old_checkpoint) && !file.exists(new_checkpoint)) {
    saveRDS(readRDS(old_checkpoint), new_checkpoint, compress = "gzip")
  }
  status$model_path[i] <- new_model
}

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

# If primary recapitulation already exists, preserve/reconstruct the completed
# numerical-stability audit instead of replacing it with a pending status during
# a later path-migration rerun.
super_path <- file.path(project_root, "11_singIST", "primary",
                        "method_faithful", "SUPERPATHWAY_RECAPITULATION.csv")
validity <- copy(status)
if (file.exists(super_path)) {
  super <- fread(super_path)[, .(
    pathway_id = pathway, reference_shift, reference_shift_status
  )]
  validity <- merge(validity, super, by = "pathway_id", all.x = TRUE)
  validity[, numeric_stability := reference_shift_status]
}
fwrite(validity, validity_path)

stopifnot(all(file.exists(status$model_path)))
cat("success_n=", sum(status$status == "SUCCESS"), "\n", sep = "")
cat("validated_n=", sum(status$reference_model_status ==
                           "REFERENCE_MODEL_VALIDATED"), "\n", sep = "")
cat("weak_n=", sum(status$reference_model_status ==
                      "REFERENCE_MODEL_WEAK"), "\n", sep = "")
cat("safe_model_files=", length(list.files(model_dir, pattern = "\\.rds$")),
    "\n", sep = "")
cat("safe_checkpoint_files=",
    length(list.files(checkpoint_dir, pattern = "\\.rds$")), "\n", sep = "")
cat("PHASE11_MODEL_PATH_MIGRATION_COMPLETE\n")
