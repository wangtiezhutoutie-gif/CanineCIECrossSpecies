project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
lib <- file.path(project_root, "11_singIST", "environment", "r_library")
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(singIST)
  library(data.table)
  library(SingleCellExperiment)
})

set.seed(20260826)
dog_dir <- file.path(project_root, "11_singIST", "dog_model")
out_dir <- file.path(dog_dir, "lodo_fc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(out_dir, "DOG_LODO_FC_STATUS.csv")
log_file <- file.path(project_root, "logs", "phase11_lodo_fc.log")
zz <- file(log_file, open = "at")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({sink(type = "message"); sink(type = "output"); close(zz)}, add = TRUE)

cat("\nPHASE11 DOG LODO FC PRECOMPUTE\n")
cat("started=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
mapping_full <- readRDS(file.path(dog_dir, "DOG_SINGIST_MAPPING_ORGANISM.rds"))
drops <- c("h_1", "h_2", "h_3", "cie_2", "cie_3", "cie_4", "cie_6")
required_fc <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
status <- if (file.exists(status_path)) fread(status_path) else data.table()

for (i in seq_along(drops)) {
  drop <- drops[[i]]
  out <- file.path(out_dir, paste0("drop_", drop, "_METHOD_FAITHFUL_FC_LIST.rds"))
  old <- if (nrow(status)) status[drop_donor == drop] else data.table()
  if (file.exists(out) && nrow(old) == 1L && old$status == "SUCCESS") {
    cat("SKIP_SUCCESS ", drop, "\n", sep = "")
    next
  }
  start <- Sys.time()
  cat("LODO_FC_START ", i, "/", length(drops), " drop=", drop, "\n", sep = "")
  row <- tryCatch({
    keep <- mapping_full$counts$donor != drop
    obj <- mapping_full
    obj$counts <- obj$counts[, keep]
    stopifnot(length(unique(obj$counts$donor)) == 6L)
    stopifnot(length(unique(obj$counts$donor[obj$counts$class == "Healthy"])) >= 2L)
    stopifnot(length(unique(obj$counts$donor[obj$counts$class == "CIE"])) >= 3L)
    mapped <- celltype_mapping(obj)
    fc <- diff_expressed(mapped, logfc.treshold = 0.25)
    stopifnot(all(vapply(fc, function(x) all(required_fc %in% colnames(x)), logical(1))))
    fc <- lapply(fc, function(x) x[, required_fc, drop = FALSE])
    saveRDS(fc, out, compress = "gzip")
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    data.table(
      drop_donor = drop, status = "SUCCESS",
      remaining_healthy_donors = length(unique(obj$counts$donor[obj$counts$class == "Healthy"])),
      remaining_CIE_donors = length(unique(obj$counts$donor[obj$counts$class == "CIE"])),
      epithelial_FDR05_n = sum(fc$Epithelial$p_val_adj <= 0.05, na.rm = TRUE),
      T_FDR05_n = sum(fc$T$p_val_adj <= 0.05, na.rm = TRUE),
      myeloid_FDR05_n = sum(fc$Myeloid$p_val_adj <= 0.05, na.rm = TRUE),
      elapsed_seconds = elapsed, output_path = out, error = NA_character_
    )
  }, error = function(e) {
    data.table(
      drop_donor = drop, status = "FAILED",
      remaining_healthy_donors = NA_integer_, remaining_CIE_donors = NA_integer_,
      epithelial_FDR05_n = NA_integer_, T_FDR05_n = NA_integer_,
      myeloid_FDR05_n = NA_integer_,
      elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = "secs")),
      output_path = out, error = conditionMessage(e)
    )
  })
  if (nrow(status)) status <- status[drop_donor != drop]
  status <- rbindlist(list(status, row), fill = TRUE)
  status[, order__ := match(drop_donor, drops)]
  setorder(status, order__)
  status[, order__ := NULL]
  fwrite(status, status_path)
  cat("LODO_FC_", row$status, " drop=", drop,
      " elapsed_s=", round(row$elapsed_seconds, 1), "\n", sep = "")
  flush.console()
}

cat("success_n=", sum(status$status == "SUCCESS"), "\n", sep = "")
cat("completed=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("PHASE11_DOG_LODO_FC_COMPLETE\n")
if (any(status$status != "SUCCESS")) quit(status = 2L)
