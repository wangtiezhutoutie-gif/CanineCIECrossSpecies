project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
lib <- file.path(project_root, "11_singIST", "environment", "r_library")
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(singIST)
  library(data.table)
  library(SingleCellExperiment)
})

set.seed(20260826)
singist_root <- file.path(project_root, "11_singIST")
dog_dir <- file.path(singist_root, "dog_model")
human_dir <- file.path(singist_root, "human_reference")
primary_dir <- file.path(singist_root, "primary", "method_faithful")
lodo_dir <- file.path(singist_root, "lodo_sensitivity")
checkpoint_dir <- file.path(lodo_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(project_root, "logs", "phase11_lodo.log")
zz <- file(log_file, open = "at")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({sink(type = "message"); sink(type = "output"); close(zz)}, add = TRUE)

cat("\nPHASE11 DOG LODO RECAPITULATION\n")
cat("started=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")

fit_status <- fread(file.path(human_dir, "HUMAN_REFERENCE_MODEL_FIT_STATUS.csv"))
stopifnot(nrow(fit_status) == 46L, all(fit_status$status == "SUCCESS"))
registered <- fread(file.path(singist_root, "superpathways",
                              "REGISTERED_PROGRAM_PANEL.csv"))
registered_ids <- unique(registered[testability == "TESTABLE", pathway_id])
fit_registered <- fit_status[pathway_id %in% registered_ids]
setorder(fit_registered, pathway_id)
mapping_full <- readRDS(file.path(dog_dir, "DOG_SINGIST_MAPPING_ORGANISM.rds"))
cache <- fread(file.path(singist_root, "method_gate",
                         "FROZEN_ENSEMBL104_ORTHOLOGY_CACHE.csv"))
drops <- c("h_1", "h_2", "h_3", "cie_2", "cie_3", "cie_4", "cie_6")

for (drop in drops) {
  cp <- file.path(checkpoint_dir, paste0("drop_", drop, "_RECAP.rds"))
  if (file.exists(cp)) {
    cat("SKIP_CHECKPOINT drop=", drop, "\n", sep = "")
    next
  }
  cat("LODO_RECAP_START drop=", drop, " pathways=", nrow(fit_registered), "\n", sep = "")
  fc <- readRDS(file.path(dog_dir, "lodo_fc",
                          paste0("drop_", drop,
                                 "_METHOD_FAITHFUL_FC_LIST.rds")))
  obj <- mapping_full
  obj$counts <- obj$counts[, obj$counts$donor != drop]
  super_rows <- list()
  cell_rows <- list()
  for (i in seq_len(nrow(fit_registered))) {
    fs <- fit_registered[i]
    model <- readRDS(fs$model_path)
    out <- singISTrecapitulations(
      obj, model, orthology_cache = cache, FC_list = fc
    )
    out$superpathway$pathway <- fs$pathway_id
    out$celltype$pathway <- fs$pathway_id
    super_rows[[i]] <- as.data.table(out$superpathway)[, `:=`(
      run_id = paste0("drop_", drop), drop_donor = drop,
      human_model_refit = FALSE
    )]
    cell_rows[[i]] <- as.data.table(out$celltype)[, `:=`(
      run_id = paste0("drop_", drop), drop_donor = drop,
      human_model_refit = FALSE
    )]
  }
  saveRDS(list(superpathway = rbindlist(super_rows),
               celltype = rbindlist(cell_rows)), cp, compress = "gzip")
  cat("LODO_RECAP_SUCCESS drop=", drop, "\n", sep = "")
  flush.console()
}

full_super <- fread(file.path(primary_dir, "SUPERPATHWAY_RECAPITULATION.csv"))[
  pathway %in% registered_ids]
full_cell <- fread(file.path(primary_dir, "CELLTYPE_RECAPITULATION.csv"))[
  pathway %in% registered_ids]
full_super[, `:=`(run_id = "Full7", drop_donor = NA_character_,
                  human_model_refit = FALSE)]
full_cell[, `:=`(run_id = "Full7", drop_donor = NA_character_,
                 human_model_refit = FALSE)]

res <- lapply(file.path(checkpoint_dir, paste0("drop_", drops, "_RECAP.rds")),
              readRDS)
lodo_super <- rbindlist(lapply(res, `[[`, "superpathway"), fill = TRUE)
lodo_cell <- rbindlist(lapply(res, `[[`, "celltype"), fill = TRUE)
super_long <- rbindlist(list(full_super, lodo_super), fill = TRUE)
cell_long <- rbindlist(list(full_cell, lodo_cell), fill = TRUE)

program_map <- registered[, .(
  registered_programs = paste(paste(lineage, program, sep = ": "),
                              collapse = "; "),
  registered_program_n = .N
), by = .(pathway = pathway_id)]
super_long <- merge(super_long, program_map, by = "pathway", all.x = TRUE, sort = FALSE)
cell_long <- merge(cell_long, program_map, by = "pathway", all.x = TRUE, sort = FALSE)

summary_rows <- lapply(registered_ids, function(pid) {
  z <- super_long[pathway == pid]
  full <- z[run_id == "Full7", recapitulation][1]
  ld <- z[run_id != "Full7"]
  delta <- ld$recapitulation - full
  influence_i <- which.max(abs(delta))
  data.table(
    pathway = pid,
    full7_recapitulation = full,
    median_LODO = median(ld$recapitulation),
    min_LODO = min(ld$recapitulation),
    max_LODO = max(ld$recapitulation),
    sign_stability = all(sign(ld$recapitulation) == sign(full)),
    same_sign_n = sum(sign(ld$recapitulation) == sign(full)),
    max_absolute_change = max(abs(delta)),
    most_influential_donor = ld$drop_donor[influence_i],
    cie_4_absolute_change = abs(ld[drop_donor == "cie_4", recapitulation][1] - full),
    cie_4_is_max_influence = ld$drop_donor[influence_i] == "cie_4",
    human_model_frozen = TRUE
  )
})
summary <- rbindlist(summary_rows)
summary <- merge(summary, program_map, by = "pathway", all.x = TRUE)

cell_summary <- cell_long[run_id != "Full7", {
  full <- cell_long[pathway == .BY$pathway & celltype == .BY$celltype &
                      run_id == "Full7", recapitulation][1]
  delta <- recapitulation - full
  i <- which.max(abs(delta))
  .(
    full7_recapitulation = full, median_LODO = median(recapitulation),
    min_LODO = min(recapitulation), max_LODO = max(recapitulation),
    sign_stability = all(sign(recapitulation) == sign(full)),
    same_sign_n = sum(sign(recapitulation) == sign(full)),
    max_absolute_change = max(abs(delta)),
    most_influential_donor = drop_donor[i],
    cie_4_absolute_change = abs(recapitulation[drop_donor == "cie_4"][1] - full),
    cie_4_is_max_influence = drop_donor[i] == "cie_4",
    human_model_frozen = TRUE
  )
}, by = .(pathway, celltype)]

fwrite(super_long, file.path(lodo_dir, "SINGIST_DOG_LODO_SUPERPATHWAY.csv"))
fwrite(cell_long, file.path(lodo_dir, "SINGIST_DOG_LODO_CELLTYPE.csv"))
fwrite(summary, file.path(lodo_dir, "SINGIST_DOG_LODO_SUMMARY.csv"))
fwrite(cell_summary, file.path(lodo_dir,
                               "SINGIST_DOG_LODO_CELLTYPE_SUMMARY.csv"))

cat("registered_program_n=", nrow(registered), "\n", sep = "")
cat("unique_registered_pathway_n=", length(registered_ids), "\n", sep = "")
cat("lodo_super_rows=", nrow(lodo_super), "\n", sep = "")
cat("cie4_max_pathway_n=", sum(summary$cie_4_is_max_influence), "\n", sep = "")
cat("completed=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("PHASE11_DOG_LODO_COMPLETE\n")
