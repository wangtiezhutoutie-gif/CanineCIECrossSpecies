options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: phase6_extract_author_rds.R <input.rds> <output_dir>")
input_rds <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)

started <- Sys.time()
obj <- readRDS(input_rds)

safe_slot <- function(x, nm) { attr(x, nm, exact = TRUE) }

meta <- safe_slot(obj, "meta.data")
if (is.null(meta) || !is.data.frame(meta)) stop("Seurat meta.data slot unavailable")
meta$author_cell_id <- rownames(meta)
meta <- meta[, c("author_cell_id", setdiff(colnames(meta), "author_cell_id")), drop = FALSE]

assays <- safe_slot(obj, "assays")
assay_names <- names(assays)
assay_lines <- character()
raw_count_available <- FALSE
raw_count_assay <- NA_character_
raw_count_nrow <- NA_integer_
raw_count_ncol <- NA_integer_

for (nm in assay_names) {
  assay <- assays[[nm]]
  sn <- names(attributes(assay))
  dims <- character()
  if ("counts" %in% sn) {
    m <- attr(assay, "counts", exact = TRUE)
    dims <- c(dims, paste0("counts=", paste(c(nrow(m), ncol(m)), collapse = "x")))
    if (!raw_count_available && length(dim(m)) == 2L && nrow(m) > 0L && ncol(m) > 0L) {
      raw_count_available <- TRUE
      raw_count_assay <- nm
      raw_count_nrow <- nrow(m)
      raw_count_ncol <- ncol(m)
    }
  }
  if ("data" %in% sn) {
    m <- attr(assay, "data", exact = TRUE)
    dims <- c(dims, paste0("data=", paste(c(nrow(m), ncol(m)), collapse = "x")))
  }
  if ("scale.data" %in% sn) {
    m <- attr(assay, "scale.data", exact = TRUE)
    dims <- c(dims, paste0("scale.data=", paste(c(nrow(m), ncol(m)), collapse = "x")))
  }
  if ("layers" %in% sn) {
    layers <- attr(assay, "layers", exact = TRUE)
    dims <- c(dims, paste0("layers=", paste(names(layers), collapse = ",")))
    if (!raw_count_available && "counts" %in% names(layers)) {
      m <- layers[["counts"]]
      raw_count_available <- TRUE
      raw_count_assay <- nm
      raw_count_nrow <- nrow(m)
      raw_count_ncol <- ncol(m)
    }
  }
  assay_lines <- c(assay_lines, paste(nm, class(assay)[1], paste(sn, collapse = ","), paste(dims, collapse = ";"), sep = "\t"))
}

active_assay <- safe_slot(obj, "active.assay")
active_ident <- safe_slot(obj, "active.ident")
reductions <- safe_slot(obj, "reductions")
commands <- safe_slot(obj, "commands")
project_name <- safe_slot(obj, "project.name")
object_version <- safe_slot(obj, "version")

write.csv(meta, file.path(output_dir, "author_metadata.csv"), row.names = FALSE, na = "")
writeLines(c(
  paste0("input_rds=", input_rds),
  paste0("object_class=", paste(class(obj), collapse = ",")),
  paste0("object_version=", paste(object_version, collapse = ",")),
  paste0("project_name=", paste(project_name, collapse = ",")),
  paste0("n_cells=", nrow(meta)),
  paste0("n_metadata_columns=", ncol(meta) - 1L),
  paste0("metadata_columns=", paste(setdiff(colnames(meta), "author_cell_id"), collapse = ",")),
  paste0("active_assay=", paste(active_assay, collapse = ",")),
  paste0("assays=", paste(assay_names, collapse = ",")),
  paste0("raw_count_available=", raw_count_available),
  paste0("raw_count_assay=", raw_count_assay),
  paste0("raw_count_dimensions=", raw_count_nrow, "x", raw_count_ncol),
  paste0("reductions=", paste(names(reductions), collapse = ",")),
  paste0("command_count=", length(commands)),
  paste0("active_ident_levels=", paste(levels(active_ident), collapse = "|")),
  paste0("started=", format(started, "%Y-%m-%dT%H:%M:%S%z")),
  paste0("finished=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "",
  "assay\tclass\tslots\tdimensions",
  assay_lines
), file.path(output_dir, "author_object_summary.txt"))

cat("AUTHOR_RDS_EXTRACTION_COMPLETE\n")
cat("N_CELLS=", nrow(meta), "\n", sep = "")
cat("RAW_COUNTS=", raw_count_available, "\n", sep = "")
