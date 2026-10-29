project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
lib <- file.path(project_root, "11_singIST", "environment", "r_library")
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(singIST)
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scuttle)
  library(S4Vectors)
})

set.seed(20260826)
dog_dir <- file.path(project_root, "11_singIST", "dog_model")
dir.create(dog_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(project_root, "logs", "phase11_dog_model.log")
zz <- file(log_file, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({sink(type = "message"); sink(type = "output"); close(zz)}, add = TRUE)

cat("PHASE11 DOG MODEL BUILD\n")
cat("started=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")

mtx_path <- file.path(dog_dir, "DOG_SINGIST_COUNTS_GENES_X_CELLS.mtx.gz")
cell_path <- file.path(dog_dir, "DOG_SINGIST_CELL_METADATA.csv.gz")
gene_path <- file.path(dog_dir, "DOG_SINGIST_GENE_METADATA.csv")

counts <- Matrix::readMM(gzfile(mtx_path))
counts <- as(counts, "dgCMatrix")
cells <- as.data.table(read.csv(gzfile(cell_path), check.names = FALSE,
                                stringsAsFactors = FALSE))
genes <- fread(gene_path)
stopifnot(nrow(counts) == nrow(genes), ncol(counts) == nrow(cells))
rownames(counts) <- genes$dog_unique_var_name
colnames(counts) <- cells$cell_id
stopifnot(!anyDuplicated(rownames(counts)), !anyDuplicated(colnames(counts)))
stopifnot(all(counts@x == round(counts@x)))

cd <- S4Vectors::DataFrame(
  sample_id = cells$sample_id,
  donor_full = cells$donor_full,
  donor = cells$donor,
  class = cells$class,
  celltype_cluster = cells$celltype_cluster,
  cl_id = cells$cl_id,
  mapping_confidence = cells$mapping_confidence,
  row.names = cells$cell_id
)
rd <- S4Vectors::DataFrame(genes, row.names = genes$dog_unique_var_name)
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = cd, rowData = rd)
# A logcounts assay is required by the official SCE branch of diff_expressed().
# This is method preprocessing only; QC, clustering, and annotations are frozen.
sce <- scuttle::logNormCounts(sce, assay.type = "counts")

mapping <- create_mapping_organism(
  organism = "Canis lupus familiaris", target_class = "CIE",
  base_class = "Healthy",
  celltype_mapping = list(
    Epithelial = "Epithelial", T = "T", Myeloid = "Myeloid"
  ),
  counts = sce
)
saveRDS(mapping, file.path(dog_dir, "DOG_SINGIST_MAPPING_ORGANISM.rds"),
        compress = "gzip")

# Freeze singIST-compatible donor x cell-type pseudobulk log-normalized profiles
# using the exact SCE route implemented by singIST::pseudobulk_log2FC().
grp <- paste0(sce$celltype_cluster, "_", sce$class, "_", sce$donor)
ids <- S4Vectors::DataFrame(
  celltype_cluster = sce$celltype_cluster,
  class = sce$class,
  donor = sce$donor,
  row.names = grp
)
sce_pb <- scuttle::aggregateAcrossCells(sce, ids = ids,
                                        use.assay.type = "counts")
sce_pb <- scuttle::computeLibraryFactors(sce_pb)
sce_pb <- scuttle::logNormCounts(
  sce_pb, assay.type = "counts",
  size.factors = SingleCellExperiment::sizeFactors(sce_pb),
  pseudo_count = 1, center_size_factors = FALSE
)
dog_pb <- t(as.matrix(assay(sce_pb, "logcounts")))
saveRDS(dog_pb, file.path(dog_dir, "DOG_SINGIST_PSEUDOBULK_LOGNORM.rds"),
        compress = "gzip")
pb_meta <- as.data.table(as.data.frame(colData(sce_pb)))[, .(
  profile_id = rownames(dog_pb), celltype = celltype_cluster,
  condition = class, donor,
  raw_library_size = as.numeric(colSums(assay(sce_pb, "counts"))),
  size_factor = as.numeric(sizeFactors(sce_pb)),
  gene_count = ncol(dog_pb),
  detected_gene_count = as.integer(colSums(assay(sce_pb, "counts") > 0)),
  normalization = paste0(
    "singIST1.0.2 SCE pseudobulk route: aggregateAcrossCells(counts); ",
    "computeLibraryFactors; logNormCounts(pseudo_count=1,center_size_factors=FALSE)"
  )
)]
fwrite(pb_meta, file.path(dog_dir, "DOG_SINGIST_PSEUDOBULK_METADATA.csv"))

# Method-faithful disease-model statistics. This is computed once and reused
# unchanged across pathways; biological_link_function would otherwise repeat
# this identical operation for every model.
mapped_for_de <- celltype_mapping(mapping)
cat("METHOD_FAITHFUL_DIFF_EXPRESSED_START\n")
method_fc <- diff_expressed(mapped_for_de, logfc.treshold = 0.25)
required_fc <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
stopifnot(identical(names(method_fc), c("Epithelial", "T", "Myeloid")))
cat("diff_expressed_return_columns=",
    paste(colnames(method_fc[[1]]), collapse = ";"), "\n", sep = "")
stopifnot(all(vapply(method_fc, function(x) all(required_fc %in% colnames(x)),
                     logical(1))))
# diff_expressed() appends avg_log2FC after the four inferential columns, while
# biological_link_function(FC_list=) requires this exact alternate order.
# Reordering is necessary for package API compatibility and changes no values.
method_fc <- lapply(method_fc, function(x) x[, required_fc, drop = FALSE])
saveRDS(method_fc, file.path(dog_dir, "DOG_METHOD_FAITHFUL_FC_LIST.rds"),
        compress = "gzip")

# Donor-aware Track B retains the exact singIST descriptive avg_log2FC and
# detection proportions, replacing only inferential fields with Phase3 Full7
# donor-level edgeR PValue and FDR. Genes not tested after Phase3 filtering are
# explicitly ineligible (p=FDR=1), never assigned cell-level significance.
donor_fc <- method_fc
audit_long <- list()
for (celltype in names(donor_fc)) {
  de_path <- file.path(project_root, "06_pseudobulk", "de_results",
                       paste0(celltype, "_CIE_vs_Healthy_edger.csv"))
  de <- fread(de_path)
  stopifnot(all(c("unique_var_name", "PValue", "FDR") %in% names(de)))
  x <- donor_fc[[celltype]]
  mi <- match(rownames(x), de$unique_var_name)
  tested <- !is.na(mi)
  x$p_val <- 1
  x$p_val_adj <- 1
  x$p_val[tested] <- de$PValue[mi[tested]]
  x$p_val_adj[tested] <- de$FDR[mi[tested]]
  x$p_val[is.na(x$p_val)] <- 1
  x$p_val_adj[is.na(x$p_val_adj)] <- 1
  donor_fc[[celltype]] <- x[, required_fc, drop = FALSE]
  audit_long[[celltype]] <- data.table(
    celltype = celltype, dog_gene = rownames(x),
    avg_log2FC = x$avg_log2FC, pct_1_CIE = x$pct.1,
    pct_2_Healthy = x$pct.2, phase3_PValue = x$p_val,
    phase3_FDR = x$p_val_adj,
    phase3_tested = tested,
    treatment_eligible_FDR05 = x$p_val_adj <= 0.05
  )
}
saveRDS(donor_fc, file.path(dog_dir, "DOG_DONOR_AWARE_FC_LIST.rds"),
        compress = "gzip")
fwrite(rbindlist(audit_long),
       file.path(dog_dir, "DOG_FC_SOURCE_AUDIT.csv.gz"), compress = "gzip")

fc_summary <- rbindlist(lapply(names(method_fc), function(celltype) {
  data.table(
    celltype = celltype,
    method_faithful_gene_n = nrow(method_fc[[celltype]]),
    method_faithful_FDR05_n = sum(method_fc[[celltype]]$p_val_adj <= 0.05,
                                  na.rm = TRUE),
    phase3_tested_n = sum(audit_long[[celltype]]$phase3_tested),
    donor_aware_FDR05_n = sum(donor_fc[[celltype]]$p_val_adj <= 0.05,
                              na.rm = TRUE),
    avg_log2FC_definition = paste0(
      "mean(CIE donor pseudobulk log-normalized expression)-",
      "mean(Healthy donor pseudobulk log-normalized expression)"
    ),
    pct_1_group = "CIE_target",
    pct_2_group = "Healthy_base",
    formal_dog_disease_inference = "Phase3 Full7 donor-level edgeR only"
  )
}))
fwrite(fc_summary, file.path(dog_dir, "DOG_FC_SOURCE_SUMMARY.csv"))

cat("cells=", ncol(sce), "\n", sep = "")
cat("genes=", nrow(sce), "\n", sep = "")
print(fc_summary)
cat("completed=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("PHASE11_DOG_MODEL_COMPLETE\n")
