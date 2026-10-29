options(timeout = 1200, warn = 1)
root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
p11 <- file.path(root, "11_singIST")
lib <- file.path(p11, "environment", "r_library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
log_file <- file.path(root, "logs", "phase11_environment.log")
zz <- file(log_file, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({sink(type = "message"); sink(type = "output"); close(zz)}, add = TRUE)

cat("PHASE11_ENVIRONMENT_START\n")
cat("R=", R.version.string, "\n", sep = "")
cat("LIB=", lib, "\n", sep = "")
options(BioC_mirror = "https://ftp.gwdg.de/pub/misc/bioconductor")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = lib,
                   type = "binary")
}
BiocManager::install(version = "3.23", ask = FALSE, update = FALSE, lib = lib)
cat("BIOCONDUCTOR=", as.character(BiocManager::version()), "\n", sep = "")

imports <- c(
  "GSEABase", "checkmate", "asmbPLS", "BiocParallel",
  "stringr", "FactoMineR", "Seurat", "SeuratObject", "biomaRt",
  "data.table", "purrr", "SingleCellExperiment", "SummarizedExperiment",
  "scran", "scuttle", "missMDA", "S4Vectors"
)
BiocManager::install(imports, ask = FALSE, update = FALSE, lib = lib,
                     type = "binary",
                     dependencies = c("Depends", "Imports", "LinkingTo"))

# msigdb is a Bioconductor ExperimentData package and has no Windows binary.
# Install the already frozen official 3.23 source-only chain in dependency order.
source_chain <- c(
  BiocBaseUtils = "BiocBaseUtils_1.14.2.tar.gz",
  AnnotationHub = "AnnotationHub_4.2.2.tar.gz",
  ExperimentHub = "ExperimentHub_3.2.2.tar.gz",
  org.Mm.eg.db = "org.Mm.eg.db_3.23.0.tar.gz",
  org.Hs.eg.db = "org.Hs.eg.db_3.23.1.tar.gz"
)
for (package_name in names(source_chain)) {
  archive <- unname(source_chain[[package_name]])
  archive_path <- file.path(p11, "reference", archive)
  stopifnot(file.exists(archive_path), file.info(archive_path)$size > 0)
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(archive_path, repos = NULL, type = "source", lib = lib)
  }
}
msigdb_tarball <- file.path(p11, "reference", "msigdb_1.20.0.tar.gz")
stopifnot(file.exists(msigdb_tarball))
install.packages(msigdb_tarball, repos = NULL, type = "source", lib = lib)
stopifnot(requireNamespace("msigdb", quietly = TRUE))
stopifnot(as.character(packageVersion("msigdb")) == "1.20.0")
cat("MSIGDB=", as.character(packageVersion("msigdb")), "\n", sep = "")

tarball <- file.path(p11, "reference", "singIST_1.0.2.tar.gz")
stopifnot(file.exists(tarball))
install.packages(tarball, repos = NULL, type = "source", lib = lib)

stopifnot(requireNamespace("singIST", quietly = TRUE))
stopifnot(as.character(packageVersion("singIST")) == "1.0.2")
cat("SINGIST=", as.character(packageVersion("singIST")), "\n", sep = "")
cat("SINGIST_PATH=", find.package("singIST"), "\n", sep = "")
cat("PHASE11_ENVIRONMENT_COMPLETE\n")
sessionInfo()
