project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
lib <- file.path(project_root, "11_singIST", "environment", "r_library")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(singIST))

cat("SINGIST_VERSION=", as.character(packageVersion("singIST")), "\n", sep = "")
funs <- c(
  "create_mapping_organism", "create_pathway", "create_superpathway",
  "create_hyperparameters", "create_superpathway_input", "fitOptimal",
  "diff_expressed", "orthology_mapping", "singIST_treat",
  "biological_link_function", "derive_contributions",
  "superpathway_recap", "celltype_recap", "gene_contrib",
  "singISTrecapitulations"
)
for (nm in funs) {
  cat("\n=== FORMALS ", nm, " ===\n", sep = "")
  print(formals(getExportedValue("singIST", nm)))
}

ext <- system.file("extdata", package = "singIST")
cat("\nEXTDATA=", ext, "\n", sep = "")
print(list.files(ext, full.names = FALSE))

load(file.path(ext, "example_superpathway_input.rda"))
load(file.path(ext, "example_superpathway_fit_model.rda"))
load(file.path(ext, "example_mapping_organism.rda"))
load(file.path(ext, "example_orthology_mapping.rda"))

for (nm in c("example_superpathway_input", "example_superpathway_fit_model",
             "example_mapping_organism", "example_orthology_mapping")) {
  cat("\n=== STR ", nm, " ===\n", sep = "")
  str(get(nm), max.level = 3)
}

cat("\n=== VIGNETTE HYPERPARAMETER LINES ===\n")
vig <- readLines(file.path(project_root, "11_singIST", "reference",
                           "singIST_1.0.2_source", "singIST", "vignettes",
                           "singISTvignette.Rmd"), warn = FALSE)
hits <- grep("quantile|fitOptimal|folds_CV|number_PLS|npermut|nsubsampling|B_accuracy|jackknife", vig)
for (i in hits) cat(sprintf("%04d: %s\n", i, vig[[i]]))
