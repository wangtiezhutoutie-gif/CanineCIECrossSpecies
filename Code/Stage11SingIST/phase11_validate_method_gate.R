project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
lib <- file.path(project_root, "11_singIST", "environment", "r_library")
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(singIST)
  library(digest)
})

set.seed(20260826)
gate_dir <- file.path(project_root, "11_singIST", "method_gate")
audit_dir <- file.path(project_root, "11_singIST", "audit")
dir.create(gate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(project_root, "logs", "phase11_method_gate.log")
zz <- file(log_file, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(zz)
}, add = TRUE)

cat("PHASE11 METHOD GATE\n")
cat("started=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("singIST=", as.character(packageVersion("singIST")), "\n", sep = "")
cat("R=", R.version.string, "\n", sep = "")

ext <- system.file("extdata", package = "singIST")
load(file.path(ext, "example_superpathway_input.rda"))
load(file.path(ext, "example_superpathway_fit_model.rda"))
load(file.path(ext, "example_mapping_organism.rda"))
load(file.path(ext, "example_orthology_mapping.rda"))

validation <- data.frame(
  check = character(), status = character(), value = character(),
  tolerance = character(), notes = character(), stringsAsFactors = FALSE
)
add_check <- function(check, status, value = "", tolerance = "", notes = "") {
  validation <<- rbind(validation, data.frame(
    check = check, status = status, value = as.character(value),
    tolerance = tolerance, notes = notes, stringsAsFactors = FALSE
  ))
}

add_check("package_load", "PASS", as.character(packageVersion("singIST")), "",
          "Frozen Bioconductor 3.23 singIST package loaded")

# Reconstruct all documented input classes through their public constructors.
pi <- example_superpathway_input$superpathway_info$pathway_info
pathway <- create_pathway(pi$standard_name, pi$dbsource, pi$collection, pi$subcollection)
sp <- create_superpathway(
  pathway,
  example_superpathway_input$superpathway_info$celltypes,
  example_superpathway_input$superpathway_info$gene_sets_celltype
)
hp0 <- example_superpathway_input$hyperparameters_info
hp <- create_hyperparameters(
  hp0$quantile_comb_table, hp0$outcome_type, as.integer(hp0$number_PLS),
  as.integer(hp0$folds_CV), as.integer(hp0$repetition_CV)
)
sp_input <- create_superpathway_input(
  sp, hp, example_superpathway_input$pseudobulk_lognorm,
  example_superpathway_input$sample_id,
  example_superpathway_input$sample_class,
  example_superpathway_input$base_class,
  example_superpathway_input$target_class
)
mapped <- create_mapping_organism(
  example_mapping_organism$organism,
  example_mapping_organism$target_class,
  example_mapping_organism$base_class,
  example_mapping_organism$celltype_mapping,
  example_mapping_organism$counts
)
add_check("create_pathway", "PASS", pathway$standard_name)
add_check("create_superpathway", "PASS", length(sp$celltypes), "",
          "Public constructor; two example cell-type blocks")
add_check("create_hyperparameters", "PASS", nrow(hp$quantile_comb_table))
add_check("create_superpathway_input", "PASS", nrow(sp_input$pseudobulk_lognorm))
add_check("create_mapping_organism", "PASS", mapped$organism)

# Actually fit the bundled example with a small gate-only permutation count.
# Formal project models use the separately preregistered npermut=1000.
fit_gate <- fitOptimal(
  sp_input, parallel = FALSE, measure = "B_accuracy",
  expected_measure_increase = 0.005, maxiter = 100,
  global_significance_full = FALSE, CIP.GIP_significance_full = FALSE,
  npermut = 10, type = "jackknife", nsubsampling = 10
)
fit_ok <- is.list(fit_gate) && !is.null(fit_gate$model_fit$`asmbPLS-DA`) &&
  is.finite(fit_gate$model_validation$pvalue_global_significance)
add_check("fitOptimal_example", if (fit_ok) "PASS" else "FAIL",
          fit_gate$model_validation$pvalue_global_significance, "",
          "Gate-only npermut=10; not a project biological result")

# Official offline wrapper route using the package-shipped cache.
wrapped <- singISTrecapitulations(
  mapped, example_superpathway_fit_model,
  orthology_cache = example_orthology_mapping
)
add_check("singISTrecapitulations_offline", "PASS",
          wrapped$superpathway$recapitulation, "",
          "Official orthology_cache API; no BioMart query")
add_check("superpathway_recapitulation", "PASS", nrow(wrapped$superpathway))
add_check("celltype_recapitulation", "PASS", nrow(wrapped$celltype))
add_check("gene_contributions", "PASS", nrow(wrapped$gene))

# The package returns diff_expressed columns in an order different from the
# exact order required by biological_link_function(FC_list=). Reorder only,
# then verify that precomputing the official DE once is numerically identical
# to letting the wrapper recompute it internally.
mapped_de <- celltype_mapping(mapped)
mapped_de$counts$test <- paste0(mapped_de$counts$celltype_cluster, "_",
                                mapped_de$counts$class)
SeuratObject::Idents(mapped_de$counts) <- "test"
precomputed_fc <- diff_expressed(mapped_de, logfc.treshold = 0.25)
fc_required <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
stopifnot(all(vapply(precomputed_fc,
                     function(x) all(fc_required %in% colnames(x)), logical(1))))
precomputed_fc <- lapply(precomputed_fc,
                         function(x) x[, fc_required, drop = FALSE])
names(precomputed_fc) <- gsub("_", " ", names(precomputed_fc))
wrapped_precomputed <- singISTrecapitulations(
  mapped, example_superpathway_fit_model,
  orthology_cache = example_orthology_mapping,
  FC_list = precomputed_fc
)

# Independent composition of the same exported official functions. This is the
# numerical-equivalence gate for the frozen-cache route. No package math is copied.
link <- biological_link_function(
  mapped, example_superpathway_fit_model,
  orthology_cache = example_orthology_mapping
)
ref <- derive_contributions(
  example_superpathway_fit_model,
  example_superpathway_fit_model$model_fit$predictor_block
)
treated <- derive_contributions(example_superpathway_fit_model, link$singIST_samples)
direct_super <- superpathway_recap(
  example_superpathway_fit_model,
  ref$superpathway_score, treated$superpathway_score
)
direct_cell_raw <- celltype_recap(
  example_superpathway_fit_model,
  ref$celltype_contribution, treated$celltype_contribution
)
direct_gene <- gene_contrib(
  example_superpathway_fit_model,
  ref$gene_contribution, treated$gene_contribution, direct_cell_raw
)
direct_cell <- direct_cell_raw[, c("pathway", "celltype", "recapitulation")]

max_abs <- function(x) if (length(x)) max(abs(x), na.rm = TRUE) else NA_real_
super_diff <- max_abs(wrapped$superpathway$recapitulation - direct_super$recapitulation)
cell_diff <- max_abs(wrapped$celltype$recapitulation - direct_cell$recapitulation)
gene_diff <- max_abs(wrapped$gene$contribution - direct_gene$contribution)
pre_super_diff <- max_abs(wrapped$superpathway$recapitulation -
                            wrapped_precomputed$superpathway$recapitulation)
pre_cell_diff <- max_abs(wrapped$celltype$recapitulation -
                           wrapped_precomputed$celltype$recapitulation)
pre_gene_diff <- max_abs(wrapped$gene$contribution -
                           wrapped_precomputed$gene$contribution)

# Explicitly compare the public orthology_mapping cache path too.
om <- orthology_mapping(
  mapped, example_superpathway_fit_model, from_species = "hsapiens",
  to_species = "mmusculus", orthology_cache = example_orthology_mapping
)
orth_identical <- identical(
  unname(lapply(om, as.data.frame)),
  unname(lapply(link$orthologs, as.data.frame))
)

# Gene contributions should add to each unrounded cell-type recapitulation.
gene_sum <- aggregate(contribution ~ pathway + celltype, direct_gene, sum)
cell_sum <- merge(
  gene_sum, direct_cell_raw[, c("pathway", "celltype", "recapitulation")],
  by = c("pathway", "celltype"), all = TRUE
)
additivity_diff <- max_abs(cell_sum$contribution - cell_sum$recapitulation)

tol <- 1e-8
eq <- data.frame(
  comparison = c(
    "official_cache_mapping_vs_biological_link_mapping",
    "wrapper_vs_exported_steps_superpathway",
    "wrapper_vs_exported_steps_celltype",
    "wrapper_vs_exported_steps_gene_contribution",
    "internal_DE_vs_precomputed_official_DE_superpathway",
    "internal_DE_vs_precomputed_official_DE_celltype",
    "internal_DE_vs_precomputed_official_DE_gene_contribution",
    "gene_contribution_sum_vs_celltype_recapitulation"
  ),
  route_A = c(
    "orthology_mapping(orthology_cache)",
    rep("singISTrecapitulations(orthology_cache)", 3),
    rep("singISTrecapitulations(internal diff_expressed)", 3),
    "sum(gene_contrib)"
  ),
  route_B = c(
    "biological_link_function(orthology_cache)$orthologs",
    rep("official exported stepwise functions", 3),
    rep("singISTrecapitulations(FC_list=precomputed official diff_expressed)", 3),
    "celltype_recap"
  ),
  max_absolute_difference = c(
    if (orth_identical) 0 else Inf, super_diff, cell_diff, gene_diff,
    pre_super_diff, pre_cell_diff, pre_gene_diff,
    additivity_diff
  ),
  tolerance = tol,
  status = c(
    if (orth_identical) "PASS" else "FAIL",
    ifelse(c(super_diff, cell_diff, gene_diff, pre_super_diff, pre_cell_diff,
             pre_gene_diff, additivity_diff) <= tol,
           "PASS", "FAIL")
  ),
  adapter_status = "ADAPTER_NOT_REQUIRED_OFFICIAL_ORTHOLOGY_CACHE_API",
  notes = c(
    "Package 1.0.2 accepts frozen input_gene/output_gene mappings and skips BioMart",
    "No custom recapitulation implementation",
    "No custom recapitulation implementation",
    "No custom contribution implementation",
    "Official diff_expressed values; columns reordered only for FC_list validator",
    "Official diff_expressed values; columns reordered only for FC_list validator",
    "Official diff_expressed values; columns reordered only for FC_list validator",
    "Unrounded official cell-type result used"
  ),
  stringsAsFactors = FALSE
)

add_check("orthology_cache_official_API", if (orth_identical) "PASS" else "FAIL",
          orth_identical, "exact",
          "Adapter not required; official API supports frozen mapping injection")
add_check("numerical_equivalence", if (all(eq$status == "PASS")) "PASS" else "FAIL",
          max(eq$max_absolute_difference), "absolute difference <= 1e-8")
add_check("gene_contribution_additivity", if (additivity_diff <= tol) "PASS" else "FAIL",
          additivity_diff, "absolute difference <= 1e-8")
add_check("precomputed_official_diff_expressed_equivalence",
          if (max(c(pre_super_diff, pre_cell_diff, pre_gene_diff)) <= tol) "PASS" else "FAIL",
          max(c(pre_super_diff, pre_cell_diff, pre_gene_diff)),
          "absolute difference <= 1e-8",
          "Enables one official DE computation reused across project pathways")

write.csv(validation,
          file.path(gate_dir, "OFFICIAL_EXAMPLE_VALIDATION.csv"), row.names = FALSE)
write.csv(eq,
          file.path(gate_dir, "FROZEN_ORTHOLOG_ADAPTER_VALIDATION.csv"),
          row.names = FALSE)

source_dir <- file.path(project_root, "11_singIST", "reference",
                        "singIST_1.0.2_source", "singIST")
src_sha <- digest(file = file.path(project_root, "11_singIST", "reference",
                                   "singIST_1.0.2.tar.gz"), algo = "sha256")
formals_text <- function(nm) {
  f <- formals(getExportedValue("singIST", nm))
  paste0("`", nm, "(", paste(names(f), collapse = ", "), ")`")
}
aud <- c(
  "# singIST 1.0.2 API and behavior audit",
  "",
  paste0("- Frozen package: Bioconductor 3.23 `singIST 1.0.2`; source SHA-256 `", src_sha, "`."),
  "- Audited sources: installed namespace, frozen source tarball, package manual, and `singISTvignette.Rmd`.",
  "- The installed package was not edited.",
  "",
  "## Verified public signatures",
  "",
  paste0("- ", vapply(c(
    "create_mapping_organism", "create_superpathway", "create_superpathway_input",
    "fitOptimal", "biological_link_function", "orthology_mapping",
    "singIST_treat", "singISTrecapitulations"
  ), formals_text, character(1))),
  "",
  "## Frozen operational behavior",
  "",
  "- `orthology_mapping()` accepts `orthology_cache`. A cache must contain `input_gene` (reference/Human identifiers) and `output_gene` (disease-model/Dog identifiers). When supplied, the function does not call BioMart and subsets the frozen cache independently for each fitted cell-type gene set.",
  "- Therefore `SINGIST_FROZEN_ORTHOLOG_ADAPTER` is not required for version 1.0.2. Phase11 uses the official cache argument with the Phase10 Ensembl 104 one-to-one table; no current Ensembl/BioMart mapping is generated.",
  "- `biological_link_function()` accepts both `FC_list` and `orthology_cache`. Its runtime validator requires each FC table, in exact order, to contain `p_val`, `avg_log2FC`, `pct.1`, `pct.2`, and `p_val_adj` (the manual prose says `p_value`, but implementation requires `p_val`).",
  "- When `FC_list` is absent, `diff_expressed()` uses Seurat `FindMarkers` or scran `findMarkers`; it then replaces the inferential method's fold-change with a donor-by-cell-type pseudobulk descriptive mean difference computed after aggregate-count library-size log normalization.",
  "- In the SingleCellExperiment route, pseudobulk is made with `scuttle::aggregateAcrossCells(use.assay.type='counts')`, `computeLibraryFactors()`, and `logNormCounts(pseudo_count=1, center_size_factors=FALSE)`; target-minus-base row means define `avg_log2FC`.",
  "- `singIST_treat()` retains `avg_log2FC` only where `p_val_adj <= 0.05`; all other mapped genes are assigned zero FC. This is why Track B uses Phase3 donor-level FDR in `p_val_adj` without relaxation.",
  "- `singISTrecapitulations()` uses the fitted Human predictor block as the original data, applies the official biological link to Human base-class samples, and derives superpathway, cell-type, and gene outputs using exported singIST functions.",
  "- Superpathway recapitulation is rounded by the package to two decimals; cell-type recapitulation and gene contributions are unrounded. Values are not clipped and can be negative or exceed 100%.",
  "- Gene contributions are contributions to cell-type recapitulation, not differential-expression statistics. Their sums are validated against the unrounded cell-type recapitulations at absolute tolerance 1e-8.",
  "- The official vignette recommends a wide quantile grid in principle and demonstrates `seq(0.1, 0.95, by=0.3)` for a small LOOCV example, maximum three PLS components, `folds_CV=1`, and `repetition_CV=1`. The Phase11 grid is frozen before formal fits.",
  "- `fitOptimal()` defaults to `measure='B_accuracy'`, `expected_measure_increase=0.005`, `maxiter=100`, `type='jackknife'`, and 100 permutations. Phase11 pre-registers 1,000 permutations and 100 resampling iterations as required by the project protocol.",
  "",
  "## Method gate conclusion",
  "",
  paste0("The official offline cache route and the same workflow composed from exported official functions had maximum absolute numerical difference `",
         format(max(eq$max_absolute_difference), scientific = TRUE),
         "` (required <= 1e-8). Status: **", if (all(eq$status == "PASS")) "PASS" else "FAIL", "**."),
  "",
  "The Phase10 Ensembl 104 mapping may therefore be injected directly through the supported `orthology_cache` API. No adapter and no package modification are used."
)
writeLines(aud, file.path(audit_dir, "SINGIST_API_BEHAVIOR_AUDIT.md"), useBytes = TRUE)
# Required copy in method_gate per delivery specification.
writeLines(aud, file.path(gate_dir, "SINGIST_API_BEHAVIOR_AUDIT.md"), useBytes = TRUE)

cat("\nOFFICIAL EXAMPLE VALIDATION\n")
print(validation)
cat("\nFROZEN ORTHOLOGY CACHE EQUIVALENCE\n")
print(eq)
cat("completed=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("PHASE11_METHOD_GATE_COMPLETE\n")

if (any(validation$status == "FAIL") || any(eq$status == "FAIL")) {
  stop("Phase11 singIST method gate failed")
}
