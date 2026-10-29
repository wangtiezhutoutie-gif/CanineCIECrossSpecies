#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

options(stringsAsFactors = FALSE)
root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
base <- file.path(root, "10_cross_species_baseline")
rank_dir <- file.path(base, "pathway_effects", "paired_ranks")
raw_dir <- file.path(base, "pathway_effects", "gsea_raw")
out_dir <- file.path(base, "pathway_effects")
audit_dir <- file.path(base, "audit")
log_path <- file.path(root, "logs", "phase10_gsea.log")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)

lineages <- c("Epithelial", "T", "Myeloid")
species_order <- c("Dog", "Human")
gene_set_dir <- file.path(root, "07_human_crohn", "GSE266616", "gsea", "gene_sets")
gene_sets <- list(
  Hallmark = file.path(gene_set_dir, "HUMAN_HALLMARK_GENE_SETS.tsv.gz"),
  Reactome = file.path(gene_set_dir, "HUMAN_REACTOME_GENE_SETS.tsv.gz"),
  GO_BP = file.path(gene_set_dir, "HUMAN_GO_BP_GENE_SETS.tsv.gz")
)
min_size <- 15L
max_size <- 500L
fgsea_eps <- 1e-10
base_seed <- 20260826L

atomic_fwrite <- function(frame, path, compress = "none") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  suffix <- if (compress == "gzip") ".csv.gz" else ".csv"
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = suffix)
  fwrite(frame, tmp, compress = compress)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Atomic rename failed: ", path)
}

read_rank <- function(lineage, species) {
  path <- file.path(rank_dir, paste0(lineage, "_", species, "_common_one2one_rank.csv"))
  rank <- fread(path)
  required <- c("human_gene_id", "human_symbol", "rank_metric", "logFC")
  if (!all(required %in% names(rank))) stop("Rank schema failure: ", path)
  if (anyDuplicated(rank$human_gene_id)) stop("Duplicate Human gene identity: ", path)
  rank <- rank[is.finite(rank_metric) & !is.na(human_gene_id) & nzchar(human_gene_id)]
  setorder(rank, human_gene_id)
  stats <- rank$rank_metric
  names(stats) <- rank$human_gene_id
  list(path = path, rank = rank, stats = stats)
}

load_membership <- function(collection, rank) {
  source <- fread(gene_sets[[collection]])
  required <- c("collection", "pathway_id", "pathway_name", "gene_symbol", "source_version", "source_species")
  if (!all(required %in% names(source))) stop("Gene-set schema failure: ", gene_sets[[collection]])
  if (any(source$source_species != "Homo sapiens")) stop("Non-Human membership source: ", collection)
  symbol_map <- unique(rank[, .(human_gene_id, gene_symbol = human_symbol)])
  membership <- merge(
    unique(source, by = c("pathway_id", "gene_symbol")),
    symbol_map,
    by = "gene_symbol",
    allow.cartesian = TRUE
  )
  membership <- unique(membership, by = c("pathway_id", "human_gene_id"))
  membership
}

run_one <- function(lineage, species, collection) {
  rank_info <- read_rank(lineage, species)
  membership <- load_membership(collection, rank_info$rank)
  pathway_info <- unique(
    membership[, .(pathway_id, pathway_name, source_version, source_species)],
    by = "pathway_id"
  )
  pathways <- split(membership$human_gene_id, membership$pathway_id)
  set.seed(
    base_seed +
      match(lineage, lineages) * 100L +
      match(collection, names(gene_sets)) * 10L +
      match(species, species_order)
  )
  fg <- as.data.table(fgseaMultilevel(
    pathways = pathways,
    stats = rank_info$stats,
    minSize = min_size,
    maxSize = max_size,
    eps = fgsea_eps,
    scoreType = "std",
    nproc = 1L
  ))
  setnames(fg, "pathway", "pathway_id")
  fg[, padj := p.adjust(pval, method = "BH")]
  result <- merge(fg, pathway_info, by = "pathway_id", all.x = TRUE, sort = FALSE)
  result[, `:=`(
    lineage = lineage,
    species = species,
    collection = collection,
    leading_edge_n = lengths(leadingEdge),
    leading_edge = vapply(
      leadingEdge,
      function(x) paste(as.character(x), collapse = ";"),
      character(1)
    ),
    direction = ifelse(
      NES > 0,
      ifelse(species == "Dog", "CIE_higher", "ActiveCrohn_higher"),
      ifelse(species == "Dog", "CIE_lower", "ActiveCrohn_lower")
    ),
    shared_universe_n = length(rank_info$stats),
    rank_key = "Human Ensembl gene ID",
    gene_set_key = "Human Ensembl gene ID reconstructed from frozen Phase7 Human symbol membership",
    rank_metric_definition = "sign(logFC)*sqrt(F)"
  )]
  result[, leadingEdge := NULL]
  setcolorder(
    result,
    c(
      "lineage", "species", "collection", "pathway_id", "pathway_name",
      "ES", "NES", "pval", "padj", "size", "leading_edge_n", "leading_edge",
      "direction", "shared_universe_n", "source_version", "source_species",
      "rank_key", "gene_set_key", "rank_metric_definition", "log2err"
    )
  )
  result[, abs_nes_sort := abs(NES)]
  setorder(result, padj, -abs_nes_sort, pathway_id)
  result[, abs_nes_sort := NULL]
  out <- file.path(raw_dir, paste0(lineage, "_", species, "_", collection, "_common_one2one_GSEA.csv"))
  atomic_fwrite(result, out)
  message(
    lineage, " / ", species, " / ", collection,
    ": tested=", nrow(result), "; FDR10=", sum(result$padj < 0.10)
  )
  result
}

jaccard <- function(left, right) {
  a <- unique(strsplit(ifelse(is.na(left), "", left), ";", fixed = TRUE)[[1]])
  b <- unique(strsplit(ifelse(is.na(right), "", right), ";", fixed = TRUE)[[1]])
  a <- a[nzchar(a)]
  b <- b[nzchar(b)]
  union_n <- length(union(a, b))
  if (!union_n) return(NA_real_)
  length(intersect(a, b)) / union_n
}

main <- function() {
  sink(log_path, split = TRUE)
  on.exit(sink(), add = TRUE)
  cat("Phase10 paired common-one2one GSEA\n")
  cat("Start UTC:", format(Sys.time(), tz = "UTC"), "\n")
  all_results <- list()
  for (lineage in lineages) {
    dog_rank <- read_rank(lineage, "Dog")$rank
    human_rank <- read_rank(lineage, "Human")$rank
    if (!identical(dog_rank$human_gene_id, human_rank$human_gene_id)) {
      stop("Paired rank universe mismatch: ", lineage)
    }
    for (collection in names(gene_sets)) {
      for (species in species_order) {
        key <- paste(lineage, species, collection, sep = "::")
        all_results[[key]] <- run_one(lineage, species, collection)
      }
    }
  }

  comparison_parts <- list()
  summary_parts <- list()
  for (lineage in lineages) {
    for (collection in names(gene_sets)) {
      dog <- copy(all_results[[paste(lineage, "Dog", collection, sep = "::")]])
      human <- copy(all_results[[paste(lineage, "Human", collection, sep = "::")]])
      if (!setequal(dog$pathway_id, human$pathway_id)) {
        stop("Dog/Human tested pathway mismatch: ", lineage, " / ", collection)
      }
      dog <- dog[, .(
        lineage, collection, pathway_id, pathway_name,
        dog_NES = NES, dog_FDR = padj, dog_pvalue = pval,
        dog_size = size, dog_leading_edge = leading_edge,
        source_version, source_species, shared_universe_n
      )]
      human <- human[, .(
        lineage, collection, pathway_id,
        human_NES = NES, human_FDR = padj, human_pvalue = pval,
        human_size = size, human_leading_edge = leading_edge
      )]
      paired <- merge(
        dog, human,
        by = c("lineage", "collection", "pathway_id"),
        all = FALSE
      )
      if (anyDuplicated(paired[, .(lineage, collection, pathway_id)])) {
        stop("Paired pathway identity duplication: ", lineage, " / ", collection)
      }
      if (any(paired$dog_size != paired$human_size)) {
        stop("Paired pathway membership size mismatch: ", lineage, " / ", collection)
      }
      paired[, `:=`(
        same_direction = sign(dog_NES) == sign(human_NES),
        delta_NES = human_NES - dog_NES,
        leading_edge_jaccard = vapply(
          seq_len(.N),
          function(i) jaccard(dog_leading_edge[[i]], human_leading_edge[[i]]),
          numeric(1)
        ),
        comparison_scope = "Traditional Ensembl104 ortholog_one2one common-universe baseline",
        region_risk = "PROVENANCE_PENDING_POST_BASELINE",
        phase10_priority = "PROVENANCE_PENDING_POST_BASELINE",
        ibdverse_replication_if_available = "PROVENANCE_PENDING_POST_BASELINE",
        notes = "Dog and Human BH adjusted separately within lineage and collection; delta_NES=Human-Dog"
      )]
      setorder(paired, pathway_id)
      comparison_parts[[paste(lineage, collection, sep = "::")]] <- paired
      summary_parts[[paste(lineage, collection, sep = "::")]] <- data.table(
        lineage = lineage,
        collection = collection,
        common_tested_pathways = nrow(paired),
        nes_spearman_rho = cor(
          paired$dog_NES, paired$human_NES,
          method = "spearman", use = "complete.obs"
        ),
        nes_pearson_r = cor(
          paired$dog_NES, paired$human_NES,
          method = "pearson", use = "complete.obs"
        ),
        direction_concordance = mean(paired$same_direction, na.rm = TRUE),
        both_FDR10_same_direction = sum(
          paired$dog_FDR < 0.10 & paired$human_FDR < 0.10 & paired$same_direction,
          na.rm = TRUE
        ),
        both_FDR10_opposite = sum(
          paired$dog_FDR < 0.10 & paired$human_FDR < 0.10 & !paired$same_direction,
          na.rm = TRUE
        ),
        dog_FDR10_human_weak_same = sum(
          paired$dog_FDR < 0.10 & paired$human_FDR >= 0.10 & paired$same_direction,
          na.rm = TRUE
        ),
        human_FDR10_dog_weak_same = sum(
          paired$human_FDR < 0.10 & paired$dog_FDR >= 0.10 & paired$same_direction,
          na.rm = TRUE
        ),
        leading_edge_median_jaccard = median(
          paired$leading_edge_jaccard,
          na.rm = TRUE
        ),
        fdr_correction = "BH separately by species, lineage, and collection",
        notes = "Same common Human Ensembl-ID universe and memberships for both species"
      )
    }
  }
  comparison <- rbindlist(comparison_parts, fill = TRUE)
  summary <- rbindlist(summary_parts, fill = TRUE)
  atomic_fwrite(
    comparison,
    file.path(out_dir, "PATHWAY_CROSS_SPECIES_BASELINE_PRE_PROVENANCE.csv.gz"),
    compress = "gzip"
  )
  atomic_fwrite(summary, file.path(out_dir, "PATHWAY_CONCORDANCE_SUMMARY.csv"))
  software <- data.table(
    software = c("R", "fgsea", "data.table"),
    version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      as.character(packageVersion("fgsea")),
      as.character(packageVersion("data.table"))
    ),
    role = c("runtime", "fgseaMultilevel", "tabular I/O")
  )
  parameters <- data.table(
    parameter = c(
      "algorithm", "minSize", "maxSize", "eps", "scoreType", "nproc",
      "rank_metric", "universe", "gene_sets", "multiple_testing"
    ),
    value = c(
      "fgseaMultilevel", min_size, max_size, fgsea_eps, "std", 1L,
      "sign(logFC)*sqrt(F)", "Lineage-specific shared Ensembl104 one2one Human IDs",
      "Frozen Phase7 Human Hallmark/Reactome/GO BP definitions",
      "BH separately by species, lineage, and collection"
    )
  )
  atomic_fwrite(software, file.path(audit_dir, "PHASE10_GSEA_SOFTWARE_VERSIONS.csv"))
  atomic_fwrite(parameters, file.path(audit_dir, "PHASE10_GSEA_PARAMETERS.csv"))
  cat("Completed UTC:", format(Sys.time(), tz = "UTC"), "\n")
  cat("Paired pathway rows:", nrow(comparison), "\n")
}

main()
