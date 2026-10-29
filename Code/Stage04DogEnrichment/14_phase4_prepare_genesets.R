#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(babelgene)
  library(msigdbr)
})

options(stringsAsFactors = FALSE)
project_root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
gsea_root <- file.path(project_root, "06_pseudobulk", "gsea")
raw_dir <- file.path(gsea_root, "gene_sets", "raw")
mapping_dir <- file.path(gsea_root, "mapping")
dir.create(mapping_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(gsea_root, "gene_sets"), recursive = TRUE, showWarnings = FALSE)

gmt_path <- file.path(raw_dir, "h.all.v2026.1.Hs.symbols.gmt")
if (!file.exists(gmt_path)) stop("Missing frozen Hallmark GMT: ", gmt_path)

supported <- msigdbr_species()
dog_rows <- supported[
  grepl("Canis", supported$species_name, ignore.case = TRUE) |
    grepl("dog", supported$species_common_name, ignore.case = TRUE),
  , drop = FALSE
]
if (nrow(dog_rows) == 0L) stop("Current msigdbr species table does not list dog")

gmt_lines <- readLines(gmt_path, warn = FALSE, encoding = "UTF-8")
gmt_parts <- strsplit(gmt_lines, "\t", fixed = TRUE)
gmt <- rbindlist(lapply(gmt_parts, function(x) {
  if (length(x) < 3L) stop("Malformed Hallmark GMT line")
  data.table(pathway_id = x[[1L]], pathway_description = x[[2L]], human_symbol = x[-c(1L, 2L)])
}))
gmt <- unique(gmt[nzchar(human_symbol)])

message("Projecting ", uniqueN(gmt$human_symbol), " human Hallmark symbols to dog with babelgene...")
orth <- as.data.table(babelgene::orthologs(unique(gmt$human_symbol), species = "dog"))
setnames(
  orth,
  c("symbol", "entrez", "ensembl"),
  c("dog_symbol", "dog_entrez", "dog_ensembl")
)
orth[, dog_entrez := as.character(dog_entrez)]
orth <- orth[!is.na(dog_entrez) & nzchar(dog_entrez)]
orth <- unique(orth)

hallmark <- merge(gmt, orth, by = "human_symbol", allow.cartesian = TRUE)
hallmark[, `:=`(
  collection = "Hallmark",
  pathway_name = gsub("^HALLMARK_", "", pathway_id),
  target_geneset_id = dog_entrez,
  target_symbol = dog_symbol,
  mapping_type = "ORTHOLOG_PROJECTED_HALLMARK_BABELGENE",
  source_version = "MSigDB 2026.1.Hs; babelgene 22.9",
  source_species = "Homo sapiens",
  target_species = "Canis lupus familiaris"
)]
hallmark <- unique(hallmark[, .(
  collection, pathway_id, pathway_name, target_geneset_id, target_symbol,
  mapping_type, source_version, source_species, target_species,
  human_symbol, human_entrez, human_ensembl, dog_ensembl, support, support_n
)])

fwrite(
  hallmark,
  file.path(gsea_root, "gene_sets", "hallmark_dog_ortholog_projected.tsv.gz"),
  sep = "\t"
)
fwrite(
  orth,
  file.path(mapping_dir, "hallmark_babelgene_ortholog_map.csv")
)
fwrite(
  as.data.table(dog_rows),
  file.path(mapping_dir, "msigdbr_supported_dog_species.csv")
)

message(
  "Hallmark complete: ", uniqueN(hallmark$pathway_id), " pathways; ",
  uniqueN(hallmark$human_symbol), " mapped human symbols; ",
  uniqueN(hallmark$target_geneset_id), " dog Entrez genes."
)
