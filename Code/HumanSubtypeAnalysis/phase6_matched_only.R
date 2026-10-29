.libPaths(c(file.path((function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT")), "envs", "phase3_r_portable", "library_user"), .libPaths()))
suppressPackageStartupMessages({library(edgeR);library(fgsea);library(data.table)})
options(stringsAsFactors=FALSE)
ROOT <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
OUT <- file.path(ROOT,"WORKING","22_DOG_AUTHOR_TCELL_CROSS_SPECIES")
WORK <- file.path(OUT,"work")
BOUT <- file.path(OUT,"02_DOG_SUBTYPE_ANALYSIS")
COUT <- file.path(OUT,"03_DOG_HUMAN_CORRESPONDENCE")
MIN_DONORS <- 3L
SEED <- 20261016L

programs <- data.table(
 program=c("Mitochondrial respiration","Mitochondrial ribosome-associated quality control","Metabolic/xenobiotic program"),
 collection=c("Reactome","Reactome","Hallmark"),
 dog_pathway_id=c("R-CFA-1428517","R-CFA-9937383","HALLMARK_XENOBIOTIC_METABOLISM"),
 human_pathway_id=c("R-HSA-1428517","R-HSA-9937383","HALLMARK_XENOBIOTIC_METABOLISM"),
 pathway_name=c("Aerobic respiration and respiratory electron transport","Mitochondrial ribosome-associated quality control","Xenobiotic Metabolism"),
 source_version=c("Reactome v97","Reactome v97","MSigDB 2026.1.Hs; babelgene 22.9")
)
slugify <- function(x) {
 x <- tolower(gsub("[^A-Za-z0-9]+","_",x))
 gsub("^_+|_+$","",x)
}
full <- fread(file.path(WORK,"native_metabolic_results.tsv"))
eligible <- unique(full$dog_author_subtype)
# Strict one-to-one matched-universe comparison for frozen supported pairs only.
matches <- fread(file.path(COUT,"01_MATCHING_DECISIONS_FROZEN.csv"))
supported <- matches[matching_classification=="SUPPORTED_EXPLORATORY_MATCH" & dog_formal_pseudobulk_eligible=="YES"]
shared_members <- fread(file.path(WORK,"shared_program_memberships.tsv"))
ortho <- fread(file.path(ROOT,"10_cross_species_baseline","mapping","DOG_HUMAN_ORTHOLOG_MASTER.csv"))
ortho <- ortho[tolower(as.character(primary_one2one))=="true",.(feature_id=dog_original_feature_id,human_symbol=toupper(human_ortholog_symbol))]
human_native <- fread(file.path(ROOT,"WORKING","21_HUMAN_TCELL_SUBTYPE_METABOLIC_ANALYSIS","04_HUMAN_SUBTYPE_METABOLIC_RESULTS.csv"))
human_files <- c("CD103+ CD4"="Subtype_CD103_CD4_DE.csv","CD103+ CD8"="Subtype_CD103_CD8_DE.csv","CD103- CD4"="Subtype_CD103_CD4_DE.csv","CD4 Tcm"="Subtype_CD4_Tcm_DE.csv","KLRG1+ CD8"="Subtype_KLRG1_CD8_DE.csv","Treg"="Subtype_Treg_DE.csv")

matched_long <- list()
if(nrow(supported)) for(i in seq_len(nrow(supported))) {
 ds <- supported$dog_author_subtype[i]; hs <- supported$human_author_subtype[i]
 dde_path <- file.path(WORK,paste0("de_",slugify(ds),".csv"))
 hde_path <- file.path(ROOT,"WORKING","21_HUMAN_TCELL_SUBTYPE_METABOLIC_ANALYSIS","work","de",human_files[[hs]])
 if(!file.exists(dde_path) || !file.exists(hde_path)) next
 dde <- fread(dde_path)
 dde <- merge(dde,ortho,by="feature_id",all=FALSE)
 dde <- dde[!is.na(human_symbol)&nzchar(human_symbol)&is.finite(rank_metric)]
 dde[,abs_rank:=abs(rank_metric)]; setorder(dde,human_symbol,-abs_rank,PValue,feature_id); dde<-dde[,.SD[1],by=human_symbol]
 hde <- fread(hde_path)
 hde[,human_symbol:=toupper(gene_symbol)]
 hde <- hde[!is.na(human_symbol)&nzchar(human_symbol)&is.finite(rank_metric)]
 hde[,abs_rank:=abs(rank_metric)]; setorder(hde,human_symbol,-abs_rank,PValue,unique_gene_id); hde<-hde[,.SD[1],by=human_symbol]
 common <- intersect(dde$human_symbol,hde$human_symbol)
 dstats <- dde[human_symbol%in%common,rank_metric]; names(dstats)<-dde[human_symbol%in%common,human_symbol]
 hstats <- hde[human_symbol%in%common,rank_metric]; names(hstats)<-hde[human_symbol%in%common,human_symbol]
 paths <- split(shared_members[human_gene_symbol%in%common,human_gene_symbol],shared_members[human_gene_symbol%in%common,dog_pathway_id])
 runfg <- function(stats,species,seedadd) {
   set.seed(SEED+10000L+i*100L+seedadd)
   if(!length(paths)) return(data.table())
   z <- as.data.table(fgseaMultilevel(pathways=paths,stats=stats,minSize=15L,maxSize=500L,eps=1e-10,scoreType="std",nproc=1L))
   setnames(z,"pathway","dog_pathway_id")
   z[,leading_edge_n:=lengths(leadingEdge)]; z[,leadingEdge:=NULL]
   z[,species:=species]
   z
 }
 dfg <- runfg(dstats,"dog",1L); hfg<-runfg(hstats,"human",2L)
 for(z in list(dfg,hfg)) if(nrow(z)) {
   z[,c("dog_author_subtype","human_author_subtype","shared_tested_gene_universe"):=
       list(ds,hs,length(common))]
   matched_long[[length(matched_long)+1L]] <- z
 }
}
ml <- if(length(matched_long)) rbindlist(matched_long,fill=TRUE) else data.table()
if(nrow(ml)) {
 ml[,matched_species_family_q:=p.adjust(pval,method="BH"),by=species]
 doglong <- ml[species=="dog"]
 humanlong <- ml[species=="human"]
 setnames(doglong,c("ES","NES","pval","matched_species_family_q","size","leading_edge_n"),paste0("dog_",c("ES","NES","P","matched_family_q","shared_pathway_size","leading_edge_n")))
 setnames(humanlong,c("ES","NES","pval","matched_species_family_q","size","leading_edge_n"),paste0("human_",c("ES","NES","P","matched_family_q","shared_pathway_size","leading_edge_n")))
 keepd<-c("dog_author_subtype","human_author_subtype","dog_pathway_id","shared_tested_gene_universe","dog_ES","dog_NES","dog_P","dog_matched_family_q","dog_shared_pathway_size","dog_leading_edge_n")
 keeph<-c("dog_author_subtype","human_author_subtype","dog_pathway_id","human_ES","human_NES","human_P","human_matched_family_q","human_shared_pathway_size","human_leading_edge_n")
 wide<-merge(doglong[,..keepd],humanlong[,..keeph],by=c("dog_author_subtype","human_author_subtype","dog_pathway_id"),all=TRUE)
 wide<-merge(wide,programs,by="dog_pathway_id",all.x=TRUE)
 wide[,effect_classification:=fifelse(is.na(dog_NES)|is.na(human_NES),"NOT_ASSESSABLE",fifelse(sign(dog_NES)==sign(human_NES),"DIRECTIONALLY_CONCORDANT","DIRECTIONALLY_DISCORDANT"))]
 wide[,comparison_boundary:="Descriptive disease-effect directions in a strict shared ortholog universe; no species-by-disease interaction was tested."]
 wide<-merge(wide,full[,.(dog_author_subtype,dog_pathway_id,dog_native_NES=NES,dog_native_P=pval,dog_native_phase6_family_q=phase6_primary_family_q)],by=c("dog_author_subtype","dog_pathway_id"),all.x=TRUE)
 wide<-merge(wide,human_native[,.(human_author_subtype=subtype,human_pathway_id=pathway_id,human_native_NES=NES,human_native_P=nominal_P,human_native_18test_q=primary_18_test_BH,human_n_healthy=n_healthy,human_n_crohn=n_crohn)],by=c("human_author_subtype","human_pathway_id"),all.x=TRUE)
} else {
 wide<-data.table(dog_author_subtype="",human_author_subtype="",dog_pathway_id="",human_pathway_id="",program="",pathway_name="",effect_classification="NOT_ASSESSABLE",comparison_boundary="No supported pair with valid donor-level estimates and testable strict-ortholog programs.")
}
fwrite(wide,file.path(WORK,"matched_state_results.tsv"),sep="	",na="")
cat("ELIGIBLE_SUBTYPES=",length(eligible)," NATIVE_VALID_P=",if(nrow(full))sum(!is.na(full$pval)) else 0L," SUPPORTED_PAIRS=",nrow(supported)," MATCHED_ROWS=",nrow(wide),"\n",sep="")
