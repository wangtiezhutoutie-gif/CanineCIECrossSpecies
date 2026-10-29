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
react <- fread(file.path(ROOT,"06_pseudobulk","gsea","gene_sets","reactome_canis_familiaris_gene_sets.tsv.gz"))
hall <- fread(file.path(ROOT,"06_pseudobulk","gsea","gene_sets","hallmark_dog_ortholog_projected_membership.tsv.gz"))
members <- unique(rbindlist(list(
 react[pathway_id %in% programs$dog_pathway_id,.(dog_pathway_id=pathway_id,target_geneset_id=as.character(target_geneset_id))],
 hall[pathway_id %in% programs$dog_pathway_id,.(dog_pathway_id=pathway_id,target_geneset_id=as.character(target_geneset_id))]
)))
mapping <- fread(file.path(ROOT,"06_pseudobulk","gsea","mapping","phase4_feature_to_geneset_id.csv"))
mapping[,target_geneset_id:=as.character(target_geneset_id)]
elig <- fread(file.path(BOUT,"01_SUBTYPE_ELIGIBILITY.csv"))
eligible <- elig[formal_pseudobulk_eligible=="YES",author_subtype]

slugify <- function(x) {
 x <- tolower(gsub("[^A-Za-z0-9]+","_",x))
 gsub("^_+|_+$","",x)
}

run_model <- function(subtype, omitted="NONE", seed_offset=0L, save_de=FALSE) {
 slug <- slugify(subtype)
 ctab <- fread(file.path(WORK,paste0("counts_",slug,".tsv")))
 meta <- fread(file.path(WORK,paste0("samples_",slug,".tsv")))
 if (omitted!="NONE") meta <- meta[donor_id!=omitted]
 hn <- meta[condition=="Healthy",uniqueN(donor_id)]
 cn <- meta[condition=="CIE",uniqueN(donor_id)]
 if (hn<MIN_DONORS || cn<MIN_DONORS)
   return(list(ok=FALSE,reason="NON_ESTIMABLE_DONOR_GATE",healthy=hn,cie=cn))
 ids <- meta$sample_id
 ann <- ctab[,.(feature_id,gene_symbol)]
 counts <- as.matrix(ctab[,..ids]); storage.mode(counts)<-"integer"; rownames(counts)<-ann$feature_id
 meta[,condition:=factor(condition,levels=c("Healthy","CIE"))]
 design <- model.matrix(~condition,data=meta)
 y <- DGEList(counts=counts,genes=ann)
 keep <- filterByExpr(y,design=design)
 if (sum(keep)<100L) return(list(ok=FALSE,reason="NON_ESTIMABLE_FILTERED_GENE_COUNT",healthy=hn,cie=cn))
 y <- y[keep,,keep.lib.sizes=FALSE]
 y <- calcNormFactors(y,method="TMM")
 y <- estimateDisp(y,design,robust=TRUE)
 fit <- glmQLFit(y,design,robust=TRUE)
 qlf <- glmQLFTest(fit,coef="conditionCIE")
 tt <- as.data.table(topTags(qlf,n=Inf,sort.by="none")$table)
 tt[,rank_metric:=sign(logFC)*sqrt(F)]
 if (save_de) fwrite(tt,file.path(WORK,paste0("de_",slug,".csv")))
 rr <- merge(tt,mapping[,.(feature_id,target_geneset_id,target_symbol)],by="feature_id",all.x=FALSE)
 rr <- rr[!is.na(target_geneset_id)&nzchar(target_geneset_id)&is.finite(rank_metric)]
 rr[,abs_rank:=abs(rank_metric)]
 setorder(rr,target_geneset_id,-abs_rank,PValue,feature_id)
 rr <- rr[,.SD[1],by=target_geneset_id]
 stats <- rr$rank_metric; names(stats)<-rr$target_geneset_id
 paths <- split(members$target_geneset_id,members$dog_pathway_id)
 set.seed(SEED+seed_offset)
 fg <- as.data.table(fgseaMultilevel(pathways=paths,stats=stats,minSize=15L,maxSize=500L,eps=1e-10,scoreType="std",nproc=1L))
 setnames(fg,"pathway","dog_pathway_id")
 fg[,leading_edge_n:=lengths(leadingEdge)]; fg[,leadingEdge:=NULL]
 z <- merge(programs,fg[,.(dog_pathway_id,ES,NES,pval,size,leading_edge_n,log2err)],by="dog_pathway_id",all.x=TRUE,sort=FALSE)
 z[,c("dog_author_subtype","omitted_donor","healthy_donor_n","cie_donor_n","tested_genes","mapped_ranked_genes","model_status","direction"):=
     list(subtype,omitted,hn,cn,sum(keep),nrow(rr),ifelse(is.na(NES),"PATHWAY_NOT_TESTABLE","ESTIMABLE"),ifelse(is.na(NES),"NOT_TESTABLE",ifelse(NES>0,"CIE_HIGHER","CIE_LOWER")))]
 list(ok=TRUE,results=z,healthy=hn,cie=cn)
}

full_list <- list(); lodo_list <- list()
for (i in seq_along(eligible)) {
 st <- eligible[[i]]
 full <- tryCatch(run_model(st,seed_offset=i*100L,save_de=TRUE),error=function(e) list(ok=FALSE,reason=paste0("MODEL_ERROR:",conditionMessage(e)),healthy=NA_integer_,cie=NA_integer_))
 if (isTRUE(full$ok)) {
   full_list[[length(full_list)+1L]] <- full$results
 } else {
   z <- copy(programs)
   z[,c("dog_author_subtype","omitted_donor","healthy_donor_n","cie_donor_n","tested_genes","mapped_ranked_genes","ES","NES","pval","size","leading_edge_n","log2err","model_status","direction"):=
       list(st,"NONE",full$healthy,full$cie,NA_integer_,NA_integer_,NA_real_,NA_real_,NA_real_,NA_integer_,NA_integer_,NA_real_,full$reason,"NOT_TESTABLE")]
   full_list[[length(full_list)+1L]] <- z
 }
 meta <- fread(file.path(WORK,paste0("samples_",slugify(st),".tsv")))
 for (j in seq_len(nrow(meta))) {
   d <- meta$donor_id[j]
   res <- tryCatch(run_model(st,d,i*100L+j),error=function(e) list(ok=FALSE,reason=paste0("MODEL_ERROR:",conditionMessage(e)),healthy=NA_integer_,cie=NA_integer_))
   if (isTRUE(res$ok)) {
     x <- res$results
   } else {
     x <- copy(programs)
     x[,c("dog_author_subtype","omitted_donor","healthy_donor_n","cie_donor_n","tested_genes","mapped_ranked_genes","ES","NES","pval","size","leading_edge_n","log2err","model_status","direction"):=
         list(st,d,res$healthy,res$cie,NA_integer_,NA_integer_,NA_real_,NA_real_,NA_real_,NA_integer_,NA_integer_,NA_real_,res$reason,"NOT_ESTIMABLE")]
   }
   x[,omitted_condition:=meta[donor_id==d,condition][1]]
   lodo_list[[length(lodo_list)+1L]] <- x
 }
}
full <- if(length(full_list)) rbindlist(full_list,fill=TRUE) else data.table()
lodo <- if(length(lodo_list)) rbindlist(lodo_list,fill=TRUE) else data.table()
if(nrow(full)) {
 full[,phase6_primary_family_q:=p.adjust(pval,method="BH")]
 full[,planned_test_family:=paste0(length(eligible)," eligible author subtypes x 3 predefined programs = ",length(eligible)*3L," planned tests")]
 full[,valid_p_values_adjusted:=sum(!is.na(pval))]
}
if(nrow(lodo)) {
 base <- full[,.(dog_author_subtype,program,full_NES=NES)]
 lodo <- merge(lodo,base,by=c("dog_author_subtype","program"),all.x=TRUE)
 lodo[,direction_match_full:=fifelse(is.na(NES),"NOT_ESTIMABLE",fifelse(sign(NES)==sign(full_NES),"YES","NO"))]
}
fwrite(full,file.path(WORK,"native_metabolic_results.tsv"),sep="	",na="")
fwrite(lodo,file.path(WORK,"native_lodo_results.tsv"),sep="	",na="")

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
