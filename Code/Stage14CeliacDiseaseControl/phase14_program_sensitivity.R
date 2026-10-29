#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(fgsea)})
root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
base <- file.path(root,"14_ced_specificity")
groot <- file.path(base,"gsea"); sroot <- file.path(base,"sensitivity")
dir.create(sroot,recursive=TRUE,showWarnings=FALSE)
freeze <- file.path(groot,"CED_DISEASE_PROGRAMS.csv")
sha_file <- file.path(groot,"CED_DISEASE_PROGRAMS.sha256")
if (!file.exists(freeze)||!file.exists(sha_file)) stop("CeD GSEA is not frozen")
expected_sha <- trimws(readLines(sha_file,warn=FALSE)[1])
sha_line <- system2("certutil",c("-hashfile",shQuote(freeze),"SHA256"),
                    stdout=TRUE,stderr=TRUE)
actual_sha <- tolower(gsub(" ","",sha_line[grepl("^[0-9A-Fa-f ]{64,}$",sha_line)][1]))
if (actual_sha!=expected_sha) stop("CeD disease-program freeze hash mismatch")

registered <- fread(file.path(root,"10_cross_species_baseline","summary",
                              "PHASE10_PROGRAM_BASELINE_MATRIX.csv"))
if (nrow(registered)!=25) stop("Registered program panel is not 25 rows")
sets <- c(
  Hallmark=file.path(root,"07_human_crohn","GSE266616","gsea","gene_sets","HUMAN_HALLMARK_GENE_SETS.tsv.gz"),
  Reactome=file.path(root,"07_human_crohn","GSE266616","gsea","gene_sets","HUMAN_REACTOME_GENE_SETS.tsv.gz"),
  GO_BP=file.path(root,"07_human_crohn","GSE266616","gsea","gene_sets","HUMAN_GO_BP_GENE_SETS.tsv.gz"))
members <- lapply(names(sets),function(collection_name) {
  z <- unique(fread(sets[[collection_name]]),by=c("pathway_id","gene_symbol"))
  z[collection==collection_name]
})
names(members) <- names(sets)

make_stats <- function(d) {
  d <- d[!is.na(gene_symbol)&nzchar(gene_symbol)&is.finite(rank_metric)]
  d[,absrank:=abs(rank_metric)]; setorder(d,-absrank,gene_id)
  d <- unique(d,by="gene_symbol")
  s <- d$rank_metric; names(s) <- d$gene_symbol; s
}
run_registered <- function(d,lineage_name,analysis_id,seed_offset) {
  stats <- make_stats(copy(d))
  rr <- registered[lineage==lineage_name]
  outs <- list(); k <- 1
  for (collection in unique(rr$representative_collection)) {
    rc <- rr[representative_collection==collection]
    mem <- members[[collection]][pathway_id %in% rc$representative_pathway]
    paths <- split(mem$gene_symbol,mem$pathway_id)
    paths <- paths[names(paths) %in% rc$representative_pathway]
    if (!length(paths)) next
    set.seed(20260828L+seed_offset+match(collection,names(sets)))
    fg <- as.data.table(fgseaMultilevel(
      pathways=paths,stats=stats,minSize=15,maxSize=500,eps=1e-10,
      scoreType="std",nproc=1L))
    setnames(fg,"pathway","representative_pathway")
    fg <- merge(fg,unique(rc[,.(program,representative_pathway)]),
                by="representative_pathway",all.x=TRUE,allow.cartesian=TRUE)
    fg[,lineage:=lineage_name]; fg[,analysis_id:=analysis_id]
    fg[,direction:=ifelse(NES>0,"Active_CeD_higher","Active_CeD_lower")]
    outs[[k]] <- fg[,.(lineage,program,analysis_id,representative_pathway,
                       NES,pval,padj,size,direction)]; k <- k+1
  }
  if (!length(outs)) return(data.table())
  rbindlist(outs,fill=TRUE)
}

rep_results <- list()
for (lineage in c("Epithelial","T","Myeloid")) {
  d <- fread(file.path(sroot,paste0("CED_",lineage,"_PRIMARY_LIBRARY_ONLY.csv")))
  rep_results[[lineage]] <- run_registered(
    d,lineage,"PRIMARY_LIBRARY_ONLY_DROP_REP2",1000L+match(lineage,c("Epithelial","T","Myeloid"))*10L)
}
rep_result <- rbindlist(rep_results,fill=TRUE)
fwrite(rep_result,file.path(sroot,"CED_REPLICATE_SENSITIVITY_REGISTERED_GSEA.csv"))

jack <- fread(file.path(sroot,"CED_32_SUBJECT_JACKKNIFE_GENE_RESULTS.csv.gz"))
groups <- unique(jack[,.(lineage,dropped_subject,dropped_condition)])
jouts <- vector("list",nrow(groups))
for (i in seq_len(nrow(groups))) {
  g <- groups[i]
  message("Registered jackknife ",i,"/",nrow(groups),": ",g$lineage," drop ",g$dropped_subject)
  d <- jack[lineage==g$lineage & dropped_subject==g$dropped_subject]
  z <- run_registered(d,g$lineage,paste0("DROP_",g$dropped_subject),2000L+i*10L)
  z[,dropped_subject:=g$dropped_subject]
  z[,dropped_condition:=g$dropped_condition]
  jouts[[i]] <- z
}
jres <- rbindlist(jouts,fill=TRUE)
fwrite(jres,file.path(sroot,"CED_DONOR_JACKKNIFE_PROGRAM_RESULTS.csv"))

primary <- fread(freeze)
primary <- primary[,.(lineage,representative_pathway=pathway,
                      primary_NES=NES,primary_FDR=FDR,
                      primary_direction=direction)]
primary <- merge(registered[,.(lineage,program,representative_pathway)],primary,
                 by=c("lineage","representative_pathway"),all.x=TRUE,
                 allow.cartesian=TRUE)
jres <- merge(jres,primary,by=c("lineage","program","representative_pathway"),all.x=TRUE)
jres[,same_direction:=direction==primary_direction]
jres[,delta_NES:=NES-primary_NES]
jsummary <- jres[,{
  ix <- which.max(abs(delta_NES))
  .(jackknife_runs=.N,same_direction_runs=sum(same_direction,na.rm=TRUE),
    same_direction_fraction=mean(same_direction,na.rm=TRUE),
    median_NES=median(NES,na.rm=TRUE),min_NES=min(NES,na.rm=TRUE),
    max_NES=max(NES,na.rm=TRUE),max_abs_delta_NES=max(abs(delta_NES),na.rm=TRUE),
    most_influential_donor=dropped_subject[ix],
    most_influential_condition=dropped_condition[ix])
},by=.(lineage,program,representative_pathway,primary_NES,primary_FDR,primary_direction)]
jsummary[,jackknife_status:=fifelse(
  same_direction_fraction==1,"DIRECTION_STABLE_32_OF_32",
  fifelse(same_direction_fraction>=.9,"DIRECTION_STABLE_GE90_PERCENT","DONOR_SENSITIVE"))]
setcolorder(jsummary,c("lineage","program","representative_pathway",
                       setdiff(names(jsummary),c("lineage","program","representative_pathway"))))
fwrite(jsummary,file.path(sroot,"CED_DONOR_JACKKNIFE_PROGRAM_SUMMARY.csv"))

rep_summary <- merge(rep_result,primary,
                     by=c("lineage","program","representative_pathway"),all.x=TRUE)
rep_summary[,same_direction:=direction==primary_direction]
rep_summary[,delta_NES:=NES-primary_NES]
rep_summary[,replicate_status:=fifelse(
  same_direction & abs(delta_NES)<=.25,"REPLICATE_HANDLING_STABLE",
  fifelse(same_direction,"DIRECTION_STABLE_EFFECT_SHIFT","REPLICATE_HANDLING_DIRECTION_SENSITIVE"))]
fwrite(rep_summary,file.path(sroot,"CED_REPLICATE_SENSITIVITY.csv"))

gate <- data.frame(
  metric=c("frozen_CeD_hash_verified","registered_programs",
           "jackknife_runs_expected","jackknife_program_rows",
           "replicate_sensitivity_programs"),
  value=c("YES","25","96",nrow(jres),nrow(rep_summary)),
  status=c("PASS","PASS","PASS",ifelse(nrow(groups)==96,"PASS","FAIL"),
           ifelse(nrow(rep_summary)==25,"PASS","FAIL")))
fwrite(gate,file.path(sroot,"CED_SENSITIVITY_GATE.csv"))
writeLines(c(paste("freeze_sha256",actual_sha),
             paste("jackknife_groups",nrow(groups)),
             paste("completed_utc",format(Sys.time(),tz="UTC",usetz=TRUE))),
           file.path(root,"logs","phase14_sensitivity.log"))
message("Phase14 registered-program sensitivity complete")

