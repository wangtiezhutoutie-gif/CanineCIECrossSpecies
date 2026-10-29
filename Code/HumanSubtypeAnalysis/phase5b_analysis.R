suppressPackageStartupMessages({library(edgeR);library(data.table);library(fgsea)})
options(stringsAsFactors=FALSE)
root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT")); out <- file.path(root,"WORKING","21_HUMAN_TCELL_SUBTYPE_METABOLIC_ANALYSIS"); work <- file.path(out,"work"); gene_set_dir <- file.path(root,"07_human_crohn","GSE266616","gsea","gene_sets")
dir.create(file.path(work,"de"),recursive=TRUE,showWarnings=FALSE); dir.create(file.path(work,"ranks"),recursive=TRUE,showWarnings=FALSE)
targets <- data.table(pathway_id=c("R-HSA-1428517","R-HSA-9937383","HALLMARK_XENOBIOTIC_METABOLISM"),collection=c("Reactome","Reactome","Hallmark"),expected_name=c("Aerobic respiration and respiratory electron transport","Mitochondrial ribosome-associated quality control","Xenobiotic Metabolism"))
membership_files <- c(Reactome=file.path(gene_set_dir,"HUMAN_REACTOME_GENE_SETS.tsv.gz"),Hallmark=file.path(gene_set_dir,"HUMAN_HALLMARK_GENE_SETS.tsv.gz"))
members <- rbindlist(lapply(names(membership_files),function(coll){dt<-fread(membership_files[[coll]]);dt[collection==coll & pathway_id %in% targets$pathway_id]}),fill=TRUE)
members <- unique(members,by=c("pathway_id","gene_symbol")); info <- unique(members[,.(pathway_id,pathway_name,collection,source_version,source_species)],by="pathway_id")
if(nrow(info)!=3L || !setequal(info$pathway_id,targets$pathway_id)) stop("Frozen target-program recovery failed")
check <- merge(targets,info,by=c("pathway_id","collection"),all=TRUE); if(any(check$expected_name!=check$pathway_name)) stop("Frozen pathway-name identity mismatch")
fwrite(members,file.path(work,"frozen_target_memberships.tsv"),sep="\t"); fwrite(info,file.path(work,"frozen_target_programs.tsv"),sep="\t")
read_matrix <- function(path){
  x<-fread(path); required<-c("unique_gene_id","gene_symbol")
  if(!all(required %in% names(x))) stop("Matrix schema failure: ",path)
  if(anyDuplicated(x$unique_gene_id)) stop("Duplicated gene IDs: ",path)
  counts<-as.matrix(x[,setdiff(names(x),required),with=FALSE]);storage.mode(counts)<-"integer";rownames(counts)<-x$unique_gene_id
  if(anyNA(counts)||any(counts<0L)) stop("Invalid count matrix: ",path)
  list(annotation=as.data.frame(x[,..required]),counts=counts)
}
fit_analysis <- function(analysis,matrix_path,meta,seed_offset=0L,write_de=TRUE){
  m<-read_matrix(matrix_path)
  if(!setequal(colnames(m$counts),meta$sample_column)) stop("Count/meta mismatch: ",analysis)
  meta<-meta[match(colnames(m$counts),meta$sample_column),,drop=FALSE]; condition<-factor(meta$condition,levels=c("Healthy","ActiveCrohn"))
  if(sum(condition=="Healthy")<5L||sum(condition=="ActiveCrohn")<5L) stop("Frozen donor gate failed: ",analysis)
  design<-model.matrix(~condition); if(qr(design)$rank!=ncol(design)) stop("Rank-deficient design: ",analysis)
  y0<-DGEList(counts=m$counts);keep<-filterByExpr(y0,design=design);if(sum(keep)<1000L) stop("Too few tested genes after Phase7 filterByExpr: ",analysis)
  y<-y0[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y,method="TMM");y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE);qlf<-glmQLFTest(fit,coef="conditionActiveCrohn")
  raw<-topTags(qlf,n=Inf,sort.by="none")$table;anno<-m$annotation[match(rownames(raw),m$annotation$unique_gene_id),,drop=FALSE]
  de<-data.table(unique_gene_id=rownames(raw),gene_symbol=anno$gene_symbol,logFC=raw$logFC,logCPM=raw$logCPM,statistic_F=raw$F,PValue=raw$PValue,FDR=raw$FDR,rank_metric=sign(raw$logFC)*sqrt(pmax(raw$F,0)))
  rank<-de[!is.na(gene_symbol)&nzchar(gene_symbol)&is.finite(rank_metric)];rank[,abs_rank:=abs(rank_metric)];setorder(rank,-abs_rank,PValue,unique_gene_id);rank<-unique(rank,by="gene_symbol");rank[,abs_rank:=NULL];setorder(rank,-rank_metric,gene_symbol)
  stats<-rank$rank_metric;names(stats)<-rank$gene_symbol;pathways<-split(members$gene_symbol,members$pathway_id);set.seed(20261015L+seed_offset)
  fg<-as.data.table(fgseaMultilevel(pathways=pathways,stats=stats,minSize=15L,maxSize=500L,eps=1e-10,scoreType="std",nproc=1L));setnames(fg,"pathway","pathway_id");fg[,adjusted_P_within_three:=p.adjust(pval,method="BH")]
  result<-merge(info,fg,by="pathway_id",all.x=TRUE,sort=FALSE)
  result$analysis<-analysis;result$n_healthy<-sum(condition=="Healthy");result$n_crohn<-sum(condition=="ActiveCrohn");result$total_cells<-sum(meta$n_cells);result$tested_gene_ids<-sum(keep);result$mapped_ranked_gene_symbols<-length(stats);result$direction<-fifelse(is.na(result$NES),"NOT_ASSESSABLE",fifelse(result$NES>0,"ActiveCrohn_higher",fifelse(result$NES<0,"ActiveCrohn_lower","No_change")));result$testability<-fifelse(is.na(result$NES),"NOT_TESTABLE","TESTABLE");result$testing_family<-"THREE_PREDEFINED_PROGRAMS_WITHIN_ANALYSIS"
  result[,leading_edge:=vapply(leadingEdge,function(x)paste(as.character(x),collapse=";"),character(1))];result[,leadingEdge:=NULL]
  safe<-gsub("[^A-Za-z0-9]+","_",analysis);if(write_de){fwrite(de,file.path(work,"de",paste0(safe,"_DE.csv")));fwrite(rank,file.path(work,"ranks",paste0(safe,"_ranked.csv")))}
  model<-data.table(analysis=analysis,n_healthy=sum(condition=="Healthy"),n_crohn=sum(condition=="ActiveCrohn"),tested_gene_ids=sum(keep),mapped_ranked_gene_symbols=length(stats),design_rank=qr(design)$rank,design_columns=ncol(design),common_dispersion=y$common.dispersion,robust_dispersion="YES",ql_fit="YES",normalization="TMM",filter="filterByExpr(design)",rank_statistic="sign(logFC)*sqrt(F)",contrast="ActiveCrohn-minus-Healthy")
  list(result=result,model=model)
}
qc<-fread(file.path(work,"pseudobulk_sample_qc.csv"));elig<-fread(file.path(work,"subtype_eligibility_preflight.csv"))
analyses<-list(Broad19_AllFrozenT=file.path(work,"Broad19_AllFrozenT_counts.csv.gz"),Broad19_VerifiedAuthorT=file.path(work,"Broad19_VerifiedAuthorT_counts.csv.gz"))
for(i in seq_len(nrow(elig))){if(elig$eligibility_status[[i]]!="ANALYZABLE")next;subtype<-elig$subtype[[i]];slug<-gsub(" ","_",gsub("-","MINUS",gsub("\\+","PLUS",subtype)));analyses[[paste0("Subtype::",subtype)]]<-file.path(work,paste0("Subtype_",slug,"_counts.csv.gz"))}
all_results<-list();all_models<-list();seed<-0L
for(analysis_name in names(analyses)){seed<-seed+1L;meta<-as.data.frame(qc[analysis==analysis_name]);fit<-fit_analysis(analysis_name,analyses[[analysis_name]],meta,seed);all_results[[analysis_name]]<-fit$result;all_models[[analysis_name]]<-fit$model;message("ANALYSIS_COMPLETE ",analysis_name," tested=",fit$model$tested_gene_ids)}
program_results<-rbindlist(all_results,fill=TRUE);model_audit<-rbindlist(all_models,fill=TRUE)
program_results[,subtype_family_adjusted_P:=NA_real_];idx<-grepl("^Subtype::",program_results$analysis)&is.finite(program_results$pval);program_results[idx,subtype_family_adjusted_P:=p.adjust(pval,method="BH")];program_results[,planned_subtype_family_tests:=18L];program_results[,valid_subtype_family_P_values:=sum(idx)]
fwrite(program_results,file.path(work,"program_results_raw.csv"));fwrite(model_audit,file.path(work,"model_audit.csv"))
lodo_rows<-list();li<-0L
for(analysis_name in names(analyses)[grepl("^Subtype::",names(analyses))]){
  primary<-program_results[analysis==analysis_name];matrix_path<-analyses[[analysis_name]];base_meta<-as.data.frame(qc[analysis==analysis_name])
  for(donor in base_meta$sample_column){
    reduced<-base_meta[base_meta$sample_column!=donor,,drop=FALSE];nh<-sum(reduced$condition=="Healthy");nc<-sum(reduced$condition=="ActiveCrohn");li<-li+1L
    if(nh<5L||nc<5L){invalid<-copy(info);invalid$analysis<-analysis_name;invalid$omitted_participant<-donor;invalid$n_healthy<-nh;invalid$n_crohn<-nc;invalid$NES<-NA_real_;invalid$pval<-NA_real_;invalid$adjusted_P_within_three<-NA_real_;invalid$direction<-"NON_ESTIMABLE";invalid$testability<-"NON_ESTIMABLE_DONOR_GATE";invalid$invalid_reason<-"Omission violates >=5 participants per condition";lodo_rows[[li]]<-invalid;next}
    matrix<-fread(matrix_path);keep_cols<-c("unique_gene_id","gene_symbol",reduced$sample_column);temp<-tempfile(pattern="phase5b_lodo_",tmpdir=work,fileext=".csv");fwrite(matrix[,..keep_cols],temp)
    fit<-tryCatch(fit_analysis(paste0(analysis_name,"::drop_",donor),temp,reduced,1000L+li,FALSE),error=function(e)e);unlink(temp)
    if(inherits(fit,"error")){invalid<-copy(info);invalid$analysis<-analysis_name;invalid$omitted_participant<-donor;invalid$n_healthy<-nh;invalid$n_crohn<-nc;invalid$NES<-NA_real_;invalid$pval<-NA_real_;invalid$adjusted_P_within_three<-NA_real_;invalid$direction<-"NON_ESTIMABLE";invalid$testability<-"NON_ESTIMABLE_MODEL_FAILURE";invalid$invalid_reason<-conditionMessage(fit);lodo_rows[[li]]<-invalid}
    else{rr<-copy(fit$result);rr$analysis<-analysis_name;rr$omitted_participant<-donor;rr$invalid_reason<-"";rr<-merge(rr,primary[,.(pathway_id,primary_NES=NES,primary_direction=direction)],by="pathway_id",all.x=TRUE);rr$NES_change<-rr$NES-rr$primary_NES;rr$sign_reversal<-fifelse(sign(rr$NES)==sign(rr$primary_NES),"NO","YES");lodo_rows[[li]]<-rr}
  }
}
fwrite(rbindlist(lodo_rows,fill=TRUE),file.path(work,"lodo_results_raw.csv"))
software<-data.table(software=c("R","edgeR","fgsea","data.table"),version=c(as.character(getRversion()),as.character(packageVersion("edgeR")),as.character(packageVersion("fgsea")),as.character(packageVersion("data.table"))));fwrite(software,file.path(work,"software_versions_R.csv"));message("PHASE5B_R_ANALYSIS_COMPLETE")

