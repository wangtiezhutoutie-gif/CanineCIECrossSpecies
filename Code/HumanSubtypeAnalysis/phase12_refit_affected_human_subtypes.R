.libPaths(c(file.path((function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT")), "envs", "phase3_r_portable", "library_user"),.libPaths()))
suppressPackageStartupMessages({library(edgeR);library(fgsea);library(data.table)})
R <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
P30 <- file.path(R,"WORKING","30_HUMAN_PSEUDOBULK_INDEX_CORRECTION","work")
O <- file.path(R,"WORKING","31_FINAL_ELIGIBILITY_AND_MANUSCRIPT_AUDIT")
W <- file.path(O,"work")
dir.create(W,recursive=TRUE,showWarnings=FALSE)
G <- file.path(R,"07_human_crohn","GSE266616","gsea","gene_sets")
P <- data.table(
 lineage=c(rep("Epithelial",3),rep("Myeloid",3)),
 frozen_program_id=c("Epithelial::Mitochondrial respiration","Epithelial::Oxidative phosphorylation","Epithelial::Xenobiotic metabolism","Myeloid::TNF / NF-kappaB","Myeloid::Inflammatory response","Myeloid::TLR4"),
 program=c("Mitochondrial respiration","Oxidative phosphorylation","Xenobiotic metabolism","TNF / NF-kappaB","Inflammatory response","TLR4"),
 collection=c("Reactome","GO_BP","Hallmark","Hallmark","GO_BP","Reactome"),
 pathway_id=c("R-HSA-1428517","GO:0006119","HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_TNFA_SIGNALING_VIA_NFKB","GO:0050727","R-HSA-166016"))
files <- c(Reactome=file.path(G,"HUMAN_REACTOME_GENE_SETS.tsv.gz"),GO_BP=file.path(G,"HUMAN_GO_BP_GENE_SETS.tsv.gz"),Hallmark=file.path(G,"HUMAN_HALLMARK_GENE_SETS.tsv.gz"))
mem <- rbindlist(lapply(names(files),function(k){x<-fread(files[[k]]);x[pathway_id%in%P$pathway_id]}))
mem <- unique(mem,by=c("pathway_id","gene_symbol"))
slug <- function(x){y<-tolower(strsplit(x,"")[[1]]);y[!grepl("[a-z0-9]",y)]<-"_";gsub("^_+|_+$","",paste0(y,collapse=""))}
empty_result <- function(name,omit,nh,nc,why){
 lin<-strsplit(name,"__",fixed=TRUE)[[1]][1];u<-copy(P[lineage==lin]);
 u[,c("analysis","omitted_donor","healthy_donors","crohn_donors","ES","NES","pval","size","leading_edge_n","tested_genes","mapped_ranked_genes","model_status","direction"):=list(name,omit,nh,nc,NA_real_,NA_real_,NA_real_,NA_integer_,NA_integer_,NA_integer_,NA_integer_,why,"NOT_ESTIMABLE")];u}
fit <- function(name,omit="NONE",untreated=FALSE,off=0){
 cts<-fread(file.path(P30,paste0("human_",slug(name),"_counts.tsv")))
 m<-fread(file.path(P30,paste0("human_",slug(name),"_meta.tsv")))
 if(anyDuplicated(m$donor))stop("DUPLICATE_DONOR")
 if(grepl("__Subtype__",name,fixed=TRUE))m<-m[n_cells>=20]
 if(omit!="NONE")m<-m[donor!=omit]
 if(untreated)m<-m[condition=="Healthy"|tolower(biologic_treatment)%in%c("no","no_biologic","none","na","")]
 if(nrow(m)&&any(m$n_cells<20))stop("LOW_CELL_SAMPLE_RETAINED")
 if(anyDuplicated(m$donor))stop("DUPLICATE_DONOR_AFTER_FILTER")
 nh<-m[condition=="Healthy",uniqueN(donor)];nc<-m[condition=="ActiveCrohn",uniqueN(donor)]
 if(nh<5||nc<5)return(list(ok=FALSE,x=empty_result(name,omit,nh,nc,"NON_ESTIMABLE_DONOR_GATE")))
 ids<-m$sample_column
 if(any(!ids%in%names(cts)))stop("COUNT_COLUMN_MISSING")
 a<-cts[,.(unique_gene_id,gene_symbol)];x<-as.matrix(cts[,..ids]);storage.mode(x)<-"integer";rownames(x)<-a$unique_gene_id
 m[,condition:=factor(condition,levels=c("Healthy","ActiveCrohn"))];d<-model.matrix(~condition,m)
 if(qr(d)$rank!=ncol(d))stop("RANK_DEFICIENT_DESIGN")
 y<-DGEList(x);k<-filterByExpr(y,design=d);y<-y[k,,keep.lib.sizes=FALSE];y<-calcNormFactors(y);y<-estimateDisp(y,d,robust=TRUE);f<-glmQLFit(y,d,robust=TRUE);q<-glmQLFTest(f,coef="conditionActiveCrohn")
 raw<-topTags(q,n=Inf,sort.by="none")$table;t<-as.data.table(raw,keep.rownames="unique_gene_id");t[,rank:=sign(logFC)*sqrt(F)];t<-merge(t,a,by="unique_gene_id");t<-t[!is.na(gene_symbol)&nzchar(gene_symbol)&is.finite(rank)];t[,ar:=abs(rank)];setorder(t,gene_symbol,-ar,PValue,unique_gene_id);t<-t[,.SD[1],by=gene_symbol]
 s<-t$rank;names(s)<-t$gene_symbol;lin<-strsplit(name,"__",fixed=TRUE)[[1]][1];pp<-P[lineage==lin];paths<-split(mem[pathway_id%in%pp$pathway_id,gene_symbol],mem[pathway_id%in%pp$pathway_id,pathway_id])
 set.seed(20261016+off);g<-as.data.table(fgseaMultilevel(paths,s,minSize=15,maxSize=500,eps=1e-10,nproc=1));setnames(g,"pathway","pathway_id");if(nrow(g)){g[,leading_edge_n:=lengths(leadingEdge)];g[,leadingEdge:=NULL]}
 u<-merge(pp,g[,.(pathway_id,ES,NES,pval,size,leading_edge_n)],by="pathway_id",all.x=TRUE)
 u[,c("analysis","omitted_donor","healthy_donors","crohn_donors","tested_genes","mapped_ranked_genes","model_status","direction"):=list(name,omit,nh,nc,sum(k),nrow(t),ifelse(is.na(NES),"PATHWAY_NOT_TESTABLE","ESTIMABLE"),ifelse(is.na(NES),"NOT_ASSESSABLE",ifelse(NES>0,"CROHN_HIGHER","CROHN_LOWER")))]
 list(ok=TRUE,x=u)}
analyses <- data.table(name=c("Epithelial__Subtype__Early Enterocyte","Epithelial__Subtype__TA"),base=c(300,600))
F<-list();L<-list();T<-list()
for(i in seq_len(nrow(analyses))){
 n<-analyses$name[i];base<-analyses$base[i]
 a<-fit(n,"NONE",FALSE,base);F[[length(F)+1]]<-a$x
 m0<-fread(file.path(P30,paste0("human_",slug(n),"_meta.tsv")));m0[,original_row:=.I];eligible<-m0[n_cells>=20]
 for(j in seq_len(nrow(eligible))){dd<-eligible$donor[j];seedoff<-base+eligible$original_row[j];b<-fit(n,dd,FALSE,seedoff);L[[length(L)+1]]<-b$x}
 tt<-fit(n,"NONE",TRUE,base+99);T[[length(T)+1]]<-tt$x
}
F<-rbindlist(F,fill=TRUE);L<-rbindlist(L,fill=TRUE);T<-rbindlist(T,fill=TRUE)
B<-F[,.(lineage,analysis,frozen_program_id,full_NES=NES)];L<-merge(L,B,by=c("lineage","analysis","frozen_program_id"),all.x=TRUE);L[,direction_match_full:=fifelse(is.na(NES),"NOT_ESTIMABLE",fifelse(sign(NES)==sign(full_NES),"YES","NO"))]
fwrite(F,file.path(W,"phase12_affected_results.tsv"),sep="\t",na="")
fwrite(L,file.path(W,"phase12_affected_lodo.tsv"),sep="\t",na="")
fwrite(T,file.path(W,"phase12_affected_treatment.tsv"),sep="\t",na="")
cat("PHASE12_AFFECTED_DONE ",nrow(F)," ",nrow(L)," ",nrow(T),"\n",sep="")
