#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(edgeR)})
root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
base <- file.path(root,"14_ced_specificity")
pb <- file.path(base,"pseudobulk")
sens <- file.path(base,"sensitivity")
fig <- file.path(base,"figures")
dir.create(sens,recursive=TRUE,showWarnings=FALSE)
dir.create(fig,recursive=TRUE,showWarnings=FALSE)
lineages <- c("Epithelial","T","Myeloid")

read_counts <- function(path) {
  d <- read.csv(gzfile(path),check.names=FALSE)
  req <- c("gene_id","gene_symbol","gene_type")
  if (!all(req %in% names(d))) stop("Missing feature columns: ",path)
  if (anyDuplicated(d$gene_id)) stop("Duplicate gene IDs: ",path)
  x <- as.matrix(d[,setdiff(names(d),req),drop=FALSE])
  storage.mode(x) <- "double"; rownames(x) <- d$gene_id
  if (anyNA(x)||any(x<0)) stop("Invalid counts: ",path)
  list(annotation=d[,req],counts=x)
}

fit_one <- function(x,meta) {
  if (!setequal(colnames(x),meta$subject_id)) stop("Subject mismatch")
  meta <- meta[match(colnames(x),meta$subject_id),,drop=FALSE]
  meta$condition <- factor(meta$condition,levels=c("Control","CeD"))
  design <- model.matrix(~condition,data=meta)
  rownames(design) <- meta$subject_id
  y0 <- DGEList(counts=x)
  keep <- filterByExpr(y0,design=design)
  if (sum(keep)<1000) stop("Fewer than 1000 genes pass filterByExpr")
  y <- y0[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y,method="TMM")
  y <- estimateDisp(y,design,robust=TRUE)
  fit <- glmQLFit(y,design,robust=TRUE)
  qlf <- glmQLFTest(fit,coef="conditionCeD")
  tab <- topTags(qlf,n=Inf,sort.by="none")$table
  list(meta=meta,design=design,keep=keep,y=y,fit=fit,qlf=qlf,tab=tab)
}

format_result <- function(z,ann,lineage,analysis) {
  tab <- z$tab
  a <- ann[match(rownames(tab),ann$gene_id),,drop=FALSE]
  out <- data.frame(
    gene_id=rownames(tab),gene_symbol=a$gene_symbol,gene_type=a$gene_type,
    logFC=tab$logFC,logCPM=tab$logCPM,F=tab$F,
    PValue=tab$PValue,FDR=tab$FDR,
    rank_metric=sign(tab$logFC)*sqrt(pmax(tab$F,0)),
    direction=ifelse(tab$logFC>0,"Active_CeD_higher",
                     ifelse(tab$logFC<0,"Active_CeD_lower","No_change")),
    lineage=lineage,analysis=analysis,
    n_control=sum(z$meta$condition=="Control"),
    n_active_CeD=sum(z$meta$condition=="CeD"),check.names=FALSE)
  out[order(out$PValue,-abs(out$logFC)),]
}

plot_qc <- function(z,res,lineage) {
  lcpm <- cpm(z$y,log=TRUE,prior.count=2)
  v <- order(apply(lcpm,1,var),decreasing=TRUE)
  pca <- prcomp(t(lcpm[v[seq_len(min(1000,length(v)))],,drop=FALSE]))
  qc <- data.frame(
    lineage=lineage,subject_id=z$meta$subject_id,condition=z$meta$condition,
    n_cells=z$meta$n_cells,raw_library_size=z$y$samples$lib.size,
    normalization_factor=z$y$samples$norm.factors,
    effective_library_size=z$y$samples$lib.size*z$y$samples$norm.factors,
    detected_genes=z$meta$detected_genes,PCA1=pca$x[,1],PCA2=pca$x[,2])
  write.csv(qc,file.path(pb,paste0("CED_",lineage,"_EDGER_SAMPLE_QC.csv")),row.names=FALSE)
  write.csv(data.frame(gene_id=rownames(lcpm),lcpm,check.names=FALSE),
            file.path(pb,paste0("CED_",lineage,"_LOGCPM.csv")),row.names=FALSE)
  col <- ifelse(qc$condition=="CeD","#B4475A","#3C78A8")
  png(file.path(fig,paste0("FIG_CED_",lineage,"_PCA.png")),1400,1100,res=160)
  plot(qc$PCA1,qc$PCA2,pch=19,col=col,cex=1.4,
       xlab="PCA1",ylab="PCA2",main=paste("CeD",lineage,"pseudobulk PCA"))
  text(qc$PCA1,qc$PCA2,labels=qc$subject_id,pos=3,cex=.65)
  legend("topright",legend=c("Control","Active CeD"),col=c("#3C78A8","#B4475A"),pch=19,bty="n")
  dev.off()
  corx <- cor(lcpm,method="spearman")
  png(file.path(fig,paste0("FIG_CED_",lineage,"_SAMPLE_CORRELATION.png")),1500,1350,res=160)
  heatmap(corx,symm=TRUE,scale="none",margins=c(7,7),
          main=paste("CeD",lineage,"sample Spearman correlation"))
  dev.off()
  sig <- res$FDR<0.05
  vc <- ifelse(sig & res$logFC>0,"#B4475A",ifelse(sig & res$logFC<0,"#3C78A8","#B9B9B9"))
  png(file.path(fig,paste0("FIG_CED_",lineage,"_VOLCANO.png")),1300,1100,res=160)
  plot(res$logFC,-log10(pmax(res$FDR,1e-300)),pch=16,cex=.55,col=vc,
       xlab="logFC (Active CeD vs Control)",ylab="-log10 FDR",
       main=paste("CeD",lineage,"donor-level edgeR QL"))
  abline(h=-log10(.05),lty=2,col="#555555"); abline(v=0,lty=3)
  dev.off()
}

meta_all <- read.csv(file.path(pb,"CED_PSEUDOBULK_SAMPLE_QC.csv"))
summaries <- list(); jack <- list(); ji <- 1
for (lineage in lineages) {
  message("Primary edgeR: ",lineage)
  inp <- read_counts(file.path(pb,paste0("CED_",lineage,"_SUBJECT_COUNTS.csv.gz")))
  meta <- meta_all[meta_all$Level1==lineage & meta_all$coverage_eligible=="YES",]
  x <- inp$counts[,meta$subject_id,drop=FALSE]
  z <- fit_one(x,meta)
  res <- format_result(z,inp$annotation,lineage,"PRIMARY_SUBJECT_AGGREGATED")
  write.csv(res,file.path(pb,paste0("CED_",lineage,"_Active_vs_Control.csv")),row.names=FALSE)
  ranked <- res[order(res$rank_metric,decreasing=TRUE),
                c("gene_id","gene_symbol","logFC","F","PValue","FDR","rank_metric")]
  write.csv(ranked,file.path(pb,paste0("CED_",lineage,"_RANKED.csv")),row.names=FALSE)
  plot_qc(z,res,lineage)
  summaries[[lineage]] <- data.frame(
    lineage=lineage,n_control=sum(z$meta$condition=="Control"),
    n_active_CeD=sum(z$meta$condition=="CeD"),tested_genes=nrow(res),
    FDR05=sum(res$FDR<.05),FDR10=sum(res$FDR<.10),
    FDR05_higher=sum(res$FDR<.05 & res$logFC>0),
    FDR05_lower=sum(res$FDR<.05 & res$logFC<0),
    rank_formula="sign(logFC)*sqrt(F)",design="~condition",
    reference="Control",robust_QL=TRUE)

  message("Primary-library-only sensitivity: ",lineage)
  pin <- read_counts(file.path(pb,paste0("CED_",lineage,"_PRIMARY_LIBRARY_ONLY_COUNTS.csv.gz")))
  pz <- fit_one(pin$counts[,meta$subject_id,drop=FALSE],meta)
  pres <- format_result(pz,pin$annotation,lineage,"PRIMARY_LIBRARY_ONLY_DROP_REP2")
  write.csv(pres,file.path(sens,paste0("CED_",lineage,"_PRIMARY_LIBRARY_ONLY.csv")),row.names=FALSE)

  message("32-subject jackknife: ",lineage)
  for (drop in meta$subject_id) {
    keep_subject <- setdiff(meta$subject_id,drop)
    jm <- meta[match(keep_subject,meta$subject_id),,drop=FALSE]
    jz <- fit_one(x[,keep_subject,drop=FALSE],jm)
    jr <- format_result(jz,inp$annotation,lineage,paste0("DROP_",drop))
    jack[[ji]] <- data.frame(
      lineage=lineage,dropped_subject=drop,
      dropped_condition=as.character(meta$condition[meta$subject_id==drop]),
      gene_id=jr$gene_id,gene_symbol=jr$gene_symbol,logFC=jr$logFC,
      F=jr$F,FDR=jr$FDR,rank_metric=jr$rank_metric)
    ji <- ji+1
  }
}
write.csv(do.call(rbind,summaries),file.path(pb,"CED_EDGER_SUMMARY.csv"),row.names=FALSE)
jout <- do.call(rbind,jack)
con <- gzfile(file.path(sens,"CED_32_SUBJECT_JACKKNIFE_GENE_RESULTS.csv.gz"),"wt")
write.csv(jout,con,row.names=FALSE); close(con)
versions <- data.frame(
  software=c("R","edgeR","limma","locfit","statmod"),
  version=c(as.character(getRversion()),as.character(packageVersion("edgeR")),
            as.character(packageVersion("limma")),as.character(packageVersion("locfit")),
            as.character(packageVersion("statmod"))),
  role=c("runtime","donor-level quasi-likelihood DE","edgeR dependency",
         "robust dispersion dependency","robust QL dependency"))
write.csv(versions,file.path(root,"00_metadata","PHASE14_EDGER_SOFTWARE_VERSIONS.csv"),row.names=FALSE)
writeLines(c(capture.output(sessionInfo()),
             paste("completed_utc",format(Sys.time(),tz="UTC",usetz=TRUE))),
           file.path(root,"logs","phase14_edger.log"))
message("Phase14 CeD edgeR complete")

