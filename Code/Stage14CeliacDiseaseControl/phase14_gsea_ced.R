#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(fgsea); library(ggplot2)})
options(stringsAsFactors=FALSE)
root <- (function(x) { if (!nzchar(x)) stop("Set CANINE_CIE_PROJECT_ROOT"); x })(Sys.getenv("CANINE_CIE_PROJECT_ROOT"))
base <- file.path(root,"14_ced_specificity")
pb <- file.path(base,"pseudobulk")
groot <- file.path(base,"gsea")
full <- file.path(groot,"full")
leading <- file.path(groot,"leading_edges")
fig <- file.path(base,"figures")
for (d in c(groot,full,leading,fig,file.path(root,"00_metadata"),file.path(root,"logs")))
  dir.create(d,recursive=TRUE,showWarnings=FALSE)
sets <- c(
  Hallmark=file.path(root,"07_human_crohn","GSE266616","gsea","gene_sets","HUMAN_HALLMARK_GENE_SETS.tsv.gz"),
  Reactome=file.path(root,"07_human_crohn","GSE266616","gsea","gene_sets","HUMAN_REACTOME_GENE_SETS.tsv.gz"),
  GO_BP=file.path(root,"07_human_crohn","GSE266616","gsea","gene_sets","HUMAN_GO_BP_GENE_SETS.tsv.gz"))
lineages <- c("Epithelial","T","Myeloid")
min_size <- 15L; max_size <- 500L; eps_value <- 1e-10; base_seed <- 20260828L

read_rank <- function(path) {
  d <- fread(path)
  req <- c("gene_id","gene_symbol","logFC","F","PValue","FDR","rank_metric")
  if (!all(req %in% names(d))) stop("Rank schema failure: ",path)
  d <- d[!is.na(gene_symbol)&nzchar(gene_symbol)&is.finite(rank_metric)]
  d[,abs_rank:=abs(rank_metric)]
  setorder(d,-abs_rank,PValue,gene_id)
  d <- unique(d,by="gene_symbol"); d[,abs_rank:=NULL]
  setorder(d,-rank_metric,gene_symbol)
  s <- d$rank_metric; names(s) <- d$gene_symbol
  list(table=d,stats=s)
}

run_one <- function(lineage,collection) {
  rank <- read_rank(file.path(pb,paste0("CED_",lineage,"_RANKED.csv")))
  mem <- unique(fread(sets[[collection]]),by=c("pathway_id","gene_symbol"))
  req <- c("collection","pathway_id","pathway_name","gene_symbol",
           "source_version","source_species")
  if (!all(req %in% names(mem))) stop("Gene-set schema failure")
  if (any(mem$source_species!="Homo sapiens")) stop("Non-Human set entered CeD GSEA")
  info <- unique(mem[,.(pathway_id,pathway_name,source_version,source_species)],
                 by="pathway_id")
  pathways <- split(mem$gene_symbol,mem$pathway_id)
  set.seed(base_seed+match(lineage,lineages)*100L+match(collection,names(sets)))
  fg <- as.data.table(fgseaMultilevel(
    pathways=pathways,stats=rank$stats,minSize=min_size,maxSize=max_size,
    eps=eps_value,scoreType="std",nproc=1L))
  setnames(fg,"pathway","pathway_id")
  fg[,padj:=p.adjust(pval,method="BH")]
  z <- merge(fg,info,by="pathway_id",all.x=TRUE,sort=FALSE)
  z[,leading_edge_n:=lengths(leadingEdge)]
  z[,leading_edge:=vapply(leadingEdge,function(x)paste(as.character(x),collapse=";"),character(1))]
  z[,leadingEdge:=NULL]
  z[,lineage:=lineage]; z[,collection:=collection]
  z[,analysis:="Primary_Active_CeD_vs_Control"]
  z[,direction:=ifelse(NES>0,"Active_CeD_higher","Active_CeD_lower")]
  z[,absNES:=abs(NES)]; setorder(z,padj,-absNES,pathway_id); z[,absNES:=NULL]
  setcolorder(z,c("lineage","collection","analysis","pathway_id","pathway_name",
                  "ES","NES","pval","padj","size","leading_edge_n","leading_edge",
                  "direction","source_version","source_species","log2err"))
  le <- rbindlist(lapply(seq_len(nrow(z)),function(i) {
    genes <- strsplit(z$leading_edge[i],";",fixed=TRUE)[[1]]
    genes <- genes[nzchar(genes)]
    if (!length(genes)) return(NULL)
    q <- rank$table[match(genes,gene_symbol)]
    q[,.(lineage=lineage,collection=collection,pathway_id=z$pathway_id[i],
         pathway_name=z$pathway_name[i],pathway_NES=z$NES[i],
         pathway_FDR=z$padj[i],gene_id,gene_symbol,logFC,F,PValue,FDR,rank_metric)]
  }),fill=TRUE)
  fwrite(z,file.path(full,paste0("CED_",lineage,"_",toupper(collection),"_GSEA.csv")))
  fwrite(le,file.path(leading,paste0("CED_",lineage,"_",toupper(collection),"_LEADING_EDGES.csv")))
  z
}

all <- list(); k <- 1
for (lineage in lineages) for (collection in names(sets)) {
  message("CeD unbiased GSEA: ",lineage," / ",collection)
  all[[k]] <- run_one(lineage,collection); k <- k+1
}
combined <- rbindlist(all,fill=TRUE)
for (collection_name in names(sets)) {
  out <- combined[collection==collection_name]
  fwrite(out,file.path(groot,paste0("CED_",toupper(collection_name),"_GSEA.csv")))
}
programs <- combined[,.(lineage,program=pathway_name,collection,pathway=pathway_id,
                        NES,FDR=padj,direction,leading_edge,
                        sensitivity_status="NOT_YET_ASSESSED",
                        primary_or_secondary="PRIMARY",
                        notes="Unbiased Human-native CeD GSEA; frozen before any Dog/Crohn/Phase10 program table was read")]
programs[,absNES:=abs(NES)]
setorder(programs,lineage,collection,FDR,-absNES,pathway)
programs[,absNES:=NULL]
freeze_path <- file.path(groot,"CED_DISEASE_PROGRAMS.csv")
fwrite(programs,freeze_path)

sha_line <- system2("certutil",c("-hashfile",shQuote(freeze_path),"SHA256"),
                    stdout=TRUE,stderr=TRUE)
sha <- tolower(gsub(" ","",sha_line[grepl("^[0-9A-Fa-f ]{64,}$",sha_line)][1]))
if (is.na(sha)||nchar(sha)!=64) stop("Could not freeze SHA-256 for CED_DISEASE_PROGRAMS.csv")
writeLines(sha,file.path(groot,"CED_DISEASE_PROGRAMS.sha256"))
freeze <- data.frame(
  artifact=freeze_path,sha256=sha,rows=nrow(programs),
  frozen_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
  posthoc_prior_results_opened="NO",
  gate="PASS_INDEPENDENT_CED_GSEA_FROZEN")
fwrite(freeze,file.path(groot,"CED_GSEA_FREEZE_AUDIT.csv"))

top <- programs[FDR<0.10]
if (!nrow(top)) top <- programs[order(FDR)][1:min(.N,60)]
top <- top[order(FDR,-abs(NES))]
top <- top[,head(.SD,12),by=.(lineage,collection,direction)]
top[,label:=paste(lineage,program,sep=" | ")]
p <- ggplot(top,aes(x=lineage,y=reorder(label,NES),fill=NES))+
  geom_tile(color="white",linewidth=.2)+
  scale_fill_gradient2(low="#356AA0",mid="white",high="#B4475A",midpoint=0)+
  facet_grid(collection~.,scales="free_y",space="free_y")+
  labs(title="CeD unbiased Human-native GSEA",x=NULL,y=NULL,fill="NES")+
  theme_bw(base_size=9)+theme(axis.text.y=element_text(size=6),panel.grid=element_blank())
ggsave(file.path(fig,"FIG_CED_GSEA_NES_HEATMAP.png"),p,width=10,height=15,dpi=180)

params <- data.frame(
  parameter=c("method","minSize","maxSize","eps","scoreType","gene_sets",
              "rank_metric","freeze_sha256"),
  value=c("fgseaMultilevel","15","500","1e-10","std",
          "Phase7 frozen Hallmark 2026.1.Hs; Reactome v97; GO BP frozen",
          "sign(logFC)*sqrt(F)",sha))
fwrite(params,file.path(groot,"CED_GSEA_PARAMETERS.csv"))
versions <- data.frame(
  software=c("R","fgsea","data.table","ggplot2"),
  version=c(as.character(getRversion()),as.character(packageVersion("fgsea")),
            as.character(packageVersion("data.table")),as.character(packageVersion("ggplot2"))),
  role=c("runtime","multilevel pre-ranked GSEA","tabular IO","figures"))
fwrite(versions,file.path(root,"00_metadata","PHASE14_GSEA_SOFTWARE_VERSIONS.csv"))
writeLines(c(capture.output(sessionInfo()),
             paste("CED_DISEASE_PROGRAMS_SHA256",sha),
             paste("completed_utc",format(Sys.time(),tz="UTC",usetz=TRUE))),
           file.path(root,"logs","phase14_gsea.log"))
message("CeD independent GSEA frozen: ",sha)

