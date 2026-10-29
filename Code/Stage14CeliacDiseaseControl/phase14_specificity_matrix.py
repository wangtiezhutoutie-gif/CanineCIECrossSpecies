#!/usr/bin/env python
"""Post-freeze Dog-Crohn-CeD specificity, concordance, provenance and figures."""
from pathlib import Path
import hashlib, json, re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

ROOT=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P14=ROOT/"14_ced_specificity"; G=P14/"gsea"; S=P14/"specificity"
PB=P14/"pseudobulk"; SEN=P14/"sensitivity"; FIG=P14/"figures"; SUM=P14/"summary"
for d in (S,FIG,SUM): d.mkdir(parents=True,exist_ok=True)

def sha256(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()

freeze=G/"CED_DISEASE_PROGRAMS.csv"
expected=(G/"CED_DISEASE_PROGRAMS.sha256").read_text().strip()
actual=sha256(freeze)
if expected!=actual: raise RuntimeError("CeD program freeze hash mismatch")

base=pd.read_csv(ROOT/"10_cross_species_baseline"/"summary"/"PHASE10_PROGRAM_BASELINE_MATRIX.csv")
if len(base)!=25: raise RuntimeError("Phase10 registered panel is not 25 rows")
ced=pd.read_csv(freeze)
p11=pd.read_csv(ROOT/"11_singIST"/"summary"/"PHASE11_VS_PHASE10_BASELINE.csv")
p12=pd.read_csv(ROOT/"12_saturn"/"summary"/"PHASE12_PROVENANCE_MATRIX.csv")
p13=pd.read_csv(ROOT/"13_concord"/"summary"/"PHASE13_PROVENANCE_MATRIX.csv")
jk=pd.read_csv(SEN/"CED_DONOR_JACKKNIFE_PROGRAM_SUMMARY.csv")
rep=pd.read_csv(SEN/"CED_REPLICATE_SENSITIVITY.csv")

key=["lineage","program"]
m=base.merge(
    ced.rename(columns={"pathway":"representative_pathway","NES":"CeD_NES",
                        "FDR":"CeD_FDR","direction":"CeD_direction",
                        "leading_edge":"CeD_leading_edge"})[
        ["lineage","representative_pathway","CeD_NES","CeD_FDR",
         "CeD_direction","CeD_leading_edge"]],
    on=["lineage","representative_pathway"],how="left")
m=m.merge(jk[key+["jackknife_status","same_direction_fraction",
                  "max_abs_delta_NES","most_influential_donor"]],
          on=key,how="left")
m=m.merge(rep[key+["replicate_status","delta_NES"]].rename(
    columns={"delta_NES":"replicate_delta_NES"}),on=key,how="left")

def sign_dir(x):
    x=str(x).lower()
    if "higher" in x: return 1
    if "lower" in x: return -1
    return 0

def relationship(r):
    if pd.isna(r.CeD_NES): return "CED_NOT_TESTABLE"
    if abs(r.CeD_NES)<1.0 and r.CeD_FDR>=0.10: return "CED_NOT_PROMINENT"
    return "CED_SHARED_DIRECTION" if sign_dir(r.human_TI_direction_frozen)==sign_dir(r.CeD_direction) else "CED_OPPOSITE_DIRECTION"

m["CeD_relationship"]=m.apply(relationship,axis=1)
m["dog_crohn_direction_conserved"]=m.apply(
    lambda r: sign_dir(r.dog_direction)==sign_dir(r.human_TI_direction_frozen)!=0,axis=1)
m["CeD_evidence"]=np.select(
    [m.CeD_NES.isna(),m.CeD_FDR<.10,m.CeD_NES.abs()>=1.3],
    ["NOT_TESTABLE","FDR_LT_0.10","DIRECTIONAL_NES_GE_1.3"],default="WEAK_NOT_PROMINENT")

def evidence_class(r):
    if not r.dog_crohn_direction_conserved: return "DOG_HUMAN_NOT_CONSERVED"
    if r.CeD_relationship=="CED_NOT_TESTABLE": return "NOT_ASSESSABLE"
    if r.region_risk=="HIGH_REGION_RISK" or "UNKNOWN" in str(r.region_risk):
        return "REGION_CONFOUNDED_OR_UNCERTAIN"
    if r.CeD_relationship=="CED_SHARED_DIRECTION" and r.CeD_FDR<.10:
        return "GENERIC_MUCOSAL_INFLAMMATORY"
    if r.CeD_relationship=="CED_OPPOSITE_DIRECTION" and r.CeD_FDR<.10:
        return "CED_DIFFERENTIATING_STRONG"
    return "CED_DIFFERENTIATING_LIMITED"

m["specificity_evidence_class"]=m.apply(evidence_class,axis=1)
def forward(r):
    c=r.specificity_evidence_class
    if c=="NOT_ASSESSABLE": return "NO"
    if c=="DOG_HUMAN_NOT_CONSERVED": return "EXPLORATORY_ONLY"
    if c=="REGION_CONFOUNDED_OR_UNCERTAIN": return "YES_WITH_REGION_CAUTION"
    if c=="GENERIC_MUCOSAL_INFLAMMATORY": return "YES_WITH_GENERIC_INFLAMMATION_CAUTION"
    return "YES_CORE"
m["forward_to_phase15"]=m.apply(forward,axis=1)
m["priority"]=np.select(
    [m.specificity_evidence_class=="CED_DIFFERENTIATING_STRONG",
     m.specificity_evidence_class=="GENERIC_MUCOSAL_INFLAMMATORY",
     m.specificity_evidence_class=="REGION_CONFOUNDED_OR_UNCERTAIN",
     m.specificity_evidence_class=="DOG_HUMAN_NOT_CONSERVED"],
    ["HIGH_CED_DIFFERENTIATING","CORE_GENERIC_INFLAMMATION",
     "REGION_CAUTION","EXPLORATORY"],default="MODERATE")
m["notes_phase14"]=np.where(
    m.CeD_relationship=="CED_NOT_PROMINENT",
    "CeD effect testable but |NES|<1 and FDR>=0.10; not evidence of absolute Crohn specificity",
    "CeD comparator interpretation is relative to this active CeD dataset only")

phase11=p11[key+["reference_model_status","interpretation"]].copy()
phase11["phase11_status"]=phase11.reference_model_status.astype(str)+";"+phase11.interpretation.astype(str)
phase12=p12[key+["SATURN_embedding_cosine","SATURN_macrogene_rho"]].copy()
phase12["phase12_status"]="cosine="+phase12.SATURN_embedding_cosine.round(4).astype(str)+";macrogene_rho="+phase12.SATURN_macrogene_rho.round(4).astype(str)
phase13=p13[key+["phase13_CONCORD_cosine","phase13_CONCORD_spearman","phase13_joint_null_p"]].copy()
phase13["phase13_status"]="cosine="+phase13.phase13_CONCORD_cosine.round(4).astype(str)+";rho="+phase13.phase13_CONCORD_spearman.round(4).astype(str)+";null_p="+phase13.phase13_joint_null_p.round(4).astype(str)
m=m.merge(phase11[key+["phase11_status"]],on=key,how="left")
m=m.merge(phase12[key+["phase12_status"]],on=key,how="left")
m=m.merge(phase13[key+["phase13_status"]],on=key,how="left")

matrix=pd.DataFrame({
    "lineage":m.lineage,"program":m.program,
    "Dog_CIE_direction":m.dog_direction,"Dog_LODO":m.dog_lodo_stability,
    "Human_Crohn_direction":m.human_TI_direction_frozen,
    "Human_Crohn_FDR":m.human_FDR,
    "Human_treatment_status":m.human_treatment_stability,
    "IBDverse_status":m.ibdverse_status,
    "CeD_direction":m.CeD_direction,"CeD_NES":m.CeD_NES,"CeD_FDR":m.CeD_FDR,
    "Dog_Crohn_baseline_status":m.baseline_status,"region_risk":m.region_risk,
    "FRE_status":m.fre_summary_status,"FRE_direction":m.FRE_direction,
    "CeD_relationship":m.CeD_relationship,
    "specificity_evidence_class":m.specificity_evidence_class,
    "CeD_jackknife":m.jackknife_status,
    "CeD_replicate_sensitivity":m.replicate_status,
    "notes":m.notes_phase14})
matrix.to_csv(S/"DOG_CROHN_CED_PROGRAM_MATRIX.csv",index=False)
matrix.specificity_evidence_class.value_counts().rename_axis(
    "specificity_evidence_class").reset_index(name="n_programs").to_csv(
        S/"CED_SPECIFICITY_CLASS_COUNTS.csv",index=False)

forward_df=pd.DataFrame({
    "lineage":m.lineage,"program":m.program,"dog_direction":m.dog_direction,
    "dog_lodo":m.dog_lodo_stability,"crohn_direction":m.human_TI_direction_frozen,
    "crohn_evidence":m.baseline_status,"ced_direction":m.CeD_direction,
    "ced_evidence":m.CeD_evidence,"ced_relationship":m.CeD_relationship,
    "specificity_evidence_class":m.specificity_evidence_class,
    "region_risk":m.region_risk,"fre_status":m.fre_summary_status,
    "human_treatment":m.human_treatment_stability,"ibdverse":m.ibdverse_status,
    "phase10_status":m.baseline_status,"phase11_status":m.phase11_status,
    "phase12_status":m.phase12_status,"phase13_status":m.phase13_status,
    "forward_to_phase15":m.forward_to_phase15,"priority":m.priority,
    "notes":m.notes_phase14})
forward_df.to_csv(S/"PHASE14_SPECIFICITY_FORWARD_LIST.csv",index=False)

# Human disease-to-disease gene concordance.
gene_rows=[]
for lineage in ["Epithelial","T","Myeloid"]:
    c=pd.read_csv(PB/f"CED_{lineage}_Active_vs_Control.csv")
    h=pd.read_csv(ROOT/"07_human_crohn"/"GSE266616"/"pseudobulk"/"de_results"/f"Human_{lineage}_ActiveCrohn_vs_Healthy.csv")
    z=h[["gene_id","logFC","FDR"]].rename(columns={"logFC":"Crohn_logFC","FDR":"Crohn_FDR"}).merge(
      c[["gene_id","logFC","FDR"]].rename(columns={"logFC":"CeD_logFC","FDR":"CeD_FDR"}),on="gene_id")
    rho=spearmanr(z.Crohn_logFC,z.CeD_logFC).statistic
    pr=pearsonr(z.Crohn_logFC,z.CeD_logFC).statistic
    gene_rows.append(dict(
      lineage=lineage,common_genes=len(z),spearman_rho=rho,pearson_r=pr,
      direction_concordance=np.mean(np.sign(z.Crohn_logFC)==np.sign(z.CeD_logFC)),
      both_FDR10_same=int(((z.Crohn_FDR<.1)&(z.CeD_FDR<.1)&
                           (np.sign(z.Crohn_logFC)==np.sign(z.CeD_logFC))).sum()),
      both_FDR10_opposite=int(((z.Crohn_FDR<.1)&(z.CeD_FDR<.1)&
                               (np.sign(z.Crohn_logFC)!=np.sign(z.CeD_logFC))).sum())))
pd.DataFrame(gene_rows).to_csv(S/"CROHN_VS_CED_GENE_CONCORDANCE.csv",index=False)

# Human pathway concordance on identical frozen collections.
path_rows=[]; path_points=[]
for lineage in ["Epithelial","T","Myeloid"]:
  for collection in ["Hallmark","Reactome","GO_BP"]:
    c=pd.read_csv(G/f"CED_{collection.upper()}_GSEA.csv")
    c=c[c.lineage==lineage]
    h=pd.read_csv(ROOT/"07_human_crohn"/"GSE266616"/"gsea"/"full"/f"Human_{lineage}_{collection}_GSEA.csv")
    z=h[["pathway_id","pathway_name","NES","padj"]].rename(columns={"NES":"Crohn_NES","padj":"Crohn_FDR"}).merge(
      c[["pathway_id","NES","padj"]].rename(columns={"NES":"CeD_NES","padj":"CeD_FDR"}),on="pathway_id")
    path_points.append(z.assign(lineage=lineage,collection=collection))
    same=np.sign(z.Crohn_NES)==np.sign(z.CeD_NES)
    both=(z.Crohn_FDR<.1)&(z.CeD_FDR<.1)
    path_rows.append(dict(
      lineage=lineage,collection=collection,common_pathways=len(z),
      spearman_rho=spearmanr(z.Crohn_NES,z.CeD_NES).statistic,
      pearson_r=pearsonr(z.Crohn_NES,z.CeD_NES).statistic,
      direction_concordance=same.mean(),
      both_FDR10_same=int((both&same).sum()),
      both_FDR10_opposite=int((both&~same).sum()),
      Crohn_FDR10_only=int(((z.Crohn_FDR<.1)&(z.CeD_FDR>=.1)).sum()),
      CeD_FDR10_only=int(((z.CeD_FDR<.1)&(z.Crohn_FDR>=.1)).sum())))
paths=pd.concat(path_points,ignore_index=True)
pd.DataFrame(path_rows).to_csv(S/"CROHN_VS_CED_PATHWAY_CONCORDANCE.csv",index=False)

# Positive controls are descriptive sanity checks only.
pc=ced[ced.program.str.contains(
  r"interferon|inflammatory|cytokine|g2m|e2f|stem|crypt|secret",
  case=False,regex=True,na=False)].sort_values(["lineage","FDR"])
pc.groupby(["lineage","collection"],as_index=False).head(10).to_csv(
  S/"CED_POSITIVE_CONTROL_SANITY.csv",index=False)

context=m[["lineage","program","dog_direction","dog_lodo_stability",
           "FRE_direction","fre_summary_status","human_TI_direction_frozen",
           "human_treatment_stability","ibdverse_status","CeD_direction",
           "CeD_NES","CeD_FDR","CeD_relationship","region_risk",
           "baseline_status","phase11_status","phase12_status","phase13_status",
           "specificity_evidence_class","forward_to_phase15"]].copy()
context.to_csv(SUM/"DISEASE_CONTEXT_MATRIX.csv",index=False)

# Ten formal figure types. Cell abundance remains descriptive.
sns.set_theme(style="whitegrid",context="paper")
qc=pd.read_csv(PB/"CED_PSEUDOBULK_SAMPLE_QC.csv")
plt.figure(figsize=(11,6))
sns.scatterplot(data=qc,x="subject_id",y="n_cells",hue="condition",
                style="Level1",s=70)
plt.xticks(rotation=90); plt.title("CeD subject × Level1 author-final cells")
plt.tight_layout(); plt.savefig(FIG/"FIG01_CED_DONOR_SAMPLE_QC.png",dpi=220); plt.close()

comp=qc.pivot_table(index=["subject_id","condition"],columns="Level1",
                    values="n_cells",aggfunc="sum").fillna(0)
prop=comp.div(comp.sum(axis=1),axis=0)
prop.plot(kind="bar",stacked=True,figsize=(12,6),
          color=["#4E79A7","#F28E2B","#59A14F"])
plt.ylabel("Author-final cell fraction (descriptive)"); plt.title("CeD Level1 composition")
plt.xticks(rotation=90); plt.tight_layout()
plt.savefig(FIG/"FIG02_CED_LEVEL1_CELL_COMPOSITION.png",dpi=220); plt.close()

fig,axs=plt.subplots(1,3,figsize=(15,4.5))
for ax,lineage in zip(axs,["Epithelial","T","Myeloid"]):
    d=pd.read_csv(PB/f"CED_{lineage}_Active_vs_Control.csv")
    sig=d.FDR<.05
    ax.scatter(d.logFC,-np.log10(np.maximum(d.FDR,1e-300)),s=5,
               c=np.where(sig,np.where(d.logFC>0,"#B4475A","#3C78A8"),"#B9B9B9"),alpha=.7)
    ax.axhline(-np.log10(.05),ls="--",c="grey",lw=.7); ax.axvline(0,ls=":",c="grey",lw=.7)
    ax.set(title=lineage,xlabel="logFC",ylabel="-log10 FDR")
fig.suptitle("Active CeD vs Control donor-level edgeR QL")
fig.tight_layout(); fig.savefig(FIG/"FIG03_CED_THREE_LINEAGE_VOLCANO.png",dpi=220); plt.close(fig)

heat=m.set_index(m.lineage+" | "+m.program)[
    ["dog_common1to1_NES","human_common1to1_NES","CeD_NES"]]
heat.columns=["Dog CIE","Human Crohn","Human CeD"]
plt.figure(figsize=(8,11)); sns.heatmap(heat,cmap="vlag",center=0,linewidths=.3)
plt.title("Registered 25-program disease context (NES)")
plt.tight_layout(); plt.savefig(FIG/"FIG05_DOG_CROHN_CED_25_PROGRAM_HEATMAP.png",dpi=220); plt.close()

g=sns.FacetGrid(paths,col="lineage",row="collection",height=2.6,margin_titles=True)
g.map_dataframe(sns.scatterplot,x="Crohn_NES",y="CeD_NES",s=16,alpha=.6)
for ax in g.axes.flat:
    lo=min(ax.get_xlim()[0],ax.get_ylim()[0]); hi=max(ax.get_xlim()[1],ax.get_ylim()[1])
    ax.plot([lo,hi],[lo,hi],ls="--",c="grey",lw=.7); ax.axhline(0,c="grey",lw=.4); ax.axvline(0,c="grey",lw=.4)
g.fig.subplots_adjust(top=.94); g.fig.suptitle("Human Crohn vs CeD pathway NES")
g.savefig(FIG/"FIG06_CROHN_VS_CED_PATHWAY_NES_SCATTER.png",dpi=220); plt.close(g.fig)

cnt=m.specificity_evidence_class.value_counts().sort_values()
cnt.plot(kind="barh",figsize=(9,4.5),color="#4E79A7")
plt.xlabel("Registered programs"); plt.title("Phase14 specificity evidence classes")
plt.tight_layout(); plt.savefig(FIG/"FIG07_SPECIFICITY_EVIDENCE_CLASSES.png",dpi=220); plt.close()

def disease_panel(rows,path,title,out):
    d=rows.set_index(rows.lineage+" | "+rows.program)[
      ["dog_common1to1_NES","human_common1to1_NES","CeD_NES"]]
    d.columns=["Dog CIE","Human Crohn","Human CeD"]
    d.plot(kind="bar",figsize=(11,5),color=["#59A14F","#4E79A7","#E15759"])
    plt.axhline(0,c="grey",lw=.7); plt.ylabel("NES"); plt.title(title)
    plt.xticks(rotation=35,ha="right"); plt.tight_layout(); plt.savefig(path/out,dpi=220); plt.close()

infl=m[m.program.str.contains("TNF|Inflammatory|Cytokine",case=False,regex=True)]
ifn=paths[(paths.collection=="Hallmark")&
          paths.pathway_name.str.contains("Interferon",case=False,na=False)]
ifn_rows=pd.DataFrame({
    "lineage":ifn.lineage,
    "program":ifn.pathway_name+" (Human IFN)",
    "dog_common1to1_NES":np.nan,
    "human_common1to1_NES":ifn.Crohn_NES,
    "CeD_NES":ifn.CeD_NES})
infl=pd.concat([infl,ifn_rows],ignore_index=True,sort=False)
disease_panel(infl,FIG,"TNF / IFN / inflammatory cytokine-JAK context",
              "FIG08_TNF_IFN_JAK_SHARED_INFLAMMATION.png")
epi=m[(m.lineage=="Epithelial")&m.program.str.contains(
    "Mitochondrial|Oxidative|Fatty|Xenobiotic",case=False,regex=True)]
disease_panel(epi,FIG,"Epithelial metabolic disease context",
              "FIG09_EPITHELIAL_METABOLIC_CONTEXT.png")
tmet=m[(m.lineage=="T")&m.program.str.contains(
    "Mitochondrial|Oxidative|Metabolic",case=False,regex=True)]
disease_panel(tmet,FIG,"T-cell metabolic contrast",
              "FIG10_T_METABOLIC_CONTRAST.png")

class_counts=m.specificity_evidence_class.value_counts().to_dict()
rel_counts=m.CeD_relationship.value_counts().to_dict()
run={"freeze_sha256_verified":actual,"registered_programs":len(m),
     "relationship_counts":rel_counts,"specificity_class_counts":class_counts,
     "gene_concordance":gene_rows,"pathway_concordance_rows":len(path_rows)}
(ROOT/"logs"/"phase14_specificity.log").write_text(json.dumps(run,indent=2),encoding="utf-8")
print(json.dumps(run,indent=2))

