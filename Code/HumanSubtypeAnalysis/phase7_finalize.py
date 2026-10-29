from pathlib import Path
import hashlib,numpy as np,pandas as pd
R=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"]);O=R/"WORKING"/"24_EPITHELIAL_MYELOID_CROSS_SPECIES";W=O/"work"
def out(rel,x):p=O/rel;p.parent.mkdir(parents=True,exist_ok=True);x.to_csv(p,index=False)
D=pd.read_csv(W/"phase7_native_program_results.tsv",sep="\t").rename(columns={"pval":"P"});DL=pd.read_csv(W/"phase7_native_program_lodo.tsv",sep="\t").rename(columns={"pval":"P"})
H=pd.read_csv(W/"phase7_human_results.tsv",sep="\t").rename(columns={"pval":"P"});HL=pd.read_csv(W/"phase7_human_lodo.tsv",sep="\t").rename(columns={"pval":"P"});HT=pd.read_csv(W/"phase7_human_treatment.tsv",sep="\t").rename(columns={"pval":"P"})
cd=["lineage","subtype","frozen_program_id","program","collection","dog_pathway_id","healthy_donors","cie_donors","tested_genes","mapped_ranked_genes","ES","NES","P","BH_q","size","leading_edge_n","direction","model_status","planned_family"]
out("02_CANINE_SUBTYPE_ANALYSIS/04_EPITHELIAL_SUBTYPE_PROGRAM_RESULTS.csv",D[D.lineage=="Epithelial"][cd]);out("02_CANINE_SUBTYPE_ANALYSIS/05_MYELOID_SUBTYPE_PROGRAM_RESULTS.csv",D[D.lineage=="Myeloid"][cd]);out("02_CANINE_SUBTYPE_ANALYSIS/06_CANINE_PROGRAM_LODO.csv",DL)
S=H[H.analysis.str.contains("__Subtype__",regex=False)].copy();S["author_subtype"]=S.analysis.str.split("__Subtype__",regex=False).str[-1];ch=["lineage","author_subtype","frozen_program_id","program","pathway_id","healthy_donors","crohn_donors","tested_genes","mapped_ranked_genes","ES","NES","P","BH_q","size","leading_edge_n","direction","model_status"]
out("04_HUMAN_SUBTYPE_ANALYSIS/02_HUMAN_EPITHELIAL_SUBTYPE_PROGRAM_RESULTS.csv",S[S.lineage=="Epithelial"][ch]);out("04_HUMAN_SUBTYPE_ANALYSIS/03_HUMAN_MYELOID_SUBTYPE_PROGRAM_RESULTS.csv",S[S.lineage=="Myeloid"][ch])
HL["author_subtype"]=HL.analysis.str.split("__Subtype__",regex=False).str[-1];HL["sensitivity_type"]="LODO";HL["omitted_or_filter"]=HL.omitted_donor
HT["author_subtype"]=HT.analysis.str.split("__Subtype__",regex=False).str[-1];HT["sensitivity_type"]="UNTREATED_CROHN_ONLY";HT["omitted_or_filter"]="EXCLUDE_BIOLOGIC_TREATED_CROHN";HT["full_NES"]=np.nan;HT["direction_match_full"]="NOT_ESTIMABLE"
cs=["lineage","author_subtype","frozen_program_id","program","sensitivity_type","omitted_or_filter","healthy_donors","crohn_donors","NES","P","direction","model_status","full_NES","direction_match_full"];out("04_HUMAN_SUBTYPE_ANALYSIS/04_HUMAN_PROGRAM_LODO_AND_TREATMENT_SENSITIVITY.csv",pd.concat([HL[cs],HT[cs]]))
PM=pd.read_csv(R/"supplementary_robustness"/"human_duodenum_ti_region_sensitivity"/"04_program_integration"/"DOG_HUMAN_REGION_25_PROGRAM_MATRIX.csv");B=H[~H.analysis.str.contains("__Subtype__",regex=False)];rows=[]
for _,p in D[["lineage","frozen_program_id","program"]].drop_duplicates().iterrows():
 q=PM[(PM.lineage==p.lineage)&(PM.frozen_program_id==p.frozen_program_id)]
 if len(q):q=q.iloc[0];rows.append([p.lineage,p.frozen_program_id,p.program,"A_ORIGINAL_27_PARTICIPANT_BROAD_LINEAGE",16,11,q.get("human_crohn_ti_NES",np.nan),np.nan,q.get("human_crohn_ti_FDR",np.nan),q.get("human_crohn_ti_direction",""),"FROZEN_READ_ONLY"])
 for key,label in [("__B_RECOVERED_PARTICIPANTS_ALL_BROAD","B_RECOVERED_PARTICIPANTS_ALL_BROAD_CELLS"),("__C_VERIFIED_AUTHOR_LABELLED_LINEAGE","C_RECOVERED_PARTICIPANTS_AUTHOR_LABELLED_LINEAGE")]:
  q=B[(B.lineage==p.lineage)&(B.frozen_program_id==p.frozen_program_id)&B.analysis.str.contains(key,regex=False)]
  if len(q):q=q.iloc[0];rows.append([p.lineage,p.frozen_program_id,p.program,label,q.healthy_donors,q.crohn_donors,q.NES,q.P,np.nan,q.direction,"ESTIMABLE"])
out("04_HUMAN_SUBTYPE_ANALYSIS/01_HUMAN_BROAD_LINEAGE_SELECTION_SENSITIVITY.csv",pd.DataFrame(rows,columns=["lineage","frozen_program_id","program","scenario","healthy_participants","crohn_participants","NES","P","FDR_or_BH_q","direction","analysis_status"]))
def match(lin,A):
 z=[]
 for _,a in A.iterrows():
  for fid in D[D.lineage==lin].frozen_program_id.unique():
   d=D[(D.lineage==lin)&(D.subtype==a.dog_subtype)&(D.frozen_program_id==fid)];h=S[(S.lineage==lin)&(S.author_subtype==a.human_subtype)&(S.frozen_program_id==fid)];ok=a.matching_classification=="SUPPORTED_EXPLORATORY_MATCH" and len(d) and len(h)
   if ok:
    d=d.iloc[0];h=h.iloc[0];cl="DIRECTIONALLY_CONCORDANT" if np.sign(d.NES)==np.sign(h.NES) else "DIRECTIONALLY_DISCORDANT";dl=DL[(DL.lineage==lin)&(DL.subtype==a.dog_subtype)&(DL.frozen_program_id==fid)&DL.NES.notna()];hl=HL[(HL.lineage==lin)&HL.analysis.str.endswith("__Subtype__"+a.human_subtype)&(HL.frozen_program_id==fid)&HL.NES.notna()];z.append([lin,a.dog_subtype,a.human_subtype,fid,d.program,a.matching_classification,d.NES,d.P,d.BH_q,h.NES,h.P,h.BH_q,cl,f"{(dl.direction_match_full=='YES').sum()}/{len(dl)}",f"{(hl.direction_match_full=='YES').sum()}/{len(hl)}","Exploratory; no species-by-disease interaction"])
   else:z.append([lin,a.dog_subtype,a.human_subtype,fid,D[(D.lineage==lin)&(D.frozen_program_id==fid)].program.iloc[0],a.matching_classification,"","","","","","","NOT_COMPARABLE","","","Frozen identity/region gate failed"])
 return pd.DataFrame(z,columns=["lineage","dog_subtype","human_subtype","frozen_program_id","program","matching_classification","dog_NES","dog_P","dog_BH_q","human_NES","human_P","human_BH_q","effect_classification","dog_lodo_direction_stability","human_lodo_direction_stability","boundary"])
ME=pd.read_csv(O/"05_CROSS_SPECIES_MATCHING"/"01_EPITHELIAL_MATCHING_DECISIONS_FROZEN.csv");MM=pd.read_csv(O/"05_CROSS_SPECIES_MATCHING"/"02_MYELOID_MATCHING_DECISIONS_FROZEN.csv");XE=match("Epithelial",ME);XM=match("Myeloid",MM);out("05_CROSS_SPECIES_MATCHING/05_EPITHELIAL_MATCHED_PROGRAM_EFFECTS.csv",XE);out("05_CROSS_SPECIES_MATCHING/06_MYELOID_MATCHED_PROGRAM_EFFECTS.csv",XM)
sigD=(D.BH_q<.05).sum();sigH=(S.BH_q<.05).sum();mac=D[(D.subtype=="Macrophage")&(D.program=="TNF / NF-kappaB")].iloc[0]
report=f"""# Phase 7 integrated epithelial and myeloid subtype report

## Gate and recovery
Canine author RDS files were not checksum-verified because the targeted Zenodo request failed at DNS resolution. Canine analyses therefore use frozen Level-2 labels and carry AUTHOR_RDS_NOT_VERIFIED. Four epithelial states (BEST4+ enterocyte, Enterocyte 1, Enterocyte 3, Goblet) and Macrophage met the 3+3 donor gate.

Exact verified-crosswalk plus participant-barcode joins recovered Human author labels for 13 Healthy and 6 Crohn participants. Early Enterocyte, Enterocyte, Ileum Goblet, TA and Macrophage met the 5+5 participant gate.

## Canine results
{sigD} of 15 state-program tests passed lineage-specific BH correction. Enterocyte 1 and Enterocyte 3 were lower for mitochondrial respiration, oxidative phosphorylation and xenobiotic metabolism. Goblet cells were higher for the two respiratory programs. BEST4+ directions were lower but did not pass the epithelial family. Macrophage TNF/NF-kappaB was higher (NES={mac.NES:.3f}); inflammatory response was also higher, while TLR4 did not pass the three-test family. Neutrophils were excluded because about 95% came from one CIE dog.

## Human results
{sigH} of 15 state-program tests passed lineage-specific BH correction. Early Enterocyte was lower for respiratory programs but higher for xenobiotic metabolism; Enterocyte was lower for all three; Ileum Goblet was higher for all three; TA was lower for respiratory programs but higher for xenobiotic metabolism. Human Macrophage TNF/NF-kappaB, inflammatory response and TLR4 were all lower and retained direction in all 19 valid LODO fits. Untreated-Crohn sensitivity was not estimable because only 2 of 6 recovered Crohn participants were untreated.

The recovered-participant broad-lineage B/C analyses were selection-sensitive. Myeloid B/C estimates were lower for all three inflammatory programs, opposite the frozen 27-participant broad-lineage directions. This blocks a clean species-biology interpretation.

## Cross-species matching
Dog Goblet versus Human Ileum Goblet was the only epithelial pair authorized for formal exploratory comparison. Mitochondrial respiration and oxidative phosphorylation were concordantly higher in disease; xenobiotic metabolism was discordant. The result remains anatomy-limited because canine samples were duodenal and Human samples terminal-ileal.

Dog versus Human Macrophage identity was supported at an exploratory level, but all three inflammatory program directions were discordant: higher in canine CIE and lower in the recovered Human subset. This is a boundary finding, not proof of an opposite conserved mechanism, because Human participant selection and untreated sensitivity remain limiting.

## Final decisions
- Canine epithelial: SUCCESS_EXPLORATORY.
- Human epithelial: SUCCESS_EXPLORATORY_SELECTION_SENSITIVE.
- Cross-species epithelial: PARTIAL_EXPLORATORY_MATCH_WITH_REGION_CAUTION.
- Canine myeloid: SUCCESS_EXPLORATORY_MACROPHAGE_ONLY.
- Human myeloid: SUCCESS_EXPLORATORY_SELECTION_AND_TREATMENT_LIMITED.
- Cross-species myeloid: DIRECTIONALLY_DISCORDANT_EXPLORATORY_BOUNDARY.

## Manuscript value
Use as Supplementary boundary evidence. The analysis supports within-lineage heterogeneity and shows that broad correspondence cannot be assigned automatically to one state. No new therapeutic node should be prioritized: Phase 7 tested predefined programs and did not perform target discovery.
"""
(O/"06_INTEGRATED_REPORT"/"01_PHASE7_INTEGRATED_REPORT.md").write_text(report,encoding="utf-8");(O/"06_INTEGRATED_REPORT"/"02_MANUSCRIPT_RECOMMENDATION.md").write_text("# Manuscript recommendation\n\nSupplementary-only, boundary-focused reporting. Goblet respiratory concordance is exploratory and anatomy-limited. Macrophage inflammatory discordance is participant-selection and treatment-sensitivity limited, not demonstrated species divergence. Do not create new therapeutic claims.\n",encoding="utf-8")
try:
 import matplotlib.pyplot as plt
 from matplotlib.colors import TwoSlopeNorm
 fig,ax=plt.subplots(1,2,figsize=(10,4.8),constrained_layout=True)
 for aa,lin in zip(ax,["Epithelial","Myeloid"]):
  d=D[D.lineage==lin].pivot(index="subtype",columns="program",values="NES");h=S[S.lineage==lin].pivot(index="author_subtype",columns="program",values="NES");d.index=["Dog "+x for x in d.index];h.index=["Human "+x for x in h.index];m=pd.concat([d,h]);v=max(3.7,np.nanmax(abs(m.values)));im=aa.imshow(m.values,aspect="auto",cmap="RdBu_r",norm=TwoSlopeNorm(vmin=-v,vcenter=0,vmax=v));aa.set_xticks(range(m.shape[1]),m.columns,rotation=35,ha="right",fontsize=8);aa.set_yticks(range(m.shape[0]),m.index,fontsize=8);aa.set_title(lin,loc="left",fontweight="bold");aa.axhline(len(d)-.5,color="#555",lw=.8)
 fig.colorbar(im,ax=ax,label="NES",shrink=.8);fig.suptitle("Predefined disease-program effects within eligible cell states",fontweight="bold");fig.savefig(O/"06_INTEGRATED_REPORT"/"MINIMAL_VALIDATED_SUBTYPE_FIGURE.pdf");fig.savefig(O/"06_INTEGRATED_REPORT"/"MINIMAL_VALIDATED_SUBTYPE_FIGURE.png",dpi=300);plt.close(fig)
except Exception as e:(O/"logs"/"figure_error.txt").write_text(str(e))
def sha(p):
 if p.stat().st_size>1500000000:return "NOT_RECOMPUTED_LARGE_FILE"
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
 return h.hexdigest()
inputs=[R/"05_annotations"/"dog_GSE254005_annotated.h5ad",R/"07_human_crohn"/"GSE266616"/"h5ad"/"human_GSE266616_TI_annotated.h5ad",R/"WORKING"/"20_HUMAN_TCELL_LABEL_FREEZE"/"source"/"scRNA_cell_metadata.csv",R/"WORKING"/"20_HUMAN_TCELL_LABEL_FREEZE"/"03_VERIFIED_GSM_SAMPLE_CROSSWALK.csv"]
rep=["# Phase 7 reproducibility log","","- Seed: 20261016","- edgeR: TMM, filterByExpr, robust dispersion, quasi-likelihood","- fgseaMultilevel: minSize 15, maxSize 500, eps 1e-10","- BH families: 12 epithelial and 3 myeloid tests, separately per species","- Matching freeze SHA256: "+(W/"matching_decisions_frozen.sha256").read_text().strip(),"","## Inputs"]+[f"- {p} | bytes={p.stat().st_size} | sha256={sha(p)}" for p in inputs]+["","## Non-estimable components","- Canine author RDS verification: DNS failure.","- Human untreated-Crohn sensitivity: only 2 untreated recovered Crohn participants.","- Independent canine single-cell validation: no dataset identified.","- GSE290695: not downloaded or analyzed."]
(O/"06_INTEGRATED_REPORT"/"03_REPRODUCIBILITY_LOG.md").write_text("\n".join(rep)+"\n",encoding="utf-8")
(O/"logs"/"PHASE7_STATUS.txt").write_text("CANINE_EPITHELIAL=SUCCESS_EXPLORATORY\nHUMAN_EPITHELIAL=SUCCESS_EXPLORATORY_SELECTION_SENSITIVE\nCROSS_SPECIES_EPITHELIAL=PARTIAL_EXPLORATORY_MATCH_WITH_REGION_CAUTION\nCANINE_MYELOID=SUCCESS_EXPLORATORY_MACROPHAGE_ONLY\nHUMAN_MYELOID=SUCCESS_EXPLORATORY_SELECTION_AND_TREATMENT_LIMITED\nCROSS_SPECIES_MYELOID=DIRECTIONALLY_DISCORDANT_EXPLORATORY_BOUNDARY\nFINAL_STATUS=PHASE7_COMPLETE_WITH_DOCUMENTED_LIMITATIONS\n")
files=[p for p in O.rglob("*") if p.is_file() and "work" not in p.parts and p.name!="OUTPUT_SHA256.tsv"];pd.DataFrame([[str(p.relative_to(O)),p.stat().st_size,sha(p)] for p in files],columns=["relative_path","bytes","sha256"]).to_csv(O/"logs"/"OUTPUT_SHA256.tsv",sep="\t",index=False)
print("FINALIZED",sigD,sigH,XE.effect_classification.value_counts().to_dict(),XM.effect_classification.value_counts().to_dict())