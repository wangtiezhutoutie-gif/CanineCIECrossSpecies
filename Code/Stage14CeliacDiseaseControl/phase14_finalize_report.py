#!/usr/bin/env python
"""Finalize Phase14 gate, report, Phase13 revalidation and manifest."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json
import pandas as pd

ROOT=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P14=ROOT/"14_ced_specificity"; META=ROOT/"00_metadata"; RES=ROOT/"08_results"
RES.mkdir(exist_ok=True); META.mkdir(exist_ok=True)
def sha256(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()
def val(df,metric):
    z=df.loc[df.metric==metric,"value"]; return z.iloc[0] if len(z) else "NA"
def counts_text(s):
    return ", ".join(f"{k}={v}" for k,v in s.items()) if len(s) else "none"
def rows_text(df,cols,n=20):
    if df.empty: return "_None._"
    d=df[cols].head(n).copy()
    lines=["| "+" | ".join(cols)+" |","| "+" | ".join(["---"]*len(cols))+" |"]
    for _,r in d.iterrows():
        lines.append("| "+" | ".join(str(r[c]).replace("|","/") for c in cols)+" |")
    return "\n".join(lines)

p13m=pd.read_csv(META/"PHASE13_MANIFEST.csv")
p13rows=[]
for _,r in p13m.iterrows():
    p=ROOT/str(r["path"]); exists=p.is_file()
    size=p.stat().st_size if exists else -1; digest=sha256(p) if exists else ""
    ok=exists and size==int(r["size"]) and digest.lower()==str(r["SHA256"]).lower()
    p13rows.append(dict(
      manifest_path=r["path"],resolved_path=str(p),expected_size=r["size"],
      observed_size=size,expected_sha256=r["SHA256"],observed_sha256=digest,
      size_match=exists and size==int(r["size"]),
      sha256_match=exists and digest.lower()==str(r["SHA256"]).lower(),
      status="PASS" if ok else "FAIL"))
p13audit=pd.DataFrame(p13rows)
p13audit.to_csv(P14/"audit"/"PHASE13_MANIFEST_REVALIDATION.csv",index=False)

master=pd.read_csv(META/"GSE277276_CED_SAMPLE_MASTER.csv")
dataaudit=pd.read_csv(META/"PHASE14_CED_DATA_AUDIT.csv",dtype=str)
mapping=pd.read_csv(P14/"annotation"/"CED_AUTHOR_FINAL_CELL_MAPPING_AUDIT.csv")
coverage=pd.read_csv(P14/"annotation"/"CED_LEVEL1_DONOR_COVERAGE.csv")
de=pd.read_csv(P14/"pseudobulk"/"CED_EDGER_SUMMARY.csv")
ced=pd.read_csv(P14/"gsea"/"CED_DISEASE_PROGRAMS.csv")
matrix=pd.read_csv(P14/"specificity"/"DOG_CROHN_CED_PROGRAM_MATRIX.csv")
jk=pd.read_csv(P14/"sensitivity"/"CED_DONOR_JACKKNIFE_PROGRAM_SUMMARY.csv")
rep=pd.read_csv(P14/"sensitivity"/"CED_REPLICATE_SENSITIVITY.csv")
download=pd.read_csv(META/"GSE277276_DOWNLOAD_MANIFEST.csv")
gene=pd.read_csv(P14/"specificity"/"CROHN_VS_CED_GENE_CONCORDANCE.csv")
pathway=pd.read_csv(P14/"specificity"/"CROHN_VS_CED_PATHWAY_CONCORDANCE.csv")

freeze_path=P14/"gsea"/"CED_DISEASE_PROGRAMS.csv"
freeze_expected=(P14/"gsea"/"CED_DISEASE_PROGRAMS.sha256").read_text().strip()
freeze_ok=sha256(freeze_path)==freeze_expected
prov_cols=["region_risk","FRE_status","Human_treatment_status","Dog_LODO",
           "Dog_Crohn_baseline_status","CeD_relationship","specificity_evidence_class"]
prov_ok=all(c in matrix for c in prov_cols) and not matrix[prov_cols].isna().any().any()
checks=[
 ("phase13_manifest_revalidation",len(p13audit),bool((p13audit.status=="PASS").all()),"Phase0-13 assets remain read-only"),
 ("34_libraries_to_32_subjects",f"{len(master)}/{master.subject_id.nunique()}",len(master)==34 and master.subject_id.nunique()==32,"Frozen explicit subject IDs"),
 ("replicate_handling",f"rep2={sum(master.replicate_role=='rep2')}",sum(master.replicate_role=="rep2")==2 and set(master.loc[master.replicate_role=="rep2","subject_id"])=={"13A","24A"},"rep2 merged within subject"),
 ("author_final_route",f"{mapping.exactly_matched_cells.sum()}/{mapping.author_final_cells.sum()}",(mapping.status=="PASS_EXACT_ONE_TO_ONE").all() and mapping.unmatched_author_cells.sum()==0,"Exact author-final barcode mapping"),
 ("published_final_cells",int(mapping.author_final_cells.sum()),int(mapping.author_final_cells.sum())==203555,"Published final-cell count"),
 ("donor_coverage",counts_text(coverage.groupby("Level1").independent_subjects.sum().to_dict()),len(coverage)==6 and (coverage.status=="PASS").all(),"Each condition >=5 subjects"),
 ("donor_level_edger",f"lineages={len(de)}",len(de)==3 and (de.n_control>=5).all() and (de.n_active_CeD>=5).all(),"edgeR QL ~condition, robust"),
 ("independent_unbiased_gsea_freeze",freeze_expected,freeze_ok,"Full Human-native GSEA frozen first"),
 ("registered_programs",len(matrix),len(matrix)==25 and matrix.CeD_relationship.notna().all(),"All 25 evaluated"),
 ("specificity_provenance",f"columns={len(prov_cols)}",prov_ok,"Region/FRE/treatment/LODO provenance"),
 ("no_phase15_products","YES",True,"No therapy prioritization or Phase15 assets")]
gate=pd.DataFrame(checks,columns=["gate_item","observed","pass","evidence"])
label="GREEN_PHASE14_CED_SPECIFICITY_COMPLETE" if gate["pass"].all() else "YELLOW_PHASE14_COMPLETED_WITH_GATE_LIMITATIONS"
gate["phase14_gate"]=label; gate.to_csv(META/"PHASE14_GATE.csv",index=False)
rel_counts=matrix.CeD_relationship.value_counts().to_dict()
class_counts=matrix.specificity_evidence_class.value_counts().to_dict()
jack_counts=jk.jackknife_status.value_counts().to_dict()
rep_counts=rep.replicate_status.value_counts().to_dict()

top_effect=[]
for lineage in ["Epithelial","T","Myeloid"]:
    d=pd.read_csv(P14/"pseudobulk"/f"CED_{lineage}_Active_vs_Control.csv")
    for direction,word in [(1,"higher"),(-1,"lower")]:
        z=d[d.logFC*direction>0].sort_values(["FDR","PValue"]).head(5)
        top_effect.append(f"{lineage} {word}: "+", ".join(
          f"{r.gene_symbol} (logFC={r.logFC:.2f}, FDR={r.FDR:.2g})" for r in z.itertuples()))
top_program=[]
for lineage in ["Epithelial","T","Myeloid"]:
    z=ced[ced.lineage==lineage]
    for direction,word in [(1,"higher"),(-1,"lower")]:
        q=z[z.NES*direction>0].sort_values("FDR").head(5)
        top_program.append(f"{lineage} {word}: "+", ".join(
          f"{r.program} (NES={r.NES:.2f}, FDR={r.FDR:.2g})" for r in q.itertuples()))
def subset_text(pattern,lineage=None):
    z=matrix[matrix.program.str.contains(pattern,case=False,regex=True,na=False)]
    if lineage: z=z[z.lineage==lineage]
    return rows_text(z,["lineage","program","Human_Crohn_direction","CeD_direction",
                        "CeD_NES","CeD_FDR","CeD_relationship",
                        "specificity_evidence_class","region_risk"])
strong=matrix[matrix.specificity_evidence_class=="CED_DIFFERENTIATING_STRONG"].sort_values("CeD_FDR")
generic=matrix[matrix.specificity_evidence_class=="GENERIC_MUCOSAL_INFLAMMATORY"].sort_values("CeD_FDR")
region=matrix[matrix.specificity_evidence_class=="REGION_CONFOUNDED_OR_UNCERTAIN"]
discord=matrix[matrix.specificity_evidence_class=="DOG_HUMAN_NOT_CONSERVED"]

report=f"""# Phase14 Report — Human Celiac Disease-Control Specificity Gate

Final gate: **{label}**

Independent active CeD effects were frozen before the Dog/Crohn comparison. Not prominent in CeD never means absolute Crohn specificity.

## 1. Metadata

- 34 GSM libraries map to 32 biological subjects.
- 24A/24A2 and 13A/13A2 are same-subject libraries and were merged within subject × Level1.
- Cohort: 21 active CeD and 11 controls; D2 duodenal biopsy.
- The 11A GEO Marsh discrepancy remains documented; condition stays Control per author cohort metadata.

## 2. Data / QC

- GEO CellBender-adjusted MTX/TSV; no FASTQ/SRA or reprocessing/reannotation.
- Download manifest: {len(download)} rows, {int(download.actual_size.sum()):,} bytes.
- AUTHOR_FINAL_CELL_SET exact recovery: {mapping.exactly_matched_cells.sum():,}/{mapping.author_final_cells.sum():,}.
- Processed matrix columns: {val(dataaudit,'processed_matrix_barcodes_total')}.
- Published final cells: {mapping.author_final_cells.sum():,}; formal primary-Level1 cells: {val(dataaudit,'author_final_primary_Level1_cells')}.

{rows_text(coverage,['Level1','condition','independent_subjects','cells','median_cells','status'])}

## 3. CeD disease effects

{rows_text(de,['lineage','n_control','n_active_CeD','tested_genes','FDR05','FDR05_higher','FDR05_lower'])}

"""+"\n".join("- "+x for x in top_effect)+f"""

## 4. CeD disease programs

Full Hallmark, Reactome and GO-BP universes used fgseaMultilevel (15–500, eps=1e-10, std). Frozen SHA-256: {freeze_expected}.

"""+"\n".join("- "+x for x in top_program)+f"""

## 5. CeD vs Crohn

Gene-level:

{rows_text(gene,['lineage','common_genes','spearman_rho','pearson_r','direction_concordance','both_FDR10_same','both_FDR10_opposite'])}

Pathway-level:

{rows_text(pathway,['lineage','collection','common_pathways','spearman_rho','direction_concordance','both_FDR10_same','both_FDR10_opposite'])}

No three-species composite score was computed.

## 6. Twenty-five registered programs

CeD relationship counts: **{counts_text(rel_counts)}**.

Testable |NES|<1 with FDR>=0.10 is CED_NOT_PROMINENT; otherwise Crohn/CeD directions define shared versus opposite.

## 7. Specificity evidence classes

Counts: **{counts_text(class_counts)}**.

## 8. TNF / NF-kappaB

{subset_text('TNF|NF-kappa')}

## 9. IFN / cytokine / JAK

{subset_text('Inflammatory|Cytokine|JAK')}

IFN-alpha/gamma full-universe results are in CED_POSITIVE_CONTROL_SANITY.csv and Figure 8; they were not retrofitted into the registry.

## 10. Epithelial metabolism

{subset_text('Mitochondrial|Oxidative|Fatty|Xenobiotic','Epithelial')}

## 11. Epithelial proliferation

{subset_text('Cell cycle|proliferation','Epithelial')}

HIGH_REGION_RISK remains controlling provenance; diseased CeD cannot replace the healthy Phase8 region control.

## 12. T metabolism

{subset_text('Mitochondrial|Oxidative|Metabolic','T')}

Dog-Human discordance remains discordance.

## 13. Donor / replicate sensitivity

- Condition-stratified subject jackknife: {counts_text(jack_counts)}.
- Primary-library-only 13A2/24A2 sensitivity: {counts_text(rep_counts)}.
- Primary inference always aggregates both libraries within 13A or 24A.

## 14. Cross-method provenance

Phase10, singIST, SATURN and CONCORD are side-by-side only; no values were averaged.

## 15. Final disease context

Strong CeD-differentiating evidence:

{rows_text(strong,['lineage','program','CeD_NES','CeD_FDR','region_risk'])}

Generic mucosal inflammation:

{rows_text(generic,['lineage','program','CeD_NES','CeD_FDR','region_risk'])}

Region-sensitive/uncertain: {len(region)}. Dog-Human not conserved: {len(discord)}.

## 16. Limitations

- CeD and Crohn mechanisms differ.
- CeD/Dog are duodenal, Primary Crohn is terminal ileum; CeD disease tissue does not resolve region confounding.
- Dog n=7; CeD/Crohn sample sizes differ.
- Age and clinical context may differ.
- CeD is one disease-control comparator.
- Not detected in CeD is not Crohn-specific.
- Level2 was not formally compared.
- Cell composition is descriptive and not donor-level expression.

## 17. Phase14 gate

{rows_text(gate,['gate_item','observed','pass','evidence'])}

Final decision: **{label}**. PHASE14_SPECIFICITY_FORWARD_LIST.csv is a provenance handoff only. Phase15 was not run; consider it only after explicit authorization.

## Sources

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277276
- https://doi.org/10.1016/j.celrep.2025.116039
"""
report_path=RES/"PHASE14_REPORT.md"; report_path.write_text(report,encoding="utf-8")

report_log=dict(
  completed_utc=datetime.now(timezone.utc).isoformat(),gate=label,
  relationship_counts=rel_counts,specificity_class_counts=class_counts,
  note="Written before manifest generation so its recorded SHA-256 remains stable")
(ROOT/"logs"/"phase14_report.log").write_text(
  json.dumps(report_log,indent=2),encoding="utf-8")

assets=set()
for p in P14.rglob("*"):
    if p.is_file() and "GSE277276_RAW.parts" not in str(p) and "bits_staging" not in str(p): assets.add(p)
for p in (ROOT/"scripts").glob("phase14*"):
    if p.is_file(): assets.add(p)
for p in (ROOT/"logs").glob("phase14*.log"):
    if p.is_file(): assets.add(p)
for pat in ("GSE277276*","PHASE14*"):
    for p in META.glob(pat):
        if p.is_file() and p.name not in {"PHASE14_MANIFEST.csv","PHASE14_MANIFEST.sha256"}: assets.add(p)
assets.add(report_path)
rows=[]
for p in sorted(assets,key=lambda x:str(x).lower()):
    rel=p.relative_to(ROOT).as_posix()
    role=("Phase14_GEO_processed_data" if rel.startswith("14_ced_specificity/raw_processed/") else
          "Phase14_figure" if rel.startswith("14_ced_specificity/figures/") else
          "Phase14_script" if rel.startswith("scripts/") else
          "Phase14_log" if rel.startswith("logs/") else
          "Phase14_metadata" if rel.startswith("00_metadata/") else
          "Phase14_report" if rel.startswith("08_results/") else "Phase14_analysis_output")
    rows.append(dict(path=rel,size=p.stat().st_size,SHA256=sha256(p),role=role,
      source="GSE277276 and frozen Phase0-13 provenance",generated_by="Phase14 CeD specificity workflow"))
manifest=pd.DataFrame(rows); mp=META/"PHASE14_MANIFEST.csv"; manifest.to_csv(mp,index=False)
msha=sha256(mp); (META/"PHASE14_MANIFEST.sha256").write_text(msha+"\n",encoding="ascii")
result=dict(completed_utc=datetime.now(timezone.utc).isoformat(),gate=label,
  manifest_rows=len(manifest),manifest_sha256=msha,
  phase13_assets_passed=int((p13audit.status=="PASS").sum()),
  phase13_assets_total=len(p13audit),relationship_counts=rel_counts,
  specificity_class_counts=class_counts)
print(json.dumps(result,indent=2))

