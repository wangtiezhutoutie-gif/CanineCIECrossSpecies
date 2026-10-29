#!/usr/bin/env python
"""Freeze Phase0-14 biological-axis evidence before opening pharmacology."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib, json
import pandas as pd

ROOT=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P15=ROOT/"15_therapeutic_translation"; OUT=P15/"project_evidence"; AUD=P15/"audit"

def sha256(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()

def uniq(x): return ";".join(sorted({str(v) for v in x if pd.notna(v) and str(v)}))
def cross_state(d,default):
    if d.empty:return default
    s=d.phase10_status.fillna("").astype(str)
    if s.str.contains("OPPOSITE").any():return "OPPOSITE"
    if s.str.contains("CONCORDANT_STRONG").any():return "MIXED" if s.str.contains("LIMITED|NOT_TESTABLE").any() else "CONSISTENT"
    if s.str.contains("CONCORDANT").any():return "MIXED"
    return default

def main():
    for d in [OUT,AUD,ROOT/"logs"]:d.mkdir(parents=True,exist_ok=True)
    mp=ROOT/"00_metadata/PHASE14_MANIFEST.csv";m=pd.read_csv(mp);vr=[]
    for r in m.itertuples(index=False):
        p=ROOT/r.path;ex=p.is_file();size=p.stat().st_size if ex else -1;h=sha256(p) if ex else "";ok=ex and size==int(r.size) and h.lower()==str(r.SHA256).lower()
        vr.append(dict(path=r.path,expected_size=r.size,observed_size=size,expected_sha256=r.SHA256,observed_sha256=h,status="PASS" if ok else "FAIL"))
    vr=pd.DataFrame(vr);vr.to_csv(AUD/"PHASE14_MANIFEST_REVALIDATION.csv",index=False)
    if not vr.status.eq("PASS").all():raise RuntimeError("Phase14 manifest revalidation failed")
    src={"forward":ROOT/"14_ced_specificity/specificity/PHASE14_SPECIFICITY_FORWARD_LIST.csv","matrix":ROOT/"14_ced_specificity/specificity/DOG_CROHN_CED_PROGRAM_MATRIX.csv","phase10":ROOT/"10_cross_species_baseline/summary/PHASE10_PROGRAM_BASELINE_MATRIX.csv","phase11":ROOT/"11_singIST/summary/PHASE11_VS_PHASE10_BASELINE.csv","phase12":ROOT/"12_saturn/summary/PHASE12_VS_PHASE10_BASELINE.csv","phase13":ROOT/"13_concord/summary/CROSS_METHOD_LINEAGE_SUMMARY.csv","phase5_axis":ROOT/"06_pseudobulk/lodo/summary/THERAPEUTIC_AXIS_LODO_ROBUSTNESS.csv","phase5_program":ROOT/"06_pseudobulk/lodo/summary/PHASE5_PROGRAM_FORWARD_LIST.csv","phase6":ROOT/"07_external_validation/FRE_2026/summary/DISCOVERY_VS_FRE_INTERPRETATION_MATRIX.csv"}
    pd.DataFrame([dict(source_id=k,path=p.relative_to(ROOT).as_posix(),size=p.stat().st_size,sha256=sha256(p)) for k,p in src.items()]).to_csv(AUD/"BIOLOGICAL_AXIS_INPUT_HASH_AUDIT.csv",index=False)
    fw=pd.read_csv(src["forward"]);p5=pd.read_csv(src["phase5_axis"]);p12=pd.read_csv(src["phase12"]).set_index("lineage");p13=pd.read_csv(src["phase13"]).set_index("lineage")
    defs=[
      ("TNF/TNFR-NF-kappaB","TNF/TNFR to NF-kappaB inflammatory signaling","TNF;TNFRSF1A;TNFRSF1B;NFKB1;RELA","DRUGGABLE_AXIS","PRIMARY_AXIS_1",["TNF / NF-kappaB"],None,"STRONG","STRONG","CONSISTENT","TNF/NF-kappaB"),
      ("JAK-STAT broad","general cytokine/interleukin JAK-STAT","JAK1;JAK2;JAK3;TYK2;STAT1;STAT3;STAT5A;STAT5B","DRUGGABLE_AXIS","PRIMARY_AXIS_2",["Cytokine / interleukin","Cytokine / interleukin / JAK-STAT"],["T","Myeloid"],"STRONG","STRONG","MIXED","JAK-STAT"),
      ("IL6-JAK-STAT3","IL6-associated JAK/STAT3 arm","IL6;IL6R;IL6ST;JAK1;JAK2;STAT3","DRUGGABLE_AXIS","PRIMARY_AXIS_3_COMPONENT_NOT_IL23",["Cytokine / interleukin / JAK-STAT"],["Myeloid"],"MODERATE","MODERATE","WEAK","IL-23/STAT3"),
      ("IL23-specific","IL23A/IL12B-IL23R/TYK2/STAT3","IL23A;IL23R;IL12B;IL12RB1;TYK2;STAT3","DRUGGABLE_AXIS","PRIMARY_AXIS_3_STRICT_SPECIFIC",[],None,"LIMITED","LIMITED","NOT_ASSESSABLE","IL-23/STAT3"),
      ("IFN-JAK/STAT","type I/II interferon signaling","IFNA1;IFNB1;IFNG;IFNAR1;IFNAR2;IFNGR1;IFNGR2;JAK1;JAK2;TYK2;STAT1","BIOLOGICAL_CONTEXT_ONLY","PHASE14_FORWARD_CONTEXT",[],None,"LIMITED","LIMITED","NOT_ASSESSABLE",""),
      ("alpha4beta7-MAdCAM1","ITGA4/ITGB7 adhesion to MADCAM1","ITGA4;ITGB7;MADCAM1","DRUGGABLE_AXIS","EXPLORATORY_ONLY",[],None,"LIMITED","LIMITED","NOT_ASSESSABLE",""),
      ("Epithelial mitochondrial metabolism","respiration and oxidative phosphorylation","PROGRAM_LEVEL_NO_SINGLE_CAUSAL_TARGET","PROGRAM_MODIFIER","PHASE14_FORWARD_PROGRAM",["Mitochondrial respiration","Oxidative phosphorylation"],["Epithelial"],"STRONG","STRONG","CONSISTENT",""),
      ("Epithelial fatty-acid metabolism","fatty-acid metabolic program","PROGRAM_LEVEL_NO_SINGLE_CAUSAL_TARGET","PROGRAM_MODIFIER","PHASE14_FORWARD_PROGRAM",["Fatty acid metabolism"],["Epithelial"],"STRONG","STRONG","MIXED",""),
      ("Epithelial xenobiotic metabolism","xenobiotic metabolic program","PROGRAM_LEVEL_NO_SINGLE_CAUSAL_TARGET","PROGRAM_MODIFIER","PHASE14_FORWARD_PROGRAM",["Xenobiotic metabolism"],["Epithelial"],"STRONG","STRONG","CONSISTENT",""),
      ("Myeloid inflammatory/cytokine program","myeloid inflammatory and cytokine structure","MULTI_TARGET_PROGRAM","PROGRAM_MODIFIER","PHASE14_FORWARD_PROGRAM",["Inflammatory response","Cytokine / interleukin / JAK-STAT"],["Myeloid"],"STRONG","STRONG","CONSISTENT",""),
      ("Myeloid TLR4","TLR4 innate inflammatory program","TLR4;MYD88;IRAK4","DRUGGABLE_AXIS","PHASE14_FORWARD_PROGRAM",["TLR4"],["Myeloid"],"STRONG","STRONG","MIXED",""),
      ("T-cell metabolic divergence","T mitochondrial/oxidative/metabolic programs","PROGRAM_LEVEL_NO_SINGLE_CAUSAL_TARGET","BIOLOGICAL_CONTEXT_ONLY","PHASE14_FORWARD_PROGRAM",["Oxidative phosphorylation","Mitochondrial respiration","Mitochondrial translation","Metabolic programs"],["T"],"STRONG","STRONG","OPPOSITE",""),
      ("Stress/UPR","cellular stress, UPR and MAPK context","MULTI_TARGET_PROGRAM","PROGRAM_MODIFIER","PHASE14_FORWARD_PROGRAM",["Stress / UPR / MAPK","Stress response"],None,"STRONG","STRONG","MIXED",""),
      ("Epithelial proliferation","cell-cycle/proliferation context","PROGRAM_LEVEL_NO_SINGLE_CAUSAL_TARGET","BIOLOGICAL_CONTEXT_ONLY","PHASE14_FORWARD_PROGRAM",["Cell cycle / proliferation"],["Epithelial"],"STRONG","STRONG","CONSISTENT","")]
    targets=sorted({g for x in defs for g in x[2].split(";") if g not in {"MULTI_TARGET_PROGRAM","PROGRAM_LEVEL_NO_SINGLE_CAUSAL_TARGET"}});der=[]
    for lin in ["Epithelial","T","Myeloid"]:
        p=ROOT/f"06_pseudobulk/de_results/{lin}_CIE_vs_Healthy_edger.csv";d=pd.read_csv(p);z=d[d.gene_symbol.astype(str).str.upper().isin(targets)].copy();z.insert(0,"source_file",p.relative_to(ROOT).as_posix());der.append(z)
    de=pd.concat(der,ignore_index=True);de.to_csv(OUT/"FROZEN_PROJECT_TARGET_GENE_EFFECTS.csv",index=False);rows=[]
    for axis,sub,tgs,role,sel,patterns,lins,assoc,cell,default,p5name in defs:
        mask=fw.program.isin(patterns)
        if lins:mask&=fw.lineage.isin(lins)
        z=fw[mask];q=p5[p5.axis.eq(p5name)] if p5name else p5.iloc[0:0]
        if lins:q=q[q.lineage.isin(lins)]
        ml=sorted(set(z.lineage)|set(q.lineage));ps12=uniq(f"{x}:cosine={p12.loc[x,'SATURN_embedding_cosine']:.4f};macrogene_rho={p12.loc[x,'SATURN_macrogene_rho']:.4f}" for x in ml if x in p12.index) or "NOT_ASSESSABLE";ps13=uniq(f"{x}:cosine={p13.loc[x,'Phase13_CONCORD_cosine']:.4f};null_p={p13.loc[x,'Phase13_CONCORD_null_p']:.4g}" for x in ml if x in p13.index) or "NOT_ASSESSABLE"
        te=de[de.gene_symbol.astype(str).str.upper().isin(tgs.split(";"))];tb=uniq(f"{r.lineage}:{r.gene_symbol}:logFC={r.logFC:.3f}:FDR={r.FDR:.3g}" for r in te.itertuples()) or "NO_MATCHED_TARGET_GENE_IN_FROZEN_DOG_DE"
        rows.append(dict(axis=axis,subaxis=sub,candidate_targets=tgs,axis_role=role,selection_source=sel,matched_phase14_program_rows=len(z),matched_programs=uniq(z.program),dog_cie_association_project_state=assoc,dog_lineages=uniq(z.lineage) or uniq(q.lineage),dog_direction=uniq(z.dog_direction) or uniq(q.full_direction),dog_lodo_stability=uniq(z.dog_lodo) or uniq(q.stability_class),target_gene_effect_basis=tb,human_crohn_direction=uniq(z.crohn_direction) or "NOT_ASSESSABLE",ibdverse_support=uniq(z.ibdverse) or "NOT_ASSESSABLE",phase10_status=uniq(z.phase10_status) or "NOT_ASSESSABLE",phase11_status="EXPLORATORY_WEAK_REFERENCE;"+(uniq(z.phase11_status) or "NOT_ASSESSABLE"),phase12_status="EXPLORATORY_FUNCTIONAL_STRUCTURE;"+ps12,phase13_status=ps13,cross_species_project_state="NOT_ASSESSABLE" if axis=="IL23-specific" else cross_state(z,default),ced_relationship=uniq(z.ced_relationship) or "NOT_ASSESSABLE",disease_specificity=uniq(z.specificity_evidence_class) or "NOT_ASSESSABLE",fre_status=uniq(z.fre_status) or "NOT_ASSESSABLE",region_risk=uniq(z.region_risk) or "NOT_ASSESSABLE",human_treatment_sensitivity=uniq(z.human_treatment) or "NOT_ASSESSABLE",cell_type_concordance_project_state=cell,project_evidence_basis=f"Phase14 programs={uniq(z.program) or 'none'};Phase5 axis rows={len(q)};target effects retained separately",major_project_caution="GENERAL_STAT3_JAK_IS_NOT_IL23_SPECIFIC" if axis=="IL23-specific" else "HUMAN_CLINICAL_AND_PHARMACOLOGY_NOT_OPENED_AT_FREEZE",literature_opened="NO"))
    out=pd.DataFrame(rows);t=datetime.now(timezone.utc).isoformat();out["project_evidence_frozen_utc"]=t;p=OUT/"BIOLOGICAL_AXIS_EVIDENCE_MATRIX.csv";out.to_csv(p,index=False);h=sha256(p);(OUT/"BIOLOGICAL_AXIS_EVIDENCE_MATRIX.sha256").write_text(h+"\n",encoding="ascii")
    pd.DataFrame([dict(artifact=p.relative_to(ROOT).as_posix(),sha256=h,rows=len(out),frozen_utc=t,human_clinical_or_pharmacology_opened="NO",gate="PASS_PROJECT_BIOLOGICAL_EVIDENCE_FROZEN_FIRST")]).to_csv(OUT/"BIOLOGICAL_AXIS_FREEZE_AUDIT.csv",index=False)
    log=dict(completed_utc=t,phase14_manifest_rows=len(m),phase14_manifest_pass=int(vr.status.eq("PASS").sum()),biological_axes=len(out),freeze_sha256=h,human_clinical_or_pharmacology_opened=False);(ROOT/"logs/phase15_project_evidence.log").write_text(json.dumps(log,indent=2),encoding="utf-8");print(json.dumps(log,indent=2))
if __name__=="__main__":main()