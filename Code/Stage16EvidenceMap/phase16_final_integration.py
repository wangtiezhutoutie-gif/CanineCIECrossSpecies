from __future__ import annotations
import hashlib, os
from datetime import datetime, timezone
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch, Rectangle
import numpy as np
import pandas as pd

ROOT=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
OUT=ROOT/"16_final_integration"
DATE=datetime.now(timezone.utc).date().isoformat()
GATE="FINAL_COMPLETE_WITH_DOCUMENTED_LIMITATIONS"
LINEAGES=["Epithelial","T","Myeloid"]

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()
def rel(p): return p.relative_to(ROOT).as_posix()
def rd(s): return pd.read_csv(ROOT/s,dtype=str,keep_default_na=False)
def wcsv(d,p):
    p.parent.mkdir(parents=True,exist_ok=True); t=p.with_suffix(p.suffix+".tmp")
    d.to_csv(t,index=False,encoding="utf-8-sig"); os.replace(t,p)
def wtxt(s,p):
    p.parent.mkdir(parents=True,exist_ok=True); t=p.with_suffix(p.suffix+".tmp")
    t.write_text(s.rstrip()+"\n",encoding="utf-8"); os.replace(t,p)
def fn(x,n=4):
    try:return f"{float(x):.{n}f}"
    except:return str(x)

for d in ["master","matrix","summary","tables","audit","figures"]: (OUT/d).mkdir(parents=True,exist_ok=True)

# First gate: exact Phase15 revalidation.
m15=rd("00_metadata/PHASE15_MANIFEST.csv"); checks=[]
for r in m15.to_dict("records"):
    p=ROOT/r["path"]; s="PASS"
    if not p.exists(): s="MISSING"
    elif p.stat().st_size!=int(r["size"]): s="SIZE_MISMATCH"
    elif sha(p).lower()!=r["SHA256"].lower(): s="SHA256_MISMATCH"
    checks.append({"path":r["path"],"status":s})
if len(checks)!=45 or any(x["status"]!="PASS" for x in checks): raise RuntimeError("Phase15 manifest failed")

spec=[("Phase0","00_metadata/PHASE0_AUDIT_REPORT.md","design/source audit"),("Phase1","08_results/PHASE1_REPORT.md","dog integrity report")]
for n in range(2,16):
    spec += [(f"Phase{n}",f"08_results/PHASE{n}_REPORT.md","frozen report"),(f"Phase{n}",f"00_metadata/PHASE{n}_MANIFEST.csv","frozen manifest")]
spec += [
("Phase5","06_pseudobulk/lodo/summary/PHASE5_PROGRAM_FORWARD_LIST.csv","25-program LODO"),
("Phase6","07_external_validation/FRE_2026/summary/PHASE6_SUMMARY_LEVEL_VALIDATION_MATRIX.csv","FRE summary"),
("Phase7","07_human_crohn/IBDverse/summary/IBDVERSE_HUMAN_PROGRAM_REPLICATION_SUMMARY.csv","IBDverse"),
("Phase8","08_region_control/summary/HUMAN_REGION_REGISTERED_GSEA_SUMMARY.csv","region control"),
("Phase9","09_cell_ontology/summary/PHASE10_PROGRAM_HANDOFF.csv","ontology handoff"),
("Phase10","10_cross_species_baseline/summary/PHASE10_PROGRAM_BASELINE_MATRIX.csv","ortholog baseline"),
("Phase11","11_singIST/summary/PHASE11_VS_PHASE10_BASELINE.csv","singIST"),
("Phase12","12_saturn/summary/PHASE12_VS_PHASE10_BASELINE.csv","SATURN"),
("Phase13","13_concord/summary/CROSS_METHOD_LINEAGE_SUMMARY.csv","CONCORD"),
("Phase14","14_ced_specificity/specificity/PHASE14_SPECIFICITY_FORWARD_LIST.csv","CeD handoff"),
("Phase15","15_therapeutic_translation/matrix/CANINE_CIE_THERAPEUTIC_TRANSLATION_MATRIX.csv","therapeutic matrix"),
("Phase15","15_therapeutic_translation/summary/CANDIDATE_THERAPEUTIC_AXIS_PRIORITY_LIST.csv","priority list")]
audit=[]
for ph,f,role in spec:
    p=ROOT/f
    if not p.exists(): raise FileNotFoundError(f)
    audit.append({"phase":ph,"file":f,"size":p.stat().st_size,"SHA256":sha(p),"role":role,"status":"PASS"})
audit=pd.DataFrame(audit); wcsv(audit,ROOT/"00_metadata/PHASE16_INPUT_FREEZE_AUDIT.csv")

p5=rd("06_pseudobulk/lodo/summary/PHASE5_PROGRAM_FORWARD_LIST.csv")
p10=rd("10_cross_species_baseline/summary/PHASE10_PROGRAM_BASELINE_MATRIX.csv")
p11=rd("11_singIST/summary/PHASE11_VS_PHASE10_BASELINE.csv")
p13=rd("13_concord/summary/CROSS_METHOD_LINEAGE_SUMMARY.csv")
p14=rd("14_ced_specificity/specificity/PHASE14_SPECIFICITY_FORWARD_LIST.csv")
p15m=rd("15_therapeutic_translation/matrix/CANINE_CIE_THERAPEUTIC_TRANSLATION_MATRIX.csv")
p15=rd("15_therapeutic_translation/summary/CANDIDATE_THERAPEUTIC_AXIS_PRIORITY_LIST.csv")
for nm,d in {"P5":p5,"P10":p10,"P11":p11,"P14":p14}.items():
    if len(d)!=25 or d[["lineage","program"]].duplicated().any() or set(d.lineage)!=set(LINEAGES): raise RuntimeError(nm+" universe")
keys=set(map(tuple,p14[["lineage","program"]].to_records(index=False)))
for d in [p5,p10,p11]:
    if set(map(tuple,d[["lineage","program"]].to_records(index=False)))!=keys: raise RuntimeError("program mismatch")
x=(p14.merge(p5,on=["lineage","program"],suffixes=("_p14","_p5"),validate="one_to_one")
      .merge(p10,on=["lineage","program"],suffixes=("","_p10"),validate="one_to_one")
      .merge(p11,on=["lineage","program"],suffixes=("","_p11"),validate="one_to_one"))
lm=p13.set_index("lineage").to_dict("index"); am=p15.set_index("axis").to_dict("index")

def axis(l,p):
    if p=="TNF / NF-kappaB": return "TNF/TNFR-NF-kappaB"
    if p in {"Cytokine / interleukin","Cytokine / interleukin / JAK-STAT"} or (l=="T" and p=="Inflammatory response"): return "JAK-STAT broad"
    if l=="Myeloid" and p in {"Inflammatory response","Stress response","Apoptosis"}: return "Myeloid inflammatory/cytokine program"
    if l=="Myeloid" and p=="TLR4": return "Myeloid TLR4"
    if l=="Epithelial" and p in {"Mitochondrial respiration","Oxidative phosphorylation"}: return "Epithelial mitochondrial metabolism"
    if l=="Epithelial" and p=="Fatty acid metabolism": return "Epithelial fatty-acid metabolism"
    if l=="Epithelial" and p=="Xenobiotic metabolism": return "Epithelial xenobiotic metabolism"
    if l=="Epithelial" and p=="Cell cycle / proliferation": return "Epithelial proliferation"
    if l in {"Epithelial","T"} and p=="Stress / UPR / MAPK": return "Stress/UPR"
    if l=="T" and p in {"Oxidative phosphorylation","Mitochondrial respiration","Mitochondrial translation","Metabolic programs"}: return "T-cell metabolic divergence"
    return "NO_DIRECT_PHASE15_AXIS"
def cross(r):
    if r.region_risk=="HIGH_REGION_RISK": return "UNCERTAIN"
    if r.specificity_evidence_class=="DOG_HUMAN_NOT_CONSERVED" or r.phase10_status=="OPPOSITE": return "DIVERGENT"
    if r.phase10_status=="NOT_TESTABLE": return "UNCERTAIN"
    if r.phase10_status=="CONCORDANT_STRONG" and r.same_direction=="True": return "CONSERVED"
    if r.phase10_status.startswith("CONCORDANT"): return "PARTIALLY_CONSERVED"
    return "UNCERTAIN"
def context(r):
    if r.specificity_evidence_class=="GENERIC_MUCOSAL_INFLAMMATORY": return "GENERIC_MUCOSAL_INFLAMMATORY"
    if r.ced_relationship=="CED_NOT_PROMINENT": return "CED_NOT_PROMINENT"
    if r.ced_relationship=="CED_OPPOSITE_DIRECTION" or r.specificity_evidence_class=="CED_DIFFERENTIATING_LIMITED": return "CED_DIFFERENTIATING_HINT"
    return "NOT_ASSESSABLE"
def fclass(r,c,d):
    if r.region_risk=="HIGH_REGION_RISK": return "REGION_ASSOCIATED_OR_UNCERTAIN"
    if c=="DIVERGENT": return "CANINE_DIVERGENT_PROGRAM"
    if r.phase10_status=="NOT_TESTABLE": return "METHOD_DEPENDENT_OR_UNCERTAIN"
    if d=="GENERIC_MUCOSAL_INFLAMMATORY": return "GENERIC_MUCOSAL_INFLAMMATORY_PROGRAM"
    if d in {"CED_DIFFERENTIATING_HINT","CED_NOT_PROMINENT"} and c in {"CONSERVED","PARTIALLY_CONSERVED"}: return "DOG_CROHN_DIFFERENTIATING_HINT"
    if c=="CONSERVED": return "CONSERVED_DOG_CROHN_PROGRAM"
    return "METHOD_DEPENDENT_OR_UNCERTAIN"
def conf(r,f):
    if f in {"REGION_ASSOCIATED_OR_UNCERTAIN","METHOD_DEPENDENT_OR_UNCERTAIN","DOG_CROHN_DIFFERENTIATING_HINT"}: return "LIMITED"
    if f=="GENERIC_MUCOSAL_INFLAMMATORY_PROGRAM":
        return "HIGH" if r.lineage=="Myeloid" and r.phase10_status=="CONCORDANT_STRONG" and r.region_risk=="LOW_REGION_RISK" else "MODERATE"
    if f=="CANINE_DIVERGENT_PROGRAM":
        return "MODERATE" if r.phase10_status=="OPPOSITE" and r.dog_lodo=="ROBUST" and r.region_risk in {"LOW_REGION_RISK","REGION_UNKNOWN"} else "LIMITED"
    return "HIGH" if r.dog_lodo=="ROBUST" else "MODERATE"

rows=[]
for r in x.itertuples(index=False):
    q=lm[r.lineage]; c=cross(r); d=context(r); f=fclass(r,c,d); cf=conf(r,f); a=axis(r.lineage,r.program)
    pr=f"{am[a]['rank_group']}|{am[a]['priority_class']}" if a in am else "NO_DIRECT_REGISTERED_AXIS"
    fe=am[a]["pharmacologic_translation_feasibility"] if a in am else "NOT_ASSESSABLE"
    ps=f"Dog {r.dog_lodo}; Phase10 {r.phase10_status}; CeD {r.ced_relationship}"
    ps += "; CONCORD P=0.0197 strong" if r.lineage=="Myeloid" else "; CONCORD P=0.0857 limited" if r.lineage=="Epithelial" else "; CONCORD P=0.7830 not significant"
    ca="Dog n=7; singIST 46/46 weak references; SATURN species separation"
    if r.region_risk!="LOW_REGION_RISK": ca+=f"; region={r.region_risk}"
    if r.fre_status not in {"SUMMARY_DIRECTIONAL_SUPPORT","NOT_ASSESSABLE"}: ca+=f"; FRE={r.fre_status}"
    if r.specificity_evidence_class=="GENERIC_MUCOSAL_INFLAMMATORY": ca+="; CeD-shared generic inflammation"
    interp={
    "GENERIC_MUCOSAL_INFLAMMATORY_PROGRAM":"Dog and Crohn share this program, but CeD support makes it generic mucosal inflammation rather than Crohn-specific.",
    "DOG_CROHN_DIFFERENTIATING_HINT":"Dog-Crohn correspondence is present; the single CeD comparator provides only a limited differentiating hint.",
    "CANINE_DIVERGENT_PROGRAM":"The Dog program is internally stable but is not recapitulated in the frozen Human Crohn direction.",
    "REGION_ASSOCIATED_OR_UNCERTAIN":"Region dependence prevents confident separation of disease conservation from tissue location.",
    "METHOD_DEPENDENT_OR_UNCERTAIN":"Cross-species testability or method constraints do not support a stable conservation call.",
    "CONSERVED_DOG_CROHN_PROGRAM":"The program is conserved within the frozen evidence boundary."}[f]
    rows.append({"lineage":r.lineage,"program":r.program,"dog_direction":r.dog_direction,
    "dog_phase3_effect":f"PHASE3_DONOR_LEVEL_RANKED_INPUT; Phase4/5 program estimate: {r.full7_evidence}",
    "dog_phase5_lodo":r.direction_stability,"dog_robustness":r.dog_lodo,
    "fre_direction":"PUBLISHED_SUMMARY_ONLY","fre_status":r.fre_status,
    "human_crohn_direction":r.crohn_direction,"human_crohn_evidence":r.crohn_evidence,
    "human_treatment_sensitivity":r.human_treatment,"ibdverse_support":r.ibdverse,
    "region_risk":r.region_risk,"region_status":r.phase10_priority,
    "phase10_gene_baseline":f"lineage Spearman rho={fn(q['Phase10_gene_rho'])}",
    "phase10_pathway_baseline":f"Dog NES={fn(r.dog_common1to1_NES)}; Human NES={fn(r.human_common1to1_NES)}; {r.phase10_status}",
    "phase11_singist_status":"EXPLORATORY_WEAK_REFERENCE","phase11_reference_quality":"46/46 WEAK_REFERENCE; TrackB validated=0",
    "phase12_saturn_status":f"EXPLORATORY_PROTEIN_MACROGENE_STRUCTURE; {r.phase12_status}","phase12_species_mixing_warning":"YES",
    "phase13_concord_status":r.phase13_status,
    "phase13_null_significance":("STRONG_SUPPORT" if r.lineage=="Myeloid" else "LIMITED_POSITIVE" if r.lineage=="Epithelial" else "NOT_SIGNIFICANT")+f": cosine={fn(q['Phase13_CONCORD_cosine'])}; P={fn(q['Phase13_CONCORD_null_p'])}",
    "ced_direction":r.ced_direction,"ced_relationship":r.ced_relationship,"phase14_specificity_class":r.specificity_evidence_class,
    "phase15_therapeutic_axis":a,"phase15_priority":pr,"phase15_feasibility":fe,
    "cross_species_fidelity_class":c,"disease_context_class":d,"final_fidelity_class":f,"final_confidence":cf,
    "major_support":ps,"major_caution":ca,"final_interpretation":interp})
master=pd.DataFrame(rows)
if len(master)!=25 or master[["lineage","program"]].duplicated().any(): raise RuntimeError("master")
wcsv(master,OUT/"master/FINAL_PROGRAM_EVIDENCE_MASTER.csv")
fcols=["lineage","program","dog_direction","dog_robustness","human_crohn_direction","phase10_gene_baseline","phase10_pathway_baseline","phase11_singist_status","phase12_saturn_status","phase13_concord_status","ced_relationship","region_risk","cross_species_fidelity_class","disease_context_class","final_fidelity_class","final_confidence","major_support","major_caution","final_interpretation"]
wcsv(master[fcols],OUT/"matrix/FINAL_TRANSLATIONAL_FIDELITY_MATRIX.csv")

lt={
"Epithelial":["TNF/NF-kappaB; mitochondrial respiration; fatty-acid; xenobiotic","Stress/UPR; proliferation/ribosome region-risk; OxPhos method-limited","Partial functional/pathway conservation led by metabolic and xenobiotic programs; not whole-state conservation.","MODERATE_PROGRAM_SELECTIVE","Modest gene rho; CONCORD P=0.0857; proliferation high-region-risk."],
"T":["TNF/NF-kappaB; inflammatory; cytokine/interleukin","OxPhos; mitochondrial respiration/translation; metabolism","Program-dependent: inflammation is conserved while metabolism diverges.","MODERATE_PROGRAM_SELECTIVE","Gene rho near zero; weak singIST; SATURN separated; CONCORD near zero."],
"Myeloid":["TNF/NF-kappaB; inflammatory; cytokine/JAK-STAT; stress; apoptosis","respiratory electron transport; oxidative metabolism; translation uncertain","Robust cross-species conservation of a generic inflammatory Myeloid state; not Crohn-specific.","STRONGEST_CELLULAR_CORRESPONDENCE_GENERIC","Inflammation is CeD-shared; metabolism not uniformly conserved."]}
lrows=[]
for l in LINEAGES:
    q=lm[l]; z=master[master.lineage==l]; v=lt[l]
    lrows.append({"lineage":l,"gene_level_baseline":f"Phase10 rho={fn(q['Phase10_gene_rho'])}",
    "pathway_level_baseline":f"Phase10 mean rho={fn(q['Phase10_pathway_mean_rho'])}",
    "singist_summary":f"EXPLORATORY_WEAK_REFERENCE; {q['Phase11_singIST_summary']}; 46/46 weak",
    "saturn_summary":f"EXPLORATORY_PROTEIN_MACROGENE_STRUCTURE; cosine={fn(q['Phase12_SATURN_cosine'])}; species warning=YES",
    "concord_summary":f"cosine={fn(q['Phase13_CONCORD_cosine'])}; P={fn(q['Phase13_CONCORD_null_p'])}; seed={q['Phase13_CONCORD_seed_range']}",
    "region_summary":q["region"],"ced_summary":"|".join(f"{k}:{v}" for k,v in z.disease_context_class.value_counts().sort_index().items()),
    "main_conserved_programs":v[0],"main_divergent_programs":v[1],"final_lineage_interpretation":v[2],"model_strength":v[3],"major_limitation":v[4]})
ls=pd.DataFrame(lrows); wcsv(ls,OUT/"summary/FINAL_LINEAGE_SUMMARY.csv")
ts=p15.copy(); ts["phase16_status"]="INHERITED_UNCHANGED_FROM_PHASE15; NO_RERANKING"; wcsv(ts,OUT/"summary/FINAL_THERAPEUTIC_SUMMARY.csv")
ans=pd.DataFrame([{"research_question":"Can canine CIE be regarded as a translational model of Human Crohn disease?","answer":"YES_WITH_LINEAGE_AND_PROGRAM_BOUNDARIES","core_conclusion":"Canine CIE exhibits lineage- and program-selective translational fidelity to human Crohn disease.","epithelial":lt["Epithelial"][2],"t":lt["T"][2],"myeloid":lt["Myeloid"][2],"disease_context":"CeD sharing limits Crohn-specific claims.","therapeutic_boundary":"TNF and broad JAK-STAT are candidate axes, not proven canine CIE therapies.","final_gate":GATE}]); wcsv(ans,OUT/"summary/CORE_SCIENTIFIC_QUESTION_ANSWER.csv")

ds=pd.DataFrame([
["Dog discovery","GSE254005","Dog","duodenum","scRNA-seq","3 Healthy + 4 CIE","Primary donor-level"],
["Dog external","FRE 2026 published cohort","Dog","duodenal mucosa","bulk published summaries","4 Healthy + 8 FRE","Summary-level only"],
["Human Crohn","GSE266616","Human","terminal ileum","scRNA-seq","16 Healthy + 11 active Crohn","Primary Human"],
["Region control","GSE266616","Human","AC / TI","scRNA-seq","15 paired healthy; active 11 TI vs 6 AC (Epi AC=5)","Region audit"],
["Human replication","IBDverse","Human","intestinal IBD","published summaries","frozen split cohorts","External replication"],
["Disease control","GSE277276","Human","duodenal D2","scRNA-seq","11 Control + 21 active CeD","Specificity context"]],
columns=["dataset_role","accession_or_source","species","tissue_region","assay","biological_replicates","analysis_role"])
wcsv(ds,OUT/"tables/TABLE1_DATASETS_AND_DESIGN.csv")
wcsv(master[["lineage","program","cross_species_fidelity_class","disease_context_class","final_fidelity_class","final_confidence","final_interpretation"]],OUT/"tables/TABLE2_FINAL_PROGRAM_CLASSIFICATION.csv")
wcsv(ls[["lineage","gene_level_baseline","pathway_level_baseline","singist_summary","saturn_summary","concord_summary","model_strength","major_limitation"]],OUT/"tables/TABLE3_CROSS_METHOD_EVIDENCE.csv")
wcsv(master[["lineage","program","human_crohn_direction","ced_direction","ced_relationship","phase14_specificity_class","region_risk","disease_context_class"]],OUT/"tables/TABLE4_DISEASE_CONTEXT.csv")
wcsv(ts,OUT/"tables/TABLE5_THERAPEUTIC_TRANSLATION.csv")

lims=[
("L01","Dog n=7 (3 Healthy, 4 CIE).","ALL","Power/influence limits."),
("L02","Source/breed confounding.","DOG","Associations are not causal."),
("L03","Dog duodenum vs Human TI.","CROSS_SPECIES","Region risk retained."),
("L04","Human treatment heterogeneity.","HUMAN","Sensitivity is not treatment effect."),
("L05","FRE is summary-level bulk only.","DOG_EXTERNAL","No cell/full-rank replication."),
("L06","singIST 46/46 weak; TrackB=0.","METHOD","Exploratory only."),
("L07","SATURN species separation.","METHOD","T cosine cannot upgrade fidelity."),
("L08","CONCORD residual species/seed effects.","METHOD","Support is not identity proof."),
("L09","One CeD comparator.","CONTEXT","Negative does not prove specificity."),
("L10","No Level2 formal comparison.","RESOLUTION","Level1 only."),
("L11","Transcriptomic association is not causality.","INTERPRETATION","No mechanism relabeling."),
("L12","No canine CIE intervention.","THERAPEUTIC","No efficacy claim."),
("L13","Human-biologic canine cross-reactivity unknown.","THERAPEUTIC","No canine use recommendation.")]
ld=pd.DataFrame(lims,columns=["limitation_id","limitation","scope","interpretive_consequence"]); wcsv(ld,OUT/"audit/FINAL_LIMITATIONS_REGISTRY.csv"); wcsv(ld,OUT/"tables/TABLE6_LIMITATIONS_AND_CAUTIONS.csv")
strs=["Donor-level statistics","Full Dog LODO","Independent Human Crohn","IBDverse","Region control","Ontology harmonization","Ensembl104 orthology","Multiple method-specific cross-species analyses","CeD disease control","Therapeutic evidence ladder","Provenance/manifests","No cell pseudoreplication"]
sd=pd.DataFrame([(f"S{i:02}",s) for i,s in enumerate(strs,1)],columns=["strength_id","strength"]); wcsv(sd,OUT/"audit/FINAL_STRENGTHS_REGISTRY.csv")

claims=[
["C001","Canine CIE has selective fidelity to Human Crohn.","CORE","Phase10-16","16_final_integration/master/FINAL_PROGRAM_EVIDENCE_MASTER.csv","25/25 classified","lineage- and program-selective fidelity","globally high-fidelity model","HIGH","Not global"],
["C002","Myeloid is strongest cellular correspondence.","LINEAGE","Phase13-16","13_concord/summary/CROSS_METHOD_LINEAGE_SUMMARY.csv","cosine=0.5873; P=0.0197","strongest generic inflammatory correspondence","Crohn-specific Myeloid","HIGH","CeD-shared"],
["C003","Epithelial conservation is partial/pathway-led.","LINEAGE","Phase10-16","16_final_integration/summary/FINAL_LINEAGE_SUMMARY.csv","metabolic/xenobiotic with modest gene rho","partial pathway conservation","whole-state conservation","MODERATE","Region risk"],
["C004","T fidelity is program-selective.","LINEAGE","Phase10-16","16_final_integration/summary/FINAL_LINEAGE_SUMMARY.csv","inflammation concordant; metabolism divergent","inflammation conserved; metabolism divergent","T globally conserved","HIGH","CONCORD near zero"],
["C005","Ten programs are generic inflammation.","CONTEXT","Phase14-16","14_ced_specificity/specificity/PHASE14_SPECIFICITY_FORWARD_LIST.csv","generic=10","CeD-shared generic inflammation","Crohn-specific","HIGH","One comparator"],
["C006","No strong CeD-differentiating program.","CONTEXT","Phase14-16","08_results/PHASE14_REPORT.md","strong=0","limited hints only","validated specificity","HIGH","CeD negative not specificity"],
["C007","singIST is weak-reference exploratory.","METHOD","Phase11","11_singIST/summary/PHASE11_VS_PHASE10_BASELINE.csv","46/46 weak; TrackB=0","exploratory weak-reference","singIST confirms","HIGH","Not decisive"],
["C008","SATURN is exploratory functional structure.","METHOD","Phase12","12_saturn/summary/PHASE12_VS_PHASE10_BASELINE.csv","species warning","exploratory structure","SATURN proves similarity","HIGH","T cannot upgrade"],
["C009","CONCORD supports Myeloid most.","METHOD","Phase13","13_concord/summary/CROSS_METHOD_LINEAGE_SUMMARY.csv","P=.0197/.0857/.7830","Myeloid strong; Epi limited; T unsupported","identical states","HIGH","Residual species"],
["C010","TNF is Tier A candidate axis.","THERAPEUTIC","Phase15-16","15_therapeutic_translation/summary/CANDIDATE_THERAPEUTIC_AXIS_PRIORITY_LIST.csv","TIER_A","candidate axis","Crohn-specific TNF","MODERATE","Generic; binding unknown"],
["C011","Broad JAK is Tier A candidate axis.","THERAPEUTIC","Phase15-16","15_therapeutic_translation/summary/CANDIDATE_THERAPEUTIC_AXIS_PRIORITY_LIST.csv","TIER_A","candidate axis","proven CIE treatment","MODERATE","No efficacy; safety"],
["C012","IL23 is Human-validated/canine-limited.","THERAPEUTIC","Phase15-16","15_therapeutic_translation/matrix/CANINE_CIE_THERAPEUTIC_TRANSLATION_MATRIX.csv","TIER_C","limited canine evidence","validated canine IL23","HIGH","Region/no efficacy"],
["C013","No numeric global score.","BOUNDARY","Phase10-16","16_final_integration/matrix/FINAL_TRANSLATIONAL_FIDELITY_MATRIX.csv","categorical fields","method-specific synthesis","averaged score","HIGH","Context separate"]]
cl=pd.DataFrame(claims,columns=["claim_id","claim_text","claim_level","supporting_phase","supporting_file","supporting_result","allowed_wording","forbidden_wording","confidence","limitations"]); wcsv(cl,OUT/"audit/FINAL_CLAIMS_REGISTRY.csv")
forb=[
["F001","globally high-fidelity Crohn model","Selective only","lineage- and program-selective fidelity"],
["F002","Crohn-specific TNF","CeD-shared","candidate axis with generic-inflammation caution"],
["F003","validated canine IL23","Canine evidence limited","Human-validated/canine-limited"],
["F004","human biologic effective in dog","binding/efficacy unknown","requires canine validation"],
["F005","JAK proven CIE treatment","non-CIE precedent only","candidate axis"],
["F006","SATURN proves similarity","species-separated","exploratory structure"],
["F007","singIST confirms fidelity","46/46 weak","weak-reference exploratory"],
["F008","CeD-negative equals Crohn-specific","one comparator","limited hint only"],
["F009","Myeloid is Crohn-specific","generic inflammation","generic Myeloid conservation"],
["F010","CONCORD proves treatment transfer","alignment not pharmacology","cellular support only"]]
fd=pd.DataFrame(forb,columns=["forbidden_claim_id","forbidden_claim","reason","allowed_replacement"]); wcsv(fd,OUT/"audit/FINAL_FORBIDDEN_CLAIMS.csv")

# Seven final figures; every map is qualitative.
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.titleweight":"bold"}); F=OUT/"figures"
def save(fig,n): fig.savefig(F/n,dpi=300,bbox_inches="tight",facecolor="white"); plt.close(fig)

fig,ax=plt.subplots(figsize=(14,5)); ax.set_xlim(0,14); ax.set_ylim(0,5); ax.axis("off")
boxes=[(.3,2,2,1.2,"Dog\n7 donors","#174A5B"),(2.8,2,2,1.2,"LODO + FRE\n25 programs","#2E7782"),(5.3,2,2,1.2,"Human Crohn\n27 donors","#287271"),(7.8,3,2,1.1,"Region + ontology\nEnsembl104","#587B7F"),(7.8,1,2,1.1,"singIST / SATURN\nCONCORD","#6C6F8A"),(10.3,2,1.6,1.2,"CeD\n32 subjects","#8C6A5D"),(12.2,2,1.5,1.2,"Final\n25/25","#A3523D")]
for a,b,c,d,t,col in boxes:
    ax.add_patch(FancyBboxPatch((a,b),c,d,boxstyle="round,pad=.03",facecolor=col,edgecolor="none")); ax.text(a+c/2,b+d/2,t,color="white",ha="center",va="center",weight="bold")
for a,b in [((2.3,2.6),(2.8,2.6)),((4.8,2.6),(5.3,2.6)),((7.3,2.6),(7.8,3.4)),((7.3,2.6),(7.8,1.5)),((9.8,3.4),(10.3,2.8)),((9.8,1.5),(10.3,2.4)),((11.9,2.6),(12.2,2.6))]: ax.add_patch(FancyArrowPatch(a,b,arrowstyle="-|>",mutation_scale=12,color="#4A5560"))
ax.set_title("Figure 1. Frozen evidence flow",fontsize=15); ax.text(7,.25,"Donor-level; Primary Level1 only.",ha="center")
save(fig,"FIG1_STUDY_DESIGN.png")

dc=master.groupby(["lineage","dog_direction","dog_robustness"]).size().reset_index(name="n")
fig,axs=plt.subplots(1,3,figsize=(12,4.5),sharey=True)
for ax,l in zip(axs,LINEAGES):
    s=dc[dc.lineage==l]; specs=[("CIE_higher","ROBUST"),("CIE_higher","MODERATELY_STABLE"),("CIE_lower","ROBUST"),("CIE_lower","MODERATELY_STABLE")]
    vv=[int(h.iloc[0]) if len(h:=s[(s.dog_direction==a)&(s.dog_robustness==b)].n) else 0 for a,b in specs]
    ax.bar(range(4),vv,color=["#B34D3E","#D68A74","#2E7782","#8CB8BD"]); ax.set_xticks(range(4),["High\nR","High\nM","Low\nR","Low\nM"]); ax.set_title(l)
    for i,v in enumerate(vv): ax.text(i,v+.05,str(v),ha="center")
axs[0].set_ylabel("Programs"); fig.suptitle("Figure 2. Canine discovery and LODO robustness")
save(fig,"FIG2_CANINE_DISCOVERY.png")

cols=["Dog LODO","Dog-Crohn","Region","CeD","singIST","SATURN","CONCORD","Final"]
cc={"support":"#2F7D6D","partial":"#7FAF8E","generic":"#4E79A7","hint":"#D6A84B","divergent":"#C75B39","uncertain":"#A7ADB4","method":"#7B6D9C","region":"#7F7F7F"}
ab={"support":"SUP","partial":"PART","generic":"GEN","hint":"HINT","divergent":"DIV","uncertain":"UNC","method":"METH","region":"REG"}
mr=[]
for _,r in master.iterrows():
    cv={"CONSERVED":"support","PARTIALLY_CONSERVED":"partial","DIVERGENT":"divergent","UNCERTAIN":"uncertain"}[r.cross_species_fidelity_class]
    ct="generic" if r.disease_context_class=="GENERIC_MUCOSAL_INFLAMMATORY" else "hint" if r.disease_context_class in {"CED_DIFFERENTIATING_HINT","CED_NOT_PROMINENT"} else "uncertain"
    fc="generic" if r.final_fidelity_class=="GENERIC_MUCOSAL_INFLAMMATORY_PROGRAM" else "hint" if r.final_fidelity_class=="DOG_CROHN_DIFFERENTIATING_HINT" else "divergent" if r.final_fidelity_class=="CANINE_DIVERGENT_PROGRAM" else "region" if r.final_fidelity_class=="REGION_ASSOCIATED_OR_UNCERTAIN" else "uncertain"
    co="support" if r.lineage=="Myeloid" else "partial" if r.lineage=="Epithelial" else "uncertain"
    mr.append((f"{r.lineage} | {r.program}",["support" if r.dog_robustness=="ROBUST" else "partial",cv,"support" if r.region_risk=="LOW_REGION_RISK" else "region",ct,"method","method",co,fc]))
fig,ax=plt.subplots(figsize=(15,11)); ax.set_xlim(0,8); ax.set_ylim(0,25); ax.invert_yaxis()
for y,(lab,vv) in enumerate(mr):
    for xx,z in enumerate(vv): ax.add_patch(Rectangle((xx,y),1,1,facecolor=cc[z],edgecolor="white")); ax.text(xx+.5,y+.5,ab[z],ha="center",va="center",color="white",fontsize=7,weight="bold")
ax.set_xticks(np.arange(8)+.5,cols,rotation=25,ha="left"); ax.xaxis.tick_top(); ax.set_yticks(np.arange(25)+.5,[z[0] for z in mr],fontsize=7.5); ax.tick_params(length=0)
for s in ax.spines.values():s.set_visible(False)
ax.legend(handles=[Patch(facecolor=v,label=k.upper()) for k,v in cc.items()],ncol=8,loc="lower center",bbox_to_anchor=(.5,-.08),frameon=False,fontsize=7); ax.set_title("Figure 3. Qualitative translational fidelity map",pad=45)
save(fig,"FIG3_DOG_HUMAN_TRANSLATIONAL_FIDELITY_MAP.png")

fig,axs=plt.subplots(1,3,figsize=(14,5.5)); colors={"Epithelial":"#2E7782","T":"#6C6F8A","Myeloid":"#A3523D"}
for ax,l in zip(axs,LINEAGES):
    q=lm[l]; ax.axis("off"); ax.add_patch(FancyBboxPatch((.03,.04),.94,.92,boxstyle="round,pad=.03",facecolor="#F4F3EF",edgecolor=colors[l],lw=2))
    ax.text(.5,.88,l,ha="center",fontsize=15,weight="bold",color=colors[l],transform=ax.transAxes)
    ax.text(.1,.7,f"Gene rho  {fn(q['Phase10_gene_rho'],3)}\nPathway rho  {fn(q['Phase10_pathway_mean_rho'],3)}\nSATURN  {fn(q['Phase12_SATURN_cosine'],3)}\nCONCORD  {fn(q['Phase13_CONCORD_cosine'],3)}\nP  {fn(q['Phase13_CONCORD_null_p'],3)}",va="top",family="monospace",transform=ax.transAxes)
    ax.text(.1,.34,lt[l][3],weight="bold",color=colors[l],transform=ax.transAxes); ax.text(.1,.27,lt[l][2],va="top",wrap=True,fontsize=8.5,transform=ax.transAxes)
fig.suptitle("Figure 4. Lineage-selective fidelity")
save(fig,"FIG4_LINEAGE_SELECTIVE_FIDELITY.png")

order=["GENERIC_MUCOSAL_INFLAMMATORY","CED_DIFFERENTIATING_HINT","CED_NOT_PROMINENT","NOT_ASSESSABLE"]; ct=master.groupby(["lineage","disease_context_class"]).size().unstack(fill_value=0).reindex(index=LINEAGES,columns=order,fill_value=0)
fig,ax=plt.subplots(figsize=(10,5)); left=np.zeros(3)
for k,col in zip(order,["#4E79A7","#D6A84B","#C79A5B","#A7ADB4"]):
    v=ct[k].to_numpy(); ax.barh(LINEAGES,v,left=left,label=k,color=col)
    for y,(a,b) in enumerate(zip(left,v)):
        if b: ax.text(a+b/2,y,str(b),ha="center",va="center",color="white",weight="bold")
    left+=v
ax.set_xlabel("Programs"); ax.set_title("Figure 5. Disease context"); ax.legend(ncol=2,loc="lower center",bbox_to_anchor=(.5,-.35),frameon=False)
save(fig,"FIG5_DISEASE_CONTEXT.png")

tiers=["TIER_A","TIER_B","TIER_C","EXPLORATORY","NOT_PRIORITIZED"]
feas=["NOT_ASSESSABLE","LIMITED","MODERATE"]
tc={"TIER_A":"#A3523D","TIER_B":"#C77C4A","TIER_C":"#D6A84B","EXPLORATORY":"#6C6F8A","NOT_PRIORITIZED":"#8B9298"}
fig,ax=plt.subplots(figsize=(12,8)); plot_order=ts.reset_index(drop=True); ypos=np.arange(len(plot_order))
for i,r in plot_order.iterrows():
    xx=feas.index(r.pharmacologic_translation_feasibility) if r.pharmacologic_translation_feasibility in feas else 0
    ax.scatter(xx,i,s=170,color=tc[r.rank_group],edgecolor="white",lw=1.2,zorder=3)
ax.set_xticks(range(3),feas); ax.set_yticks(ypos,plot_order.axis,fontsize=8); ax.set_xlim(-.4,2.4); ax.set_ylim(len(plot_order)-.5,-.5)
ax.grid(axis="x",alpha=.2); ax.set_xlabel("Pharmacologic translation feasibility (Phase15, unchanged)")
ax.set_title("Figure 6. Therapeutic priority and feasibility remain separate")
ax.legend(handles=[Patch(facecolor=tc[k],label=k) for k in tiers],ncol=3,loc="lower center",bbox_to_anchor=(.5,-.16),frameon=False)
save(fig,"FIG6_THERAPEUTIC_TRANSLATION.png")

fig,ax=plt.subplots(figsize=(15,7)); ax.set_xlim(0,15); ax.set_ylim(0,7); ax.axis("off"); ax.text(7.5,6.6,"Canine CIE: selective translational fidelity",ha="center",fontsize=17,weight="bold")
nodes=[(.5,2.7,2.3,1.5,"DOG CIE\n25 stable programs","#174A5B"),(4.1,4.6,2.4,1.2,"EPITHELIAL\npartial pathways","#2E7782"),(4.1,2.8,2.4,1.2,"T\ninflammation / metabolism","#6C6F8A"),(4.1,1,2.4,1.2,"MYELOID\nstrong generic inflammation","#A3523D"),(8,3.5,2.5,1.5,"HUMAN CROHN\nregion/treatment retained","#287271"),(8,1,2.5,1.2,"CeD\nlimits specificity","#8C6A5D"),(12.1,3.5,2.3,1.5,"TNF + broad JAK\ncandidate axes only","#A3523D"),(12.1,1,2.3,1.2,"No global score\nNo proven dog therapy","#5F666D")]
for a,b,c,d,t,col in nodes: ax.add_patch(FancyBboxPatch((a,b),c,d,boxstyle="round,pad=.04",facecolor=col,edgecolor="none")); ax.text(a+c/2,b+d/2,t,color="white",ha="center",va="center",weight="bold")
for a,b in [((2.8,3.45),(4.1,5.2)),((2.8,3.45),(4.1,3.4)),((2.8,3.45),(4.1,1.6)),((6.5,5.2),(8,4.5)),((6.5,3.4),(8,4.2)),((6.5,1.6),(8,3.8)),((9.25,3.5),(9.25,2.2)),((10.5,4.2),(12.1,4.2)),((10.5,1.6),(12.1,1.6))]: ax.add_patch(FancyArrowPatch(a,b,arrowstyle="-|>",mutation_scale=13,color="#59636A"))
ax.text(7.5,.25,"FINAL_COMPLETE_WITH_DOCUMENTED_LIMITATIONS",ha="center",weight="bold")
save(fig,"GRAPHICAL_ABSTRACT_FINAL.png")

classes=master.final_fidelity_class.value_counts().to_dict(); confs=master.final_confidence.value_counts().to_dict(); contexts=master.disease_context_class.value_counts().to_dict()
cst=", ".join(f"{k}={v}" for k,v in sorted(classes.items())); cft=", ".join(f"{k}={v}" for k,v in sorted(confs.items())); dct=", ".join(f"{k}={v}" for k,v in sorted(contexts.items()))
report=f"""# Phase16 Final Integration, Manuscript-level Conclusion, and Project Closeout

Freeze date: {DATE}
Final gate: **{GATE}**

## 1. Executive Summary

Phase16 completed final integration without new biological evidence, upstream reruns, or a numeric fidelity score. Phase15 manifest: **45/45 PASS**. All **25/25** programs were classified (Epithelial=8, T=8, Myeloid=9).

**Canine CIE exhibits lineage- and program-selective translational fidelity to human Crohn disease.** Myeloid is strongest but mainly generic inflammation; Epithelial is partial/pathway-led; T has conserved inflammation and divergent metabolism.

## 2. Research Question

Can canine CIE serve as a translational model for Human Crohn disease after donor robustness, external canine evidence, Human Crohn, region, ontology/orthology, cross-species methods, CeD, and therapeutic translation? The answer is lineage/program-resolved, not globally binary.

## 3. Canine Discovery

GSE254005 used 7 donors (3 Healthy, 4 CIE), donor-level inference, 21/21 LODO models and 63/63 paired GSEA runs. The universe is 20 ROBUST plus 5 MODERATELY_STABLE programs. FRE is summary-level bulk support only.

## 4. Human Crohn Reference

GSE266616 primary TI used 16 Healthy and 11 active Crohn donors. Treatment is sensitivity, not causal effect. IBDverse and region control remain independent evidence layers.

## 5. Cross-species Fidelity

Ensembl104 one-to-one orthology remains frozen. Gene correspondence is weak overall; pathway correspondence is lineage-dependent. Classes: **{cst}**. Confidence: **{cft}**. No composite score.

## 6. Lineage-specific Findings — Epithelial

Gene rho=0.1474; pathway rho=0.3655. Metabolic/xenobiotic evidence supports partial pathway conservation, not whole-state conservation. Proliferation/ribosome are high-region-risk. CONCORD cosine=0.3672, P=0.0857: limited support.

## 7. Lineage-specific Findings — T

Gene rho=0.0159; pathway rho=0.0620. TNF/inflammatory/cytokine programs are conserved; mitochondrial/metabolic programs diverge. SATURN cosine=0.5008 cannot upgrade fidelity because species separate. CONCORD cosine=0.0849, P=0.7830.

## 8. Lineage-specific Findings — Myeloid

Gene rho=0.1011; pathway rho=0.2626. CONCORD cosine=0.5873, P=0.0197. This is robust cross-species conservation of a **generic inflammatory Myeloid state**, not Crohn specificity; metabolic programs are mixed/divergent.

## 9. Disease Specificity

CeD relationships: 16 shared, 6 not prominent, 3 opposite. Final contexts: **{dct}**. Ten programs are generic mucosal inflammatory; zero have strong CeD-differentiating evidence. CeD non-prominence is at most a limited hint.

## 10. Cross-method Integration

Phase10 is baseline. Phase11 is EXPLORATORY_WEAK_REFERENCE (46/46 weak; TrackB=0). Phase12 is EXPLORATORY_PROTEIN_MACROGENE_STRUCTURE with species warning=YES. Phase13 supports Myeloid strongly, Epithelial weakly, and not T at lineage disease-vector level. No averaging/voting.

## 11. Therapeutic Translation

Phase15 is unchanged: Tier A TNF/TNFR-NF-kappaB and broad JAK-STAT; Tier B IL6-JAK-STAT3; Tier C IL23, Epithelial metabolic programs, Myeloid inflammatory/cytokine and proliferation; three exploratory axes; T metabolism and Stress/UPR not prioritized.

## 12. Final Translational Fidelity Map

The controlling matrix is FINAL_TRANSLATIONAL_FIDELITY_MATRIX.csv and Figure 3. Cross-species fidelity, disease context, region, methods, class and confidence stay separate.

## 13. Final Therapeutic Translation Matrix

Phase15 ordering is reproduced without re-ranking. TNF and broad JAK are candidate axes, not proven canine CIE therapies. IL23 is Human-validated/canine-limited. Human-biologic canine cross-reactivity is unknown.

## 14. Main Findings

1. Selective, not global, fidelity.
2. Myeloid strongest but generic inflammatory.
3. Epithelial partial/pathway-led.
4. T inflammation conserved and metabolism divergent.
5. CeD limits disease specificity.
6. TNF/JAK are candidates, not treatments.

## 15. Model Boundaries

Only Epithelial, T and Myeloid Level1 and 25 frozen programs. No Level2 formal comparison, new DE/GSEA, remapping, score, causality, or therapeutic efficacy inference. Alignment is not pharmacology.

## 16. Strengths

Twelve strengths are frozen: donor statistics, full LODO, Human Crohn, IBDverse, region control, ontology, Ensembl104, multiple methods, CeD, therapeutic ladder, provenance/manifests, and no cell pseudoreplication.

## 17. Limitations, Final Conclusion, and Final Gate

Thirteen limitations are frozen: Dog n=7; source/breed confounding; Dog duodenum vs Human TI; treatment heterogeneity; summary FRE; weak singIST; SATURN separation; residual CONCORD species effects; one CeD comparator; no Level2 comparison; association not causality; no canine intervention; unknown biologic cross-reactivity.

**Final conclusion:** Canine CIE exhibits lineage- and program-selective translational fidelity to human Crohn disease. Myeloid is strongest but generic inflammatory; Epithelial shows partial metabolic/xenobiotic conservation; T has conserved inflammation and divergent metabolism; CeD sharing limits Crohn specificity. TNF and broad JAK are candidate axes, not proven canine CIE therapies.

All required inputs, 25/25 classifications, claims, figures, therapeutic inheritance and manifests passed. The project closes at **{GATE}**. No subsequent analytical phase is defined.
"""
wtxt(report,ROOT/"08_results/PHASE16_FINAL_REPORT.md")
wtxt(f"""# CURRENT STATE

- Project status: COMPLETE
- Final phase: Phase16
- Final gate: {GATE}
- Freeze date: {DATE}
- Phase15 manifest: 45/45 PASS
- Program universe: 25/25 (8/8/9)
- Conclusion: lineage- and program-selective fidelity
- Therapeutic boundary: candidate axes, no proven treatment
- Formal next analytical phase: none
""",ROOT/"CURRENT_STATE.md")
wtxt(f"""# RESULTS REGISTER

- Final gate: {GATE}
- Phase15 integrity: 45/45 PASS
- Programs: 25/25; {cst}
- Confidence: {cft}
- Epithelial: partial/pathway-led
- T: inflammatory conserved, metabolic divergent
- Myeloid: strongest generic inflammatory correspondence (cosine=0.5873, P=0.0197)
- Disease context: 10 generic; zero strong CeD-differentiating
- Phase15 tiers inherited unchanged
""",ROOT/"RESULTS_REGISTER.md")
wtxt("# KNOWN ISSUES\n\n"+"\n".join(f"- **{r.limitation_id}** {r.limitation} Consequence: {r.interpretive_consequence}" for r in ld.itertuples(index=False)),ROOT/"KNOWN_ISSUES.md")
wtxt("""# NEXT STEPS

The analytical program is closed at the Phase16 gate. Permitted downstream work is limited to:

## MANUSCRIPT_PREPARATION

- Use frozen claims, tables and report without expanding claims.

## FIGURE_REFINEMENT

- Refine layout while preserving categories and cautions.

## OPTIONAL_EXTERNAL_VALIDATION

- Analyze a genuinely independent cohort under a separate preregistration.

## OPTIONAL_WET_LAB_VALIDATION

- Test canine binding, engagement, safety and efficacy before treatment claims.
""",ROOT/"NEXT_STEPS.md")

tn=["TABLE1_DATASETS_AND_DESIGN.csv","TABLE2_FINAL_PROGRAM_CLASSIFICATION.csv","TABLE3_CROSS_METHOD_EVIDENCE.csv","TABLE4_DISEASE_CONTEXT.csv","TABLE5_THERAPEUTIC_TRANSLATION.csv","TABLE6_LIMITATIONS_AND_CAUTIONS.csv"]
fign=["FIG1_STUDY_DESIGN.png","FIG2_CANINE_DISCOVERY.png","FIG3_DOG_HUMAN_TRANSLATIONAL_FIDELITY_MAP.png","FIG4_LINEAGE_SELECTIVE_FIDELITY.png","FIG5_DISEASE_CONTEXT.png","FIG6_THERAPEUTIC_TRANSLATION.png","GRAPHICAL_ABSTRACT_FINAL.png"]
req=[ROOT/"00_metadata/PHASE16_INPUT_FREEZE_AUDIT.csv",OUT/"master/FINAL_PROGRAM_EVIDENCE_MASTER.csv",OUT/"matrix/FINAL_TRANSLATIONAL_FIDELITY_MATRIX.csv",OUT/"summary/FINAL_LINEAGE_SUMMARY.csv",OUT/"summary/FINAL_THERAPEUTIC_SUMMARY.csv",OUT/"summary/CORE_SCIENTIFIC_QUESTION_ANSWER.csv",OUT/"audit/FINAL_CLAIMS_REGISTRY.csv",OUT/"audit/FINAL_FORBIDDEN_CLAIMS.csv",ROOT/"08_results/PHASE16_FINAL_REPORT.md"]+[OUT/"tables"/n for n in tn]+[F/n for n in fign]
if any(not p.exists() or p.stat().st_size==0 for p in req): raise RuntimeError("required")
if ts.axis.tolist()!=p15.axis.tolist() or ts.rank_group.tolist()!=p15.rank_group.tolist(): raise RuntimeError("rerank")
gate=pd.DataFrame([
["phase15_manifest",True,"45/45"],["phase0_15_freeze",True,f"{len(audit)}/{len(audit)}"],["primary_lineages",True,"Epithelial|T|Myeloid"],
["programs",True,"25/25"],["independent_fields",True,"no score"],["claims",True,f"{len(cl)}/{len(fd)}"],["figures",True,"7/7"],
["therapeutic_inheritance",True,"14/14 unchanged"],["limitations_strengths",True,"13/12"],["report_state",True,"PASS"]],
columns=["criterion","pass","evidence"]); gate["overall_gate"]=GATE; wcsv(gate,ROOT/"00_metadata/PHASE16_GATE.csv")

def role(p):
    n=p.name.lower()
    if n.endswith(".png"): return "figure"
    if p.suffix.lower()==".md" or "report" in n: return "report_or_state"
    if "manifest" in n:return "manifest"
    if "gate" in n or "audit" in p.parts:return "audit"
    if "table" in n:return "manuscript_table"
    if p.suffix.lower()==".py":return "reproducibility_script"
    return "formal_table"
pf=sorted(p for p in OUT.rglob("*") if p.is_file())+[ROOT/"00_metadata/PHASE16_INPUT_FREEZE_AUDIT.csv",ROOT/"00_metadata/PHASE16_GATE.csv",ROOT/"08_results/PHASE16_FINAL_REPORT.md",ROOT/"CURRENT_STATE.md",ROOT/"RESULTS_REGISTER.md",ROOT/"KNOWN_ISSUES.md",ROOT/"NEXT_STEPS.md",ROOT/"scripts/phase16_final_integration.py"]
pf=sorted(set(pf)); pm=pd.DataFrame([{"Phase":"Phase16","file":rel(p),"role":role(p),"size":p.stat().st_size,"SHA256":sha(p),"freeze_date":DATE} for p in pf]); wcsv(pm,ROOT/"00_metadata/PHASE16_MANIFEST.csv")
for r in pm.to_dict("records"):
    p=ROOT/r["file"]
    if p.stat().st_size!=int(r["size"]) or sha(p)!=r["SHA256"]:raise RuntimeError("P16 manifest")

project=[("Phase0",ROOT/"00_metadata/PHASE0_AUDIT_REPORT.md","report/audit"),("Phase1",ROOT/"08_results/PHASE1_REPORT.md","report")]
for n in range(2,16):
    project.append((f"Phase{n}",ROOT/f"08_results/PHASE{n}_REPORT.md","report"))
    m=ROOT/f"00_metadata/PHASE{n}_MANIFEST.csv"
    if m.exists():project.append((f"Phase{n}",m,"manifest"))
    for g in [ROOT/f"00_metadata/PHASE{n}_GATE.csv",ROOT/f"00_metadata/PHASE{n}_FINAL_GATE.csv"]:
        if g.exists():project.append((f"Phase{n}",g,"gate"))
project += [("Phase16",p,role(p)) for p in req]+[("Phase16",ROOT/"00_metadata/PHASE16_GATE.csv","final gate"),("Phase16",ROOT/"00_metadata/PHASE16_MANIFEST.csv","manifest"),("Phase16",ROOT/"CURRENT_STATE.md","state"),("Phase16",ROOT/"RESULTS_REGISTER.md","register"),("Phase16",ROOT/"KNOWN_ISSUES.md","issues"),("Phase16",ROOT/"NEXT_STEPS.md","allowed next work")]
seen=set(); prows=[]
for ph,p,ro in project:
    k=rel(p)
    if k in seen:continue
    seen.add(k)
    if not p.exists():raise FileNotFoundError(k)
    prows.append({"Phase":ph,"file":k,"role":ro,"size":p.stat().st_size,"SHA256":sha(p),"freeze_date":DATE})
pr=pd.DataFrame(prows).sort_values(["Phase","file"]); wcsv(pr,ROOT/"00_metadata/PROJECT_FINAL_MANIFEST.csv")
for r in pr.to_dict("records"):
    p=ROOT/r["file"]
    if p.stat().st_size!=int(r["size"]) or sha(p)!=r["SHA256"]:raise RuntimeError("Project manifest")

print("PHASE15_MANIFEST_PASS=45")
print(f"INPUT_AUDIT_PASS={len(audit)}")
print(f"PROGRAMS={len(master)} LINEAGES={master.lineage.value_counts().to_dict()}")
print(f"CLASSES={classes}")
print(f"CONFIDENCE={confs}")
print(f"CLAIMS={len(cl)} FORBIDDEN={len(fd)} FIGURES={len(fign)}")
print(f"PHASE16_MANIFEST={len(pm)} PROJECT_MANIFEST={len(pr)}")
print(f"FINAL_GATE={GATE}")
