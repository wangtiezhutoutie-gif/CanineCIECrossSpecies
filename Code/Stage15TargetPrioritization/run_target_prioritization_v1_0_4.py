from __future__ import annotations
import csv, hashlib, json
from collections import Counter
from datetime import datetime
from pathlib import Path
import pandas as pd

ROOT=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
BASE=ROOT/"TARGET_PRIORITIZATION_V1_0"
RUN_ID="TP_V1_0_20260916_005"
DESIGN_VERSION="Target_Prioritization_v1.0"
WORKFLOW_VERSION="Target_Prioritization_Workflow_v1.0.4"
CUTOFF="2026-09-16"
LINEAGES=["Epithelial","T","Myeloid"]
EXPECTED={
 "TRACK_A_FROZEN_CANDIDATES.tsv":"5f6d0979642824bd16bc665e38ef18e977b266bf7b34f2e33c584a8f4eb579ee",
 "TARGET_PROGRAM_LINEAGE_ANCHOR_SNAPSHOT.tsv":"29332ed03b2f8027a3149d101f450abeb66a2b187b2a4d60905eab683419db4a",
 "TRACK_B_FROZEN_CANDIDATES.tsv":"17336de3a0ab3aa95da13c2b8f5240af5db8bc6f4ceb6fab7686ddb8f193bab2"}

def sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for b in iter(lambda:f.read(1048576),b""): h.update(b)
 return h.hexdigest()
def rd(p,sep=None):
 return pd.read_csv(p,sep=sep or ("\t" if str(p).lower().endswith(".tsv") else ","),dtype=str,keep_default_na=False)
def wt(p,df):
 Path(p).parent.mkdir(parents=True,exist_ok=True); df.to_csv(p,sep="\t",index=False,lineterminator="\n")
def split(v): return [x.strip() for x in str(v).split(";") if x.strip()]
def fnum(v):
 try:return f"{float(v):.4g}"
 except:return "NA" if not v else str(v)
def direction(v,up,down):
 try:
  x=float(v); return up if x>0 else down if x<0 else "NO_CLEAR_DIRECTION"
 except:return "NOT_ASSESSABLE"
def passfail(v):
 return {"SUPPORTED":"PASS","PARTIALLY_SUPPORTED":"PARTIAL","NO_POSITIVE_SUPPORT":"FAIL_NO_POSITIVE_SUPPORT","DISCORDANT":"FAIL_DISCORDANT","CONTRADICTED":"FAIL_DISCORDANT","NOT_ASSESSABLE":"NOT_ASSESSABLE"}.get(v,v)
def trans(v):
 return {"SUPPORTED":"YES","PARTIALLY_SUPPORTED":"PARTIAL","NO_POSITIVE_SUPPORT":"NO","CONTRADICTED":"NO","NOT_ASSESSABLE":"NOT_ASSESSABLE"}.get(v,"NOT_ASSESSABLE")

for d in ["00_DESIGN","01_INPUT_MANIFEST","02_FROZEN_CANDIDATES","03_EVIDENCE","04_DECISION_AUDIT","05_FINAL_RESULTS","06_FIGURE_TABLES","07_LOGS"]:(BASE/d).mkdir(parents=True,exist_ok=True)
snap=BASE/"02_FROZEN_CANDIDATES"
for n,h in EXPECTED.items():
 if sha(snap/n)!=h: raise RuntimeError("FROZEN_SNAPSHOT_HASH_MISMATCH:"+n)
A=rd(snap/"TRACK_A_FROZEN_CANDIDATES.tsv")
ANCH=rd(snap/"TARGET_PROGRAM_LINEAGE_ANCHOR_SNAPSHOT.tsv").set_index("TARGET_ENTITY_ID")
B=rd(snap/"TRACK_B_FROZEN_CANDIDATES.tsv")
CFG=rd(BASE/"00_DESIGN"/"DECISION_CLASSIFICATION_V1_0_FROZEN.tsv")
if (len(A),len(B),len(CFG))!=(9,24,33):raise RuntimeError("FROZEN_POOL_ROW_COUNT_MISMATCH")
keys=set(zip(CFG.ENTRY_TRACK,CFG.TARGET_ENTITY_ID))
expected={("TRACK_A",x) for x in A.TARGET_ENTITY_ID}|{("TRACK_B",x) for x in B.TARGET_ENTITY_ID}
if keys!=expected:raise RuntimeError("DECISION_CONFIG_KEY_MISMATCH")
CFG=CFG.set_index(["ENTRY_TRACK","TARGET_ENTITY_ID"])

p_program=ROOT/"08_results"/"supplementary_robustness"/"human_duodenum_ti_region_sensitivity"/"04_program_integration"/"DOG_HUMAN_REGION_25_PROGRAM_MATRIX.csv"
P=rd(p_program)
if len(P)!=25 or P.frozen_program_id.nunique()!=25:raise RuntimeError("FROZEN_PROGRAM_ID_NOT_UNIQUE")
P=P.set_index("frozen_program_id")
components=sorted(set(sum([split(x) for x in pd.concat([A.TARGET_ENTITY_COMPONENTS,B.TARGET_ENTITY_COMPONENTS])],[])))
programs=set()
for _,r in A.iterrows():programs.update(split(ANCH.loc[r.TARGET_ENTITY_ID,"ALLOWED_FROZEN_PROGRAMS"]))
for _,r in B.iterrows():programs.update(split(r.SOURCE_PROGRAMS))
if programs-set(P.index):raise RuntimeError("FROZEN_PROGRAM_NOT_UNIQUELY_IDENTIFIED:"+";".join(sorted(programs-set(P.index))))

# Recompute all candidate gene effects from six frozen DE tables.
ge=[]; gem={}; de_paths=[]
for lin in LINEAGES:
 dp=ROOT/"06_pseudobulk"/"de_results"/f"{lin}_CIE_vs_Healthy_edger.csv"
 hp=ROOT/"07_human_crohn"/"GSE266616"/"pseudobulk"/"de_results"/f"Human_{lin}_ActiveCrohn_vs_Healthy.csv"
 de_paths += [dp,hp]
 d=rd(dp); h=rd(hp)
 d=d[d.gene_symbol.isin(components)].drop_duplicates("gene_symbol").set_index("gene_symbol")
 h=h[h.gene_symbol.isin(components)].drop_duplicates("gene_symbol").set_index("gene_symbol")
 for g in components:
  dv=d.loc[g] if g in d.index else None; hv=h.loc[g] if g in h.index else None
  dd=direction(dv.logFC if dv is not None else "","CIE_higher","CIE_lower")
  hd=direction(hv.logFC if hv is not None else "","ActiveCrohn_higher","ActiveCrohn_lower")
  if dv is None or hv is None:pair="NOT_ASSESSABLE"
  elif "NO_CLEAR_DIRECTION" in (dd,hd):pair="NO_CLEAR_DIRECTION"
  else:pair="SAME_DIRECTION" if (float(dv.logFC)>0)==(float(hv.logFC)>0) else "OPPOSITE_DIRECTION"
  q={"gene":g,"lineage":lin,"dog_tested":str(dv is not None).upper(),"dog_logFC":"" if dv is None else dv.logFC,"dog_FDR":"" if dv is None else dv.FDR,"dog_direction":dd,"human_tested":str(hv is not None).upper(),"human_logFC":"" if hv is None else hv.logFC,"human_FDR":"" if hv is None else hv.FDR,"human_direction":hd,"gene_effect_direction_class":pair}
  ge.append(q);gem[(g,lin)]=q
GE=pd.DataFrame(ge);wt(BASE/"03_EVIDENCE"/"CANDIDATE_GENE_EFFECTS_RUN004.tsv",GE)

# Recompute leading-edge membership from all 18 final exports.
term={(P.loc[x,"lineage"],P.loc[x,"term_id"]):x for x in programs}
leh=[];le_paths=[]
for sp in ["Dog","Human"]:
 for lin in LINEAGES:
  for filecol in ["Hallmark","Reactome","GOBP"]:
   if sp=="Dog":p=ROOT/"06_pseudobulk"/"gsea"/"leading_edges"/f"{lin}_{filecol}_leading_edges.csv"
   else:p=ROOT/"07_human_crohn"/"GSE266616"/"gsea"/"leading_edges"/f"Human_{lin}_{'GO_BP' if filecol=='GOBP' else filecol}_leading_edges.csv"
   le_paths.append(p);x=rd(p);x=x[x.gene_symbol.isin(components)]
   for _,z in x.iterrows():
    tid=z.pathway_id.replace("R-CFA-","R-HSA-") if sp=="Dog" else z.pathway_id
    pid=term.get((lin,tid))
    if pid:leh.append({"source_species":sp,"lineage":lin,"frozen_program_id":pid,"term_id":tid,"gene":z.gene_symbol,"pathway_direction":z.pathway_direction,"pathway_NES":z.pathway_NES,"pathway_FDR":z.pathway_padj,"gene_logFC":z.logFC,"gene_FDR":z.FDR})
LE=pd.DataFrame(leh).drop_duplicates().sort_values(["source_species","lineage","frozen_program_id","gene"])
wt(BASE/"03_EVIDENCE"/"CANDIDATE_LEADING_EDGE_RUN004.tsv",LE)

# Resolve components against the existing frozen ortholog table.
p_ortho=ROOT/"10_cross_species_baseline"/"mapping"/"DOG_HUMAN_ORTHOLOG_MASTER.csv"
O=rd(p_ortho);OM={}
for g in components:
 x=O[(O.dog_original_symbol==g)&(O.human_ortholog_symbol==g)&(O.primary_one2one.str.lower()=="true")]
 OM[g]="HIGH_CONFIDENCE_ONE2ONE" if len(x)>=1 else "NOT_ASSESSABLE"
prov_path=BASE/"03_EVIDENCE"/"evidence_provenance.tsv";PROV=rd(prov_path)
req={"target_entity_id","drug_if_applicable","direct_target_entity","mechanism_of_action","disease","evidence_level","study_type","publication_year","PMID","DOI","source_title","evidence_cutoff_date","short_evidence_statement","supports_or_limits","source_species"}
if not req.issubset(PROV.columns):raise RuntimeError("EVIDENCE_PROVENANCE_FIELDS_MISSING")
def gene_text(comps,sp):
 out=[]
 for g in comps:
  for lin in LINEAGES:
   r=gem[(g,lin)];dire=r["dog_direction" if sp=="Dog" else "human_direction"];fdr=r["dog_FDR" if sp=="Dog" else "human_FDR"]
   out.append(f"{g}/{lin}:{dire};FDR={fnum(fdr)};pair={r['gene_effect_direction_class']}")
 return " | ".join(out)
def le_state(comps,ps,sp):
 x=LE[(LE.source_species==sp)&LE.gene.isin(comps)&LE.frozen_program_id.isin(ps)]
 if not ps:return "NOT_ASSESSABLE","No frozen source program"
 if x.empty:return "NO","No candidate component in frozen leading edges"
 detail=";".join(f"{r.gene}@{r.frozen_program_id}({r.pathway_direction})" for _,r in x.iterrows())
 return "YES",f"HITS={len(x)}:{detail}"

def build(track,base):
 cid=base.TARGET_ENTITY_ID;cfg=CFG.loc[(track,cid)]
 if track=="TRACK_A":
  ps=split(ANCH.loc[cid,"ALLOWED_FROZEN_PROGRAMS"]);lins=ANCH.loc[cid,"ALLOWED_LINEAGES"]
  reason="Human-first direct therapeutic target within predefined discovery universe";disc=base.CANDIDATE_DISCOVERY_SOURCE;human=base.HUMAN_THERAPEUTIC_EVIDENCE
  direct="SUPPORTED" if base.DIRECT_TARGET_STATUS=="DIRECT" else "NOT_ASSESSABLE";hgate="SUPPORTED"
 else:
  ps=split(base.SOURCE_PROGRAMS);lins=base.SOURCE_LINEAGES
  reason="Canonical regulatory node from frozen 25-program universe";disc="FROZEN_25_PROGRAMS";human="NONE";direct="NOT_ASSESSABLE";hgate="NOT_ASSESSABLE"
 comps=split(base.TARGET_ENTITY_COMPONENTS)
 pe=[];re=[];mt=[];ce=[];lc=[]
 for pid in ps:
  p=P.loc[pid]
  pe.append(f"{pid}[Dog={p.dog_disease_direction};Human={p.phase16_human_crohn_direction};Phase10={p.phase16_human_crohn_evidence};robust={p.phase16_dog_robustness};region={p.phase16_region_risk};CeD={p.phase16_ced_relationship}]")
  re.append(f"{pid}={p.region_interpretation};phase16={p.phase16_region_risk}")
  mt.append(f"{pid}={p.dog_effect_qc_robustness}")
  cm={"CED_SHARED_DIRECTION":"CROHN_CED_SHARED","CED_OPPOSITE_DIRECTION":"CROHN_CED_DISCORDANT","CED_NOT_PROMINENT":"CED_NOT_PROMINENT"}
  ce.append(f"{pid}={cm.get(p.phase16_ced_relationship,'CED_NOT_ASSESSABLE')}")
  lc.append(f"{pid}:CONCORD={p.phase16_phase13_concord_status};SATURN={p.phase16_phase12_saturn_status};singIST={p.phase16_phase11_singist_status}")
 dl,dd=le_state(comps,ps,"Dog");hl,hd=le_state(comps,ps,"Human")
 cross="NOT_ASSESSABLE" if "NOT_ASSESSABLE" in (dl,hl) else "BOTH" if (dl,hl)==("YES","YES") else "DOG_ONLY" if dl=="YES" else "HUMAN_ONLY" if hl=="YES" else "NONE"
 od=";".join(f"{g}={OM[g]}" for g in comps);entity="COMPLEX" if len(comps)>1 else "ONE_TO_ONE" if OM[comps[0]]!="NOT_ASSESSABLE" else "NOT_ASSESSABLE"
 cr=PROV[((PROV.target_entity_id==cid)|PROV.target_component_if_applicable.isin(comps))&PROV.source_species.str.lower().str.contains("dog|canine",regex=True)]
 independent=" | ".join(cr.short_evidence_statement) if len(cr) else "No independent canine target-level evidence identified in the frozen evidence audit"
 kind509_present=bool(((PROV.target_entity_id=="TNF")&(PROV.drug_if_applicable.str.upper()=="KIND-509")&(PROV.evidence_level=="CANINE_INTERVENTIONAL_PILOT_SIGNAL")).any())
 support_text=str(cfg.SUPPORT);limit_text=str(cfg.LIMIT)
 tt=trans(cfg.TARGET_TRANSLATION) if track=="TRACK_A" else "NOT_ASSESSABLE"
 sd="UNKNOWN" if track=="TRACK_A" and cfg.SPECIFIC_DRUG in ["NOT_ASSESSABLE","UNKNOWN"] else trans(cfg.SPECIFIC_DRUG)
 if cid=="TNF" and kind509_present:
  pharm="EMERGING";efficacy="TESTED"
  interventional="PILOT_POSITIVE_SIGNAL_IN_HISTORICAL_CANINE_IBD"
  evidence_limitation="small pilot study;historical canine IBD definition;not established standard therapy for contemporary CIE"
  support_text=support_text+"; preliminary canine KIND-509 randomized blinded placebo-controlled pilot intervention signal"
  limit_text="Dog target-gene evidence limited; generic inflammation; Epithelial mt-QC sensitivity; small pilot study under a historical canine IBD definition; not established standard therapy for contemporary CIE; specific human-biologic cross-reactivity remains unknown"
 else:
  pharm="EMERGING" if cfg.PHARM=="EMERGING" else "UNKNOWN" if track=="TRACK_A" or cid=="TNF" else "NOT_APPLICABLE"
  efficacy="UNTESTED" if track=="TRACK_A" else "NOT_ASSESSABLE"
  interventional="NOT_APPLICABLE";evidence_limitation="NOT_APPLICABLE"
 return {"ANALYSIS_RUN_ID":RUN_ID,"TARGET_ENTITY_ID":cid,"TARGET_ENTITY_NAME":base.TARGET_ENTITY_NAME,"TARGET_ENTITY_TYPE":base.TARGET_ENTITY_TYPE,"TARGET_ENTITY_COMPONENTS":base.TARGET_ENTITY_COMPONENTS,"ENTRY_TRACK":track,"ENTRY_REASON":reason,"CANDIDATE_DISCOVERY_SOURCE":disc,"SOURCE_VERSION_OR_DATE":base.get("SOURCE_VERSION_OR_DATE","Frozen program assets through 2026-09-16"),"SEARCH_DATE":base.get("SEARCH_DATE",CUTOFF),"SOURCE_PROGRAM":";".join(ps),"SOURCE_LINEAGE":lins,"ORTHOLOGY_STATUS":od,"ENTITY_ORTHOLOGY_STATUS":entity,"HUMAN_THERAPEUTIC_GATE_STATUS":hgate,"DIRECT_TARGET_GATE_STATUS":direct,"PROGRAM_CONCORDANCE_GATE_STATUS":cfg.PROGRAM_GATE,"LINEAGE_GATE_STATUS":cfg.LINEAGE_GATE,"TARGET_LEVEL_GATE_STATUS":cfg.TARGET_GATE,"THERAPEUTIC_DIRECTION_GATE_STATUS":cfg.THERAPEUTIC_GATE,"ROBUSTNESS_GATE_STATUS":cfg.ROBUSTNESS_GATE,"SUPPORTING_EVIDENCE":support_text,"CONTRADICTORY_OR_LIMITING_EVIDENCE":limit_text,"DOWNGRADE_REASON":cfg.DOWNGRADE,"FINAL_ROLE":cfg.FINAL_ROLE,"PRIMARY_SUPPORTING_LINEAGE":cfg.PRIMARY_LINEAGE,"LINEAGE_CONTEXT_STATUS":cfg.LINEAGE_CONTEXT,"HUMAN_IBD_THERAPEUTIC_EVIDENCE":human,"CANINE_PHARMACOLOGICAL_FEASIBILITY":pharm,"CANINE_CIE_EFFICACY":efficacy,"CANINE_INTERVENTIONAL_EVIDENCE":interventional,"EVIDENCE_LIMITATION":evidence_limitation,"TARGET_TRANSLATION_SUPPORTED":tt,"SPECIFIC_DRUG_TRANSLATION_SUPPORTED":sd,"DOG_LEADING_EDGE":dl,"HUMAN_LEADING_EDGE":hl,"CROSS_SPECIES_LEADING_EDGE_SUPPORT":cross,"DOG_LEADING_EDGE_DETAIL":dd,"HUMAN_LEADING_EDGE_DETAIL":hd,"DOG_GENE_EVIDENCE":gene_text(comps,"Dog"),"HUMAN_GENE_EVIDENCE":gene_text(comps,"Human"),"PROGRAM_LEVEL_EVIDENCE":" | ".join(pe),"REGION_EVIDENCE":" | ".join(re),"MT_QC_EVIDENCE":" | ".join(mt),"CED_EVIDENCE":" | ".join(ce),"LINEAGE_CONTEXT_EVIDENCE":" | ".join(lc),"INDEPENDENT_CANINE_EVIDENCE":independent,"DOG_TRANSCRIPTOMIC_EXTERNAL_REPLICATION":"NOT_AVAILABLE","UC_ICI_DISEASE_CONTEXT":"NOT_ASSESSABLE","EVIDENCE_COUNTING_CLASS":"PRIMARY_CROSS_SPECIES_DISEASE_EVIDENCE plus non-duplicated robustness, lineage-context, independent canine biology, Human therapeutic, and canine pharmacology layers"}

rows=[]
for _,r in A.iterrows():rows.append(build("TRACK_A",r))
for _,r in B.iterrows():rows.append(build("TRACK_B",r))
F=pd.DataFrame(rows)
if len(F)!=33:raise RuntimeError("FINAL_ROW_COUNT_MISMATCH")
wt(BASE/"05_FINAL_RESULTS"/"ALL_TARGETS_FINAL.tsv",F)
wt(BASE/"05_FINAL_RESULTS"/"TRACK_A_THERAPEUTIC_TRANSLATION.tsv",F[F.ENTRY_TRACK=="TRACK_A"])
wt(BASE/"05_FINAL_RESULTS"/"TRACK_B_MOLECULAR_NODES.tsv",F[F.ENTRY_TRACK=="TRACK_B"])
PRI=F[F.FINAL_ROLE.isin(["PRIORITIZED_THERAPEUTIC_TARGET","PRIORITIZED_MOLECULAR_NODE"])]
wt(BASE/"05_FINAL_RESULTS"/"FINAL_PRIORITIZED_TARGETS.tsv",PRI)

AU=[]
for _,r in F.iterrows():
 AU.append({"ANALYSIS_RUN_ID":RUN_ID,"target":r.TARGET_ENTITY_ID,"ENTRY_TRACK":r.ENTRY_TRACK,"ENTRY_REASON":r.ENTRY_REASON,"SOURCE_PROGRAM":r.SOURCE_PROGRAM,"SOURCE_LINEAGE":r.SOURCE_LINEAGE,"ORTHOLOGY_STATUS":r.ENTITY_ORTHOLOGY_STATUS,"PASS_FAIL_HUMAN_THERAPEUTIC_GATE":passfail(r.HUMAN_THERAPEUTIC_GATE_STATUS),"PASS_FAIL_DIRECT_TARGET_GATE":passfail(r.DIRECT_TARGET_GATE_STATUS),"PASS_FAIL_PROGRAM_CONCORDANCE_GATE":passfail(r.PROGRAM_CONCORDANCE_GATE_STATUS),"PASS_FAIL_LINEAGE_GATE":passfail(r.LINEAGE_GATE_STATUS),"PASS_FAIL_TARGET_LEVEL_GATE":passfail(r.TARGET_LEVEL_GATE_STATUS),"PASS_FAIL_THERAPEUTIC_DIRECTION_GATE":passfail(r.THERAPEUTIC_DIRECTION_GATE_STATUS),"PASS_FAIL_ROBUSTNESS_GATE":passfail(r.ROBUSTNESS_GATE_STATUS),"SUPPORTING_EVIDENCE":r.SUPPORTING_EVIDENCE,"CONTRADICTORY_OR_LIMITING_EVIDENCE":r.CONTRADICTORY_OR_LIMITING_EVIDENCE,"DOWNGRADE_REASON":r.DOWNGRADE_REASON,"FINAL_ROLE":r.FINAL_ROLE})
AUD=pd.DataFrame(AU);wt(BASE/"04_DECISION_AUDIT"/"target_decision_audit.tsv",AUD)
fig=PRI[["TARGET_ENTITY_ID","FINAL_ROLE","ENTRY_TRACK","PRIMARY_SUPPORTING_LINEAGE","SOURCE_PROGRAM","HUMAN_IBD_THERAPEUTIC_EVIDENCE","PROGRAM_CONCORDANCE_GATE_STATUS","TARGET_LEVEL_GATE_STATUS","INDEPENDENT_CANINE_EVIDENCE","CONTRADICTORY_OR_LIMITING_EVIDENCE","TARGET_TRANSLATION_SUPPORTED","SPECIFIC_DRUG_TRANSLATION_SUPPORTED"]].copy()
fig.columns=["target","final_role","track","primary_lineage","conserved_program","human_therapeutic_level","program_concordance","target_level_support","canine_independent_support","main_limitation","target_translation","specific_drug_translation"]
wt(BASE/"06_FIGURE_TABLES"/"Figure_Target_Prioritization.tsv",fig)# Concise evidence cards and biological-language handoff.
cards=["# Target Evidence Cards","",f"Run: `{RUN_ID}`",""]
for _,r in PRI.iterrows():
 cards += [f"## {r.TARGET_ENTITY_ID} — {r.ENTRY_TRACK}","",f"- Final role: {r.FINAL_ROLE}",f"- Primary lineage: {r.PRIMARY_SUPPORTING_LINEAGE}",f"- Anchored program: {r.SOURCE_PROGRAM}",f"- Human therapeutic evidence: {r.HUMAN_IBD_THERAPEUTIC_EVIDENCE}",f"- Dog–Human biological evidence: {r.PROGRAM_CONCORDANCE_GATE_STATUS}; leading edge={r.CROSS_SPECIES_LEADING_EDGE_SUPPORT}",f"- Target-level evidence: {r.TARGET_LEVEL_GATE_STATUS}",f"- Robustness: {r.ROBUSTNESS_GATE_STATUS}; mt-QC={r.MT_QC_EVIDENCE}; region={r.REGION_EVIDENCE}",f"- Independent canine evidence: {r.INDEPENDENT_CANINE_EVIDENCE}",f"- Canine interventional evidence: {r.CANINE_INTERVENTIONAL_EVIDENCE}",f"- Evidence limitation: {r.EVIDENCE_LIMITATION}",f"- Main supporting evidence: {r.SUPPORTING_EVIDENCE}",f"- Main limitation: {r.CONTRADICTORY_OR_LIMITING_EVIDENCE}",f"- Target translation status: {r.TARGET_TRANSLATION_SUPPORTED}",f"- Specific-drug translation status: {r.SPECIFIC_DRUG_TRANSLATION_SUPPORTED}","- Future validation interpretation: Candidate for prospective canine validation; not a validated canine treatment.",""]
(BASE/"05_FINAL_RESULTS"/"TARGET_EVIDENCE_CARDS.md").write_text("\n".join(cards)+"\n",encoding="utf-8")
summary=f'''# Final Biological Summary

Run `{RUN_ID}` applied the frozen v1.0 rules to all 9 Track A and 24 Track B entries without adding candidates or rerunning upstream analyses.

## Human Crohn targets most suitable for canine validation

TNF and JAK1 met the frozen criteria for prioritized therapeutic targets. TNF was recovered as a prioritized translational target and is additionally supported by preliminary canine interventional evidence from a randomized, blinded, placebo-controlled pilot study conducted under a historical canine IBD framework. This supports TNF as a positive translational benchmark rather than a newly proposed canine CIE target. JAK1 remains supported chiefly by concordant Myeloid and T-cell cytokine/JAK-STAT activity.

## Main limitations

The KIND-509 evidence is a small sponsor-reported pilot in 10 dogs diagnosed under a historical canine IBD definition. It does not establish anti-TNF as a mature standard treatment for contemporary canine CIE and does not establish the dose, safety, pharmacokinetics, immunogenicity, target engagement, or cross-reactivity of any specific human anti-TNF biologic in dogs. TNF also retains limited direct canine gene-level evidence, Epithelial mt-QC sensitivity, and CeD sharing. JAK1 remains program-anchored, with incomplete direct target-level evidence and untested canine CIE efficacy.

## Conserved mechanistic nodes

TNF, JAK1, STAT3, and the Epithelial mitochondrial component NDUFS1 met the frozen criteria for prioritized molecular nodes. STAT3 is supported principally by the conserved immune JAK-STAT program. NDUFS1 anchors the robust Epithelial mitochondrial-respiration decrease, but it is a mechanistic node rather than a mature therapeutic recommendation.

## Human targets not promoted for canine translation

The p40, p19, and TL1A targets remain exploratory because target-specific canine evidence is limited or the frozen IL-23-specific axis is not assessable. Alpha4beta7, alpha4, S1PR1/5, and ITGB7 were not prioritized because the frozen canine lineage/program evidence did not establish target-level translational conservation.

## Boundaries

No independent canine transcriptomic replication dataset was available, and formal UC/ICI disease-context analysis was not assessable. Region, mt-QC, CeD, CONCORD, SATURN, and singIST were retained as distinct evidence layers and were not counted as independent replications. No result validates a canine therapeutic, supports direct use of a human biologic in dogs, or establishes a universal intestinal inflammatory mechanism.
'''
(BASE/"05_FINAL_RESULTS"/"FINAL_BIOLOGICAL_SUMMARY.md").write_text(summary,encoding="utf-8")

# Hash only files read or explicitly checked by this run.
support=[snap/"TRACK_A_FROZEN_CANDIDATES.tsv",snap/"TRACK_A_FROZEN_CANDIDATES.sha256",snap/"TARGET_PROGRAM_LINEAGE_ANCHOR_SNAPSHOT.tsv",snap/"TARGET_PROGRAM_LINEAGE_ANCHOR_SNAPSHOT.sha256",snap/"TRACK_B_FROZEN_CANDIDATES.tsv",snap/"TRACK_B_FROZEN_CANDIDATES.sha256",BASE/"00_DESIGN"/"DECISION_CLASSIFICATION_V1_0_FROZEN.tsv",BASE/"00_DESIGN"/"DECISION_CLASSIFICATION_V1_0_FROZEN.sha256",Path(__file__),prov_path,p_program,ROOT/"16_final_integration"/"master"/"FINAL_PROGRAM_EVIDENCE_MASTER.csv",ROOT/"08_results"/"supplementary_robustness"/"epithelial_mt_qc_sensitivity"/"10_final_report"/"EPITHELIAL_MT_QC_SENSITIVITY_SUMMARY.csv",ROOT/"14_ced_specificity"/"specificity"/"DOG_CROHN_CED_PROGRAM_MATRIX.csv",ROOT/"13_concord"/"summary"/"CROSS_METHOD_LINEAGE_SUMMARY.csv",ROOT/"12_saturn"/"macrogenes"/"MACROGENE_CONCORDANCE_SUMMARY.csv",ROOT/"11_singIST"/"summary"/"SINGIST_CORE_GENE_DRIVERS.csv",ROOT/"08_results"/"supplementary_robustness"/"dog_external_dataset_availability_audit"/"01_DOG_EXTERNAL_CANDIDATE_AUDIT.csv",ROOT/"CURRENT_STATE.md",ROOT/"RESULTS_REGISTER.md",p_ortho]+de_paths+le_paths
man=[]
for p in support:
 if not p.exists():raise RuntimeError("DATA_FILE_MISSING:"+str(p))
 s=p.stat();low=str(p).lower()
 phase="TARGET_PRIORITIZATION_V1_0" if "target_prioritization_v1_0" in low else "Phase16" if "16_final_integration" in low else "Phase14" if "14_ced" in low else "Phase13" if "13_concord" in low else "Phase12" if "12_saturn" in low else "Phase11" if "11_singist" in low else "Phase10" if "10_cross" in low else "Phase7" if "07_human_crohn" in low else "Phase3-5" if "06_pseudobulk" in low else "Supplementary/frozen project state"
 man.append({"file_path":str(p),"file_name":p.name,"file_size":s.st_size,"modification_time":datetime.fromtimestamp(s.st_mtime).astimezone().isoformat(),"checksum_sha256":sha(p),"source_phase":phase,"source_status":"FROZEN_OR_FINAL_READ_ONLY","version_if_available":"project-frozen/current approved","purpose":"Candidate identity, decision configuration, gene/program/robustness/context evidence, provenance, or workflow verification"})
M=pd.DataFrame(man);wt(BASE/"01_INPUT_MANIFEST"/"input_manifest.tsv",M)
(BASE/"01_INPUT_MANIFEST"/"input_manifest.json").write_text(json.dumps(man,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")

roles=Counter(F.FINAL_ROLE)
allowed={"PRIORITIZED_THERAPEUTIC_TARGET","PRIORITIZED_MOLECULAR_NODE","MECHANISTIC_CANDIDATE","EXPLORATORY_TARGET","NOT_PRIORITIZED_FOR_CANINE_TRANSLATION"}
checks={"snapshot_hashes_match":True,"track_a_rows":int((F.ENTRY_TRACK=="TRACK_A").sum()),"track_b_rows":int((F.ENTRY_TRACK=="TRACK_B").sum()),"all_rows":len(F),"audit_rows":len(AUD),"prioritized_rows":len(PRI),"roles_permitted":set(roles).issubset(allowed),"support_and_limit_complete":bool((F.SUPPORTING_EVIDENCE.str.len()>0).all() and (F.CONTRADICTORY_OR_LIMITING_EVIDENCE.str.len()>0).all()),"leading_edge_vocab_valid":bool(F.DOG_LEADING_EDGE.isin(["YES","NO","NOT_ASSESSABLE"]).all() and F.HUMAN_LEADING_EDGE.isin(["YES","NO","NOT_ASSESSABLE"]).all() and F.CROSS_SPECIES_LEADING_EDGE_SUPPORT.isin(["BOTH","DOG_ONLY","HUMAN_ONLY","NONE","NOT_ASSESSABLE"]).all()),"canine_pharmacology_vocab_valid":bool(F.CANINE_PHARMACOLOGICAL_FEASIBILITY.isin(["ESTABLISHED","EMERGING","UNKNOWN","NOT_APPLICABLE"]).all()),"target_translation_vocab_valid":bool(F.TARGET_TRANSLATION_SUPPORTED.isin(["YES","PARTIAL","NO","NOT_ASSESSABLE"]).all()),"specific_drug_vocab_valid":bool(F.SPECIFIC_DRUG_TRANSLATION_SUPPORTED.isin(["YES","NO","UNKNOWN","NOT_ASSESSABLE"]).all()),"external_replication_fixed":bool((F.DOG_TRANSCRIPTOMIC_EXTERNAL_REPLICATION=="NOT_AVAILABLE").all()),"uc_ici_fixed":bool((F.UC_ICI_DISEASE_CONTEXT=="NOT_ASSESSABLE").all())}
tnf_check=F[F.TARGET_ENTITY_ID=="TNF"]
checks["kind509_provenance_present"]=bool(((PROV.target_entity_id=="TNF")&(PROV.drug_if_applicable.str.upper()=="KIND-509")).any())
checks["tnf_efficacy_tested"]=bool((tnf_check.CANINE_CIE_EFFICACY=="TESTED").all())
checks["tnf_interventional_evidence_limited"]=bool((tnf_check.CANINE_INTERVENTIONAL_EVIDENCE=="PILOT_POSITIVE_SIGNAL_IN_HISTORICAL_CANINE_IBD").all())
checks["pass"]=checks["track_a_rows"]==9 and checks["track_b_rows"]==24 and checks["all_rows"]==33 and checks["audit_rows"]==33 and checks["prioritized_rows"]==6 and all(v for k,v in checks.items() if isinstance(v,bool) and k!="pass")
if not checks["pass"]:raise RuntimeError("FINAL_INTEGRITY_CHECK_FAILED:"+json.dumps(checks))
(BASE/"07_LOGS"/"FINAL_INTEGRITY_CHECK_RUN004.json").write_text(json.dumps({"ANALYSIS_RUN_ID":RUN_ID,"checks":checks,"role_counts":dict(roles)},indent=2)+"\n",encoding="utf-8")
meta_path=BASE/"00_DESIGN"/"RUN_METADATA.json";meta=json.loads(meta_path.read_text(encoding="utf-8-sig"));meta.update({"DOG_TARGET_RESULTS_UNLOCKED":True,"RUN_END_TIME":datetime.now().astimezone().isoformat(),"RUN_STATUS":"COMPLETE_PENDING_CACHE_CLEANUP","TRACK_A_CANDIDATES":9,"TRACK_B_CANDIDATES":24,"PRIORITIZED_THERAPEUTIC_TARGETS":roles["PRIORITIZED_THERAPEUTIC_TARGET"],"PRIORITIZED_MOLECULAR_NODES":roles["PRIORITIZED_MOLECULAR_NODE"],"MECHANISTIC_CANDIDATES":roles["MECHANISTIC_CANDIDATE"],"EXPLORATORY_TARGETS":roles["EXPLORATORY_TARGET"],"NOT_PRIORITIZED_FOR_CANINE_TRANSLATION":roles["NOT_PRIORITIZED_FOR_CANINE_TRANSLATION"],"NOT_ASSESSABLE_PROJECT_LEVEL_ITEMS":"UC_ICI_DISEASE_CONTEXT;DOG_TRANSCRIPTOMIC_EXTERNAL_REPLICATION_UNAVAILABLE","FINAL_STATUS":"PENDING_CACHE_CLEANUP_AND_FINAL_FREEZE"});meta_path.write_text(json.dumps(meta,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
(BASE/"07_LOGS"/"RUN005_EXECUTION_LOG.txt").write_text(f"ANALYSIS_RUN_ID={RUN_ID}\nDESIGN_VERSION={DESIGN_VERSION}\nSCRIPT_WORKFLOW_VERSION={WORKFLOW_VERSION}\nTRACK_A_CANDIDATES=9\nTRACK_B_CANDIDATES=24\nALL_TARGET_ROWS=33\nPRIORITIZED_ROWS={len(PRI)}\nSTATUS=FORMAL_OUTPUTS_COMPLETE_PENDING_CACHE_CLEANUP_AND_FINAL_FREEZE\n",encoding="utf-8")
print(json.dumps({"run_id":RUN_ID,"rows":len(F),"prioritized":len(PRI),"roles":dict(roles),"integrity":checks["pass"]}))