from __future__ import annotations
import json
from datetime import datetime
from pathlib import Path
import pandas as pd

BASE=(Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"]) / "TARGET_PRIORITIZATION_V1_0")
PROVENANCE=BASE/"03_EVIDENCE"/"evidence_provenance.tsv"
OLD_SCRIPT=BASE/"00_DESIGN"/"run_target_prioritization_v1_0_3.py"
NEW_SCRIPT=BASE/"00_DESIGN"/"run_target_prioritization_v1_0_4.py"
NEW_RUN="TP_V1_0_20260916_005"

def update_provenance():
 df=pd.read_csv(PROVENANCE,sep="\t",dtype=str,keep_default_na=False)
 for col in ["relationship_to_modern_CIE","evidence_interpretation"]:
  if col not in df.columns:df[col]=""
 if (df["drug_if_applicable"].str.upper()=="KIND-509").any():raise RuntimeError("KIND509_ALREADY_PRESENT_BEFORE_CORRECTION")
 row={c:"" for c in df.columns}
 row.update({"target_entity_id":"TNF","target_entity_name":"Tumor necrosis factor","target_entity_type":"LIGAND","target_component_if_applicable":"TNF","drug_if_applicable":"KIND-509","direct_target_entity":"canine TNF-alpha","target_selectivity":"caninized monoclonal anti-TNF-alpha antibody","mechanism_of_action":"neutralization of canine TNF-alpha","disease":"historical canine IBD","evidence_level":"CANINE_INTERVENTIONAL_PILOT_SIGNAL","study_type":"randomized blinded placebo-controlled pilot field study","publication_year":"2021","source_title":"Kindred Biosciences 2020 Form 10-K KIND-509 results; TNF-alpha and NGF antibodies for veterinary use patent","source_species":"Dog","source_version_or_date":"SEC Form 10-K FY2020, accession 0001561743-21-000024; US20240101715A1","access_date":"2026-09-16","evidence_cutoff_date":"2026-09-16","short_evidence_statement":"A randomized, blinded, placebo-controlled pilot effectiveness study enrolled 10 dogs diagnosed under a historical canine IBD framework and reported a preliminary efficacy signal for the caninized anti-TNF-alpha antibody KIND-509 over four weeks.","supports_or_limits":"SUPPORTS: preliminary canine interventional efficacy signal; LIMITS: small n=10 sponsor-reported pilot, historical canine IBD definition, not an established standard therapy for contemporary CIE, and does not establish direct usability of a human anti-TNF biologic in dogs","source_url_or_identifier":"https://www.sec.gov/Archives/edgar/data/1561743/000156174321000024/kin-20201231.htm;https://patents.google.com/patent/US20240101715A1/en","relationship_to_modern_CIE":"overlapping/compatible spectrum, but not equivalent to all contemporary CIE","evidence_interpretation":"preliminary canine interventional efficacy signal"})
 out=pd.concat([df,pd.DataFrame([row],columns=df.columns)],ignore_index=True)
 out.to_csv(PROVENANCE,sep="\t",index=False,lineterminator="\n")

def update_workflow_script():
 text=OLD_SCRIPT.read_text(encoding="utf-8")
 text=text.replace('RUN_ID="TP_V1_0_20260916_004"',f'RUN_ID="{NEW_RUN}"')
 text=text.replace('WORKFLOW_VERSION="Target_Prioritization_Workflow_v1.0.3"','WORKFLOW_VERSION="Target_Prioritization_Workflow_v1.0.4"')
 text=text.replace('RUN004_EXECUTION_LOG.txt','RUN005_EXECUTION_LOG.txt')
 old=''' cr=PROV[((PROV.target_entity_id==cid)|PROV.target_component_if_applicable.isin(comps))&PROV.source_species.str.lower().str.contains("dog|canine",regex=True)]
 independent=" | ".join(cr.short_evidence_statement) if len(cr) else "No independent canine target-level evidence identified in the frozen evidence audit"
 tt=trans(cfg.TARGET_TRANSLATION) if track=="TRACK_A" else "NOT_ASSESSABLE"
 sd="UNKNOWN" if track=="TRACK_A" and cfg.SPECIFIC_DRUG in ["NOT_ASSESSABLE","UNKNOWN"] else trans(cfg.SPECIFIC_DRUG)
 pharm="EMERGING" if cfg.PHARM=="EMERGING" else "UNKNOWN" if track=="TRACK_A" or cid=="TNF" else "NOT_APPLICABLE"
 efficacy="UNTESTED" if track=="TRACK_A" else "NOT_ASSESSABLE"
'''
 new=''' cr=PROV[((PROV.target_entity_id==cid)|PROV.target_component_if_applicable.isin(comps))&PROV.source_species.str.lower().str.contains("dog|canine",regex=True)]
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
'''
 if old not in text:raise RuntimeError("WORKFLOW_EVIDENCE_BLOCK_NOT_FOUND")
 text=text.replace(old,new)
 text=text.replace('"SUPPORTING_EVIDENCE":cfg.SUPPORT,"CONTRADICTORY_OR_LIMITING_EVIDENCE":cfg.LIMIT','"SUPPORTING_EVIDENCE":support_text,"CONTRADICTORY_OR_LIMITING_EVIDENCE":limit_text')
 text=text.replace('"CANINE_CIE_EFFICACY":efficacy,"TARGET_TRANSLATION_SUPPORTED"','"CANINE_CIE_EFFICACY":efficacy,"CANINE_INTERVENTIONAL_EVIDENCE":interventional,"EVIDENCE_LIMITATION":evidence_limitation,"TARGET_TRANSLATION_SUPPORTED"')
 text=text.replace('f"- Independent canine evidence: {r.INDEPENDENT_CANINE_EVIDENCE}",f"- Main supporting evidence:','f"- Independent canine evidence: {r.INDEPENDENT_CANINE_EVIDENCE}",f"- Canine interventional evidence: {r.CANINE_INTERVENTIONAL_EVIDENCE}",f"- Evidence limitation: {r.EVIDENCE_LIMITATION}",f"- Main supporting evidence:')
 text=text.replace("TNF and JAK1 met the frozen criteria for prioritized therapeutic targets. Their strongest support comes from concordant immune-lineage programs, especially Myeloid and T-cell TNF/NF-kappaB and cytokine/JAK-STAT activity. This supports testing the biological targets in canine CIE; it does not establish efficacy, dosing, safety, or cross-reactivity of a specific human drug in dogs.","TNF and JAK1 met the frozen criteria for prioritized therapeutic targets. TNF was recovered as a prioritized translational target and is additionally supported by preliminary canine interventional evidence from a randomized, blinded, placebo-controlled pilot study conducted under a historical canine IBD framework. This supports TNF as a positive translational benchmark rather than a newly proposed canine CIE target. JAK1 remains supported chiefly by concordant Myeloid and T-cell cytokine/JAK-STAT activity.")
 text=text.replace("TNF has limited direct canine gene-level evidence, its Epithelial TNF/NF-kappaB result is mt-QC sensitive, the inflammatory signal is shared with CeD, and canine anti-TNF efficacy/cross-reactivity is unknown. JAK1 is program-anchored and has emerging canine pharmacological feasibility, but direct target-level transcriptional support is incomplete and canine CIE efficacy remains untested.","The KIND-509 evidence is a small sponsor-reported pilot in 10 dogs diagnosed under a historical canine IBD definition. It does not establish anti-TNF as a mature standard treatment for contemporary canine CIE and does not establish the dose, safety, pharmacokinetics, immunogenicity, target engagement, or cross-reactivity of any specific human anti-TNF biologic in dogs. TNF also retains limited direct canine gene-level evidence, Epithelial mt-QC sensitivity, and CeD sharing. JAK1 remains program-anchored, with incomplete direct target-level evidence and untested canine CIE efficacy.")
 marker='checks["pass"]='
 insert='''tnf_check=F[F.TARGET_ENTITY_ID=="TNF"]
checks["kind509_provenance_present"]=bool(((PROV.target_entity_id=="TNF")&(PROV.drug_if_applicable.str.upper()=="KIND-509")).any())
checks["tnf_efficacy_tested"]=bool((tnf_check.CANINE_CIE_EFFICACY=="TESTED").all())
checks["tnf_interventional_evidence_limited"]=bool((tnf_check.CANINE_INTERVENTIONAL_EVIDENCE=="PILOT_POSITIVE_SIGNAL_IN_HISTORICAL_CANINE_IBD").all())
'''
 if marker not in text:raise RuntimeError("WORKFLOW_CHECK_MARKER_NOT_FOUND")
 text=text.replace(marker,insert+marker,1)
 NEW_SCRIPT.write_text(text,encoding="utf-8")

def initialize_new_run():
 metadata={"ANALYSIS_RUN_ID":NEW_RUN,"DESIGN_VERSION":"Target_Prioritization_v1.0","DESIGN_STATUS":"FROZEN","EXECUTION_MODE":"EVIDENCE_OMISSION_CORRECTION_RERUN","SCRIPT_WORKFLOW_VERSION":"Target_Prioritization_Workflow_v1.0.4","RUN_START_TIME":datetime.now().astimezone().isoformat(),"THERAPEUTIC_EVIDENCE_CUTOFF_DATE":"2026-09-16","PARENT_INVALIDATED_RUN_ID":"TP_V1_0_20260916_004","INVALIDATED_RUN_HISTORY":["TP_V1_0_20260916_001","TP_V1_0_20260916_002","TP_V1_0_20260916_003","TP_V1_0_20260916_004"],"RESTART_BASIS":"VERIFIED_EVIDENCE_INPUT_OMISSION_CORRECTION_FROM_FROZEN_SNAPSHOTS","ERROR_CORRECTED":"KIND-509_CANINE_ANTI_TNF_EVIDENCE_OMISSION","DOG_TARGET_RESULTS_UNLOCKED":True,"TRACK_A_CANDIDATE_POOL_STATUS":"FROZEN_UNCHANGED","TRACK_A_CANDIDATES":9,"TARGET_PROGRAM_LINEAGE_ANCHOR_STATUS":"FROZEN_UNCHANGED","TRACK_B_CANDIDATE_POOL_STATUS":"FROZEN_UNCHANGED","TRACK_B_CANDIDATES":24,"FROZEN_SPEC_SHA256":"0fd9bd5b4831ccb3530ea3ea5050e32bf0d79be79caa2ffa396118b986bb3962","FROZEN_SNAPSHOT_REVERIFICATION_STATUS":"PENDING_SCRIPT_REVERIFICATION","RUN_STATUS":"READY_FOR_FULL_CLASSIFICATION_RERUN","FINAL_STATUS":"PENDING_EXECUTION"}
 (BASE/"00_DESIGN"/"RUN_METADATA.json").write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")

if __name__=="__main__":
 update_provenance();update_workflow_script();initialize_new_run()
 print(json.dumps({"new_run":NEW_RUN,"workflow":"Target_Prioritization_Workflow_v1.0.4","provenance_rows":len(pd.read_csv(PROVENANCE,sep="\t"))}))