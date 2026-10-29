from pathlib import Path
import json
import pandas as pd

BASE=(Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"]) / "TARGET_PRIORITIZATION_V1_0")
OLD=BASE/"07_LOGS"/"INVALIDATED_RUN_004"/"05_FINAL_RESULTS"/"ALL_TARGETS_FINAL.tsv"
NEW=BASE/"05_FINAL_RESULTS"/"ALL_TARGETS_FINAL.tsv"
OUT=BASE/"07_LOGS"/"EVIDENCE_OMISSION_RERUN_COMPARISON.tsv"
old=pd.read_csv(OLD,sep="\t",dtype=str,keep_default_na=False)
new=pd.read_csv(NEW,sep="\t",dtype=str,keep_default_na=False)
for df in (old,new):df["candidate_key"]=df["ENTRY_TRACK"]+"::"+df["TARGET_ENTITY_ID"]
if len(old)!=33 or len(new)!=33 or set(old.candidate_key)!=set(new.candidate_key):raise RuntimeError("OLD_NEW_CANDIDATE_SET_MISMATCH")
o=old.set_index("candidate_key");n=new.set_index("candidate_key")
comparison=[]
for key in old.candidate_key:
 old_role=o.loc[key,"FINAL_ROLE"];new_role=n.loc[key,"FINAL_ROLE"];changed="YES" if old_role!=new_role else "NO"
 if key in ("TRACK_A::TNF","TRACK_B::TNF"):
  reason="KIND-509 evidence updated canine intervention fields; frozen gates and final role unchanged" if changed=="NO" else "KIND-509 evidence changed role under unchanged frozen rules"
 else:reason="No role change; frozen candidate, gate, and classification criteria unchanged" if changed=="NO" else "UNEXPLAINED_ROLE_CHANGE"
 comparison.append({"candidate":key,"old_final_role":old_role,"new_final_role":new_role,"changed_yes_no":changed,"change_reason":reason})
pd.DataFrame(comparison).to_csv(OUT,sep="\t",index=False,lineterminator="\n")
expected_tnf={"ANALYSIS_RUN_ID","SUPPORTING_EVIDENCE","CONTRADICTORY_OR_LIMITING_EVIDENCE","CANINE_PHARMACOLOGICAL_FEASIBILITY","CANINE_CIE_EFFICACY","INDEPENDENT_CANINE_EVIDENCE"}
shared=[c for c in old.columns if c in new.columns and c!="candidate_key"]
unexplained=[]
for key in old.candidate_key:
 for col in shared:
  if o.loc[key,col]!=n.loc[key,col]:
   allowed=col=="ANALYSIS_RUN_ID" or (key.endswith("::TNF") and col in expected_tnf)
   if not allowed:unexplained.append({"candidate":key,"column":col,"old":o.loc[key,col],"new":n.loc[key,col]})
role_changes=sum(x["changed_yes_no"]=="YES" for x in comparison)
other_changes=sum(x["changed_yes_no"]=="YES" and not x["candidate"].endswith("::TNF") for x in comparison)
tnf_a=n.loc["TRACK_A::TNF"]
tnf_b=n.loc["TRACK_B::TNF"]
checks={"old_rows":len(old),"new_rows":len(new),"candidate_set_identical":set(old.candidate_key)==set(new.candidate_key),"role_changes":role_changes,"other_role_changes":other_changes,"unexplained_non_evidence_changes":unexplained,"tnf_track_a_role":tnf_a.FINAL_ROLE,"tnf_track_b_role":tnf_b.FINAL_ROLE,"tnf_efficacy_all_tested":tnf_a.CANINE_CIE_EFFICACY=="TESTED" and tnf_b.CANINE_CIE_EFFICACY=="TESTED","tnf_interventional_status_all_limited":tnf_a.CANINE_INTERVENTIONAL_EVIDENCE=="PILOT_POSITIVE_SIGNAL_IN_HISTORICAL_CANINE_IBD" and tnf_b.CANINE_INTERVENTIONAL_EVIDENCE=="PILOT_POSITIVE_SIGNAL_IN_HISTORICAL_CANINE_IBD"}
checks["pass"]=checks["candidate_set_identical"] and role_changes==0 and other_changes==0 and not unexplained and checks["tnf_track_a_role"]=="PRIORITIZED_THERAPEUTIC_TARGET" and checks["tnf_track_b_role"]=="PRIORITIZED_MOLECULAR_NODE" and checks["tnf_efficacy_all_tested"] and checks["tnf_interventional_status_all_limited"]
(BASE/"07_LOGS"/"EVIDENCE_OMISSION_RERUN_VALIDATION.json").write_text(json.dumps(checks,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
if not checks["pass"]:raise RuntimeError("EVIDENCE_CORRECTION_COMPARISON_FAILED:"+json.dumps(checks))
print(json.dumps(checks))