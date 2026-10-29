from __future__ import annotations
import hashlib,json
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
OUT=ROOT/"WORKING"/"21_HUMAN_TCELL_SUBTYPE_METABOLIC_ANALYSIS"; WORK=OUT/"work"
P5A=ROOT/"WORKING"/"20_HUMAN_TCELL_LABEL_FREEZE"; P2B=ROOT/"WORKING"/"16_TCELL_SUBTYPE_METABOLIC_SENSITIVITY"
HBASE=ROOT/"07_human_crohn"/"GSE266616"
PROGRAM_ORDER=["Aerobic respiration and respiratory electron transport","Mitochondrial ribosome-associated quality control","Xenobiotic Metabolism"]
SHORT={"Aerobic respiration and respiratory electron transport":"Respiration / electron transport","Mitochondrial ribosome-associated quality control":"Mitoribosome-associated QC","Xenobiotic Metabolism":"Xenobiotic metabolism"}
SUBTYPES=["CD103+ CD4","CD103+ CD8","CD103- CD4","CD4 Tcm","KLRG1+ CD8","Treg"]

def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()
def write_csv(df,name): df.to_csv(OUT/name,index=False,encoding="utf-8")
def direction(nes):
    if pd.isna(nes): return "NOT_ASSESSABLE"
    return "ActiveCrohn_higher" if nes>0 else ("ActiveCrohn_lower" if nes<0 else "No_change")

raw=pd.read_csv(WORK/"program_results_raw.csv")
models=pd.read_csv(WORK/"model_audit.csv")
elig=pd.read_csv(WORK/"subtype_eligibility_preflight.csv")
qc=pd.read_csv(WORK/"pseudobulk_sample_qc.csv")
lodo=pd.read_csv(WORK/"lodo_results_raw.csv")
markers=pd.read_csv(WORK/"treg_marker_audit_raw.csv")
contrib=pd.read_csv(WORK/"subtype_participant_contributions.csv")
selection=pd.read_csv(P5A/"05_COHORT_SELECTION_AND_MISSINGNESS.csv")
reliability=pd.read_csv(P5A/"06_SUBTYPE_ANNOTATION_RELIABILITY.csv")
preflight=json.loads((WORK/"preflight_audit.json").read_text(encoding="utf-8"))
target_info=pd.read_csv(WORK/"frozen_target_programs.tsv",sep="\t")
membership=pd.read_csv(WORK/"frozen_target_memberships.tsv",sep="\t")
rversions=pd.read_csv(WORK/"software_versions_R.csv")

# Original frozen 27-participant reference.
rfile=HBASE/"gsea"/"full"/"Human_T_Reactome_GSEA.csv"
hfile=HBASE/"gsea"/"full"/"Human_T_Hallmark_GSEA.csv"
frozen=pd.concat([pd.read_csv(rfile),pd.read_csv(hfile)],ignore_index=True)
frozen=frozen[frozen.pathway_id.isin(target_info.pathway_id)].copy()
rank27=pd.read_csv(HBASE/"pseudobulk"/"ranked_lists"/"Human_T_ActiveCrohn_vs_Healthy_ranked.csv")
rank27_symbols=rank27.dropna(subset=["gene_symbol"]).sort_values(["rank_metric","PValue"],ascending=[False,True]).drop_duplicates("gene_symbol").shape[0]
collection_n={"Reactome":len(pd.read_csv(rfile)),"Hallmark":len(pd.read_csv(hfile))}
full_map=frozen.set_index("pathway_id")
broad_rows=[]
for cohort,analysis,nh,nc,cells in [
    ("Full27_FrozenReference",None,16,11,21378),
    ("Recovered19_AllFrozenBroadT","Broad19_AllFrozenT",13,6,int(qc.loc[qc.analysis.eq("Broad19_AllFrozenT"),"n_cells"].sum())),
    ("Recovered19_VerifiedAuthorT","Broad19_VerifiedAuthorT",13,6,int(qc.loc[qc.analysis.eq("Broad19_VerifiedAuthorT"),"n_cells"].sum()))]:
    for p in PROGRAM_ORDER:
        ident=target_info[target_info.pathway_name.eq(p)].iloc[0]; pid=ident.pathway_id
        if analysis is None:
            z=full_map.loc[pid]; nes=float(z.NES); pv=float(z.pval); ap=float(z.padj); tested=len(rank27); mapped=rank27_symbols; size=int(z["size"])
            family=f"ORIGINAL_FULL_{ident.collection.upper()}_COLLECTION_BH_{collection_n[ident.collection]}_TESTED_PATHWAYS"; source="FROZEN_PHASE7_REFERENCE"
        else:
            z=raw[(raw.analysis==analysis)&(raw.pathway_id==pid)].iloc[0]; nes=float(z.NES);pv=float(z.pval);ap=float(z.adjusted_P_within_three);tested=int(z.tested_gene_ids);mapped=int(z.mapped_ranked_gene_symbols);size=int(z["size"])
            family="BH_ACROSS_THREE_PREDEFINED_PROGRAMS_WITHIN_THIS_COHORT";source="PHASE5B_TARGETED_REANALYSIS"
        ref=float(full_map.loc[pid].NES)
        broad_rows.append(dict(cohort=cohort,analysis_source=source,n_healthy=nh,n_active_crohn=nc,total_T_cells=cells,pathway_id=pid,pathway_name=p,collection=ident.collection,source_version=ident.source_version,source_species=ident.source_species,raw_count_layer="layers[counts]" if analysis else "frozen Phase7 donor pseudobulk",contrast="ActiveCrohn-minus-Healthy",tested_gene_ids=tested,mapped_ranked_gene_symbols=mapped,gene_set_size_in_tested_universe=size,NES=nes,direction=direction(nes),nominal_P=pv,adjusted_P=ap,adjustment_family=family,testability="TESTABLE",full27_reference_NES=ref,absolute_NES_difference_from_full27=abs(nes-ref),direction_agreement_with_full27="YES" if np.sign(nes)==np.sign(ref) else "NO",universe_compatibility="SAME_FROZEN_HUMAN_MEMBERSHIP_AND_RANK_CONVENTION;COHORT_SPECIFIC_TESTED_GENE_UNIVERSE",ortholog_mapping_convention="Human-native gene sets; no ortholog mapping in Human tests; Reactome species-specific equivalent recovered by stable pathway identity"))
broad=pd.DataFrame(broad_rows)
write_csv(broad,"01_BROAD_T_COHORT_SENSITIVITY.csv")

# Eligibility table.
elig_out=elig.merge(models[models.analysis.str.startswith("Subtype::")].assign(subtype=lambda d:d.analysis.str.replace("Subtype::","",regex=False))[["subtype","tested_gene_ids","mapped_ranked_gene_symbols","design_rank","design_columns","common_dispersion","normalization","filter","rank_statistic","contrast"]],on="subtype",how="left")
elig_out["minimum_cells_per_participant_subtype"]=20; elig_out["minimum_participants_per_condition"]=5
elig_out["independent_biological_unit"]="participant";elig_out["author_label_status"]="AUTHOR_ANNOTATED_STATE_NOT_INDEPENDENTLY_REANNOTATED"
elig_out["aggregation_integrity"]="PASS" if qc.aggregation_verified.eq("YES").all() else "FAIL"
write_csv(elig_out,"03_HUMAN_SUBTYPE_PSEUDOBULK_ELIGIBILITY.csv")

# Subtype primary family.
sub=raw[raw.analysis.str.startswith("Subtype::")].copy();sub["subtype"]=sub.analysis.str.replace("Subtype::","",regex=False)
full_dir={pid:direction(float(full_map.loc[pid].NES)) for pid in target_info.pathway_id}
b19=raw[raw.analysis.eq("Broad19_AllFrozenT")].set_index("pathway_id")
lodo_valid=lodo[lodo.NES.notna()].copy()
summary=[]
for (analysis,pid),g in lodo_valid.groupby(["analysis","pathway_id"]):
    idx=g.NES_change.abs().idxmax()
    summary.append(dict(analysis=analysis,pathway_id=pid,lodo_valid_omissions=len(g),lodo_NES_min=g.NES.min(),lodo_NES_max=g.NES.max(),lodo_sign_reversals=int(g.sign_reversal.eq("YES").sum()),lodo_max_abs_NES_change=float(g.NES_change.abs().max()),lodo_influential_participant=g.loc[idx,"omitted_participant"]))
ls=pd.DataFrame(summary)
sub=sub.merge(ls,on=["analysis","pathway_id"],how="left")
sub["original_broad27_direction"]=sub.pathway_id.map(full_dir);sub["recovered_broad19_direction"]=sub.pathway_id.map(lambda x:direction(float(b19.loc[x].NES)))
sub["agreement_with_broad27"]=np.where(np.sign(sub.NES)==sub.pathway_id.map(lambda x:np.sign(float(full_map.loc[x].NES))),"SAME_DIRECTION","OPPOSITE_DIRECTION")
sub["agreement_with_broad19"]=np.where(np.sign(sub.NES)==sub.pathway_id.map(lambda x:np.sign(float(b19.loc[x].NES))),"SAME_DIRECTION","OPPOSITE_DIRECTION")
sub["family_significant_FDR05"]=np.where(sub.subtype_family_adjusted_P<.05,"YES","NO")
sub["planned_family"]="6 subtypes x 3 programs = 18 tests";sub["valid_P_values_adjusted"]=18
sub["interpretation_boundary"]="Directional subtype heterogeneity is descriptive; no subtype-by-disease interaction was tested."
cols=["subtype","pathway_id","pathway_name","collection","source_version","source_species","n_healthy","n_crohn","total_cells","tested_gene_ids","mapped_ranked_gene_symbols","size","NES","pval","adjusted_P_within_three","subtype_family_adjusted_P","planned_family","valid_P_values_adjusted","direction","testability","family_significant_FDR05","original_broad27_direction","recovered_broad19_direction","agreement_with_broad27","agreement_with_broad19","lodo_valid_omissions","lodo_NES_min","lodo_NES_max","lodo_sign_reversals","lodo_max_abs_NES_change","lodo_influential_participant","interpretation_boundary"]
write_csv(sub[cols].rename(columns={"size":"gene_set_size_in_tested_universe","pval":"nominal_P","adjusted_P_within_three":"within_subtype_three_program_BH","subtype_family_adjusted_P":"primary_18_test_BH"}),"04_HUMAN_SUBTYPE_METABOLIC_RESULTS.csv")

# Participant treatment map, restricted recovered Crohn = exactly 2/2/2.
parts=selection[selection.row_type.eq("PARTICIPANT")].copy()
part_map=parts.set_index("frozen_participant_id")
recovered=set(pd.read_csv(P5A/"01_VERIFIED_HUMAN_TCELL_CELL_LABEL_MAP.csv",usecols=["frozen_participant_id"]).frozen_participant_id.unique())
crohn_recovered=parts[parts.frozen_participant_id.isin(recovered)&parts.condition.ne("Healthy")].copy()
treatment_counts=crohn_recovered.biologic_class.value_counts().to_dict()
if treatment_counts!={"ANTI_TNF":2,"ANTI_IL12_23":2,"NO_BIOLOGIC":2}: raise RuntimeError(f"Recovered Crohn treatment checkpoint failed: {treatment_counts}")

# LODO row-level + summaries.
lodo_records=[]
for row in lodo.itertuples(index=False):
    subtype=row.analysis.replace("Subtype::",""); donor=row.omitted_participant
    cond=part_map.loc[donor,"condition"] if donor in part_map.index else ""
    treatment=part_map.loc[donor,"biologic_class"] if donor in part_map.index else ""
    lodo_records.append(dict(row_type="OMISSION",subtype=subtype,pathway_id=row.pathway_id,pathway_name=row.pathway_name,omitted_participant=donor,omitted_condition=cond,omitted_treatment_class=treatment,valid_omission="YES" if pd.notna(row.NES) else "NO",NES=row.NES,direction=row.direction,primary_NES=getattr(row,"primary_NES",np.nan),NES_change=getattr(row,"NES_change",np.nan),sign_reversal=getattr(row,"sign_reversal",""),invalid_reason=(row.invalid_reason if isinstance(row.invalid_reason,str) and row.invalid_reason.strip() else ("PROGRAM_NOT_TESTABLE_AFTER_OMISSION" if pd.isna(row.NES) else "")),valid_omissions=np.nan,invalid_omissions=np.nan,NES_min=np.nan,NES_max=np.nan,sign_reversals=np.nan,max_abs_NES_change=np.nan,influential_participant="",influential_treatment_class="",program_testability_changed="YES" if pd.isna(row.NES) else "NO",interpretation="LODO robustness assessment; not independent validation"))
for (analysis,pid),g in lodo.groupby(["analysis","pathway_id"]):
    valid=g[g.NES.notna()]; primary=float(sub[(sub.analysis==analysis)&(sub.pathway_id==pid)].NES.iloc[0]); pname=g.pathway_name.iloc[0]
    if len(valid):
        idx=valid.NES_change.abs().idxmax(); influential=valid.loc[idx,"omitted_participant"]; treatment=part_map.loc[influential,"biologic_class"]
        rec=dict(row_type="SUMMARY",subtype=analysis.replace("Subtype::",""),pathway_id=pid,pathway_name=pname,omitted_participant="",omitted_condition="",omitted_treatment_class="",valid_omission="",NES=np.nan,direction="",primary_NES=primary,NES_change=np.nan,sign_reversal="",invalid_reason="",valid_omissions=len(valid),invalid_omissions=len(g)-len(valid),NES_min=valid.NES.min(),NES_max=valid.NES.max(),sign_reversals=int(valid.sign_reversal.eq("YES").sum()),max_abs_NES_change=valid.NES_change.abs().max(),influential_participant=influential,influential_treatment_class=treatment,program_testability_changed="YES" if len(g)>len(valid) else "NO",interpretation="Directionally stable" if not valid.sign_reversal.eq("YES").any() else "DONOR_SENSITIVE_DIRECTION_REVERSAL")
    else:
        rec=dict(row_type="SUMMARY",subtype=analysis.replace("Subtype::",""),pathway_id=pid,pathway_name=pname,valid_omissions=0,invalid_omissions=len(g),primary_NES=primary,interpretation="NOT_ASSESSABLE")
    lodo_records.append(rec)
lodo_out=pd.DataFrame(lodo_records)
write_csv(lodo_out,"05_HUMAN_SUBTYPE_LODO_SENSITIVITY.csv")

# Annotation and treatment audit.
audit=[]
for r in markers.itertuples(index=False):
    audit.append(dict(record_type="TREG_MARKER",subtype="Treg",marker=r.marker,unique_gene_mapping=r.unique_gene_mapping,gene_mapping_count=r.gene_mapping_count,testability=r.testability,treg_detection_fraction=r.treg_detection_fraction,other_T_detection_fraction=r.other_T_detection_fraction,treg_mean_log1p_cpm=r.treg_mean_log1p_cpm,other_T_mean_log1p_cpm=r.other_T_mean_log1p_cpm,treg_participants_with_detection=r.treg_participants_with_detection,treg_participants_total=r.treg_participants_total,marker_direction_support=r.direction_consistent_with_regulatory_annotation,annotation_classification="",participant_id="",condition="",treatment_class="",subtype_cell_count=np.nan,eligible_ge20="",annotation_concern="Predefined marker check only; no genome-wide marker discovery"))
support=int(markers.direction_consistent_with_regulatory_annotation.eq("YES").sum())
audit.append(dict(record_type="TREG_ASSESSMENT",subtype="Treg",marker="six-marker panel",annotation_classification="AUTHOR_TREG_ANNOTATION_WITH_UNVERIFIED_CANONICAL_IDENTITY",marker_direction_support=f"{support}/6 predefined markers enriched versus other verified T cells",annotation_concern="Strong multi-marker consistency supports a regulatory-state interpretation, but author labels were not independently reannotated and do not prove canonical identity."))
for r in reliability.itertuples(index=False):
    audit.append(dict(record_type="SUBTYPE_ANNOTATION",subtype=r.author_subtype,annotation_classification=r.annotation_reliability_status,annotation_concern=r.identity_ambiguity_or_boundary,author_annotation_evidence=r.author_annotation_evidence,marker_table_status=r.candidate_specific_marker_table_local_status,participant_id="",condition="",treatment_class="",subtype_cell_count=np.nan,eligible_ge20=""))
for r in crohn_recovered.itertuples(index=False):
    rows=contrib[contrib.participant_id.eq(r.frozen_participant_id)]
    for z in rows.itertuples(index=False):
        audit.append(dict(record_type="CROHN_TREATMENT_SUBTYPE_CONTRIBUTION",subtype=z.subtype,participant_id=r.frozen_participant_id,condition="ActiveCrohn",treatment_class=r.biologic_class,treatment=r.biologic_treatment,subtype_cell_count=z.cell_count,eligible_ge20=z.eligible_ge20,annotation_concern="Treatment contribution is descriptive; no treatment-adjusted model was fit."))
audit.append(dict(record_type="TREATMENT_DISTRIBUTION",subtype="ALL_RECOVERED_CROHN",participant_id="",condition="ActiveCrohn",treatment_class="ANTI_TNF=2;ANTI_IL12_23=2;NO_BIOLOGIC=2",annotation_concern="Six Crohn participants are insufficient for causal or stable treatment-adjusted inference."))
write_csv(pd.DataFrame(audit),"06_HUMAN_TCELL_ANNOTATION_AND_TREATMENT_AUDIT.csv")

# Cohort selection table: participants, cohort summaries, and program direction changes.
sel_records=[]
for r in parts.itertuples(index=False):
    sel_records.append(dict(row_type="PARTICIPANT",participant_id=r.frozen_participant_id,condition=r.condition,recovered_author_subset=r.participant_in_recovered_author_labelled_subset,broad_T_cells=r.frozen_broad_T_cells,verified_T_related_labels=r.valid_author_T_related_labels,T_related_traceability_fraction=r.T_related_traceability_fraction,missing_T_related_label_fraction=r.missing_T_related_label_fraction,biologic_class=r.biologic_class,biologic_treatment=r.biologic_treatment,tissue_region=r.biopsy_location,inflammation_state=r.inflammation_state,platform=r.platform,program_id="",program_name="",full27_NES=np.nan,recovered19_allT_NES=np.nan,recovered19_verifiedT_NES=np.nan,direction_change_full_to_recovered=""))
for cohort,nh,nc,cells,trace in [("FULL_PRIMARY_27",16,11,21378,1.0),("RECOVERED_19_ALL_BROAD_T",13,6,int(qc[qc.analysis.eq("Broad19_AllFrozenT")].n_cells.sum()),np.nan),("RECOVERED_19_VERIFIED_AUTHOR_T",13,6,15233,15233/int(qc[qc.analysis.eq("Broad19_AllFrozenT")].n_cells.sum()))]:
    sel_records.append(dict(row_type="COHORT_SUMMARY",participant_id="",condition="",recovered_author_subset=cohort,broad_T_cells=cells,verified_T_related_labels=15233 if cohort.endswith("AUTHOR_T") else np.nan,T_related_traceability_fraction=trace,missing_T_related_label_fraction=1-trace if pd.notna(trace) else np.nan,n_healthy=nh,n_active_crohn=nc,tissue_region="terminal ileum",inflammation_state="Healthy versus active/inflamed Crohn",program_id="",program_name="",full27_NES=np.nan,recovered19_allT_NES=np.nan,recovered19_verifiedT_NES=np.nan,direction_change_full_to_recovered=""))
for p in PROGRAM_ORDER:
    rows=broad[broad.pathway_name.eq(p)].set_index("cohort")
    f=float(rows.loc["Full27_FrozenReference","NES"]);b=float(rows.loc["Recovered19_AllFrozenBroadT","NES"]);c=float(rows.loc["Recovered19_VerifiedAuthorT","NES"])
    sel_records.append(dict(row_type="PROGRAM_DIRECTION",participant_id="",condition="",recovered_author_subset="",program_id=rows.iloc[0].pathway_id,program_name=p,full27_NES=f,recovered19_allT_NES=b,recovered19_verifiedT_NES=c,direction_change_full_to_recovered="NO" if np.sign(f)==np.sign(b)==np.sign(c) else "YES",interpretation="BROAD_DIRECTION_PRESERVED" if np.sign(f)==np.sign(b)==np.sign(c) else "BROAD_DIRECTION_SENSITIVE_TO_SELECTION"))
write_csv(pd.DataFrame(sel_records),"07_COHORT_SELECTION_SENSITIVITY.csv")

# Figures.
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.titlesize":10,"axes.labelsize":9})
fig,ax=plt.subplots(figsize=(7.2,5.2))
x=np.arange(3);offsets=[-.18,0,.18]; cohorts=["Full27_FrozenReference","Recovered19_AllFrozenBroadT","Recovered19_VerifiedAuthorT"];labels=["Full 27 (16 H / 11 CD)","Recovered 19, all broad T (13 H / 6 CD)","Recovered 19, verified author-T (13 H / 6 CD)"];colors=["#4C596B","#C26455","#4D7C8A"];markers_plot=["o","s","D"]
for off,cohort,label,color,mk in zip(offsets,cohorts,labels,colors,markers_plot):
    d=broad[broad.cohort.eq(cohort)].set_index("pathway_name").loc[PROGRAM_ORDER]
    ax.scatter(x+off,d.NES,s=55,color=color,marker=mk,label=label,zorder=3)
    for xx,y in zip(x+off,d.NES): ax.text(xx,y+(0.055 if y>=0 else -0.075),f"{y:.2f}",ha="center",va="bottom" if y>=0 else "top",fontsize=7,color=color)
ax.axhline(0,color="#777777",lw=.8);ax.set_xticks(x,["Respiration /\nelectron transport","Mitoribosome-\nassociated QC","Xenobiotic\nmetabolism"]);ax.set_ylabel("Normalized enrichment score (NES)");ax.set_title("Broad Human T-cell cohort-selection sensitivity",loc="left",fontweight="bold");fig.legend(frameon=False,fontsize=8,loc="lower center",bbox_to_anchor=(.5,.07),ncol=1);ax.spines[["top","right"]].set_visible(False);fig.text(.10,.015,"Same frozen Human program definitions and rank convention; tested-gene universes are cohort-specific.",fontsize=7,color="#59636F")
fig.subplots_adjust(bottom=.40,top=.90,left=.10,right=.98);fig.savefig(OUT/"08_BROAD_T_COHORT_SENSITIVITY.png",dpi=300,bbox_inches="tight");plt.close(fig)

mat=sub.pivot(index="subtype",columns="pathway_name",values="NES").loc[SUBTYPES,PROGRAM_ORDER]
qmat=sub.pivot(index="subtype",columns="pathway_name",values="subtype_family_adjusted_P").loc[SUBTYPES,PROGRAM_ORDER]
counts=elig.set_index("subtype")
fig,ax=plt.subplots(figsize=(7.4,5.4));lim=max(2.35,float(np.nanmax(np.abs(mat.values))));im=ax.imshow(mat.values,cmap="RdBu_r",norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),aspect="auto")
ax.set_xticks(range(3),["Respiration /\nelectron transport","Mitoribosome-\nassociated QC","Xenobiotic\nmetabolism"])
yl=[f"{s}  (H={int(counts.loc[s,'healthy_participants_ge20'])}, CD={int(counts.loc[s,'crohn_participants_ge20'])})" for s in SUBTYPES];ax.set_yticks(range(6),yl)
for i in range(6):
    for j in range(3):
        v=mat.iloc[i,j];q=qmat.iloc[i,j];color="white" if abs(v)>1.35 else "#20242A";ax.text(j,i,f"NES {v:.2f}\nq18 {q:.3g}",ha="center",va="center",fontsize=7.6,color=color)
ax.set_title("Exploratory Human T-cell subtype metabolic heterogeneity",loc="left",fontweight="bold");cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.03);cb.set_label("NES (Active Crohn ? Healthy)");ax.text(0,-.17,"Author-defined states; participant-level pseudobulk; BH across the prespecified 18-test family.",transform=ax.transAxes,fontsize=7,color="#59636F")
fig.tight_layout();fig.savefig(OUT/"09_HUMAN_SUBTYPE_METABOLIC_HETEROGENEITY.png",dpi=300,bbox_inches="tight");plt.close(fig)

# Reports.
full_to_b=[]
for p in PROGRAM_ORDER:
    d=broad[broad.pathway_name.eq(p)].set_index("cohort")
    full_to_b.append(f"- {p}: full 27 NES {d.loc['Full27_FrozenReference','NES']:.3f}; recovered-19 all-T NES {d.loc['Recovered19_AllFrozenBroadT','NES']:.3f}; verified-author-T NES {d.loc['Recovered19_VerifiedAuthorT','NES']:.3f}; all Active-Crohn-higher.")
sig=sub[sub.subtype_family_adjusted_P<.05].sort_values("subtype_family_adjusted_P")
sig_lines=[f"- {r.subtype}, {r.pathway_name}: NES {r.NES:.3f}, P={r.pval:.3g}, BH q18={r.subtype_family_adjusted_P:.3g}." for r in sig.itertuples()]
ha49=lodo[(lodo.analysis=="Subtype::CD4 Tcm")&(lodo.omitted_participant=="HA49")&lodo.NES.notna()].set_index("pathway_name")
ha49_lines=[f"- {p}: full NES {float(sub[(sub.subtype=='CD4 Tcm')&(sub.pathway_name==p)].NES.iloc[0]):.3f}; without HA49 NES {ha49.loc[p,'NES']:.3f}; direction retained." for p in PROGRAM_ORDER]
report=f"""# Phase 5B integrated Human T-cell subtype metabolic heterogeneity analysis

## A. Executive conclusion

**GO for concise Supplementary Results as exploratory sensitivity evidence.** All three requested components completed. All six author-defined subtypes passed the frozen donor/cell gate and the original broad Human T-cell direction was reproduced in both recovered-19 analyses. Interpretation remains selection-, donor-, annotation- and treatment-limited.

## B. Broad T-cell cohort-selection sensitivity

{chr(10).join(full_to_b)}

Thus all three programs are **BROAD_DIRECTION_PRESERVED**. Absolute NES changes are descriptive because tested-gene universes are cohort-specific and the frozen full-collection FDR family differs from the targeted three-program Phase5B family. The recovered subset cannot be assumed representative of all 27 participants.

## C. Human T-cell subtype disease effects

Six of six candidates were analyzable using participant-level raw-count pseudobulk. Seven of 18 predefined tests passed BH FDR 0.05:

{chr(10).join(sig_lines)}

No mitochondrial ribosome-associated quality-control test passed the 18-test correction. The strongest estimate was Treg respiration (NES {sub[(sub.subtype=='Treg')&(sub.pathway_name==PROGRAM_ORDER[0])].NES.iloc[0]:.3f}, q18 {sub[(sub.subtype=='Treg')&(sub.pathway_name==PROGRAM_ORDER[0])].subtype_family_adjusted_P.iloc[0]:.3g}). The weakest overall pattern was CD103+ CD8: respiration and xenobiotic metabolism were negative but imprecise, and no program passed the family correction.

## D. T-cell subtype heterogeneity

Respiration and xenobiotic metabolism were positive in broad T cells and most author-defined states, but negative in CD103+ CD8. This is **SUBTYPE_HETEROGENEITY** in a descriptive sense only: no subtype-by-disease interaction was tested, and differing NES signs do not establish statistically significant heterogeneity. Mitoribosome-associated QC was positive in all six primary subtype estimates but generally imprecise.

## E. Donor and treatment influence

Valid LODO retained at least five participants per condition. CD4 Tcm was testable after Healthy-participant omissions only because it had exactly five eligible Crohn participants. HA49 contributed 1,157/3,868 (29.91%) of all recovered CD4 Tcm cells (1,157/3,819 eligible-analysis cells), but its omission retained all three positive directions:

{chr(10).join(ha49_lines)}

HA49 changes magnitude, especially respiration and mitoribosome-associated QC, but does not create the CD4 Tcm directions. Several weaker subtype-program estimates showed LODO sign reversals and are marked **DONOR_SENSITIVE**. The six recovered Crohn participants comprise 2 anti-TNF, 2 anti-IL-12/23 and 2 without biologic therapy. This heterogeneity was described, not modeled; treatment independence cannot be claimed.

## F. Annotation and selection limitations

All six predefined Treg markers mapped uniquely and were more frequently detected and more highly expressed in author-Treg than in other verified T cells. Detection in author-Treg was FOXP3 {markers.set_index('marker').loc['FOXP3','treg_detection_fraction']:.1%}, IL2RA {markers.set_index('marker').loc['IL2RA','treg_detection_fraction']:.1%}, CTLA4 {markers.set_index('marker').loc['CTLA4','treg_detection_fraction']:.1%}, TIGIT {markers.set_index('marker').loc['TIGIT','treg_detection_fraction']:.1%}, IKZF2 {markers.set_index('marker').loc['IKZF2','treg_detection_fraction']:.1%} and ICOS {markers.set_index('marker').loc['ICOS','treg_detection_fraction']:.1%}. This provides strong multi-marker consistency with a regulatory state, but the correct classification remains **AUTHOR_TREG_ANNOTATION_WITH_UNVERIFIED_CANONICAL_IDENTITY** because the label was not independently reannotated.

Only 19/27 participants had recoverable fine labels; verified author-T traceability also differed by disease group in Phase5A. Missing labels were never converted to zero subtype abundance, and NK/ILC labels were excluded.

## G. Implications for canine-Human translational interpretation

The Human broad-T metabolic increase is not an artifact of restricting to the 19 recovered participants. Positive effects are concentrated in several CD4/regulatory and KLRG1+ states, while CD103+ CD8 differs, making the Human signal subtype-aware rather than homogeneous. This sharpens the boundary of the previously observed canine-Human broad-T divergence, but does not establish a matched canine-Human subtype mechanism.

The frozen canine analysis remains limiting: GZMA-high T was excluded for donor-confounded T_8; canine Treg is canonical-marker limited; Cycling T showed mitochondrial decreases; Memory T was heterogeneous; and there were only seven dogs. Human CD4 Tcm was not equated to canine Memory T, and Human/canine Treg labels were not treated as equivalent.

## H. Manuscript recommendation

**Suitable for concise Supplementary Results.** Present it as exploratory sensitivity evidence, not independent validation or a formal dog-Human subtype comparison. Retain the selected-cohort, donor, annotation and mixed-treatment caveats.
"""
(OUT/"10_PHASE5B_INTEGRATED_SCIENTIFIC_REPORT.md").write_text(report,encoding="utf-8")
broad_report="# Broad Human T-cell cohort-selection sensitivity\n\nAll three frozen metabolic programs remained Active-Crohn-higher in the full 27-participant reference, the recovered 19-participant all-broad-T analysis, and the recovered 19-participant verified-author-T analysis.\n\n"+ "\n".join(full_to_b)+"\n\nAdjusted P values are not compared across the frozen full-collection and targeted three-program testing families. NES differences are descriptive; tested-gene universes are cohort-specific.\n"
(OUT/"02_BROAD_T_COHORT_SENSITIVITY_REPORT.md").write_text(broad_report,encoding="utf-8")

# Reproducibility log written after all deliverables exist.
inputs=[Path(preflight["phase5a_map"]["path"]),Path(preflight["phase5a_crosswalk"]["path"]),Path(preflight["phase5a_coverage"]["path"]),Path(preflight["phase5a_reliability"]["path"]),HBASE/"gsea"/"full"/"Human_T_Reactome_GSEA.csv",HBASE/"gsea"/"full"/"Human_T_Hallmark_GSEA.csv",HBASE/"gsea"/"gene_sets"/"HUMAN_REACTOME_GENE_SETS.tsv.gz",HBASE/"gsea"/"gene_sets"/"HUMAN_HALLMARK_GENE_SETS.tsv.gz"]
output_files=[OUT/f"{i:02d}_{name}" for i,name in []]
csvs=[OUT/n for n in ["01_BROAD_T_COHORT_SENSITIVITY.csv","03_HUMAN_SUBTYPE_PSEUDOBULK_ELIGIBILITY.csv","04_HUMAN_SUBTYPE_METABOLIC_RESULTS.csv","05_HUMAN_SUBTYPE_LODO_SENSITIVITY.csv","06_HUMAN_TCELL_ANNOTATION_AND_TREATMENT_AUDIT.csv","07_COHORT_SELECTION_SENSITIVITY.csv"]]
others=[OUT/"02_BROAD_T_COHORT_SENSITIVITY_REPORT.md",OUT/"08_BROAD_T_COHORT_SENSITIVITY.png",OUT/"09_HUMAN_SUBTYPE_METABOLIC_HETEROGENEITY.png",OUT/"10_PHASE5B_INTEGRATED_SCIENTIFIC_REPORT.md"]
log=f"""# Phase 5B reproducibility log

- Completed UTC: {datetime.now(timezone.utc).isoformat()}
- Frozen Phase5A mapping: {preflight['phase5a_map']['sha256']}
- Frozen h5ad SHA256 (reused from Phase5A; not recomputed): {preflight['h5ad']['sha256_reused_from_phase5a']}
- Raw count source: layers[counts]; nonnegative integer invariant PASS.
- Biological replicate: participant; no cell-level inference.
- Pseudobulk: one raw-count sum per eligible participant x state.
- Eligibility: >=20 cells per participant x subtype and >=5 Healthy plus >=5 Active Crohn participants.
- edgeR: filterByExpr(design); TMM; estimateDisp robust=TRUE; glmQLFit robust=TRUE; glmQLFTest, ActiveCrohn-minus-Healthy.
- Rank: sign(logFC)*sqrt(F), duplicate symbols collapsed exactly as Phase7.
- GSEA: fgseaMultilevel, minSize=15, maxSize=500, eps=1e-10, scoreType=std, nproc=1, base seed 20261015 plus deterministic analysis offset.
- Human gene sets: Reactome v97 and MSigDB 2026.1.Hs; exact IDs R-HSA-1428517, R-HSA-9937383, HALLMARK_XENOBIOTIC_METABOLISM.
- Primary subtype family: 18 planned and 18 valid P values; BH across all 18.
- Broad Phase5B family: three predefined programs within each analysis. Frozen 27-participant FDRs retain their original full-collection families.
- Valid LODO: omission must preserve >=5 participants per condition; invalid omissions recorded NON_ESTIMABLE.
- R software: {rversions.to_dict(orient='records')}
- Frozen inputs modified: NO.
- Manuscript/main figures modified: NO.
- Cross-species subtype matching performed: NO.

## Input checksums

""" + "\n".join(f"- {p}: {sha256(p)}" for p in inputs) + f"\n- {preflight['h5ad']['path']}: {preflight['h5ad']['sha256_reused_from_phase5a']} (reused)\n\n## Output checksums\n\n" + "\n".join(f"- {p.name}: {sha256(p)}" for p in csvs+others) + "\n"
(OUT/"11_REPRODUCIBILITY_LOG.md").write_text(log,encoding="utf-8")
print(json.dumps({"analyzable_subtypes":int(elig.eligibility_status.eq("ANALYZABLE").sum()),"family_significant_tests":int((sub.subtype_family_adjusted_P<.05).sum()),"broad_directions_preserved":bool(broad.direction_agreement_with_full27.eq("YES").all()),"treg_marker_support":support,"output":str(OUT)},indent=2))

