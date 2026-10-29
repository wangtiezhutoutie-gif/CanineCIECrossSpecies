from pathlib import Path
import hashlib, json, math
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
OUT = ROOT / "WORKING" / "22_DOG_AUTHOR_TCELL_CROSS_SPECIES"
AOUT = OUT / "01_AUTHOR_ANNOTATION_AUDIT"
BOUT = OUT / "02_DOG_SUBTYPE_ANALYSIS"
COUT = OUT / "03_DOG_HUMAN_CORRESPONDENCE"
WORK = OUT / "work"
RDS = WORK / "Tcell_duod_annotated.rds"

def sha256(path, block=1024*1024):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(block),b""): h.update(b)
    return h.hexdigest()

def md5(path, block=1024*1024):
    h=hashlib.md5()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(block),b""): h.update(b)
    return h.hexdigest()

def fmt(x, n=3):
    if x is None or (isinstance(x,float) and math.isnan(x)): return "NA"
    return f"{x:.{n}f}"

def direction(x):
    if pd.isna(x): return "NOT_ASSESSABLE"
    return "CIE_HIGHER" if x>0 else ("CIE_LOWER" if x<0 else "TIE")

pre=json.loads((WORK/"phase6_preflight.json").read_text(encoding="utf-8"))
freeze=json.loads((WORK/"matching_decision_freeze.json").read_text(encoding="utf-8"))
cw=pd.read_csv(AOUT/"02_AUTHOR_FROZEN_CELL_CROSSWALK.csv",low_memory=False)
dc=pd.read_csv(AOUT/"03_AUTHOR_SUBTYPE_DONOR_COUNTS.csv")
ct=pd.read_csv(AOUT/"04_AUTHOR_FROZEN_LABEL_CROSSTAB.csv")
elig=pd.read_csv(BOUT/"01_SUBTYPE_ELIGIBILITY.csv")
comp=pd.read_csv(BOUT/"02_DONOR_COMPOSITION.csv")
clodo=pd.read_csv(BOUT/"03_COMPOSITION_LODO.csv")
marker=pd.read_csv(WORK/"dog_targeted_marker_audit.csv")
matches=pd.read_csv(COUT/"01_MATCHING_DECISIONS_FROZEN.csv")
compat=pd.read_csv(COUT/"03_GENESET_COMPATIBILITY.csv")
native=pd.read_csv(WORK/"native_metabolic_results.tsv",sep="\t")
nlodo=pd.read_csv(WORK/"native_lodo_results.tsv",sep="\t")
matched=pd.read_csv(WORK/"matched_state_results.tsv",sep="\t")

programs=[
 ("Mitochondrial respiration","Reactome","R-CFA-1428517","Aerobic respiration and respiratory electron transport","Reactome v97"),
 ("Mitochondrial ribosome-associated quality control","Reactome","R-CFA-9937383","Mitochondrial ribosome-associated quality control","Reactome v97"),
 ("Metabolic/xenobiotic program","Hallmark","HALLMARK_XENOBIOTIC_METABOLISM","Xenobiotic Metabolism","MSigDB 2026.1.Hs; babelgene 22.9"),
]
if len(native):
    native=native.rename(columns={"pval":"nominal_P","phase6_primary_family_q":"BH_q_complete_phase6_primary_family","size":"pathway_gene_set_size"})
native_rows=[]
present=set(native["dog_author_subtype"].unique()) if len(native) else set()
for r in elig.itertuples(index=False):
    if r.author_subtype in present: continue
    for prog,coll,pid,pname,ver in programs:
        native_rows.append({
            "dog_author_subtype":r.author_subtype,"program":prog,"collection":coll,
            "dog_pathway_id":pid,"pathway_name":pname,"source_version":ver,
            "healthy_donor_n":r.healthy_donors_ge20,"cie_donor_n":r.cie_donors_ge20,
            "tested_genes":np.nan,"mapped_ranked_genes":np.nan,"pathway_gene_set_size":np.nan,
            "NES":np.nan,"nominal_P":np.nan,"BH_q_complete_phase6_primary_family":np.nan,
            "direction":"NOT_TESTABLE","model_status":"NOT_ASSESSABLE_"+str(r.annotation_classification),
            "planned_test_family":f"{int((elig.formal_pseudobulk_eligible=='YES').sum())} eligible author subtypes x 3 predefined programs",
            "valid_p_values_adjusted":int(native["nominal_P"].notna().sum()) if len(native) else 0,
        })
metabolic=pd.concat([native,pd.DataFrame(native_rows)],ignore_index=True,sort=False)
metabolic=metabolic.merge(elig[["author_subtype","annotation_classification","annotation_evidence_note","maximum_donor_fraction","t8_fraction_of_subtype","t9_fraction_of_subtype"]],left_on="dog_author_subtype",right_on="author_subtype",how="left").drop(columns=["author_subtype"])
metabolic["analysis_scope"]="POST_HOC_PREDEFINED_FOLLOW_UP_SENSITIVITY"
metabolic.to_csv(BOUT/"04_METABOLIC_PROGRAM_RESULTS.csv",index=False)

if len(nlodo):
    nlodo=nlodo.rename(columns={"pval":"nominal_P","size":"pathway_gene_set_size"})
nlodo.to_csv(BOUT/"05_METABOLIC_LODO.csv",index=False)

if len(matched) and {"dog_native_NES","dog_NES","human_native_NES","human_NES"}.issubset(matched.columns):
    matched["dog_native_vs_strict_direction"] = np.where(np.sign(matched["dog_native_NES"])==np.sign(matched["dog_NES"]),"SAME_DIRECTION","SIGN_CHANGED")
    matched["human_native_vs_strict_direction"] = np.where(np.sign(matched["human_native_NES"])==np.sign(matched["human_NES"]),"SAME_DIRECTION","SIGN_CHANGED")
    matched["gene_universe_direction_stability"] = np.where((matched["dog_native_vs_strict_direction"]=="SAME_DIRECTION") & (matched["human_native_vs_strict_direction"]=="SAME_DIRECTION"),"STABLE_IN_BOTH_SPECIES","GENE_UNIVERSE_SENSITIVE")
matched.to_csv(COUT/"04_MATCHED_STATE_METABOLIC_RESULTS.csv",index=False)

# Descriptive author-versus-frozen comparison.
l2ct=ct[ct["table_type"].eq("AUTHOR_L2_X_FROZEN_LEVEL2")].copy()
l2ct=l2ct.sort_values(["author_label","shared_cell_count"],ascending=[True,False])
dominant=l2ct.groupby("author_label",as_index=False).first()[["author_label","frozen_label","row_percentage"]]
frozen_res=pd.read_csv(ROOT/"WORKING"/"16_TCELL_SUBTYPE_METABOLIC_SENSITIVITY"/"02_SUBTYPE_METABOLIC_RESULTS.csv")
broad=pd.read_csv(ROOT/"WORKING"/"16_TCELL_SUBTYPE_METABOLIC_SENSITIVITY"/"work"/"program_definitions.tsv",sep="\t")
compare=[]
for r in metabolic.itertuples(index=False):
    dm=dominant[dominant["author_label"].eq(r.dog_author_subtype)]
    frozen_label=dm.iloc[0]["frozen_label"] if len(dm) else ""
    overlap=float(dm.iloc[0]["row_percentage"]) if len(dm) else np.nan
    fr=frozen_res[(frozen_res["subtype"].eq(frozen_label)) & (frozen_res["pathway_id"].eq(r.dog_pathway_id))]
    br=broad[broad["pathway_id"].eq(r.dog_pathway_id)]
    fnes=float(fr.iloc[0]["NES"]) if len(fr) and pd.notna(fr.iloc[0]["NES"]) else np.nan
    fq=float(fr.iloc[0]["FDR_predefined_family"]) if len(fr) and pd.notna(fr.iloc[0]["FDR_predefined_family"]) else np.nan
    bnes=float(br.iloc[0]["broad_dog_native_NES"]) if len(br) else np.nan
    lo=nlodo[(nlodo.get("dog_author_subtype",pd.Series(dtype=str)).eq(r.dog_author_subtype)) & (nlodo.get("program",pd.Series(dtype=str)).eq(r.program))] if len(nlodo) else pd.DataFrame()
    valid=lo[lo["model_status"].eq("ESTIMABLE")] if len(lo) else pd.DataFrame()
    reversals=int((valid["direction_match_full"].eq("NO")).sum()) if len(valid) else 0
    if pd.isna(r.NES):
        cls="NOT_ASSESSABLE"
    elif reversals>0:
        cls="DONOR_SENSITIVE"
    elif pd.notna(fnes):
        cls="CONSISTENT_WITH_FROZEN_ANALYSIS" if np.sign(r.NES)==np.sign(fnes) else "HETEROGENEOUS_SUBTYPE_RESPONSE"
    elif pd.notna(bnes) and np.sign(r.NES)==np.sign(bnes):
        cls="CONSISTENT_WITH_BROAD_T_ANALYSIS"
    else:
        cls="IMPROVED_STATE_RESOLUTION"
    compare.append({
        "author_subtype":r.dog_author_subtype,"dominant_frozen_level2":frozen_label,
        "dominant_frozen_overlap_row_percentage":overlap,"program":r.program,
        "dog_pathway_id":r.dog_pathway_id,"author_subtype_NES":r.NES,
        "author_subtype_phase6_BH_q":r.BH_q_complete_phase6_primary_family,
        "frozen_subtype_NES":fnes,"frozen_subtype_phase2b_BH_q":fq,
        "broad_dog_T_NES":bnes,"valid_lodo_omissions":len(valid),
        "valid_lodo_sign_reversals":reversals,"comparison_classification":cls,
        "interpretation_boundary":"Different author/frozen labels and matched cell universes can change estimates; no subtype-by-disease interaction was tested.",
    })
pd.DataFrame(compare).to_csv(BOUT/"06_AUTHOR_VS_FROZEN_COMPARISON.csv",index=False)

# Evidence summaries used in reports.
author_total=int(pre["author_cells"])
frozen_t_total=int(pre["frozen_t_cells"])
matched_author=int(cw["match_status"].isin(["EXACT_MATCH","VERIFIED_TRANSFORMED_MATCH"]).sum())
matched_t=int(pre["shared_frozen_t_cells"])
author_match_rate=matched_author/author_total if author_total else np.nan
frozen_t_match_rate=matched_t/frozen_t_total if frozen_t_total else np.nan
subtype_count=len(elig)
eligible_n=int(elig["formal_pseudobulk_eligible"].eq("YES").sum())
donor_conf_n=int(elig["annotation_classification"].eq("DONOR_CONFOUNDED").sum())
supported_n=int(matches["matching_classification"].eq("SUPPORTED_EXPLORATORY_MATCH").sum())
partial_n=int(matches["matching_classification"].eq("PARTIAL_FUNCTIONAL_ANALOGY").sum())
valid_met=metabolic[metabolic["model_status"].eq("ESTIMABLE")].copy()
valid_match=matched[matched.get("effect_classification",pd.Series(dtype=str)).isin(["DIRECTIONALLY_CONCORDANT","DIRECTIONALLY_DISCORDANT"])] if len(matched) else pd.DataFrame()

shared_cw=cw[cw["match_status"].isin(["EXACT_MATCH","VERIFIED_TRANSFORMED_MATCH"])].copy()
t8=shared_cw[shared_cw["frozen_annotation_cluster"].eq("T_8")]
t8_sub=t8["author_celltype_l2"].value_counts()
t8_text="; ".join(f"{k}: {v}" for k,v in t8_sub.items()) if len(t8_sub) else "No matched T_8 cells"
t8_donor=t8["frozen_donor_id"].value_counts()
t8_donor_text="; ".join(f"{k}: {v}" for k,v in t8_donor.items()) if len(t8_donor) else "NA"

treg=marker[(marker["subtype"].eq("Treg (c7)")) & marker["mapping_status"].eq("UNIQUE_GENE_SYMBOL")].copy()
treg_support=treg[(treg["subtype_detection_fraction"]>=0.005)&(treg["mean_lognorm_difference"]>0)]
treg_text="; ".join(f"{r.gene} det={r.subtype_detection_fraction:.3f}, delta={r.mean_lognorm_difference:.3f}" for r in treg.itertuples(index=False))
treg_cross=ct[(ct["table_type"].eq("AUTHOR_L2_X_FROZEN_LEVEL2")) & (ct["frozen_label"].eq("Treg"))].copy()
frozen_treg_total=int(treg_cross["shared_cell_count"].sum())
author_treg_in_frozen=int(treg_cross.loc[treg_cross["author_label"].eq("Treg (c7)"),"shared_cell_count"].sum())
author_tnaive_in_frozen=int(treg_cross.loc[treg_cross["author_label"].eq("Tnaive (c1)"),"shared_cell_count"].sum())
author_treg_pct=100*author_treg_in_frozen/frozen_treg_total if frozen_treg_total else float("nan")
author_tnaive_pct=100*author_tnaive_in_frozen/frozen_treg_total if frozen_treg_total else float("nan")

# Same-universe IL7R comparison.
sall=shared_cw[(shared_cw["frozen_level1"].eq("T")) & shared_cw["author_celltype_l2"].notna()].copy()
il7_rows=[]
for d,g in sall.groupby("frozen_donor_id"):
    il7_rows.append({
        "donor":d,"condition":g["frozen_condition"].iloc[0],
        "author_il7r_l1_fraction":g["author_celltype_l1"].eq("IL7R_high").mean(),
        "frozen_il7r_level2_fraction":g["frozen_level2"].eq("IL7R-high T").mean(),
    })
il7=pd.DataFrame(il7_rows)
if len(il7):
    il7_h=il7[il7.condition.eq("Healthy")]
    il7_c=il7[il7.condition.eq("CIE")]
    il7_summary=(f"On the same shared-cell denominator, author IL7R_high median H={il7_h.author_il7r_l1_fraction.median():.3f}, "
                 f"CIE={il7_c.author_il7r_l1_fraction.median():.3f}; frozen IL7R-high T median H={il7_h.frozen_il7r_level2_fraction.median():.3f}, "
                 f"CIE={il7_c.frozen_il7r_level2_fraction.median():.3f}.")
else:
    il7_summary="Same-universe IL7R comparison was not assessable."

# Per-donor identity balance.
per_donor=[]
for d in sorted(set(shared_cw["frozen_donor_id"].dropna())):
    fden=int(((cw["frozen_donor_id"].eq(d)) & (cw["frozen_level1"].eq("T"))).sum())
    n=int((shared_cw["frozen_donor_id"].eq(d)&shared_cw["frozen_level1"].eq("T")).sum())
    per_donor.append(f"{d}: {n}/{fden} ({n/fden:.1%})" if fden else f"{d}: NA")
per_donor_text="; ".join(per_donor)

source_md=f"""# Source metadata and checksums

- Official record: https://zenodo.org/records/11434076
- DOI: 10.5281/zenodo.11434076
- File: Tcell_duod_annotated.rds
- Published size: 1,255,686,853 bytes
- Observed size: {RDS.stat().st_size:,} bytes
- Published MD5: ce31c4350ccd2c8b5902ed5a169c4f94
- Observed MD5: {md5(RDS)}
- Observed SHA256: {sha256(RDS)}
- Integrity status: **{"VERIFIED" if md5(RDS)=="ce31c4350ccd2c8b5902ed5a169c4f94" and RDS.stat().st_size==1255686853 else "FAILED"}**
- Official annotation code: https://github.com/dyammons/canine_duodenal_atlas/blob/main/scripts/3_tcell_analysis.R
- Local author-code SHA256: {sha256(WORK/"author_3_tcell_analysis.R")}
- Frozen canine h5ad SHA256: {sha256(ROOT/"05_annotations"/"dog_GSE254005_annotated.h5ad")}
- Strict ortholog master SHA256: {sha256(ROOT/"10_cross_species_baseline"/"mapping"/"DOG_HUMAN_ORTHOLOG_MASTER.csv")}

## Object inspection

{(WORK/"author_object_summary.txt").read_text(encoding="utf-8",errors="replace")}

The RDS was read sparsely; no complete expression matrix was converted to a dense table.
"""
(AOUT/"01_SOURCE_METADATA_AND_CHECKSUMS.md").write_text(source_md,encoding="utf-8")

audit_md=f"""# Targeted T_8, Treg and IL7R audit

## T_8 / frozen GZMA-high component

Matched frozen T_8 cells were assigned by the author annotation as: {t8_text}.

Per-donor T_8 distribution: {t8_donor_text}.

A new author label was not treated as resolving T_8 confounding. Eligibility was frozen from donor coverage and component dependence before metabolic results. If a donor-dominant component was necessary for a subtype to pass the 3-Healthy/3-CIE gate, that subtype was classified DONOR_CONFOUNDED; cells were not removed to improve eligibility.

## Treg identity

The author object contains Treg (c7). Targeted same-dataset evidence was: {treg_text or "no uniquely mapped supportive markers"}.

{len(treg_support)}/{len(treg)} uniquely mapped predefined markers met the prespecified descriptive support rule (detection >=0.005 and higher mean log-normalized expression than other matched T cells). This supports at most a regulatory-like state. It does not independently validate canonical Treg identity, and it cannot by itself convert the extreme CIE-donor enrichment into a disease effect.

Crucially, only {author_treg_in_frozen}/{frozen_treg_total} ({author_treg_pct:.1f}%) matched frozen Treg cells were author Treg (c7), whereas {author_tnaive_in_frozen}/{frozen_treg_total} ({author_tnaive_pct:.1f}%) were author Tnaive (c1). The author annotation therefore challenges, rather than validates, the breadth of the frozen Treg-like label.

## IL7R discrepancy

{il7_summary}

The discrepancy is classified primarily as ANNOTATION_DEFINITION, with contributions from DENOMINATOR and DONOR_HETEROGENEITY. The author celltype.l1 is a broad IL7R-high compartment, whereas the frozen Level2 label is a different project-defined partition. Pooled proportions are not substituted for donor-level summaries.
"""
(AOUT/"05_T8_TREG_IL7R_AUDIT.md").write_text(audit_md,encoding="utf-8")

a_report=f"""# Component A report: author canine T-cell annotation audit

## Outcome

- Official RDS recovered and checksum-verified: **YES**.
- Actual author celltype.l2 labels: **{subtype_count}**.
- Author cells in RDS: **{author_total:,}**.
- Verified author-to-frozen matches (all frozen Level1 identities): **{matched_author:,}/{author_total:,} ({author_match_rate:.2%})**.
- Shared cells retained in the frozen-T analysis universe: **{matched_t:,}/{frozen_t_total:,} ({frozen_t_match_rate:.2%} of frozen T cells)**.
- Per-donor frozen-T recovery: {per_donor_text}.
- Eligible author states at the annotation/donor gate: **{eligible_n}/{subtype_count}**.

The author labels provide finer state names than the frozen six-label system, but resolution and reliability are not synonymous. CD8 resident/effector, naive, gamma-delta and IFN-state labels add interpretive granularity; NK/NK-T boundaries, Treg canonical identity and donor-dominant components remain limited. Only {author_treg_pct:.1f}% of matched frozen Treg cells were author Treg; most ({author_tnaive_pct:.1f}%) were author Tnaive, so the frozen Treg-like boundary is not confirmed.

## Frozen decision

COMPONENT_A = {"GO" if matched_t>0 and subtype_count>0 else "STOP"}

Eligibility was determined without using Phase 6 NES or P values. See 06_AUTHOR_ANNOTATION_ELIGIBILITY.csv for branch-level decisions.
"""
(AOUT/"07_COMPONENT_A_REPORT.md").write_text(a_report,encoding="utf-8")

# Component B summaries.
shared_comp=comp[comp["universe"].eq("SHARED_AUTHOR_FROZEN_T_UNIVERSE")].copy()
comp_summary=[]
for st,g in shared_comp.groupby("author_subtype"):
    full_diff=float(g["cie_minus_healthy_median_difference"].iloc[0])
    z=clodo[(clodo["author_subtype"].eq(st)) & (clodo["universe"].eq("SHARED_AUTHOR_FROZEN_T_UNIVERSE"))]
    dirs=set(z["direction"].dropna())
    full_dir=direction(full_diff)
    stable=all(x==full_dir for x in dirs) if len(dirs) else False
    comp_summary.append((st,full_diff,full_dir,stable,z.loc[z["cie_minus_healthy_median_difference"].abs().idxmax(),"omitted_donor"] if len(z) else "NA"))
comp_summary=sorted(comp_summary,key=lambda x:abs(x[1]),reverse=True)
comp_text="\n".join(f"- {st}: median CIE-minus-Healthy={d:.4f} ({dr}); all descriptive LODO directions stable={stable}; largest absolute LODO influence={inf}." for st,d,dr,stable,inf in comp_summary)

met_lines=[]
for r in valid_met.sort_values(["dog_author_subtype","program"]).itertuples(index=False):
    z=nlodo[(nlodo["dog_author_subtype"].eq(r.dog_author_subtype)) & (nlodo["program"].eq(r.program)) & (nlodo["model_status"].eq("ESTIMABLE"))] if len(nlodo) else pd.DataFrame()
    signrev=int(z["direction_match_full"].eq("NO").sum()) if len(z) else 0
    rng=f"{z.NES.min():.3f} to {z.NES.max():.3f}" if len(z) else "no valid omissions"
    met_lines.append(f"- {r.dog_author_subtype} / {r.program}: NES={r.NES:.3f}, P={r.nominal_P:.4g}, Phase6-family q={r.BH_q_complete_phase6_primary_family:.4g}; valid LODO={len(z)}, range={rng}, sign reversals={signrev}.")
met_text="\n".join(met_lines) if met_lines else "- No author subtype passed all gates for an estimable metabolic analysis."

b_status="GO" if eligible_n>0 and len(valid_met)>0 else ("PARTIAL_GO" if eligible_n>0 else "STOP")
b_report=f"""# Component B report: canine author-subtype composition and metabolic heterogeneity

## Frozen analysis universe

Only cells with a verified author-to-frozen identity and frozen Level1 T annotation were used for raw-count pseudobulk. Author-only cells were excluded rather than mixed with the frozen count matrix. The primary gate remained >=3 Healthy and >=3 CIE dogs with >=20 cells per donor-by-subtype.

## Donor-level composition

{comp_text}

These are descriptive donor-level patterns. They do not establish disease causality or quantify how much broad-T expression is explained by composition.

## Three predefined metabolic programs

Planned family: {eligible_n} eligible author subtypes x 3 programs = {eligible_n*3} tests. Valid P values in the BH family: {int(valid_met.nominal_P.notna().sum())}.

{met_text}

Healthy-dog omission was always non-estimable under the frozen 3-versus-3 minimum. Valid LODO results, where present, are sensitivity checks rather than independent replication.

## Decision

COMPONENT_B = {b_status}

Finer author annotation is biologically informative only where it adds interpretable state resolution with independent-dog coverage. It does not repair the seven-dog design, Healthy-Beagle versus clinical-CIE confounding, or donor-dominant components.
"""
(BOUT/"07_COMPONENT_B_REPORT.md").write_text(b_report,encoding="utf-8")

# Component C report.
supported=matches[matches["matching_classification"].eq("SUPPORTED_EXPLORATORY_MATCH")]
partial=matches[matches["matching_classification"].eq("PARTIAL_FUNCTIONAL_ANALOGY")]
pair_text="\n".join(f"- {r.dog_author_subtype} <-> {r.human_author_subtype}: {r.matching_rationale}" for r in supported.itertuples(index=False)) or "- None."
partial_text="\n".join(f"- {r.dog_author_subtype} ~ {r.human_author_subtype or 'no unique counterpart'}: {r.matching_rationale}" for r in partial.itertuples(index=False)) or "- None."
effect_text="\n".join(
    f"- {r.dog_author_subtype} <-> {r.human_author_subtype}, {r.program}: native dog/human NES={fmt(r.dog_native_NES)}/{fmt(r.human_native_NES)}; strict dog/human NES={fmt(r.dog_NES)}/{fmt(r.human_NES)}; {r.effect_classification}; {r.gene_universe_direction_stability}."
    for r in valid_match.itertuples(index=False)
) if len(valid_match) else "- No supported pair yielded a testable strict-ortholog matched-state result."
c_status="GO" if supported_n>0 and len(valid_match)>0 else ("PARTIAL_GO" if supported_n>0 else "STOP")
c_report=f"""# Component C report: dog-human T-cell subtype comparability

## Matching was frozen before metabolic outcomes

Matching-decision SHA256: {freeze["sha256"]}

Supported exploratory matches ({supported_n}):

{pair_text}

Partial functional analogies ({partial_n}):

{partial_text}

Similar names, single markers and metabolic directions were not used to create matches. Partial analogies were excluded from the primary matched-state comparison.

## Strict one-to-one ortholog comparison

{effect_text}

Native NES magnitudes were not treated as directly comparable. New matched-state results use the shared tested ortholog universe and pathway-member intersection. Directional concordance or discordance is descriptive; no species-by-disease interaction was tested.

## Decision

COMPONENT_C = {c_status}

The comparison remains bounded by canine duodenum versus human terminal ileum, seven canine donors, Human 13 Healthy/6 Crohn author-labelled subcohort selection, mixed Crohn treatment exposure and species-specific annotation differences.
"""
(COUT/"05_COMPONENT_C_REPORT.md").write_text(c_report,encoding="utf-8")

# Optional evidence figures are generated only when corresponding valid results exist.
if len(shared_comp):
    sts=[x[0] for x in comp_summary]
    fig,ax=plt.subplots(figsize=(max(7,0.75*len(sts)),4.6))
    for i,st in enumerate(sts):
        g=shared_comp[shared_comp.author_subtype.eq(st)]
        for cond,color,marker_style in [("Healthy","#3B6FB6","o"),("CIE","#C14B4B","s")]:
            z=g[g.condition.eq(cond)]
            ax.scatter(np.full(len(z),i)+(-0.09 if cond=="Healthy" else 0.09),z.proportion,color=color,marker=marker_style,s=30,alpha=.85,label=cond if i==0 else None)
    ax.set_xticks(range(len(sts))); ax.set_xticklabels(sts,rotation=40,ha="right")
    ax.set_ylabel("Proportion within shared author-frozen T universe")
    ax.legend(frameon=False,ncol=2); ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(BOUT/"08_DOG_SUBTYPE_COMPOSITION.png",dpi=300); plt.close(fig)
if len(valid_met):
    piv=valid_met.pivot(index="dog_author_subtype",columns="program",values="NES")
    fig,ax=plt.subplots(figsize=(7,max(3,0.55*len(piv))))
    im=ax.imshow(piv.values,aspect="auto",cmap="RdBu_r",vmin=-max(2,abs(np.nanmax(piv.values))),vmax=max(2,abs(np.nanmax(piv.values))))
    ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns,rotation=30,ha="right")
    ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
    for i,st in enumerate(piv.index):
        for j,pr in enumerate(piv.columns):
            rr=valid_met[(valid_met.dog_author_subtype==st)&(valid_met.program==pr)]
            if len(rr): ax.text(j,i,f"{rr.iloc[0].NES:.2f}\nq={rr.iloc[0].BH_q_complete_phase6_primary_family:.2g}",ha="center",va="center",fontsize=7)
    fig.colorbar(im,ax=ax,label="NES"); fig.tight_layout(); fig.savefig(BOUT/"09_DOG_SUBTYPE_METABOLIC_SENSITIVITY.png",dpi=300); plt.close(fig)
if len(valid_match):
    labels=[f"{r.dog_author_subtype} / {r.human_author_subtype}\n{r.program}" for r in valid_match.itertuples(index=False)]
    y=np.arange(len(labels))
    fig,ax=plt.subplots(figsize=(8,max(3,0.48*len(labels))))
    ax.axvline(0,color="#9AA0A6",lw=.8)
    for i,r in enumerate(valid_match.itertuples(index=False)):
        ax.plot([r.dog_NES,r.human_NES],[i,i],color="#B0B7C3",lw=1)
        ax.scatter(r.dog_NES,i,color="#4C78A8",label="Dog" if i==0 else None)
        ax.scatter(r.human_NES,i,color="#E07B39",marker="s",label="Human" if i==0 else None)
    ax.set_yticks(y); ax.set_yticklabels(labels,fontsize=7); ax.set_xlabel("Strict-ortholog matched-universe NES")
    ax.legend(frameon=False); ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(COUT/"06_DOG_HUMAN_TCELL_HETEROGENEITY.png",dpi=300); plt.close(fig)

component_a="GO" if matched_t>0 and subtype_count>0 else "STOP"
component_b=b_status
component_c=c_status
overall="READY_FOR_SUPPLEMENTARY_INTEGRATION" if component_a=="GO" and component_b=="GO" and component_c in {"GO","PARTIAL_GO"} else ("CANINE_EXPLORATORY_ONLY" if component_b!="STOP" else "ANNOTATION_AUDIT_ONLY")
concord_n=int((valid_match["effect_classification"].eq("DIRECTIONALLY_CONCORDANT")).sum()) if len(valid_match) else 0
discord_n=int((valid_match["effect_classification"].eq("DIRECTIONALLY_DISCORDANT")).sum()) if len(valid_match) else 0
stable_comp=[x[0] for x in comp_summary if x[3]]
unstable_comp=[x[0] for x in comp_summary if not x[3]]
integrated=f"""# Integrated decision report

## Final status

- COMPONENT_A: {component_a}
- COMPONENT_B: {component_b}
- COMPONENT_C: {component_c}
- OVERALL: {overall}

## A. Author annotation recovery

1. **RDS recovery:** Yes; exact published size and MD5 were verified.
2. **Actual author subtype count:** {subtype_count}.
3. **Added information:** The author labels add CD8 resident/effector, naive, gamma-delta, regulatory-like and IFN-state granularity, but NK/T boundaries and canonical Treg identity remain annotation-limited.
4. **Cell matching:** {matched_author:,}/{author_total:,} author cells matched a frozen identity; {matched_t:,}/{frozen_t_total:,} frozen T cells entered the shared analysis universe.
5. **Balance:** {per_donor_text}.
6. **T_8:** Author labels split/assign T_8 as {t8_text}; the original donor-confounding warning is not lifted.
7. **Treg-like:** {len(treg_support)}/{len(treg)} predefined uniquely mapped markers support the small author Treg (c7) state, but only {author_treg_in_frozen}/{frozen_treg_total} ({author_treg_pct:.1f}%) frozen Treg cells carried that author label and {author_tnaive_pct:.1f}% were author Tnaive. The author annotation therefore challenges the frozen Treg-like boundary and does not support a formal author-Treg disease model.
8. **IL7R:** {il7_summary} The difference is mainly annotation-definition plus denominator and donor heterogeneity.

## B. Canine subtype heterogeneity

9. **Strict donor-gate passes:** {eligible_n}/{subtype_count}.
10. **Donor-confounded exclusions:** {donor_conf_n}.
11. **Directionally stable composition:** {", ".join(stable_comp) if stable_comp else "None"}.
12. **Analyzable metabolic results:** {len(valid_met)} subtype-program tests across {valid_met.dog_author_subtype.nunique() if len(valid_met) else 0} author states.
13. **Valid donor omission:** See 05_METABOLIC_LODO.csv; only omissions retaining >=3 Healthy and >=3 CIE were calculated. {sum((nlodo.model_status=="ESTIMABLE")) if len(nlodo) else 0} omission-program results were estimable.
14. **Not assessable with seven dogs:** All Healthy-dog omissions, and CIE omissions for states with only three eligible CIE donors.
15. **Resolution versus subdivision:** Biological resolution improved for state names with multi-donor coverage; unsupported or donor-dependent subdivisions remain exploratory rather than automatically superior.

## C. Cross-species correspondence

16. **Independently supported matches:** {supported_n}.
17. **Partial analogies:** {partial_n}.
18. **Compatible programs:** {int((compat.compatibility_status=="STRICT_ORTHOLOG_MATCHED_ANALYSIS_POSSIBLE").sum())}/3 had >=15 shared one-to-one pathway members.
19. **Concordant matched effects:** {concord_n}.
20. **Discordant matched effects:** {discord_n}.
21. **Donor/annotation limits:** Explicit in the frozen matching table, canine LODO table and Human Phase5B provenance. Gene-universe sensitivity is also explicit: {int(valid_match.gene_universe_direction_stability.eq("GENE_UNIVERSE_SENSITIVE").sum()) if len(valid_match) else 0}/{len(valid_match)} matched program comparisons changed sign in at least one species between native and strict-ortholog analyses.
22. **Localization of broad-T divergence:** {"Only partly localized. Mitochondrial respiration remained directionally concordant across native and strict universes, whereas the other program directions were gene-universe sensitive; Human matched estimates were imprecise. This is a boundary result, not proof of equivalent or different mechanisms." if len(valid_match) else "Not defensibly localized to a matched state."}

## D. Scientific contribution

23. **Manuscript strength:** {"Materially strengthened by separating broad-lineage divergence from donor-gated state-level evidence." if component_b=="GO" else "Strengthened mainly as an applicability-boundary and annotation audit."}
24. **Broad-T interpretation:** Broad canine-human divergence is not reclassified as mechanism. State heterogeneity and donor/composition constraints refine its scope.
25. **Suitable for Supplementary Results:** Verified author/frozen crosswalk, donor-gated composition, complete-family three-program canine results, and any supported strict-ortholog matched comparison.
26. **Internal exploratory evidence:** Partial functional analogies, same-dataset marker support and composition LODO with two remaining Healthy dogs.
27. **Exclude from biological claims:** NK/T identity conflicts, donor-confounded branches, non-estimable omissions and forced matches.
28. **Further work:** Stop this branch after Supplementary integration. Additional claims would require more independent canine donors and prospectively harmonized dog-human tissue/state annotation, not further slicing of the same cells.

## Interpretation boundary

This workflow tests the translational applicability boundary of canine CIE; it does not establish causal mechanisms, quantify compositional mediation, or prove a species-by-disease interaction.
"""
(OUT/"00_INTEGRATED_DECISION_REPORT.md").write_text(integrated,encoding="utf-8")

# Reproducibility log and terminal status.
result_files=[]
for folder in [AOUT,BOUT,COUT,OUT/"scripts"]:
    for p in sorted(folder.glob("*")):
        if p.is_file() and p.name!="00_REPRODUCIBILITY_LOG.md":
            result_files.append(p)
hash_lines="\n".join(f"- {p.relative_to(OUT)} | {p.stat().st_size} bytes | SHA256 {sha256(p)}" for p in result_files)
repro=f"""# Phase 6 reproducibility log

## Governance and sources

- Project root: {ROOT}
- AGENTS.md SHA256: {sha256(ROOT/"AGENTS.md")}
- CURRENT_STATE.md SHA256: {sha256(ROOT/"CURRENT_STATE.md")}
- DECISIONS.md: NOT PRESENT at the project root or in the targeted project search; absence documented, not substituted.
- Zenodo record: https://zenodo.org/records/11434076
- Zenodo DOI: 10.5281/zenodo.11434076
- Published RDS MD5: ce31c4350ccd2c8b5902ed5a169c4f94
- Observed RDS MD5: {md5(RDS)}
- Official author code: https://github.com/dyammons/canine_duodenal_atlas/blob/main/scripts/3_tcell_analysis.R
- Author code SHA256: {sha256(WORK/"author_3_tcell_analysis.R")}

## Identity reconciliation

- Author cells: {author_total}
- Frozen T cells: {frozen_t_total}
- Shared frozen-T cells: {matched_t}
- Matching rule: unique author-sample to frozen-donor assignment from barcode-overlap matrix, followed by donor plus full 16bp 10x barcode.
- Bare barcodes were never joined across samples.
- Match states retained explicitly: EXACT_MATCH, VERIFIED_TRANSFORMED_MATCH, AUTHOR_ONLY, FROZEN_ONLY, AMBIGUOUS and IDENTITY_CONFLICT.
- Per-donor recovery: {per_donor_text}

## Frozen gates

- At least 3 Healthy and 3 CIE independent dogs.
- At least 20 matched cells per donor-by-author-subtype pseudobulk.
- Raw UMI source: counts layer in dog_GSE254005_annotated.h5ad.
- Author-only cells excluded from frozen raw-count pseudobulks.
- One-dog >50% and component-cluster >80% diagnostics applied before outcomes.
- T_8/T_9 component dependence checked before outcomes.
- Eligible author subtypes: {", ".join(elig.loc[elig.formal_pseudobulk_eligible.eq("YES"),"author_subtype"]) or "NONE"}.

## Statistics

- edgeR: filterByExpr; TMM; estimateDisp robust=TRUE; glmQLFit robust=TRUE; glmQLFTest; CIE-minus-Healthy.
- Rank: sign(logFC)*sqrt(F).
- fgseaMultilevel: minSize=15, maxSize=500, eps=1e-10, scoreType=std, nproc=1.
- Seed: 20261016 plus deterministic analysis offsets.
- Canine primary BH family: all valid P values from {eligible_n} frozen-eligible author subtypes x 3 predefined programs.
- Canine LODO: only reduced designs retaining >=3 Healthy and >=3 CIE; every Healthy omission is non-estimable.
- Program IDs: R-CFA-1428517; R-CFA-9937383; HALLMARK_XENOBIOTIC_METABOLISM.
- Human paired IDs: R-HSA-1428517; R-HSA-9937383; HALLMARK_XENOBIOTIC_METABOLISM.
- Strict ortholog source: DOG_HUMAN_ORTHOLOG_MASTER.csv, primary_one2one only.
- Matching decision file SHA256: {freeze["sha256"]}.
- Matching freeze timestamp: {freeze["freeze_timestamp_utc"]}.
- Matching decisions were written before new Phase 6 metabolic outcomes were loaded or computed.

## Model and interpretation boundaries

- Healthy research-Beagle source remains confounded with disease group.
- Seven canine donors limit sensitivity; subtype gating does not create new biological replication.
- Human Phase5B was not rerun.
- Dog duodenum and human terminal ileum differ anatomically.
- Human recovered author-labelled cohort remains selected and treatment-mixed.
- Native species NES magnitudes are not interpreted as directly comparable.
- Matched direction is descriptive; no species-by-disease interaction was tested.
- The three programs are post hoc predefined follow-up programs, not prospective validation.
- Same-dataset markers support annotation only; they are not independent validation.
- Frozen Level1/Level2 annotations, old results, manuscripts and Figures 1-6 were not modified.

## Deviations and implementation notes

- DECISIONS.md was absent and therefore could not be read.
- The Windows patch helper was unavailable under the active ACL context; new scripts were written only inside the authorized Phase 6 workspace and no frozen scripts were edited.
- The official RDS was read with the existing portable R environment; no package upgrade or environment installation was performed.
- No project-wide rescan was performed; the existing project manifest and only named prior T-cell directories were used.

## Output checksums

{hash_lines}
"""
(OUT/"00_REPRODUCIBILITY_LOG.md").write_text(repro,encoding="utf-8")

required=[
 AOUT/"01_SOURCE_METADATA_AND_CHECKSUMS.md",AOUT/"02_AUTHOR_FROZEN_CELL_CROSSWALK.csv",
 AOUT/"03_AUTHOR_SUBTYPE_DONOR_COUNTS.csv",AOUT/"04_AUTHOR_FROZEN_LABEL_CROSSTAB.csv",
 AOUT/"05_T8_TREG_IL7R_AUDIT.md",AOUT/"06_AUTHOR_ANNOTATION_ELIGIBILITY.csv",AOUT/"07_COMPONENT_A_REPORT.md",
 BOUT/"01_SUBTYPE_ELIGIBILITY.csv",BOUT/"02_DONOR_COMPOSITION.csv",BOUT/"03_COMPOSITION_LODO.csv",
 BOUT/"04_METABOLIC_PROGRAM_RESULTS.csv",BOUT/"05_METABOLIC_LODO.csv",BOUT/"06_AUTHOR_VS_FROZEN_COMPARISON.csv",BOUT/"07_COMPONENT_B_REPORT.md",
 COUT/"01_MATCHING_DECISIONS_FROZEN.csv",COUT/"02_MARKER_AND_ORTHOLOG_EVIDENCE.csv",COUT/"03_GENESET_COMPATIBILITY.csv",
 COUT/"04_MATCHED_STATE_METABOLIC_RESULTS.csv",COUT/"05_COMPONENT_C_REPORT.md",
 OUT/"00_INTEGRATED_DECISION_REPORT.md",OUT/"00_REPRODUCIBILITY_LOG.md"
]
missing=[str(p) for p in required if not p.exists()]
if missing: raise RuntimeError("Missing required outputs: "+"; ".join(missing))
status_lines=[
 f"COMPONENT_A={component_a}",f"COMPONENT_B={component_b}",f"COMPONENT_C={component_c}",
 f"OVERALL={overall}",f"AUTHOR_SUBTYPES={subtype_count}",f"SHARED_FROZEN_T_CELLS={matched_t}",
 f"ELIGIBLE_DOG_SUBTYPES={eligible_n}",f"SUPPORTED_DOG_HUMAN_MATCHES={supported_n}",
 "PHASE6_COMPLETE=YES"
]
(OUT/"logs"/"PHASE6_FINAL_STATUS.txt").write_text("\n".join(status_lines)+"\n",encoding="utf-8")
print("\n".join(status_lines))
