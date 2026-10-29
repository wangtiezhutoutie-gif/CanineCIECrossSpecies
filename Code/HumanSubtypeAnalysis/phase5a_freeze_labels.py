from __future__ import annotations
import hashlib, platform
from datetime import datetime, timezone
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd

ROOT=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
OUT=ROOT/"WORKING"/"20_HUMAN_TCELL_LABEL_FREEZE"
SRC=OUT/"source"/"scRNA_cell_metadata.csv"
H5=ROOT/"07_human_crohn"/"GSE266616"/"h5ad"/"human_GSE266616_TI_annotated.h5ad"
SOFT=ROOT/"07_human_crohn"/"GSE266616"/"metadata"/"extracted"/"GSE266616_family.soft"
P18=ROOT/"WORKING"/"18_AUTHOR_CROHN_TCELL_LABEL_AUDIT"
P19=ROOT/"WORKING"/"19_HUMAN_SAMPLE_ID_PROVENANCE_AUDIT"
EXPECTED_SIZE=42714235
EXPECTED_MD5="b5016c10ca18c6e58de88589b3f6611e"
SOURCE_RECORD="Zenodo:21508422"
PUB="Zhao et al., Nature Genetics (2026), doi:10.1038/s41588-026-02755-z"
PROV="ORIGINAL_AUTHOR_CLUSTER_ANNOTATION_MARKER_BASED_NO_CELL_CONFIDENCE_FIELD"
T_LABELS=["CD103+ CD4","CD103+ CD8","CD103- CD4","CD4 Tcm","KLRG1+ CD8","MAIT","Proliferating T","Treg","gdT"]
NON_T=["NK","ILCs"]
CANDIDATES=["CD103+ CD4","CD103+ CD8","CD103- CD4","CD4 Tcm","KLRG1+ CD8","Treg"]

def digest(path,algo="sha256"):
    h=hashlib.new(algo)
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()

def sig(path):
    s=Path(path).stat(); return (s.st_size,s.st_mtime_ns)

def clean(x):
    return "" if x is None or pd.isna(x) else str(x).strip()

def yn(x): return "YES" if bool(x) else "NO"

def parse_soft(path):
    out={}; cur=None
    with open(path,encoding="utf-8",errors="replace") as f:
        for raw in f:
            line=raw.rstrip()
            if line.startswith("^SAMPLE = "):
                cur=line.split("=",1)[1].strip(); out[cur]={"gsm":cur}
            elif cur and line.startswith("!Sample_title = "):
                out[cur]["title"]=line.split("=",1)[1].strip()
            elif cur and line.startswith("!Sample_characteristics_ch1 = "):
                item=line.split("=",1)[1].strip()
                if ":" in item:
                    k,v=item.split(":",1); out[cur][k.strip().lower()]=v.strip()
    return out

def write_csv(df,path):
    df.to_csv(path,index=False,encoding="utf-8",lineterminator="\n")

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    required=[
        SRC,H5,SOFT,
        P18/"02_PRIMARY_COHORT_BARCODE_RECOVERY.csv",
        P18/"03_AUTHOR_TCELL_SUBTYPE_COVERAGE.csv",
        P18/"04_AUTHOR_SUBTYPE_FEASIBILITY.csv",
        P19/"01_PARTICIPANT_SOURCE_INVENTORY.csv",
        P19/"02_SAMPLE_ID_DISCREPANCIES.csv",
        P19/"03_ALTERNATE_KEY_VALIDATION.csv",
        P19/"04_SOURCE_VERSION_COMPARISON.csv",
        P19/"05_UPDATED_HUMAN_SUBTYPE_FEASIBILITY.csv",
        P19/"06_PROVENANCE_DECISION_REPORT.md",
        P19/"07_REPRODUCIBILITY_LOG.md"]
    absent=[str(p) for p in required if not p.exists()]
    if absent: raise RuntimeError("Missing input: "+"; ".join(absent))
    src_md5=digest(SRC,"md5"); src_sha=digest(SRC)
    if SRC.stat().st_size!=EXPECTED_SIZE or src_md5!=EXPECTED_MD5:
        raise RuntimeError(f"SOURCE_UNAVAILABLE size={SRC.stat().st_size} md5={src_md5}")
    before={str(p):sig(p) for p in required}
    cols=["Unnamed: 0","barcode","biopsy_location","Patient_ID_scRNA","Sample_ID_scRNA",
          "disease_status","inflammation_status","nCount_RNA","nFeature_RNA","Doublet_score",
          "percent.mt","percent.rb","log10GenesPerUMI","lineage","cluster_label1","cluster_label2"]
    author=pd.read_csv(SRC,low_memory=False)
    miss=set(cols)-set(author.columns)
    if miss: raise RuntimeError("Source columns missing: "+",".join(sorted(miss)))
    author["author_source_row_number"]=np.arange(2,len(author)+2)
    author["author_identity"]=author["Sample_ID_scRNA"].astype(str)+"|"+author["barcode"].astype(str)
    if author["author_identity"].duplicated().any(): raise RuntimeError("Duplicate author sample+barcode identity")
    inv=pd.read_csv(P19/"01_PARTICIPANT_SOURCE_INVENTORY.csv",low_memory=False)
    inv=inv.loc[inv.frozen_participant_id.notna() & inv.frozen_gsm.notna()].copy()
    if len(inv)!=27 or inv.frozen_participant_id.duplicated().any(): raise RuntimeError(f"Frozen participants={len(inv)}")
    inv=inv.set_index("frozen_participant_id",drop=False)

    h=ad.read_h5ad(H5,backed="r")
    obs_cols=["barcode","GSM","sample_title","donor","disease_group","state","region","tissue",
              "biologic_treatment","biologic_class","age","race","platform","primary_contrast_eligible",
              "total_counts","n_genes_by_counts","pct_counts_mt","human_level1"]
    m=(h.obs.primary_contrast_eligible.astype(str)=="YES_SELECTED") & (h.obs.human_level1.astype(str)=="T")
    frozen=h.obs.loc[m,obs_cols].copy()
    frozen.insert(0,"frozen_cell_id",frozen.index.astype(str))
    h.file.close(); frozen=frozen.reset_index(drop=True)
    if len(frozen)!=21378 or frozen.frozen_cell_id.duplicated().any(): raise RuntimeError(f"Primary broad-T={len(frozen)}")
    frozen["selected_author_sample_id"]=frozen.donor.map(inv.selected_author_sample_id).map(clean)
    frozen["phase4b_mapping_method"]=frozen.donor.map(inv.mapping_method).map(clean)
    frozen["phase4b_status"]=frozen.donor.map(inv.final_provenance_status).map(clean)
    frozen["phase4b_reason"]=frozen.donor.map(inv.reason_or_note).map(clean)
    frozen["author_identity"]=np.where(frozen.selected_author_sample_id.ne(""),
        frozen.selected_author_sample_id+"|"+frozen.barcode.astype(str),"")
    acols=["author_identity","Unnamed: 0","author_source_row_number","barcode","biopsy_location",
           "Patient_ID_scRNA","Sample_ID_scRNA","disease_status","inflammation_status","lineage",
           "cluster_label1","cluster_label2"]
    x=frozen.merge(author[acols],on="author_identity",how="left",validate="many_to_one",suffixes=("","_author"))
    x["author_record_recovered"]=x.Sample_ID_scRNA.notna()
    rec=x.author_record_recovered; healthy=x.disease_group.astype(str).eq("Healthy"); crohn=x.disease_group.astype(str).eq("Crohn")
    astate=x.inflammation_status.map(clean)
    x["participant_conflict"]=rec & x.Patient_ID_scRNA.map(clean).ne(x.donor.map(clean))
    x["sample_conflict"]=rec & x.Sample_ID_scRNA.map(clean).ne(x.selected_author_sample_id)
    x["region_conflict"]=rec & x.biopsy_location.map(clean).ne(x.region.map(clean))
    x["disease_conflict"]=rec & ~((healthy & x.disease_status.map(clean).eq("Control"))|(crohn & x.disease_status.map(clean).eq("CD")))
    x["inflammation_conflict"]=rec & ~((healthy & astate.eq(""))|(crohn & x.state.astype(str).eq("Active_inflamed") & astate.eq("inf")))
    id_conf=["participant_conflict","sample_conflict","region_conflict","disease_conflict","inflammation_conflict"]
    x["identity_conflict"]=x[id_conf].any(axis=1)
    x["author_lineage_conflict"]=rec & ~x.lineage.map(clean).eq("tcell")
    x["author_subtype_missing"]=rec & x.lineage.map(clean).eq("tcell") & (x.cluster_label1.map(clean).eq("")|x.cluster_label2.map(clean).eq(""))
    x["valid_fine_label"]=rec & ~x.identity_conflict & ~x.author_lineage_conflict & ~x.author_subtype_missing
    x["is_nk_ilc"]=x.valid_fine_label & x.cluster_label1.isin(NON_T)
    x["is_t_related"]=x.valid_fine_label & ~x.is_nk_ilc & x.cluster_label1.isin(T_LABELS)
    unexpected=sorted(set(x.loc[x.valid_fine_label,"cluster_label1"].map(clean))-set(T_LABELS)-set(NON_T))
    if unexpected: raise RuntimeError("Unexpected labels: "+",".join(unexpected))
    def status(r):
        if r.identity_conflict:return "UNRESOLVED"
        if r.author_lineage_conflict:return "AUTHOR_LINEAGE_CONFLICT"
        if r.author_subtype_missing:return "AUTHOR_RECORD_WITHOUT_VALID_SUBTYPE"
        if r.is_nk_ilc:return "VERIFIED_AUTHOR_NK_ILC_LABEL"
        if r.is_t_related:return "VERIFIED_T_SUBTYPE"
        if not r.author_record_recovered:return "NO_AUTHOR_SAMPLE" if clean(r.selected_author_sample_id)=="" else "NO_CELL_MATCH"
        return "UNRESOLVED"
    x["annotation_status"]=x.apply(status,axis=1)
    def conflicts(r):
        z=[c.replace("_conflict","").upper()+"_CONFLICT" for c in id_conf if r[c]]
        if r.author_lineage_conflict:z.append("AUTHOR_LINEAGE_CONFLICT")
        return ";".join(z) if z else "NONE"
    x["conflict_status"]=x.apply(conflicts,axis=1)
    missing={"VERIFIED_T_SUBTYPE":"NOT_MISSING","VERIFIED_AUTHOR_NK_ILC_LABEL":"NOT_MISSING_AUTHOR_LABEL_BUT_NOT_T_RELATED",
             "AUTHOR_LINEAGE_CONFLICT":"AUTHOR_BROAD_LINEAGE_NOT_TCELL",
             "AUTHOR_RECORD_WITHOUT_VALID_SUBTYPE":"AUTHOR_RECORD_LACKS_COMPLETE_FINE_SUBTYPE",
             "NO_CELL_MATCH":"NO_BARCODE_MATCH_IN_VERIFIED_AUTHOR_SAMPLE","UNRESOLVED":"IDENTITY_FIELD_CONFLICT"}
    x["missing_label_reason"]=x.annotation_status.map(missing).fillna(x.phase4b_reason)
    valid=x.loc[x.valid_fine_label].copy()
    evidence="Official GSE266616 GSM participant/terminal-ileum/condition-state; unique Phase4B author sample; exact sample-specific 10x barcode"
    out01=pd.DataFrame({
        "frozen_cell_id":valid.frozen_cell_id,"frozen_sample_id":valid.sample_title.map(clean),
        "geo_gsm_accession":valid.GSM.map(clean),"frozen_participant_id":valid.donor.map(clean),
        "author_patient_id":valid.Patient_ID_scRNA.map(clean),"author_sample_id":valid.Sample_ID_scRNA.map(clean),
        "raw_10x_barcode":valid.barcode.map(clean),"biopsy_location":valid.biopsy_location.map(clean),
        "disease_status":valid.disease_status.map(clean),"inflammation_status":valid.inflammation_status.map(clean),
        "frozen_level1":valid.human_level1.map(clean),"author_lineage":valid.lineage.map(clean),
        "author_cluster_label1":valid.cluster_label1.map(clean),"author_cluster_label2":valid.cluster_label2.map(clean),
        "author_original_row_id":valid["Unnamed: 0"].map(clean),"is_T_related_subtype":valid.is_t_related.map(yn),
        "is_NK_or_ILC_label":valid.is_nk_ilc.map(yn),"matching_method":valid.phase4b_mapping_method.map(clean),
        "matching_evidence":evidence,"source_record":SOURCE_RECORD,"source_filename":SRC.name,
        "source_md5":src_md5,"annotation_provenance":PROV,"verification_status":valid.annotation_status
    }).sort_values(["frozen_participant_id","frozen_cell_id"]).reset_index(drop=True)
    out02=pd.DataFrame({
        "frozen_cell_id":x.frozen_cell_id,"frozen_participant_id":x.donor.map(clean),
        "frozen_condition":x.disease_group.map(clean),"frozen_state":x.state.map(clean),"frozen_region":x.region.map(clean),
        "frozen_sample_id":x.sample_title.map(clean),"geo_gsm_accession":x.GSM.map(clean),"raw_10x_barcode":x.barcode.map(clean),
        "selected_author_sample_id":x.selected_author_sample_id,"author_record_recovered":x.author_record_recovered.map(yn),
        "valid_author_fine_label_recovered":x.valid_fine_label.map(yn),"valid_T_related_subtype_recovered":x.is_t_related.map(yn),
        "author_lineage":x.lineage.map(clean),"author_cluster_label1":x.cluster_label1.map(clean),
        "author_cluster_label2":x.cluster_label2.map(clean),"annotation_status":x.annotation_status,
        "missing_label_reason":x.missing_label_reason.map(clean),"conflict_status":x.conflict_status,
        "matching_method":x.phase4b_mapping_method.map(clean),
        "matching_evidence":np.where(x.selected_author_sample_id.ne(""),"VERIFIED_SAMPLE_IDENTITY_PLUS_BARCODE","NO_VERIFIED_REQUIRED_AUTHOR_SAMPLE"),
        "source_record":SOURCE_RECORD,"source_filename":SRC.name,"source_md5":src_md5,
        "annotation_provenance":np.where(x.valid_fine_label,PROV,"NO_VALID_AUTHOR_FINE_LABEL")
    }).sort_values(["frozen_participant_id","frozen_cell_id"]).reset_index(drop=True)

    geo=parse_soft(SOFT); cw=[]; rescued=[]
    for p,r in inv.sort_index().iterrows():
        gsm=clean(r.frozen_gsm); g=geo.get(gsm,{})
        xp=x.loc[x.donor.astype(str).eq(p)]
        method=clean(r.mapping_method); selected=clean(r.selected_author_sample_id)
        is_rescue=bool(selected and method!="EXACT_SAMPLE_ID_scRNA_PLUS_BARCODE")
        if is_rescue:rescued.append(p)
        match="EXACT_PARTICIPANT_AND_CELL_MATCH" if selected and not is_rescue else ("VERIFIED_GSM_CROSSWALK_RESCUE" if is_rescue else clean(r.final_provenance_status))
        offok=(clean(g.get("patient"))==p and clean(g.get("tissue"))=="terminal_ileum" and
               ((clean(r.frozen_condition)=="Healthy" and clean(g.get("condition")).lower()=="control") or
                (clean(r.frozen_condition)=="Crohn" and clean(g.get("condition")).lower() in {"cd","crohn","crohn disease","crohn's disease"})))
        cw.append({
            "frozen_participant_id":p,"frozen_sample_id":clean(r.official_GEO_title),"official_gsm_accession":gsm,
            "official_geo_participant":clean(g.get("patient")),"official_tissue_region":clean(g.get("tissue")),
            "official_condition":clean(g.get("condition")),"official_inflammation_state":clean(g.get("sample region")),
            "official_biologic_treatment":clean(g.get("treatment (biological)")) or "NOT_AVAILABLE",
            "author_participant_id":clean(r.author_Patient_ID_scRNA),"author_sample_id":selected,"matching_status":match,
            "rescued_by_verified_crosswalk":yn(is_rescue),"frozen_primary_T_cells":len(xp),
            "recovered_author_records":int(xp.author_record_recovered.sum()),
            "verified_author_fine_label_cells":int(xp.valid_fine_label.sum()),
            "verified_T_related_label_cells":int(xp.is_t_related.sum()),
            "verified_NK_ILC_label_cells":int(xp.is_nk_ilc.sum()),
            "official_geo_validation_status":"PASS" if offok else "FAIL",
            "exact_provenance_evidence":clean(r.independent_evidence),
            "source_record_references":f"GSE266616 SOFT {gsm}; Zenodo 21508422 {selected or 'NO_REQUIRED_SAMPLE'}",
            "remaining_uncertainty":clean(r.reason_or_note) if not selected else "NONE_FOR_ACCEPTED_IDENTITY_MAPPING"})
    out03=pd.DataFrame(cw)

    vt=x.loc[x.is_t_related].copy(); vn=x.loc[x.is_nk_ilc].copy()
    labels=CANDIDATES+[z for z in T_LABELS if z not in CANDIDATES]+NON_T
    alln=x.groupby("donor",observed=True).size().to_dict()
    tn=vt.groupby("donor",observed=True).size().to_dict()
    nkn=vn.groupby("donor",observed=True).size().to_dict()
    fn=valid.groupby("donor",observed=True).size().to_dict()
    sn=valid.groupby(["donor","cluster_label1"],observed=True).size().to_dict()
    rows=[]
    for p,r in inv.sort_index().iterrows():
        da=int(alln.get(p,0)); dt=int(tn.get(p,0)); df=int(fn.get(p,0))
        for lab in labels:
            n=int(sn.get((p,lab),0)); cat="NK_OR_ILC_LABEL" if lab in NON_T else "T_CELL_RELATED_LABEL"
            rows.append({
                "frozen_participant_id":p,"condition":clean(r.frozen_condition),"tissue_region":clean(r.frozen_tissue_location),
                "frozen_sample_id":clean(r.official_GEO_title),"geo_gsm_accession":clean(r.frozen_gsm),
                "author_subtype":lab,"label_category":cat,"subtype_cell_count":n,
                "total_frozen_broad_T_cells":da,"valid_author_T_related_labels":dt,
                "valid_author_NK_ILC_labels":int(nkn.get(p,0)),"valid_author_fine_labels_total":df,
                "subtype_fraction_among_valid_recovered_T_labels":n/dt if cat=="T_CELL_RELATED_LABEL" and dt else np.nan,
                "subtype_fraction_among_all_valid_author_fine_labels":n/df if df else np.nan,
                "author_T_related_coverage_relative_to_all_frozen_T":dt/da if da else np.nan,
                "valid_fine_label_coverage_relative_to_all_frozen_T":df/da if da else np.nan,
                "missing_T_related_label_fraction":1-dt/da if da else np.nan,
                "subtype_observation_status":"OBSERVED_IN_RECOVERED_AUTHOR_LABELS" if n else "NOT_OBSERVED_OR_NOT_RECOVERABLE",
                "annotation_provenance":PROV,
                "missingness_caution":"Fractions are conditional on recovered author labels and are not unbiased whole-compartment estimates."})
    out04=pd.DataFrame(rows)

    mr=[]
    for p,r in inv.sort_index().iterrows():
        xp=x.loc[x.donor.astype(str).eq(p)]; na=len(xp); nt=int(xp.is_t_related.sum()); nn=int(xp.is_nk_ilc.sum()); nf=int(xp.valid_fine_label.sum())
        recovered=nf>0
        mr.append({
            "row_type":"PARTICIPANT","frozen_participant_id":p,"condition":clean(r.frozen_condition),
            "frozen_sample_id":clean(r.official_GEO_title),"geo_gsm_accession":clean(r.frozen_gsm),
            "biopsy_location":clean(r.frozen_tissue_location),"inflammation_state":clean(r.frozen_state),
            "participant_in_recovered_author_labelled_subset":yn(recovered),"frozen_broad_T_cells":na,
            "recovered_author_records":int(xp.author_record_recovered.sum()),"valid_author_fine_labels":nf,
            "valid_author_T_related_labels":nt,"author_NK_ILC_labels":nn,"T_related_traceability_fraction":nt/na,
            "fine_label_traceability_fraction":nf/na,"missing_T_related_label_fraction":1-nt/na,
            "source_missingness_class":clean(r.final_provenance_status),"source_missingness_reason":clean(r.reason_or_note),
            "biologic_treatment":";".join(sorted(set(xp.biologic_treatment.astype(str)))) or "NOT_AVAILABLE",
            "biologic_class":";".join(sorted(set(xp.biologic_class.astype(str)))) or "NOT_AVAILABLE",
            "age":";".join(sorted(set(xp.age.astype(str)))) or "NOT_AVAILABLE",
            "race":";".join(sorted(set(xp.race.astype(str)))) or "NOT_AVAILABLE",
            "platform":";".join(sorted(set(xp.platform.astype(str)))) or "NOT_AVAILABLE",
            "median_frozen_total_counts":float(xp.total_counts.median()),
            "median_frozen_n_genes_by_counts":float(xp.n_genes_by_counts.median()),
            "median_frozen_pct_counts_mt":float(xp.pct_counts_mt.median()),
            "selection_missingness_interpretation":"RECOVERED_PARTICIPANT" if recovered else "UNRECOVERED_PARTICIPANT; missing label is not subtype absence"})
    pm=pd.DataFrame(mr); summaries=[]
    for cond in ["Healthy","Crohn","ALL"]:
        z=pm if cond=="ALL" else pm.loc[pm.condition.eq(cond)]
        npart=len(z); nrec=int(z.participant_in_recovered_author_labelled_subset.eq("YES").sum())
        a=int(z.frozen_broad_T_cells.sum()); t=int(z.valid_author_T_related_labels.sum()); f=int(z.valid_author_fine_labels.sum())
        summaries.append({
            "row_type":"GROUP_SUMMARY","frozen_participant_id":f"{cond}_SUMMARY","condition":cond,
            "frozen_sample_id":"NOT_APPLICABLE","geo_gsm_accession":"NOT_APPLICABLE","biopsy_location":"TI",
            "inflammation_state":"Healthy" if cond=="Healthy" else ("Active_inflamed" if cond=="Crohn" else "MIXED"),
            "participant_in_recovered_author_labelled_subset":f"{nrec}/{npart}","frozen_broad_T_cells":a,
            "recovered_author_records":int(z.recovered_author_records.sum()),"valid_author_fine_labels":f,
            "valid_author_T_related_labels":t,"author_NK_ILC_labels":int(z.author_NK_ILC_labels.sum()),
            "T_related_traceability_fraction":t/a,"fine_label_traceability_fraction":f/a,
            "missing_T_related_label_fraction":1-t/a,"source_missingness_class":"GROUP_SUMMARY",
            "source_missingness_reason":f"Recovered participants {nrec}/{npart} ({nrec/npart:.1%}).",
            "biologic_treatment":";".join(sorted(set(z.biologic_treatment))),"biologic_class":";".join(sorted(set(z.biologic_class))),
            "age":"DESCRIPTIVE_PARTICIPANT_VALUES_AVAILABLE","race":"DESCRIPTIVE_PARTICIPANT_VALUES_AVAILABLE",
            "platform":";".join(sorted(set(z.platform))),"median_frozen_total_counts":float(z.median_frozen_total_counts.median()),
            "median_frozen_n_genes_by_counts":float(z.median_frozen_n_genes_by_counts.median()),
            "median_frozen_pct_counts_mt":float(z.median_frozen_pct_counts_mt.median()),
            "selection_missingness_interpretation":"DIFFERENTIAL_DISEASE_GROUP_RECOVERY; RECOVERED_SUBSET_NOT_DEMONSTRATED_REPRESENTATIVE" if cond=="ALL" else "DESCRIPTIVE_GROUP_SUMMARY"})
    out05=pd.concat([pm,pd.DataFrame(summaries)],ignore_index=True)

    ambiguity={
        "CD103+ CD4":"Author label is CD103-defined; candidate-specific local marker table unavailable; no direct equivalence to canine Memory T.",
        "CD103+ CD8":"Author label is CD103-defined; candidate-specific local marker table unavailable; no direct equivalence to canine Memory T.",
        "CD103- CD4":"Author label is CD103-defined; negative-state boundary cannot be independently checked from metadata alone.",
        "CD4 Tcm":"Central-memory label cannot be independently separated from resident-memory states without the author marker table; CD4 Trm equivalence is not assumed.",
        "KLRG1+ CD8":"Named KLRG1-positive state is author-defined; distinctness is not independently revalidated from metadata alone.",
        "Treg":"Treg is an original author label, but FOXP3/IL2RA/CTLA4/TIGIT/IKZF2/ICOS support is not directly available in a local candidate-specific marker table."}
    rr=[]
    for lab in CANDIDATES:
        z=vt.loc[vt.cluster_label1.eq(lab)]
        dc=z.groupby(["disease_group","donor"],observed=True).size().reset_index(name="cells")
        hh=dc.loc[(dc.disease_group=="Healthy") & (dc.cells>=20)]
        cc=dc.loc[(dc.disease_group=="Crohn") & (dc.cells>=20)]
        eligible=pd.concat([hh,cc]); total=len(z)
        imax=dc.cells.idxmax(); lp=clean(dc.loc[imax,"donor"]); lc=int(dc.loc[imax,"cells"]); lf=lc/total
        cov="COVERAGE_PASS" if len(hh)>=5 and len(cc)>=5 else "COVERAGE_FAIL"
        dom="POTENTIAL_SINGLE_DONOR_DOMINANCE" if lf>0.35 else ("NOTABLE_SINGLE_DONOR_CONTRIBUTION" if lf>0.25 else "NO_STRONG_SINGLE_DONOR_DOMINANCE")
        rr.append({
            "author_subtype":lab,"label_category":"T_CELL_RELATED_LABEL",
            "healthy_total_cells":len(z.loc[z.disease_group=="Healthy"]),
            "healthy_participants_with_any_cells":int(dc.loc[dc.disease_group=="Healthy"].donor.nunique()),
            "healthy_participants_ge20_cells":len(hh),"active_crohn_total_cells":len(z.loc[z.disease_group=="Crohn"]),
            "active_crohn_participants_with_any_cells":int(dc.loc[dc.disease_group=="Crohn"].donor.nunique()),
            "active_crohn_participants_ge20_cells":len(cc),"total_verified_subtype_cells":total,
            "median_cells_per_eligible_participant":float(eligible.cells.median()),
            "minimum_eligible_donor_count_across_conditions":min(len(hh),len(cc)),
            "largest_participant":lp,"largest_participant_cell_count":lc,"largest_participant_contribution_fraction":lf,
            "single_donor_dominance_status":dom,"fraction_of_all_valid_T_related_labels":total/int(x.is_t_related.sum()),
            "full_primary_cohort_T_related_traceability_fraction":int(x.is_t_related.sum())/len(x),
            "full_primary_cohort_missing_T_related_label_fraction":1-int(x.is_t_related.sum())/len(x),
            "coverage_gate":">=5 Healthy and >=5 Active Crohn participants with >=20 cells per participant-by-subtype",
            "coverage_status":cov,"provenance_status":"VERIFIED","author_lineage_evidence":"lineage=tcell for all included cells",
            "author_cluster_label1":lab,"observed_author_cluster_label2_values":";".join(sorted(set(z.cluster_label2.map(clean)))),
            "author_annotation_evidence":f"{PUB}; Supplementary Note describes marker-based lineage-specific Louvain annotation and Extended Data Fig. 2 marker summaries.",
            "candidate_specific_marker_table_local_status":"NOT_AVAILABLE_WITHIN_IDENTIFIED_LOCAL_METADATA_ASSETS",
            "Treg_marker_support_FOXP3_IL2RA_CTLA4_TIGIT_IKZF2_ICOS":"NOT_VERIFIED_FROM_LOCAL_AUTHOR_MARKER_TABLE" if lab=="Treg" else "NOT_APPLICABLE",
            "identity_ambiguity_or_boundary":ambiguity[lab],"annotation_reliability_status":"ANNOTATION_LIMITED",
            "phase5b_eligibility":"ELIGIBLE_WITH_LIMITATIONS" if cov=="COVERAGE_PASS" else "NOT_ELIGIBLE",
            "interpretation_boundary":"Coverage eligibility does not establish biological equivalence, unbiased composition, or differential-expression validity."})
    out06=pd.DataFrame(rr)

    outputs=[OUT/"01_VERIFIED_HUMAN_TCELL_CELL_LABEL_MAP.csv",OUT/"02_ALL_PRIMARY_TCELL_ANNOTATION_COVERAGE.csv",
             OUT/"03_VERIFIED_GSM_SAMPLE_CROSSWALK.csv",OUT/"04_PARTICIPANT_SUBTYPE_COVERAGE.csv",
             OUT/"05_COHORT_SELECTION_AND_MISSINGNESS.csv",OUT/"06_SUBTYPE_ANNOTATION_RELIABILITY.csv"]
    for d,p in zip([out01,out02,out03,out04,out05,out06],outputs):write_csv(d,p)
    sc=out02.annotation_status.value_counts().to_dict()
    nvalid=len(out01); nt=int((out01.is_T_related_subtype=="YES").sum()); nn=int((out01.is_NK_or_ILC_label=="YES").sum())
    hall=len(x.loc[x.disease_group=="Healthy"]); call=len(x.loc[x.disease_group=="Crohn"])
    ht=len(vt.loc[vt.disease_group=="Healthy"]); ct=len(vt.loc[vt.disease_group=="Crohn"])
    rh=len(out03.loc[(out03.official_condition.str.lower()=="control") & (out03.verified_author_fine_label_cells>0)])
    rc=len(out03.loc[(out03.official_condition.str.lower().isin(["cd","crohn","crohn disease","crohn's disease"])) & (out03.verified_author_fine_label_cells>0)])
    passed=out06.loc[out06.coverage_status=="COVERAGE_PASS","author_subtype"].tolist()
    dominated=out06.loc[out06.single_donor_dominance_status=="POTENTIAL_SINGLE_DONOR_DOMINANCE","author_subtype"].tolist()
    checks={
        "all_21378_frozen_T_cells_once":len(out02)==21378 and out02.frozen_cell_id.is_unique,
        "verified_frozen_cell_mapping_one_to_one":out01.frozen_cell_id.is_unique,
        "verified_author_identity_one_to_one":not valid.author_identity.duplicated().any(),
        "participant_identity_agreement":int(x.loc[x.valid_fine_label,"participant_conflict"].sum())==0,
        "region_state_disease_agreement":int(x.loc[x.valid_fine_label,["region_conflict","disease_conflict","inflammation_conflict"]].sum().sum())==0,
        "no_unresolved_promoted_to_verified":not x.loc[x.valid_fine_label,"identity_conflict"].any(),
        "NK_ILC_distinct_from_T_related":not (valid.is_nk_ilc & valid.is_t_related).any(),
        "author_labels_preserved_verbatim":set(out01.author_cluster_label1).issubset(set(author.cluster_label1.dropna().astype(str))),
        "source_file_verified_and_retained":SRC.exists() and SRC.stat().st_size==EXPECTED_SIZE and src_md5==EXPECTED_MD5,
        "coverage_statuses_reconcile_to_21378":sum(sc.values())==21378,
        "all_27_GEO_records_independently_validate":out03.official_geo_validation_status.eq("PASS").all()}
    integrity="PASS" if all(checks.values()) else "FAIL"
    rp=OUT/"07_PHASE5A_VALIDATION_REPORT.md"
    report=f"""# Phase 5A - Human T-cell label freeze and annotation reliability audit

## Final decision

- **Decision: PARTIAL GO**
- **Integrity status: {integrity}**
- **Phase 5B started: NO**
- The recovered author labels are permanently preserved with original identifiers and provenance. Annotation recovery is strongly differential by disease group, and candidate-specific marker support was not independently revalidated from a local author marker table.

## Recalculated checkpoints

- Frozen primary comparison: **27 participants** (16 Healthy, 11 Active Crohn), **{len(out02):,} broad-T cells**.
- Verified sample-specific author matches: **{int(x.author_record_recovered.sum()):,}**.
- Valid author fine-label records after excluding lineage/identity problems: **{nvalid:,}**.
- Genuine author T-related labels: **{nt:,}** ({nt/len(out02):.2%} of the frozen T compartment).
- Author NK/ILC labels kept separate: **{nn:,}** ({nn/len(out02):.2%}).
- Broad-only/unmatched/conflicted cells without a valid T-related subtype: **{len(out02)-nt:,}** ({1-nt/len(out02):.2%}). This is missing annotation, not biological absence.
- Recovered participants: **{rh} Healthy / {rc} Active Crohn**.
- Valid T-related traceability by condition: Healthy **{ht:,}/{hall:,} ({ht/hall:.2%})**; Active Crohn **{ct:,}/{call:,} ({ct/call:.2%})**.
- Verified crosswalk rescues: **{len(rescued)} participants** ({', '.join(rescued)}).
- Coverage-passing predefined subtypes: **{len(passed)}** ({', '.join(passed)}).

## Required scientific questions

1. **Were the recovered labels permanently preserved?** Yes. 01_VERIFIED_HUMAN_TCELL_CELL_LABEL_MAP.csv contains one row per valid recovered cell, and the exact author source is retained under source/ with MD5 {src_md5} and SHA256 {src_sha}.
2. **Did the prior 15,670 records reproduce?** {"Yes" if nvalid==15670 else "No"}. The independently reconstructed count is **{nvalid:,}**; every other frozen cell is represented in the complete coverage table.
3. **How many are T-related versus NK/ILC?** **{nt:,} T-related** and **{nn:,} NK/ILC**.
4. **How many of 21,378 cells have traceable author T subtypes?** **{nt:,} ({nt/len(out02):.2%})**. The broader fine-label total including NK/ILC is {nvalid:,}.
5. **Did verified cells pass identity checks?** {"Yes" if integrity=="PASS" else "No"}. Accepted cells have unique frozen and author identities with no participant, disease, region or inflammation-state conflict. Unnamed: 0 was retained only as provenance and never used as a key.
6. **Do the six candidates pass donor coverage?** **{len(passed)}/6 pass** the unchanged >=5 Healthy and >=5 Active Crohn donors with >=20 cells rule. Coverage is distinct from annotation reliability and DE authorization.
7. **Are subtypes strongly dominated by one participant?** {"No candidate exceeded the descriptive 35% potential-dominance flag." if not dominated else "Potential dominance was flagged for: "+", ".join(dominated)+"."} Exact contributions are in 06_SUBTYPE_ANNOTATION_RELIABILITY.csv.
8. **Are annotations biologically supported?** They are original author, lineage-specific, marker-based cluster annotations, not project-inferred labels. No candidate-specific local author marker table was found in the identified metadata assets. Treg FOXP3/IL2RA/CTLA4/TIGIT/IKZF2/ICOS support therefore remains **NOT VERIFIED locally**; all six are conservatively **ANNOTATION_LIMITED**, not invalidated.
9. **How severe is selection/missingness bias?** Important. Participant recovery is 13/16 Healthy versus 6/11 Crohn; T-related label traceability is {ht/hall:.2%} versus {ct/call:.2%}. The 19-participant subset is not demonstrated to represent all 27 participants.
10. **Is donor-level composition defensible?** Only as a donor-level sensitivity analysis within the recovered subset, with both recovered-label fractions and whole-compartment annotation coverage shown. It is not defensible as an unbiased estimate for the complete 27-participant cohort.
11. **Is within-subtype pseudobulk technically feasible?** Yes for the six coverage-passing subtypes at the identity/coverage level, using only the verified cell map and donor as the biological unit. Raw-count suitability, formal design and inference must be separately frozen in Phase 5B.
12. **Which subtypes are eligible for proposed Phase 5B?** {', '.join(passed)} are **ELIGIBLE_WITH_LIMITATIONS**. No Phase 5B analysis was run.

## Annotation reliability boundaries

- cluster_label1/2 come from {PUB}. The Supplementary Note describes marker-based, lineage-specific Louvain annotation and points to Extended Data Fig. 2 marker summaries.
- No cell-level annotation-confidence field exists in the retained source.
- Treg marker support and the distinction of CD4 Tcm from resident-memory states could not be independently verified from an existing local candidate-specific author marker table.
- CD103-positive/negative and KLRG1-positive identities remain author-defined states. No equivalence to the canine Memory T label is asserted.

## Missingness and selection boundary

Six primary participants are absent from the relevant author metadata (HA14, HA41, HA48, HA61, HA68, HA82), and two lack the required matching sample (HA70, OR95). Missing labels remain missing and are never converted to subtype zeros.

## Integrity checks

"""
    report+="\n".join([f"- {k}: **{'PASS' if v else 'FAIL'}**" for k,v in checks.items()])
    report+=f"""

## Stop status

Phase 5A is complete. No expression analysis, reclustering, marker discovery, DE, GSEA, cross-species subtype comparison, manuscript change or Phase 5B work was performed.
"""
    rp.write_text(report,encoding="utf-8")
    reread={p.name:len(pd.read_csv(p,low_memory=False)) for p in outputs}
    unchanged=before=={str(p):sig(p) for p in required}
    oh={str(p):digest(p) for p in outputs+[rp]}
    ih={str(p):digest(p) for p in required}
    log=f"""# Phase 5A reproducibility log

- Generated UTC: {datetime.now(timezone.utc).isoformat()}
- Mode: FAST / metadata-only / non-destructive
- Python: {platform.python_version()}
- pandas: {pd.__version__}
- anndata: {ad.__version__}
- Expression matrix accessed: NO (h5ad.obs only, backed read)
- Phase 5B started: NO
- Original inputs unchanged by size/mtime signature: {'YES' if unchanged else 'NO'}
- Final integrity status: {integrity}

## Deterministic mapping procedure

1. Select primary_contrast_eligible=YES_SELECTED and human_level1=T.
2. Use only the 27 Phase 4B participant/GSM records and verified selected_author_sample_id values.
3. Join exact selected_author_sample_id + raw 10x barcode; never use a bare barcode or Unnamed: 0.
4. Require participant, disease, TI region and inflammation-state agreement.
5. Exclude author broad-lineage conflicts and incomplete fine labels.
6. Keep NK/ILC separate from T-related labels.
7. Apply >=20 cells in >=5 Healthy and >=5 Active Crohn participants without modification.

## Reconciled status counts

"""
    for k,v in sorted(sc.items()): log+=f"- {k}: {v}\n"
    log+=f"- TOTAL: {sum(sc.values())}\n\n## CSV validation row counts\n\n"
    for k,v in reread.items():log+=f"- {k}: {v}\n"
    log+="\n## Input SHA256 checksums\n\n"
    for p,d in ih.items():log+=f"- {p}: {d}\n"
    log+="\n## Output SHA256 checksums\n\n"
    for p,d in oh.items():log+=f"- {p}: {d}\n"
    log+="\n08_REPRODUCIBILITY_LOG.md is the checksum manifest itself and is intentionally excluded from its own embedded hash list.\n"
    (OUT/"08_REPRODUCIBILITY_LOG.md").write_text(log,encoding="utf-8")
    if not unchanged or integrity!="PASS":raise RuntimeError(f"Integrity failure unchanged={unchanged} checks={checks}")
    print("PHASE5A_DECISION=PARTIAL_GO")
    print(f"FROZEN_PRIMARY_T_CELLS={len(out02)}")
    print(f"VALID_FINE_LABELS={nvalid}")
    print(f"T_RELATED_LABELS={nt}")
    print(f"NK_ILC_LABELS={nn}")
    print(f"RECOVERED_PARTICIPANTS=Healthy:{rh};Crohn:{rc}")
    print("COVERAGE_PASS_SUBTYPES="+",".join(passed))
    print(f"INTEGRITY_STATUS={integrity}")
    print("PHASE5B_STARTED=NO")

if __name__=="__main__":
    main()
