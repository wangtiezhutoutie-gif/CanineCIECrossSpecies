from __future__ import annotations
import hashlib,json,time
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
ROOT=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"]); OUT=ROOT/"WORKING"/"21_HUMAN_TCELL_SUBTYPE_METABOLIC_ANALYSIS"; WORK=OUT/"work"
P5A=ROOT/"WORKING"/"20_HUMAN_TCELL_LABEL_FREEZE"; H5AD=ROOT/"07_human_crohn"/"GSE266616"/"h5ad"/"human_GSE266616_TI_annotated.h5ad"
MAP=P5A/"01_VERIFIED_HUMAN_TCELL_CELL_LABEL_MAP.csv"; CROSSWALK=P5A/"03_VERIFIED_GSM_SAMPLE_CROSSWALK.csv"; COVERAGE=P5A/"04_PARTICIPANT_SUBTYPE_COVERAGE.csv"; SELECTION=P5A/"05_COHORT_SELECTION_AND_MISSINGNESS.csv"; RELIABILITY=P5A/"06_SUBTYPE_ANNOTATION_RELIABILITY.csv"
SUBTYPES=["CD103+ CD4","CD103+ CD8","CD103- CD4","CD4 Tcm","KLRG1+ CD8","Treg"]; MARKERS=["FOXP3","IL2RA","CTLA4","TIGIT","IKZF2","ICOS"]; MIN_CELLS=20

def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()
def write_csv(df,path): path.parent.mkdir(parents=True,exist_ok=True); df.to_csv(path,index=False,encoding="utf-8")
def write_matrix(gene,vectors,path):
    if not vectors: raise RuntimeError(f"No pseudobulk vectors for {path.name}")
    frame=gene.copy()
    for sample,vec in vectors.items():
        if len(vec)!=len(frame): raise RuntimeError(f"Vector length mismatch for {sample}")
        frame[sample]=vec.astype(np.int64,copy=False)
    frame.to_csv(path,index=False,compression="gzip",encoding="utf-8")
def aggregate(layer,obs_index,cells,group_col,analysis,gene):
    if cells.index.duplicated().any(): raise RuntimeError(f"Duplicated selected cells in {analysis}")
    positions=obs_index.get_indexer(cells.index)
    if (positions<0).any(): raise RuntimeError(f"Selected cells absent from h5ad in {analysis}")
    cells=cells.copy(); cells["_position"]=positions; vectors={}; rows=[]
    for group,sub in cells.groupby(group_col,sort=True):
        idx=np.sort(sub["_position"].to_numpy(dtype=np.int64)); selected=layer[idx,:]
        if not sparse.issparse(selected): selected=sparse.csr_matrix(selected)
        data=selected.data
        if data.size and (data.min()<0 or not np.allclose(data,np.rint(data))): raise RuntimeError(f"Raw-count invariant failed in {analysis}/{group}")
        vec=np.asarray(selected.sum(axis=0)).ravel().astype(np.int64); vectors[str(group)]=vec
        conditions=sorted(set(sub["condition"].astype(str)))
        if len(conditions)!=1: raise RuntimeError(f"Multiple conditions in {analysis}/{group}")
        rows.append(dict(analysis=analysis,sample_column=str(group),participant_id=str(group),condition=conditions[0],n_cells=int(len(sub)),library_size=int(vec.sum()),detected_genes=int((vec>0).sum()),aggregation_verified="YES" if int(vec.sum())==int(selected.sum()) else "NO",matrix_gene_count=int(len(gene))))
    return vectors,rows

def main():
    started=time.time(); WORK.mkdir(parents=True,exist_ok=True)
    for p in (MAP,CROSSWALK,COVERAGE,SELECTION,RELIABILITY,H5AD):
        if not p.is_file(): raise FileNotFoundError(p)
    label_map=pd.read_csv(MAP,dtype=str,keep_default_na=False); crosswalk=pd.read_csv(CROSSWALK,dtype=str,keep_default_na=False); coverage=pd.read_csv(COVERAGE,dtype=str,keep_default_na=False); selection=pd.read_csv(SELECTION,dtype=str,keep_default_na=False); reliability=pd.read_csv(RELIABILITY,dtype=str,keep_default_na=False)
    if len(label_map)!=15670 or label_map["frozen_cell_id"].duplicated().any(): raise RuntimeError("Phase5A verified map checkpoint failed")
    verified_t=label_map[label_map["is_T_related_subtype"].eq("YES")].copy(); nk_ilc=label_map[label_map["is_NK_or_ILC_label"].eq("YES")].copy()
    if len(verified_t)!=15233 or len(nk_ilc)!=437: raise RuntimeError(f"Phase5A T/NK-ILC checkpoint failed: T={len(verified_t)}, NK/ILC={len(nk_ilc)}")
    rc=label_map[["frozen_participant_id","disease_status"]].drop_duplicates(); rc["condition"]=rc["disease_status"].map({"Control":"Healthy","CD":"ActiveCrohn"})
    if rc["condition"].isna().any() or rc["frozen_participant_id"].duplicated().any(): raise RuntimeError("Recovered participant-condition mapping is not unique")
    recovered=sorted(rc["frozen_participant_id"].unique()); cond_lookup=rc.set_index("frozen_participant_id")["condition"].to_dict()
    if len(recovered)!=19 or rc["condition"].value_counts().to_dict()!={"Healthy":13,"ActiveCrohn":6}: raise RuntimeError("Recovered cohort checkpoint failed")
    a=ad.read_h5ad(H5AD,backed="r")
    if "counts" not in a.layers: raise RuntimeError("h5ad lacks layers[counts]")
    layer=a.layers["counts"]; probe=layer[:100,:]
    if not sparse.issparse(probe): probe=sparse.csr_matrix(probe)
    if probe.data.size and (probe.data.min()<0 or not np.allclose(probe.data,np.rint(probe.data))): raise RuntimeError("counts layer is not raw UMI data")
    obs=a.obs.copy(); obs_index=a.obs_names.copy(); gene=pd.DataFrame({"unique_gene_id":a.var["gene_id"].astype(str).to_numpy(),"gene_symbol":a.var["gene_symbol"].astype(str).to_numpy()})
    if obs_index.duplicated().any() or gene["unique_gene_id"].duplicated().any(): raise RuntimeError("Duplicated h5ad identity")
    mask=obs["primary_contrast_eligible"].astype(str).eq("YES_SELECTED") & obs["region"].astype(str).eq("TI") & obs["human_level1"].astype(str).eq("T")
    full_t=obs.loc[mask,["donor","disease_group"]].copy(); full_t["condition"]=full_t["disease_group"].astype(str).map({"Healthy":"Healthy","Crohn":"ActiveCrohn"})
    if len(full_t)!=21378 or full_t["donor"].nunique()!=27: raise RuntimeError("Frozen broad-T checkpoint failed")
    broad19=full_t[full_t["donor"].astype(str).isin(recovered)].copy(); broad19["participant_id"]=broad19["donor"].astype(str); broad19.index=broad19.index.astype(str)
    vt=verified_t.set_index("frozen_cell_id").copy(); vt["participant_id"]=vt["frozen_participant_id"].astype(str); vt["condition"]=vt["participant_id"].map(cond_lookup)
    all_qc=[]
    vectors,rows=aggregate(layer,obs_index,broad19,"participant_id","Broad19_AllFrozenT",gene); write_matrix(gene,vectors,WORK/"Broad19_AllFrozenT_counts.csv.gz"); all_qc.extend(rows)
    vectors,rows=aggregate(layer,obs_index,vt,"participant_id","Broad19_VerifiedAuthorT",gene); write_matrix(gene,vectors,WORK/"Broad19_VerifiedAuthorT_counts.csv.gz"); all_qc.extend(rows)
    eligibility=[]; contributions=[]
    for subtype in SUBTYPES:
        sub=vt[vt["author_cluster_label2"].eq(subtype)].copy(); per=sub.groupby(["participant_id","condition"]).size().rename("n_cells").reset_index(); eligible=per[per["n_cells"]>=MIN_CELLS]
        nh=int((eligible["condition"]=="Healthy").sum()); nc=int((eligible["condition"]=="ActiveCrohn").sum()); analyzable=nh>=5 and nc>=5
        largest=per.sort_values(["n_cells","participant_id"],ascending=[False,True]).iloc[0]; eligible_ids=set(eligible["participant_id"])
        eligibility.append(dict(subtype=subtype,total_verified_cells=int(len(sub)),participants_with_any_cells=int(per["participant_id"].nunique()),healthy_participants_ge20=nh,crohn_participants_ge20=nc,eligible_participants_total=int(len(eligible)),eligible_cells_total=int(sub[sub["participant_id"].isin(eligible_ids)].shape[0]),largest_participant=str(largest["participant_id"]),largest_participant_cells=int(largest["n_cells"]),largest_participant_fraction=float(largest["n_cells"]/len(sub)),raw_count_status="PASS",design_estimable_preflight="YES" if analyzable else "NO",eligibility_status="ANALYZABLE" if analyzable else "NOT_ESTIMABLE",eligibility_reason=">=20 cells and >=5 independent participants per group" if analyzable else "DONOR_COVERAGE_BELOW_FROZEN_GATE"))
        for row in per.itertuples(index=False): contributions.append(dict(subtype=subtype,participant_id=row.participant_id,condition=row.condition,cell_count=int(row.n_cells),eligible_ge20="YES" if int(row.n_cells)>=20 else "NO"))
        if analyzable:
            selected=sub[sub["participant_id"].isin(eligible_ids)].copy(); vectors,rows=aggregate(layer,obs_index,selected,"participant_id",f"Subtype::{subtype}",gene)
            slug=subtype.replace("+","PLUS").replace("-","MINUS").replace(" ","_"); write_matrix(gene,vectors,WORK/f"Subtype_{slug}_counts.csv.gz"); all_qc.extend(rows)
    symbol_to_idx={}
    for i,symbol in enumerate(gene["gene_symbol"]): symbol_to_idx.setdefault(str(symbol),[]).append(i)
    positions=obs_index.get_indexer(vt.index); is_treg=vt["author_cluster_label2"].eq("Treg").to_numpy(); participants=vt["participant_id"].to_numpy(); totals=obs.loc[vt.index,"total_counts"].astype(float).to_numpy(); marker_rows=[]
    for marker in MARKERS:
        hits=symbol_to_idx.get(marker,[]); base={"marker":marker,"unique_gene_mapping":"YES" if len(hits)==1 else "NO","gene_mapping_count":len(hits)}
        if len(hits)!=1:
            marker_rows.append(base|{"testability":"NOT_ASSESSABLE_NONUNIQUE_GENE","treg_detection_fraction":np.nan,"other_T_detection_fraction":np.nan,"treg_mean_log1p_cpm":np.nan,"other_T_mean_log1p_cpm":np.nan,"treg_participants_with_detection":np.nan,"treg_participants_total":np.nan,"direction_consistent_with_regulatory_annotation":"NOT_ASSESSABLE"}); continue
        values=layer[positions,hits[0]]; values=values.toarray().ravel() if sparse.issparse(values) else np.asarray(values).ravel(); detected=values>0
        norm=np.log1p(np.divide(values*1e4,totals,out=np.zeros_like(values,dtype=float),where=totals>0)); tp=sorted(set(participants[is_treg])); nd=sum(bool(detected[is_treg & (participants==p)].any()) for p in tp)
        td,od=float(detected[is_treg].mean()),float(detected[~is_treg].mean()); te,oe=float(norm[is_treg].mean()),float(norm[~is_treg].mean())
        marker_rows.append(base|{"testability":"TESTABLE","treg_detection_fraction":td,"other_T_detection_fraction":od,"treg_mean_log1p_cpm":te,"other_T_mean_log1p_cpm":oe,"treg_participants_with_detection":int(nd),"treg_participants_total":int(len(tp)),"direction_consistent_with_regulatory_annotation":"YES" if td>od and te>oe else "NO"})
    write_csv(pd.DataFrame(all_qc),WORK/"pseudobulk_sample_qc.csv"); write_csv(pd.DataFrame(eligibility),WORK/"subtype_eligibility_preflight.csv"); write_csv(pd.DataFrame(contributions),WORK/"subtype_participant_contributions.csv"); write_csv(pd.DataFrame(marker_rows),WORK/"treg_marker_audit_raw.csv"); write_csv(selection,WORK/"phase5a_cohort_selection_copy.csv"); write_csv(reliability,WORK/"phase5a_subtype_reliability_copy.csv")
    audit={"phase5a_map":{"path":str(MAP),"sha256":sha256(MAP),"rows":len(label_map)},"phase5a_crosswalk":{"path":str(CROSSWALK),"sha256":sha256(CROSSWALK),"rows":len(crosswalk)},"phase5a_coverage":{"path":str(COVERAGE),"sha256":sha256(COVERAGE),"rows":len(coverage)},"phase5a_reliability":{"path":str(RELIABILITY),"sha256":sha256(RELIABILITY),"rows":len(reliability)},"h5ad":{"path":str(H5AD),"sha256_reused_from_phase5a":"13a5759fecbfa1fb4a773a2b98494bef02669b051238f496cc7c1d2e829d5eaf","shape":list(a.shape),"raw_layer":"layers[counts]"},"checkpoints":{"full_T_cells":len(full_t),"full_participants":int(full_t["donor"].nunique()),"recovered_participants":len(recovered),"verified_T_cells":len(verified_t),"excluded_NK_ILC":len(nk_ilc)},"runtime_seconds":round(time.time()-started,3)}
    (WORK/"preflight_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8"); a.file.close(); print(json.dumps(audit["checkpoints"],indent=2))
if __name__=="__main__": main()

