from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
R=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"]);O=R/"WORKING"/"24_EPITHELIAL_MYELOID_CROSS_SPECIES";W=O/"work"
H=ad.read_h5ad(R/"07_human_crohn"/"GSE266616"/"h5ad"/"human_GSE266616_TI_annotated.h5ad",backed="r")
obs=H.obs.copy();obs["cell_id"]=obs.index.astype(str)
for c in ["donor","disease_group","human_level1","primary_contrast_eligible","biologic_treatment"]:obs[c]=obs[c].astype(str)
obs=obs[obs.primary_contrast_eligible=="YES_SELECTED"].copy()
uid=H.var_names.astype(str);sym=H.var["gene_symbol"].astype(str).values if "gene_symbol" in H.var else uid
X=H.layers["counts"] if "counts" in H.layers else H.X
def slug(s):return "".join(c.lower() if c.isalnum() else "_" for c in s).strip("_")
def aggregate(name,ids,meta_extra=None):
 idx=obs.index.get_indexer(pd.Index(ids));idx=idx[idx>=0];z=obs.loc[obs.index[idx]].copy();cols=[];mats=[];meta=[]
 for donor,g in z.groupby("donor"):
  pos=obs.index.get_indexer(g.index);a=np.asarray(X[pos,:].sum(axis=0)).ravel().astype(np.int64);cond="Healthy" if g.disease_group.iloc[0]=="Healthy" else "ActiveCrohn";cols.append(donor);mats.append(a);meta.append([donor,donor,cond,g.biologic_treatment.iloc[0],len(g),int(a.sum()),name])
 q=pd.DataFrame(np.column_stack(mats),columns=cols);q.insert(0,"gene_symbol",sym);q.insert(0,"unique_gene_id",uid);q.to_csv(W/f"human_{slug(name)}_counts.tsv",sep="\t",index=False);pd.DataFrame(meta,columns=["sample_column","donor","condition","biologic_treatment","n_cells","total_UMI","analysis"]).to_csv(W/f"human_{slug(name)}_meta.tsv",sep="\t",index=False)
for lin in ["Epithelial","Myeloid"]:
 broad=obs[obs.human_level1==lin];mapped=pd.read_csv(O/f"03_HUMAN_ANNOTATION/{'01_HUMAN_EPITHELIAL_LABEL_MAP.csv' if lin=='Epithelial' else '02_HUMAN_MYELOID_LABEL_MAP.csv'}")
 recovered=set(mapped.donor.astype(str));aggregate(f"{lin}__B_RECOVERED_PARTICIPANTS_ALL_BROAD",broad[broad.donor.isin(recovered)].index)
 aggregate(f"{lin}__C_VERIFIED_AUTHOR_LABELLED_LINEAGE",mapped.h5ad_cell_id.astype(str))
 elig=pd.read_csv(O/"03_HUMAN_ANNOTATION"/"04_HUMAN_SUBTYPE_ELIGIBILITY.csv")
 for _,e in elig[(elig.lineage==lin)&(elig.formal_pseudobulk_eligible=="YES")].iterrows():aggregate(f"{lin}__Subtype__{e.author_subtype}",mapped[mapped.author_subtype==e.author_subtype].h5ad_cell_id.astype(str))
print("HUMAN_PSEUDOBULK_READY")