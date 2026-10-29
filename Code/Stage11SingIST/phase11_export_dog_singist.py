from __future__ import annotations

import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
SOURCE = ROOT / "09_cell_ontology" / "dog" / "dog_GSE254005_ontology.h5ad"
OUT = ROOT / "11_singIST" / "dog_model"
LOG = ROOT / "logs" / "phase11_dog_object.log"
OUT.mkdir(parents=True, exist_ok=True)

lineages = ["Epithelial", "T", "Myeloid"]
adata = ad.read_h5ad(SOURCE, backed="r")
obs_all = adata.obs
mask = obs_all["harmonized_level1"].astype(str).isin(lineages).to_numpy()
idx = np.flatnonzero(mask)

counts = adata.layers["counts"][idx, :]
if not sparse.issparse(counts):
    counts = sparse.csr_matrix(counts)
else:
    counts = counts.tocsr()
if counts.data.size and not np.allclose(counts.data, np.rint(counts.data), atol=0, rtol=0):
    raise RuntimeError("Dog frozen counts layer is not integer-valued")
counts.data = np.rint(counts.data).astype(np.int32, copy=False)
counts.eliminate_zeros()

obs = obs_all.iloc[idx].copy()
obs_out = pd.DataFrame(
    {
        "cell_id": obs.index.astype(str),
        "sample_id": obs["sample_id"].astype(str).to_numpy(),
        "donor_full": obs["donor_id"].astype(str).to_numpy(),
        "donor": obs["sample_prefix"].astype(str).to_numpy(),
        "class": obs["condition"].astype(str).to_numpy(),
        "celltype_cluster": obs["harmonized_level1"].astype(str).to_numpy(),
        "cl_id": obs["harmonized_cl_id"].astype(str).to_numpy(),
        "mapping_confidence": obs["mapping_confidence"].astype(str).to_numpy(),
    }
)

var = adata.var.copy()
var_out = pd.DataFrame(
    {
        "dog_unique_var_name": var.index.astype(str),
        "dog_feature_id": var["feature_id"].astype(str).to_numpy(),
        "dog_original_symbol": var["gene_symbol"].astype(str).to_numpy(),
        "feature_type": var["feature_type"].astype(str).to_numpy(),
    }
)

if counts.shape != (len(obs_out), len(var_out)):
    raise RuntimeError(f"Shape mismatch: counts={counts.shape}")
if obs_out["donor"].nunique() != 7:
    raise RuntimeError("Dog donor count is not 7")
if set(obs_out["class"]) != {"Healthy", "CIE"}:
    raise RuntimeError("Dog class labels are not Healthy/CIE")
if set(obs_out["celltype_cluster"]) != set(lineages):
    raise RuntimeError("Dog Level1 scope mismatch")

mtx_path = OUT / "DOG_SINGIST_COUNTS_GENES_X_CELLS.mtx.gz"
with gzip.open(mtx_path, "wb", compresslevel=6) as handle:
    mmwrite(handle, counts.T.tocoo(), field="integer", symmetry="general")

obs_path = OUT / "DOG_SINGIST_CELL_METADATA.csv.gz"
obs_out.to_csv(obs_path, index=False, compression={"method": "gzip", "compresslevel": 6})
var_path = OUT / "DOG_SINGIST_GENE_METADATA.csv"
var_out.to_csv(var_path, index=False)

audit = (
    obs_out.groupby(["donor", "donor_full", "class", "celltype_cluster"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
    .sort_values(["celltype_cluster", "class", "donor"], kind="stable")
)
audit["source_object"] = str(SOURCE)
audit["counts_layer"] = "counts"
audit["annotation_field"] = "harmonized_level1"
audit["level_scope"] = "PRIMARY_LEVEL1_ONLY"
audit.to_csv(OUT / "DOG_SINGIST_OBJECT_AUDIT.csv", index=False)

with LOG.open("w", encoding="utf-8", newline="\n") as handle:
    handle.write("PHASE11 DOG OBJECT EXPORT\n")
    handle.write(f"source={SOURCE}\n")
    handle.write(f"selected_cells={counts.shape[0]}\n")
    handle.write(f"genes={counts.shape[1]}\n")
    handle.write(f"nonzero_counts={counts.nnz}\n")
    handle.write(f"donors={obs_out['donor'].nunique()}\n")
    handle.write(f"healthy_donors={obs_out.loc[obs_out['class'].eq('Healthy'), 'donor'].nunique()}\n")
    handle.write(f"cie_donors={obs_out.loc[obs_out['class'].eq('CIE'), 'donor'].nunique()}\n")
    for lineage in lineages:
        handle.write(f"{lineage}_cells={int(obs_out['celltype_cluster'].eq(lineage).sum())}\n")
    handle.write("qc_reannotation=NO\n")
    handle.write("batch_integration=NO\n")
    handle.write("PHASE11_DOG_OBJECT_EXPORT_COMPLETE\n")

adata.file.close()
print(f"DOG_EXPORT_COMPLETE cells={counts.shape[0]} genes={counts.shape[1]} nnz={counts.nnz}")
