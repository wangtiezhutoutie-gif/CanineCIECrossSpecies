from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


SOURCES = {
    "dog": (Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"]) / "09_cell_ontology" / "dog" / "dog_GSE254005_ontology.h5ad"),
    "human": (Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"]) / "09_cell_ontology" / "human" / "human_GSE266616_TI_ontology.h5ad"),
}

PRIMARY = ["Epithelial", "T", "Myeloid"]


for species, path in SOURCES.items():
    a = ad.read_h5ad(path, backed="r")
    donor_col = "donor_id" if species == "dog" else "donor"
    condition_col = "condition" if species == "dog" else "disease_group"
    obs = a.obs[["harmonized_level1", donor_col, condition_col]].copy()
    obs.columns = ["lineage", "donor", "condition"]
    obs["species"] = species
    obs = obs.loc[obs["lineage"].isin(PRIMARY)]
    print(f"\n{species.upper()} SHAPE {a.shape}; PRIMARY {len(obs)}")
    print(obs.groupby(["lineage", "donor", "condition"], observed=True).size().to_string())
    print("LINEAGE TOTALS")
    print(obs.groupby("lineage", observed=True).size().to_string())
    var = a.var.copy()
    stable_col = "feature_id" if species == "dog" else "gene_id"
    ids = var[stable_col].astype(str).str.replace(r"\.\d+$", "", regex=True)
    print(
        "GENES",
        len(var),
        "STABLE_NONEMPTY",
        int(ids.ne("").sum()),
        "STABLE_UNIQUE",
        int(ids.nunique()),
        "ENSEMBL_PREFIX",
        int(ids.str.startswith("ENS").sum()),
    )
    print(var[[stable_col, "gene_symbol"]].head(10).to_string())
    a.file.close()
