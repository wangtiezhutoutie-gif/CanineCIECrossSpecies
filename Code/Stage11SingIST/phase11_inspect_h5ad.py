from pathlib import Path
import anndata as ad
import pandas as pd

root = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
paths = [
    root / "09_cell_ontology" / "dog" / "dog_GSE254005_ontology.h5ad",
    root / "09_cell_ontology" / "human" / "human_GSE266616_TI_ontology.h5ad",
]
for path in paths:
    print(f"FILE={path}")
    x = ad.read_h5ad(path, backed="r")
    print("shape", x.shape)
    print("X", type(x.X), getattr(x.X, "dtype", None))
    print("layers", list(x.layers.keys()))
    print("obs_columns", list(x.obs.columns))
    print("var_columns", list(x.var.columns))
    print("obs_dtypes")
    print(x.obs.dtypes.to_string())
    print("obs_head")
    print(x.obs.head(3).to_string())
    for col in x.obs.columns:
        if any(token in col.lower() for token in ("donor", "condition", "level1", "celltype", "ontology")):
            vals = x.obs[col].astype(str).value_counts(dropna=False)
            print(f"VALUE_COUNTS={col}")
            print(vals.head(30).to_string())
    print("var_head")
    print(x.var.head(5).to_string())
    print()
    x.file.close()
