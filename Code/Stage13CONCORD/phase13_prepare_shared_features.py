#!/usr/bin/env python
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P13 = ROOT / "13_concord"
MAPPING = ROOT / "10_cross_species_baseline" / "mapping" / "DOG_HUMAN_ORTHOLOG_MASTER.csv"
SAMPLING = ROOT / "12_saturn" / "sampling" / "SATURN_PRIMARY_CELL_SAMPLING.csv"
DOG_INPUT = ROOT / "12_saturn" / "input" / "dog_saturn_input.h5ad"
HUMAN_INPUT = ROOT / "12_saturn" / "input" / "human_saturn_input.h5ad"
CONDITION = ROOT / "12_saturn" / "input" / "CELL_CONDITION_SIDECAR.csv"
OUT = P13 / "input"
AUDIT = P13 / "audit"
LOG = ROOT / "logs" / "phase13_feature_mapping.log"
SEED = 20260827
N_HVG = 10_000
ALLOWED_OBS = ["cell_id", "species", "donor_id", "harmonized_level1", "dataset_domain"]


def log(msg: str) -> None:
    stamp = datetime.now(timezone.utc).isoformat()
    line = f"[{stamp}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("", encoding="utf-8")

    log("Reading frozen Phase10 ortholog mapping and Phase12 cell cohort")
    mapping = pd.read_csv(MAPPING, low_memory=False)
    primary = mapping["primary_one2one"].astype(str).str.lower().eq("true")
    pairs = mapping.loc[primary].copy()
    require(len(pairs) == 14_945, f"Expected 14945 primary one-to-one pairs, found {len(pairs)}")
    dog_col = "dog_original_feature_id"
    human_col = "human_ortholog_release104_id"
    require(not pairs[dog_col].isna().any(), "Missing dog IDs in primary mapping")
    require(not pairs[human_col].isna().any(), "Missing human IDs in primary mapping")
    require(not pairs[dog_col].duplicated().any(), "Duplicate dog IDs in primary mapping")
    require(not pairs[human_col].duplicated().any(), "Duplicate human IDs in primary mapping")
    require((pairs["mapping_status"] == "ONE2ONE_PRIMARY").all(), "Primary mapping status changed")

    sampling = pd.read_csv(SAMPLING, low_memory=False)
    selected = sampling.loc[sampling["selected"].eq("YES")].copy()
    require(len(selected) == 28_046, f"Expected 28046 selected cells, found {len(selected)}")
    require(not selected["cell_id"].duplicated().any(), "Duplicate selected cell IDs")
    expected_counts = {
        ("dog", "Epithelial"): 6000,
        ("dog", "T"): 6000,
        ("dog", "Myeloid"): 2023,
        ("human", "Epithelial"): 6000,
        ("human", "T"): 6000,
        ("human", "Myeloid"): 2023,
    }
    observed_counts = selected.groupby(["species", "lineage"]).size().to_dict()
    require(observed_counts == expected_counts, f"Frozen sampling counts changed: {observed_counts}")

    log("Loading Phase12 raw-count input objects")
    dog = ad.read_h5ad(DOG_INPUT)
    human = ad.read_h5ad(HUMAN_INPUT)
    require(dog.n_obs == 14_023 and human.n_obs == 14_023, "Species cell totals differ from Phase12")
    require(dog.n_vars == 20_257 and human.n_vars == 36_601, "Unexpected source feature dimensions")
    for name, a in [("dog", dog), ("human", human)]:
        require(a.obs_names.is_unique, f"{name} obs_names are not unique")
        require(a.var_names.is_unique, f"{name} var_names are not unique")
        require(a.obs["cell_id"].astype(str).tolist() == a.obs_names.astype(str).tolist(), f"{name} cell_id != obs_names")
        require(sparse.issparse(a.X), f"{name} counts are not sparse")
        require(np.isfinite(a.X.data).all(), f"{name} counts contain non-finite values")
        require((a.X.data >= 0).all(), f"{name} counts contain negative values")
        require(np.allclose(a.X.data, np.rint(a.X.data)), f"{name} X is not integer raw counts")

    selected_ids = set(selected["cell_id"].astype(str))
    object_ids = set(dog.obs_names.astype(str)) | set(human.obs_names.astype(str))
    require(selected_ids == object_ids, f"Phase12 cell ID match is not 100%: selected={len(selected_ids)} objects={len(object_ids)} overlap={len(selected_ids & object_ids)}")

    log("Resolving the frozen 1:1 coordinate system without disease-effect fields")
    dog_present = pairs[dog_col].isin(dog.var_names)
    human_present = pairs[human_col].isin(human.var_names)
    pairs["dog_id_present_in_input"] = dog_present
    pairs["human_id_present_in_input"] = human_present
    both_id = dog_present & human_present

    dog_index = dog.var_names.get_indexer(pairs.loc[both_id, dog_col])
    human_index = human.var_names.get_indexer(pairs.loc[both_id, human_col])
    require((dog_index >= 0).all() and (human_index >= 0).all(), "Internal feature index resolution failed")
    dog_nnz = np.asarray(dog.X[:, dog_index].getnnz(axis=0)).ravel()
    human_nnz = np.asarray(human.X[:, human_index].getnnz(axis=0)).ravel()
    pairs["dog_detected_cells"] = 0
    pairs["human_detected_cells"] = 0
    pairs.loc[both_id, "dog_detected_cells"] = dog_nnz
    pairs.loc[both_id, "human_detected_cells"] = human_nnz
    pairs["eligible_shared_expression"] = both_id & pairs["dog_detected_cells"].gt(0) & pairs["human_detected_cells"].gt(0)
    eligible = pairs.loc[pairs["eligible_shared_expression"]].copy()
    require(len(eligible) >= N_HVG, f"Only {len(eligible)} eligible shared expressed genes")
    require(not eligible[dog_col].duplicated().any(), "Eligible dog feature collision")
    require(not eligible[human_col].duplicated().any(), "Eligible human feature collision")
    log(f"Frozen pairs=14945; IDs in both inputs={both_id.sum()}; expressed in both species={len(eligible)}")

    dog_idx = dog.var_names.get_indexer(eligible[dog_col])
    human_idx = human.var_names.get_indexer(eligible[human_col])
    dog_x = dog.X[:, dog_idx].tocsr().astype(np.float32)
    human_x = human.X[:, human_idx].tocsr().astype(np.float32)
    combined_x = sparse.vstack([dog_x, human_x], format="csr")

    dog_obs = dog.obs[["cell_id", "species", "donor_id", "harmonized_level1"]].copy()
    human_obs = human.obs[["cell_id", "species", "donor_id", "harmonized_level1"]].copy()
    obs = pd.concat([dog_obs, human_obs], axis=0)
    obs["species"] = obs["species"].astype(str)
    obs["dataset_domain"] = obs["species"]
    obs = obs[ALLOWED_OBS]
    require(list(obs.columns) == ALLOWED_OBS, "Training obs schema mismatch")
    require(set(obs["species"]) == {"dog", "human"}, "Species labels changed")
    require(set(obs["harmonized_level1"]) == {"Epithelial", "T", "Myeloid"}, "Unauthorized Level1 labels")
    forbidden_columns = {"condition", "state", "treatment", "treatment_class", "region", "healthy", "cie", "crohn", "phase10_status", "phase11_status", "saturn_status"}
    require(not (set(map(str.lower, obs.columns)) & forbidden_columns), "Condition leakage in obs columns")

    var = pd.DataFrame(index=pd.Index(eligible[human_col].astype(str), name="human_ortholog_gene_id"))
    var["dog_feature_id"] = eligible[dog_col].astype(str).to_numpy()
    var["human_feature_id"] = eligible[human_col].astype(str).to_numpy()
    var["dog_symbol"] = eligible["dog_original_symbol"].astype(str).to_numpy()
    var["human_symbol"] = eligible["human_ortholog_symbol"].astype(str).to_numpy()
    var["homology_type"] = eligible["homology_type_release104"].astype(str).to_numpy()
    var["goc_score"] = eligible["goc_score_if_available"].to_numpy()
    var["mapping_confidence"] = eligible["mapping_confidence"].astype(str).to_numpy()

    combined = ad.AnnData(X=combined_x, obs=obs, var=var)
    combined.obs_names = combined.obs["cell_id"].astype(str)
    combined.obs_names.name = "cell_id"
    require(combined.obs_names.is_unique and combined.var_names.is_unique, "Combined object identities are not unique")

    log(f"Running condition-blind combined HVG selection on {combined.n_obs} cells x {combined.n_vars} eligible genes")
    sc.pp.highly_variable_genes(
        combined,
        n_top_genes=N_HVG,
        flavor="seurat_v3",
        check_values=True,
        span=0.3,
    )
    require(int(combined.var["highly_variable"].sum()) == N_HVG, "HVG count is not 10000")
    selected_hvg = combined.var["highly_variable"].to_numpy()
    hvg = combined[:, selected_hvg].copy()
    hvg.layers["counts"] = hvg.X.copy()

    mapping_out = pairs.rename(columns={
        dog_col: "dog_feature_id",
        human_col: "human_feature_id",
        "dog_original_symbol": "dog_symbol",
        "human_ortholog_symbol": "human_symbol",
        "homology_type_release104": "homology_type",
        "goc_score_if_available": "goc_score",
    })
    keep_mapping = [
        "dog_feature_id", "human_feature_id", "dog_symbol", "human_symbol",
        "homology_type", "goc_score", "mapping_confidence", "orthology_confidence",
        "dog_id_present_in_input", "human_id_present_in_input",
        "dog_detected_cells", "human_detected_cells", "eligible_shared_expression",
    ]
    selected_ids_hvg = set(hvg.var_names.astype(str))
    mapping_out["selected_hvg"] = mapping_out["human_feature_id"].astype(str).isin(selected_ids_hvg)
    mapping_out[keep_mapping + ["selected_hvg"]].to_csv(OUT / "CONCORD_SHARED_FEATURE_MAPPING.csv", index=False)

    log("Computing species-bias audit and normalizing/log-transforming training input")
    bias = pd.DataFrame({
        "human_feature_id": hvg.var_names.astype(str),
        "dog_feature_id": hvg.var["dog_feature_id"].astype(str).to_numpy(),
        "dog_symbol": hvg.var["dog_symbol"].astype(str).to_numpy(),
        "human_symbol": hvg.var["human_symbol"].astype(str).to_numpy(),
        "hvg_rank": hvg.var["highly_variable_rank"].to_numpy(),
        "dog_mean_raw": np.asarray(hvg.X[:dog.n_obs].mean(axis=0)).ravel(),
        "human_mean_raw": np.asarray(hvg.X[dog.n_obs:].mean(axis=0)).ravel(),
        "dog_detection_fraction": np.asarray((hvg.X[:dog.n_obs] > 0).mean(axis=0)).ravel(),
        "human_detection_fraction": np.asarray((hvg.X[dog.n_obs:] > 0).mean(axis=0)).ravel(),
    })
    bias["mean_ratio_dog_over_human"] = (bias["dog_mean_raw"] + 1e-8) / (bias["human_mean_raw"] + 1e-8)
    bias.to_csv(OUT / "CONCORD_HVG_SPECIES_BIAS_AUDIT.csv", index=False)

    sc.pp.normalize_total(hvg, target_sum=1e4)
    sc.pp.log1p(hvg)
    require(np.isfinite(hvg.X.data).all(), "Normalized training matrix contains non-finite values")
    require(list(hvg.obs.columns) == ALLOWED_OBS, "Training obs acquired unauthorized columns")
    hvg.uns["phase13_preprocessing"] = {
        "orthology": "Phase10 frozen Ensembl release104 one-to-one",
        "hvg": "combined balanced 28046 cells; seurat_v3; condition-blind",
        "n_frozen_pairs": 14945,
        "n_eligible_shared_expressed": int(len(eligible)),
        "n_hvg": N_HVG,
        "normalize_total_target_sum": 10000,
        "log1p": True,
        "selection_seed": SEED,
    }
    h5ad_path = OUT / "concord_primary_input.h5ad"
    hvg.write_h5ad(h5ad_path, compression="gzip")

    sidecar = pd.read_csv(CONDITION, low_memory=False)
    require(not sidecar["cell_id"].duplicated().any(), "Condition sidecar cell IDs are duplicated")
    sidecar = sidecar.set_index("cell_id").loc[hvg.obs_names].reset_index()
    require(len(sidecar) == 28_046 and not sidecar.isna().all(axis=1).any(), "Condition sidecar join coverage failed")
    sidecar.to_csv(OUT / "CONCORD_CELL_CONDITION_SIDECAR.csv", index=False)

    audit_rows = [
        {"check": "phase12_cell_id_match", "value": "28046/28046", "status": "PASS"},
        {"check": "dog_cells", "value": dog.n_obs, "status": "PASS"},
        {"check": "human_cells", "value": human.n_obs, "status": "PASS"},
        {"check": "frozen_one2one_pairs", "value": len(pairs), "status": "PASS"},
        {"check": "ids_present_both_inputs", "value": int(both_id.sum()), "status": "PASS"},
        {"check": "eligible_expressed_both_species", "value": len(eligible), "status": "PASS"},
        {"check": "selected_hvg", "value": hvg.n_vars, "status": "PASS"},
        {"check": "mapping_collision", "value": 0, "status": "PASS"},
        {"check": "condition_columns_in_training_obs", "value": 0, "status": "PASS"},
        {"check": "level2_columns_in_training_obs", "value": 0, "status": "PASS"},
        {"check": "nonfinite_training_values", "value": 0, "status": "PASS"},
        {"check": "resource_fallback", "value": "NO", "status": "PASS"},
    ]
    pd.DataFrame(audit_rows).to_csv(AUDIT / "CONCORD_INPUT_CONTRACT_AUDIT.csv", index=False)
    with (OUT / "CONCORD_INPUT_CHECKSUMS.json").open("w", encoding="utf-8") as fh:
        json.dump({
            "concord_primary_input.h5ad": sha256(h5ad_path),
            "CONCORD_SHARED_FEATURE_MAPPING.csv": sha256(OUT / "CONCORD_SHARED_FEATURE_MAPPING.csv"),
            "CONCORD_CELL_CONDITION_SIDECAR.csv": sha256(OUT / "CONCORD_CELL_CONDITION_SIDECAR.csv"),
            "CONCORD_HVG_SPECIES_BIAS_AUDIT.csv": sha256(OUT / "CONCORD_HVG_SPECIES_BIAS_AUDIT.csv"),
        }, fh, indent=2)

    log(f"PASS: wrote {h5ad_path.name}, shape={hvg.shape}, condition leakage=0, resource fallback=NO")


if __name__ == "__main__":
    main()

