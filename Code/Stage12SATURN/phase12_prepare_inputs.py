"""Prepare condition-blind, donor-balanced SATURN inputs from frozen Phase9 H5ADs."""

from __future__ import annotations

import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
OUT = ROOT / "12_saturn"
SEED = 20260827
CAP = 6000
LINEAGES = ("Epithelial", "T", "Myeloid")
SOURCES = {
    "dog": ROOT / "09_cell_ontology" / "dog" / "dog_GSE254005_ontology.h5ad",
    "human": ROOT / "09_cell_ontology" / "human" / "human_GSE266616_TI_ontology.h5ad",
}


def donor_quotas(counts: pd.Series, target: int) -> dict[str, int]:
    remaining = int(target)
    available = {str(k): int(v) for k, v in counts.sort_index().items()}
    allocation = {k: 0 for k in available}
    while remaining:
        active = [k for k in sorted(available) if available[k] > allocation[k]]
        if not active:
            raise RuntimeError("Quota redistribution exhausted cells before target")
        base, rem = divmod(remaining, len(active))
        proposal = {k: base + (i < rem) for i, k in enumerate(active)}
        if base == 0:
            proposal = {k: int(i < rem) for i, k in enumerate(active)}
        used = 0
        for donor in active:
            room = available[donor] - allocation[donor]
            take = min(room, int(proposal[donor]))
            allocation[donor] += take
            used += take
        if used == 0:
            raise RuntimeError("Quota algorithm made no progress")
        remaining -= used
    return allocation


def matrix_audit(x: sparse.spmatrix) -> dict[str, float | int | str]:
    if not sparse.issparse(x):
        x = sparse.csr_matrix(x)
    data = x.data
    return {
        "n_cells": int(x.shape[0]),
        "n_genes": int(x.shape[1]),
        "nnz": int(x.nnz),
        "minimum": float(data.min()) if data.size else 0.0,
        "maximum": float(data.max()) if data.size else 0.0,
        "negative_n": int(np.sum(data < 0)),
        "noninteger_n": int(np.sum(~np.isclose(data, np.rint(data), atol=1e-8))),
        "finite": str(bool(np.isfinite(data).all())).upper(),
    }


def stable_seed(species: str, lineage: str) -> int:
    suffix = int(hashlib.sha256(f"{species}|{lineage}".encode()).hexdigest()[:8], 16)
    return (SEED + suffix) % (2**32)


def main() -> None:
    (OUT / "input").mkdir(parents=True, exist_ok=True)
    (OUT / "sampling").mkdir(parents=True, exist_ok=True)
    candidates: dict[str, pd.DataFrame] = {}
    target_by_lineage: dict[str, int] = {}

    for species, source in SOURCES.items():
        a = ad.read_h5ad(source, backed="r")
        donor_col = "donor_id" if species == "dog" else "donor"
        obs = pd.DataFrame(
            {
                "cell_id": a.obs_names.astype(str),
                "species": species,
                "lineage": a.obs["harmonized_level1"].astype(str).to_numpy(),
                "donor": a.obs[donor_col].astype(str).to_numpy(),
                "source_position": np.arange(a.n_obs, dtype=np.int64),
            }
        )
        keep = obs["lineage"].isin(LINEAGES)
        if species == "human":
            keep &= a.obs["primary_contrast_eligible"].astype(str).to_numpy() == "YES_SELECTED"
        candidates[species] = obs.loc[keep].reset_index(drop=True)
        a.file.close()

    for lineage in LINEAGES:
        target_by_lineage[lineage] = min(
            CAP,
            *[int((candidates[s]["lineage"] == lineage).sum()) for s in ("dog", "human")],
        )

    manifests = []
    selected_positions: dict[str, list[int]] = {"dog": [], "human": []}
    for species in ("dog", "human"):
        frame = candidates[species].copy()
        frame["selected"] = False
        frame["assigned_quota"] = 0
        frame["selection_seed"] = SEED
        frame["quota_rule"] = "iterative_equal_donor_quota_no_upsample"
        for lineage in LINEAGES:
            idx = frame.index[frame["lineage"] == lineage]
            group = frame.loc[idx]
            target = target_by_lineage[lineage]
            quotas = donor_quotas(group.groupby("donor").size(), target)
            rng = np.random.default_rng(stable_seed(species, lineage))
            for donor, quota in quotas.items():
                donor_idx = group.index[group["donor"] == donor].to_numpy()
                chosen = rng.choice(donor_idx, size=quota, replace=False)
                frame.loc[donor_idx, "assigned_quota"] = quota
                frame.loc[chosen, "selected"] = True
        if int(frame["selected"].sum()) != sum(target_by_lineage.values()):
            raise RuntimeError(f"Selected total mismatch for {species}")
        selected_positions[species] = sorted(
            frame.loc[frame["selected"], "source_position"].astype(int).tolist()
        )
        manifests.append(frame)

    manifest = pd.concat(manifests, ignore_index=True)
    manifest["selected"] = manifest["selected"].map({True: "YES", False: "NO"})
    manifest.to_csv(OUT / "sampling" / "SATURN_PRIMARY_CELL_SAMPLING.csv", index=False)

    sampling_summary = (
        manifest.groupby(["species", "lineage"], observed=True)
        .agg(original_cells=("cell_id", "size"), selected_cells=("selected", lambda x: int((x == "YES").sum())), donor_coverage=("donor", "nunique"))
        .reset_index()
    )
    sampling_summary["target_cap"] = CAP
    sampling_summary["resource_fallback"] = "NO"
    sampling_summary.to_csv(OUT / "sampling" / "SATURN_PRIMARY_SAMPLING_SUMMARY.csv", index=False)

    sidecars = []
    matrix_rows = []
    allowed_obs = ["cell_id", "species", "harmonized_level1", "donor_id"]
    for species, source in SOURCES.items():
        a = ad.read_h5ad(source, backed="r")
        positions = np.asarray(selected_positions[species], dtype=np.int64)
        counts = a.layers["counts"][positions, :]
        counts = sparse.csr_matrix(counts, dtype=np.float32)
        source_obs = a.obs.iloc[positions].copy()
        cell_ids = a.obs_names[positions].astype(str)
        donor_col = "donor_id" if species == "dog" else "donor"
        train_obs = pd.DataFrame(
            {
                "cell_id": cell_ids,
                "species": species,
                "harmonized_level1": source_obs["harmonized_level1"].astype(str).to_numpy(),
                "donor_id": source_obs[donor_col].astype(str).to_numpy(),
            },
            index=cell_ids,
        )
        train_obs.index.name = "cell_id_index"
        stable_col = "feature_id" if species == "dog" else "gene_id"
        stable_ids = a.var[stable_col].astype(str).str.replace(r"\.\d+$", "", regex=True)
        if stable_ids.duplicated().any():
            raise RuntimeError(f"Duplicate stable IDs in {species}")
        train_var = pd.DataFrame(
            {"gene_symbol": a.var["gene_symbol"].astype(str).to_numpy()},
            index=stable_ids.to_numpy(),
        )
        train_var.index.name = "ensembl_gene_id"
        out_adata = ad.AnnData(X=counts, obs=train_obs[allowed_obs], var=train_var)
        out_adata.uns["phase12_input"] = {
            "source": str(source),
            "matrix": "source layers[counts]",
            "condition_blind": "YES",
            "sampling_seed": SEED,
            "lineages": list(LINEAGES),
            "human_primary_contrast_only": "YES" if species == "human" else "NOT_APPLICABLE",
        }
        out_path = OUT / "input" / f"{species}_saturn_input.h5ad"
        out_adata.write_h5ad(out_path, compression="gzip")
        matrix_rows.append({"species": species, "path": str(out_path), **matrix_audit(counts)})

        if species == "dog":
            side = pd.DataFrame(
                {
                    "cell_id": cell_ids,
                    "species": species,
                    "donor_id": source_obs["donor_id"].astype(str).to_numpy(),
                    "condition": source_obs["condition"].astype(str).to_numpy(),
                    "state": source_obs["condition"].astype(str).to_numpy(),
                    "treatment": "NOT_AVAILABLE",
                    "treatment_class": "NOT_AVAILABLE",
                    "region": "duodenum",
                    "lineage": source_obs["harmonized_level1"].astype(str).to_numpy(),
                    "primary_contrast_eligible": "YES_SELECTED",
                }
            )
        else:
            side = pd.DataFrame(
                {
                    "cell_id": cell_ids,
                    "species": species,
                    "donor_id": source_obs["donor"].astype(str).to_numpy(),
                    "condition": source_obs["disease_group"].astype(str).to_numpy(),
                    "state": source_obs["state"].astype(str).to_numpy(),
                    "treatment": source_obs["biologic_treatment"].astype(str).to_numpy(),
                    "treatment_class": source_obs["biologic_class"].astype(str).to_numpy(),
                    "region": source_obs["region"].astype(str).to_numpy(),
                    "lineage": source_obs["harmonized_level1"].astype(str).to_numpy(),
                    "primary_contrast_eligible": source_obs["primary_contrast_eligible"].astype(str).to_numpy(),
                }
            )
        sidecars.append(side)
        a.file.close()

    pd.concat(sidecars, ignore_index=True).to_csv(OUT / "input" / "CELL_CONDITION_SIDECAR.csv", index=False)
    pd.DataFrame(matrix_rows).to_csv(OUT / "audit" / "SATURN_RAW_COUNT_AUDIT.csv", index=False)

    forbidden = [
        "condition", "Healthy", "CIE", "Crohn", "Active Crohn", "treatment", "region",
        "Phase10 status", "Phase11 result",
    ]
    feature_rows = []
    for species in ("dog", "human"):
        for field in allowed_obs:
            feature_rows.append({"species": species, "field": field, "read_by_training": "YES", "forbidden": "NO"})
        for field in forbidden:
            feature_rows.append({"species": species, "field": field, "read_by_training": "NO", "forbidden": "YES"})
    pd.DataFrame(feature_rows).to_csv(OUT / "audit" / "SATURN_TRAINING_FEATURE_AUDIT.csv", index=False)
    print(sampling_summary.to_string(index=False))
    print(pd.DataFrame(matrix_rows).to_string(index=False))


if __name__ == "__main__":
    main()
