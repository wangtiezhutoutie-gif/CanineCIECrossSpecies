"""Build donor x lineage raw-count pseudobulk matrices for Phase 3."""

from __future__ import annotations

import hashlib
import logging
import os
import platform
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

CODE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_ROOT))
PROJECT_ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"] )

from config.paths import LOG_DIR  # noqa: E402
from config.phase3_config import (  # noqa: E402
    ANNOTATED_H5AD,
    GENE_MAPPING_MASTER,
    MIN_STATE_CELLS_PER_DONOR,
    PB_MATRICES_DIR,
    PB_METADATA_DIR,
    PHASE2_MANIFEST,
    PHASE2_REPORT,
    PRIMARY_LINEAGES,
    SECONDARY_LINEAGES,
    ensure_phase3_directories,
)
from Stage02DogAnnotation.phase2_utils import configure_logger, peak_ram_mb, validate_phase1_clean  # noqa: E402

SCRIPT_VERSION = "1.0.0"


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def verify_phase2_manifest() -> None:
    manifest = pd.read_csv(PHASE2_MANIFEST)
    missing: list[str] = []
    mismatch: list[str] = []
    for row in manifest.itertuples(index=False):
        path = PROJECT_ROOT / row.path_relative_to_project
        if not path.exists():
            missing.append(str(path))
        elif sha256(path).lower() != str(row.sha256).lower():
            mismatch.append(str(path))
    if missing or mismatch:
        raise RuntimeError(f"Phase2 manifest validation failed missing={missing} mismatch={mismatch}")


def atomic_csv(frame: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temp, index=False, compression=compression)
    os.replace(temp, path)


def donor_order(obs: pd.DataFrame) -> list[str]:
    donor_condition = obs[["donor_id", "condition"]].drop_duplicates().copy()
    donor_condition["condition_order"] = donor_condition["condition"].map({"Healthy": 0, "CIE": 1})
    return donor_condition.sort_values(["condition_order", "donor_id"])["donor_id"].astype(str).tolist()


def validate_input(adata: ad.AnnData, mapping: pd.DataFrame) -> None:
    validate_phase1_clean(adata)
    if adata.shape != (38_186, 20_257):
        raise RuntimeError(f"Unexpected annotated shape {adata.shape}")
    required_obs = {"donor_id", "sample_id", "condition", "level1", "level2", "annotation_cluster"}
    required_var = {"feature_id", "gene_symbol"}
    if required_obs - set(adata.obs):
        raise RuntimeError(f"Missing obs fields {required_obs - set(adata.obs)}")
    if required_var - set(adata.var):
        raise RuntimeError(f"Missing var fields {required_var - set(adata.var)}")
    if not adata.var_names.is_unique:
        raise RuntimeError("unique_var_name index is not unique")
    if adata.obs["donor_id"].nunique() != 7:
        raise RuntimeError("Expected seven donors")
    donor_conditions = adata.obs[["donor_id", "condition"]].drop_duplicates()
    if donor_conditions.groupby("condition")["donor_id"].nunique().to_dict() != {"CIE": 4, "Healthy": 3}:
        raise RuntimeError("Condition/donor counts differ from frozen design")
    counts = adata.layers["counts"]
    if not sparse.issparse(counts):
        raise RuntimeError("counts layer must remain sparse")
    values = counts.data
    if values.size and (values.min() < 0 or not np.allclose(values, np.rint(values))):
        raise RuntimeError("counts layer is not nonnegative integer-valued")
    if len(mapping) != adata.n_vars or mapping["feature_id"].duplicated().any():
        raise RuntimeError("Gene mapping master does not provide one stable row per feature_id")
    if not np.array_equal(mapping["feature_id"].astype(str).to_numpy(), adata.var["feature_id"].astype(str).to_numpy()):
        raise RuntimeError("Gene mapping row order differs from annotated h5ad")


def aggregate_lineage(
    adata: ad.AnnData,
    lineage: str,
    donors: list[str],
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    counts = adata.layers["counts"].tocsr()
    obs_donor = adata.obs["donor_id"].astype(str).to_numpy()
    obs_lineage = adata.obs["level1"].astype(str).to_numpy()
    columns: dict[str, np.ndarray] = {}
    metadata: list[dict[str, object]] = []
    for donor in donors:
        mask = (obs_donor == donor) & (obs_lineage == lineage)
        n_cells = int(mask.sum())
        summed = np.asarray(counts[mask].sum(axis=0)).ravel().astype(np.int64)
        columns[donor] = summed
        donor_obs = adata.obs.loc[mask]
        metadata.append(
            {
                "donor_id": donor,
                "sample_id": ";".join(sorted(donor_obs["sample_id"].astype(str).unique())),
                "condition": str(donor_obs["condition"].iloc[0]) if n_cells else "",
                "lineage": lineage,
                "analysis_tier": "Primary",
                "n_cells": n_cells,
                "library_size": int(summed.sum()),
                "n_detected_features": int((summed > 0).sum()),
                "eligible": n_cells >= MIN_STATE_CELLS_PER_DONOR,
                "exclusion_reason": "" if n_cells >= MIN_STATE_CELLS_PER_DONOR else f"n_cells<{MIN_STATE_CELLS_PER_DONOR}",
            }
        )
    matrix = pd.DataFrame(columns)
    matrix.insert(0, "feature_id", adata.var["feature_id"].astype(str).to_numpy())
    return matrix, metadata


def level1_eligibility(adata: ad.AnnData) -> pd.DataFrame:
    counts = pd.crosstab(adata.obs["level1"].astype(str), adata.obs["donor_id"].astype(str))
    donor_condition = adata.obs[["donor_id", "condition"]].drop_duplicates().set_index("donor_id")["condition"].astype(str)
    healthy = donor_condition[donor_condition == "Healthy"].index
    cie = donor_condition[donor_condition == "CIE"].index
    rows: list[dict[str, object]] = []
    for lineage in SECONDARY_LINEAGES:
        row = counts.loc[lineage] if lineage in counts.index else pd.Series(0, index=counts.columns)
        h_ge20 = int((row.reindex(healthy, fill_value=0) >= MIN_STATE_CELLS_PER_DONOR).sum())
        c_ge20 = int((row.reindex(cie, fill_value=0) >= MIN_STATE_CELLS_PER_DONOR).sum())
        eligible = h_ge20 == 3 and c_ge20 == 4
        reasons = []
        if h_ge20 < 3:
            reasons.append("Healthy donor(s) <20 cells")
        if c_ge20 < 4:
            reasons.append("CIE donor(s) <20 cells")
        rows.append(
            {
                "level1": lineage,
                "analysis_tier": "Secondary/Exploratory",
                "healthy_donors_present": int((row.reindex(healthy, fill_value=0) > 0).sum()),
                "cie_donors_present": int((row.reindex(cie, fill_value=0) > 0).sum()),
                "healthy_donors_ge20cells": h_ge20,
                "cie_donors_ge20cells": c_ge20,
                "total_cells": int(row.sum()),
                "eligible": eligible,
                "reason": "ELIGIBLE_SECONDARY" if eligible else "; ".join(reasons),
            }
        )
    return pd.DataFrame(rows)


def level2_eligibility(adata: ad.AnnData) -> pd.DataFrame:
    obs = adata.obs.copy()
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["condition"] = obs["condition"].astype(str)
    obs["level1"] = obs["level1"].astype(str)
    obs["level2"] = obs["level2"].astype(str)
    obs["annotation_cluster"] = obs["annotation_cluster"].astype(str)
    donor_condition = obs[["donor_id", "condition"]].drop_duplicates().set_index("donor_id")["condition"]
    healthy = donor_condition[donor_condition == "Healthy"].index
    cie = donor_condition[donor_condition == "CIE"].index
    cluster_dominance = obs["annotation_cluster"].value_counts()
    cluster_donor = pd.crosstab(obs["annotation_cluster"], obs["donor_id"])
    dominant_clusters = set((cluster_donor.max(axis=1) / cluster_donor.sum(axis=1)).loc[lambda x: x >= 0.80].index)
    rows: list[dict[str, object]] = []
    for (level1, level2), group in obs.groupby(["level1", "level2"], observed=True):
        counts = group["donor_id"].value_counts()
        h_counts = counts.reindex(healthy, fill_value=0)
        c_counts = counts.reindex(cie, fill_value=0)
        component_clusters = set(group["annotation_cluster"].unique())
        confounded_components = sorted(component_clusters & dominant_clusters)
        dominant_fraction = float(counts.max() / len(group))
        h_ge20 = int((h_counts >= MIN_STATE_CELLS_PER_DONOR).sum())
        c_ge20 = int((c_counts >= MIN_STATE_CELLS_PER_DONOR).sum())
        donor_confounded = dominant_fraction >= 0.80 or bool(confounded_components)
        eligible = h_ge20 >= 3 and c_ge20 >= 3 and not donor_confounded
        reasons: list[str] = []
        if h_ge20 < 3:
            reasons.append("Healthy donors >=20 cells <3")
        if c_ge20 < 3:
            reasons.append("CIE donors >=20 cells <3")
        if donor_confounded:
            reasons.append("DONOR_CONFOUNDED")
        rows.append(
            {
                "level2": level2,
                "level1": level1,
                "healthy_donors_present": int((h_counts > 0).sum()),
                "cie_donors_present": int((c_counts > 0).sum()),
                "healthy_donors_ge20cells": h_ge20,
                "cie_donors_ge20cells": c_ge20,
                "total_cells": len(group),
                "dominant_donor": str(counts.index[0]),
                "dominant_donor_fraction": dominant_fraction,
                "contains_donor_confounded_annotation_cluster": bool(confounded_components),
                "confounded_annotation_clusters": ";".join(confounded_components),
                "eligible": eligible,
                "status": "ELIGIBLE_EXPLORATORY" if eligible else "NOT_ELIGIBLE_FOR_STATE_PSEUDOBULK",
                "reason": "Pre-specified criteria passed" if eligible else "; ".join(reasons),
            }
        )
    return pd.DataFrame(rows).sort_values(["eligible", "level1", "level2"], ascending=[False, True, True])


def main() -> None:
    start = time.time()
    ensure_phase3_directories()
    logger = configure_logger(LOG_DIR / "phase3_build_pseudobulk.log", "phase3_build_pseudobulk")
    logger.info("START script_version=%s python=%s input=%s counts_source=layers[counts]", SCRIPT_VERSION, platform.python_version(), ANNOTATED_H5AD)
    try:
        if not PHASE2_REPORT.exists():
            raise RuntimeError(f"Missing Phase2 report: {PHASE2_REPORT}")
        verify_phase2_manifest()
        logger.info("PHASE2_MANIFEST status=pass rows=%d", len(pd.read_csv(PHASE2_MANIFEST)))
        adata = ad.read_h5ad(ANNOTATED_H5AD)
        mapping = pd.read_csv(GENE_MAPPING_MASTER)
        validate_input(adata, mapping)
        donors = donor_order(adata.obs)
        annotation = pd.DataFrame(
            {
                "feature_id": adata.var["feature_id"].astype(str).to_numpy(),
                "original_gene_symbol": adata.var["gene_symbol"].astype(str).to_numpy(),
                "unique_var_name": adata.var_names.astype(str),
            }
        )
        atomic_csv(annotation, PB_METADATA_DIR / "pseudobulk_gene_annotation.csv")
        sample_rows: list[dict[str, object]] = []
        for lineage in PRIMARY_LINEAGES:
            matrix, rows = aggregate_lineage(adata, lineage, donors)
            if len(rows) != 7 or not all(bool(row["eligible"]) for row in rows):
                raise RuntimeError(f"Primary lineage eligibility failed: {lineage}")
            path = PB_MATRICES_DIR / f"pseudobulk_{lineage}_counts.csv.gz"
            atomic_csv(matrix, path, compression="gzip")
            sample_rows.extend(rows)
            logger.info(
                "PSEUDOBULK lineage=%s cells=%d samples=%d features=%d library_min=%d library_max=%d output=%s",
                lineage,
                int(sum(row["n_cells"] for row in rows)),
                len(rows),
                adata.n_vars,
                min(int(row["library_size"]) for row in rows),
                max(int(row["library_size"]) for row in rows),
                path,
            )
        sample_metadata = pd.DataFrame(sample_rows)
        atomic_csv(sample_metadata, PB_METADATA_DIR / "pseudobulk_sample_metadata.csv")
        level1_table = level1_eligibility(adata)
        level2_table = level2_eligibility(adata)
        atomic_csv(level1_table, PB_METADATA_DIR / "secondary_level1_pseudobulk_eligibility.csv")
        atomic_csv(level2_table, PB_METADATA_DIR / "level2_pseudobulk_eligibility.csv")
        logger.info(
            "ELIGIBILITY primary=3 secondary_eligible=%d/%d level2_eligible=%d/%d",
            int(level1_table["eligible"].sum()),
            len(level1_table),
            int(level2_table["eligible"].sum()),
            len(level2_table),
        )
        logger.info("END status=success elapsed_seconds=%.3f peak_ram_mb=%.1f", time.time() - start, peak_ram_mb())
    except Exception:
        logger.exception("END status=error elapsed_seconds=%.3f peak_ram_mb=%.1f", time.time() - start, peak_ram_mb())
        raise


if __name__ == "__main__":
    main()
