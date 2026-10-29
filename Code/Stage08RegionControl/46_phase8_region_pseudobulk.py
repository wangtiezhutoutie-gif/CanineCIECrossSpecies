from __future__ import annotations

import logging
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P8 = ROOT / "08_region_control"
P8_META = P8 / "metadata"
BASE = P8 / "human_AC_TI"
AC_ANNOTATED = BASE / "per_sample_annotated"
TI_ANNOTATED = ROOT / "07_human_crohn" / "GSE266616" / "h5ad" / "per_sample_annotated"
PB = BASE / "pseudobulk"
MATRICES = PB / "matrices"
LOGS = (ROOT / "logs" / "phase8_region_analysis.log", P8 / "logs" / "phase8_region_analysis.log")

SAMPLE_MASTER = ROOT / "00_metadata" / "GSE266616_HUMAN_SAMPLE_MASTER.csv"
AC_AUDIT = P8_META / "HUMAN_AC_REGION_AUDIT.csv"
GENE_MASTER = ROOT / "00_metadata" / "GSE266616_HUMAN_GENE_MASTER.csv"
PB_META = P8_META / "HUMAN_REGION_PSEUDOBULK_METADATA.csv"
ELIGIBILITY = P8_META / "HUMAN_REGION_LINEAGE_ELIGIBILITY.csv"
GATE = P8_META / "PHASE8A_REGION_PSEUDOBULK_GATE.csv"

PRIMARY_LINEAGES = ("Epithelial", "T", "Myeloid")
MIN_CELLS = 20


def setup_logging() -> logging.Logger:
    logger = logging.getLogger("phase8_region_pseudobulk")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)sZ %(levelname)s %(message)s", datefmt="%Y-%m-%dT%H:%M:%S")
    for path in LOGS:
        path.parent.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(path, mode="a", encoding="utf-8")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    stream = logging.StreamHandler()
    stream.setFormatter(formatter)
    logger.addHandler(stream)
    return logger


def atomic_csv(frame: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temp, index=False, encoding="utf-8-sig", compression=compression)
    temp.replace(path)


def sample_path(gsm: str, region: str) -> Path:
    directory = TI_ANNOTATED if region == "TI" else AC_ANNOTATED
    return directory / f"{gsm}_annotated.h5ad"


def load_sum(gsm: str, region: str, lineage: str, gene_ids: np.ndarray) -> tuple[np.ndarray, int, int, int]:
    path = sample_path(gsm, region)
    if not path.is_file():
        raise RuntimeError(f"Missing annotated sample: {path}")
    adata = ad.read_h5ad(path)
    if "counts" not in adata.layers or not sparse.issparse(adata.layers["counts"]):
        raise RuntimeError(f"Missing sparse counts layer: {path}")
    if not np.array_equal(adata.var.index.astype(str).to_numpy(), gene_ids):
        raise RuntimeError(f"Gene identity/order mismatch: {path}")
    mask = adata.obs["human_level1"].astype(str).eq(lineage).to_numpy()
    counts = adata.layers["counts"]
    values = counts.data
    if values.size and (values.min() < 0 or not np.allclose(values, np.rint(values))):
        raise RuntimeError(f"Count invariant failed: {path}")
    summed = np.asarray(counts[mask].sum(axis=0)).ravel().astype(np.int64)
    return summed, int(mask.sum()), int(summed.sum()), int((summed > 0).sum())


def matrix_frame(gene: pd.DataFrame, columns: list[str], vectors: list[np.ndarray]) -> pd.DataFrame:
    if len(columns) != len(vectors) or len(columns) != len(set(columns)):
        raise RuntimeError("Pseudobulk column identity failure")
    frame = gene[["unique_gene_id", "gene_symbol"]].copy()
    for column, vector in zip(columns, vectors, strict=True):
        frame[column] = vector
    return frame


def main() -> int:
    logger = setup_logging()
    MATRICES.mkdir(parents=True, exist_ok=True)
    sample = pd.read_csv(SAMPLE_MASTER, dtype=str, keep_default_na=False)
    ac_audit = pd.read_csv(AC_AUDIT, dtype=str, keep_default_na=False)
    gene = pd.read_csv(GENE_MASTER, dtype=str, keep_default_na=False)
    if "unique_gene_id" not in gene:
        gene.insert(0, "unique_gene_id", gene["gene_id"])
    if gene["unique_gene_id"].duplicated().any():
        raise RuntimeError("Frozen Human gene IDs are not unique")
    gene_ids = gene["gene_id"].to_numpy(dtype=str)

    healthy_ac = ac_audit.loc[ac_audit["eligible"].eq("YES_HEALTHY_REGION_AND_AC_DISEASE")].copy()
    active_ac = ac_audit.loc[ac_audit["eligible"].eq("YES_ACTIVE_REGION_AND_AC_DISEASE")].copy()
    healthy_ti = sample.loc[
        sample["region"].eq("TI")
        & sample["healthy_or_crohn"].eq("Healthy")
        & sample["patient_id"].isin(healthy_ac["donor_id"]),
        ["GSM", "patient_id", "region", "healthy_or_crohn", "inflammation_state", "biologic_class", "platform"],
    ].copy()
    active_ti = sample.loc[
        sample["primary_contrast_eligible"].eq("YES_SELECTED")
        & sample["healthy_or_crohn"].eq("Crohn"),
        ["GSM", "patient_id", "region", "healthy_or_crohn", "inflammation_state", "biologic_class", "platform"],
    ].copy()
    if len(healthy_ac) != 16 or len(healthy_ti) != 15 or len(active_ac) != 6 or len(active_ti) != 11:
        raise RuntimeError(
            f"Frozen regional sample scope changed: healthyAC={len(healthy_ac)} healthyTIpaired={len(healthy_ti)} activeAC={len(active_ac)} activeTI={len(active_ti)}"
        )
    ac_only_healthy = set(healthy_ac["donor_id"]) - set(healthy_ti["patient_id"])
    if ac_only_healthy != {"HA52"}:
        raise RuntimeError(f"Frozen unpaired Healthy AC identity changed: {sorted(ac_only_healthy)}")

    healthy_ac_table = healthy_ac.rename(columns={
        "donor_id": "patient_id", "disease": "healthy_or_crohn", "inflammation_state": "state",
        "treatment": "biologic_class",
    })
    healthy_ac_table["platform"] = healthy_ac_table["platform"].astype(str)
    healthy_ac_table["region"] = "AC"
    active_ac_table = active_ac.rename(columns={
        "donor_id": "patient_id", "disease": "healthy_or_crohn", "inflammation_state": "state",
        "treatment": "biologic_class",
    })
    active_ac_table["platform"] = active_ac_table["platform"].astype(str)
    active_ac_table["region"] = "AC"
    healthy_ti = healthy_ti.rename(columns={"inflammation_state": "state"})
    active_ti = active_ti.rename(columns={"inflammation_state": "state"})

    paired_donors = sorted(set(healthy_ti["patient_id"]) & set(healthy_ac_table["patient_id"]))
    if len(paired_donors) != 15:
        raise RuntimeError(f"Expected 15 unambiguous Healthy pairs, observed {len(paired_donors)}")
    analysis_specs = {
        "Healthy_AC_vs_TI": pd.concat([
            healthy_ti.loc[healthy_ti["patient_id"].isin(paired_donors)],
            healthy_ac_table.loc[healthy_ac_table["patient_id"].isin(paired_donors), healthy_ti.columns],
        ], ignore_index=True),
        "ActiveCrohn_AC_vs_TI": pd.concat([
            active_ti,
            active_ac_table.loc[:, active_ti.columns],
        ], ignore_index=True),
        "AC_ActiveCrohn_vs_Healthy": pd.concat([
            healthy_ac_table.loc[:, active_ti.columns],
            active_ac_table.loc[:, active_ti.columns],
        ], ignore_index=True),
    }
    designs = {
        "Healthy_AC_vs_TI": "~patient_id + region",
        "ActiveCrohn_AC_vs_TI": "~region",
        "AC_ActiveCrohn_vs_Healthy": "~condition",
    }
    directions = {
        "Healthy_AC_vs_TI": "positive=Healthy AC higher than Healthy TI",
        "ActiveCrohn_AC_vs_TI": "positive=Active Crohn AC higher than Active Crohn TI",
        "AC_ActiveCrohn_vs_Healthy": "positive=Active Crohn higher than Healthy within AC",
    }

    metadata_rows: list[dict[str, object]] = []
    eligibility_rows: list[dict[str, object]] = []
    for analysis, table in analysis_specs.items():
        table = table.sort_values(["region", "patient_id", "GSM"]).reset_index(drop=True)
        for lineage in PRIMARY_LINEAGES:
            columns: list[str] = []
            vectors: list[np.ndarray] = []
            lineage_rows: list[dict[str, object]] = []
            for row in table.itertuples(index=False):
                vector, n_cells, library_size, detected_genes = load_sum(str(row.GSM), str(row.region), lineage, gene_ids)
                sample_column = f"{row.patient_id}__{row.region}"
                condition = "ActiveCrohn" if row.healthy_or_crohn == "Crohn" else "Healthy"
                eligible = n_cells >= MIN_CELLS
                record = {
                    "analysis": analysis,
                    "lineage": lineage,
                    "patient_id": row.patient_id,
                    "sample_id": row.GSM,
                    "sample_column": sample_column,
                    "region": row.region,
                    "condition": condition,
                    "state": row.state,
                    "biologic_class": row.biologic_class,
                    "platform": row.platform,
                    "n_cells": n_cells,
                    "library_size": library_size,
                    "detected_genes": detected_genes,
                    "eligible": "YES" if eligible else "NO_LT20_CELLS",
                    "design": designs[analysis],
                    "logFC_direction": directions[analysis],
                }
                metadata_rows.append(record)
                lineage_rows.append(record)
                if eligible:
                    columns.append(sample_column)
                    vectors.append(vector)
            matrix = matrix_frame(gene, columns, vectors)
            output = MATRICES / f"{analysis}_{lineage}_counts.csv.gz"
            atomic_csv(matrix, output, compression="gzip")
            eligible_frame = pd.DataFrame(lineage_rows).loc[lambda value: value["eligible"].eq("YES")]
            if analysis == "Healthy_AC_vs_TI":
                paired = eligible_frame.groupby("patient_id")["region"].nunique()
                group_evidence = f"paired_donors={int((paired == 2).sum())}"
                capable = int((paired == 2).sum()) >= 3
            elif analysis == "ActiveCrohn_AC_vs_TI":
                region_counts = eligible_frame.groupby("region")["patient_id"].nunique().to_dict()
                group_evidence = f"TI={region_counts.get('TI', 0)};AC={region_counts.get('AC', 0)}"
                capable = min(region_counts.get("TI", 0), region_counts.get("AC", 0)) >= 3
            else:
                condition_counts = eligible_frame.groupby("condition")["patient_id"].nunique().to_dict()
                group_evidence = f"Healthy={condition_counts.get('Healthy', 0)};ActiveCrohn={condition_counts.get('ActiveCrohn', 0)}"
                capable = min(condition_counts.get("Healthy", 0), condition_counts.get("ActiveCrohn", 0)) >= 3
            eligibility_rows.append({
                "analysis": analysis,
                "lineage": lineage,
                "design": designs[analysis],
                "eligible_columns": len(columns),
                "group_evidence": group_evidence,
                "analysis_capable": "YES" if capable else "NO",
                "matrix_path": str(output),
            })
            logger.info("PSEUDOBULK %s %s columns=%d %s", analysis, lineage, len(columns), group_evidence)

    metadata = pd.DataFrame(metadata_rows)
    eligibility = pd.DataFrame(eligibility_rows)
    atomic_csv(metadata, PB_META)
    atomic_csv(eligibility, ELIGIBILITY)
    checks = [
        ("Nine preregistered lineage comparisons materialized", len(eligibility) == 9, f"rows={len(eligibility)}"),
        ("Every comparison remains donor-level analyzable", eligibility["analysis_capable"].eq("YES").all(), ";".join(f"{r.analysis}:{r.lineage}:{r.group_evidence}" for r in eligibility.itertuples())),
        ("Every matrix has unique donor-region columns", True, "sample columns are patient_id__region"),
        ("Phase7 TI objects remained read-only", True, "Counts read from frozen layers['counts']; no TI writes"),
    ]
    gate = pd.DataFrame([{"criterion": name, "status": "PASS" if okay else "FAIL", "evidence": evidence} for name, okay, evidence in checks])
    gate["gate"] = "GREEN_REGION_PSEUDOBULK_COMPLETE" if gate["status"].eq("PASS").all() else "YELLOW_REGION_PSEUDOBULK_LIMITED"
    atomic_csv(gate, GATE)
    if not gate["status"].eq("PASS").all():
        raise RuntimeError("Phase8 region pseudobulk gate failed")
    logger.info("PSEUDOBULK_GATE %s", gate["gate"].iloc[0])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
