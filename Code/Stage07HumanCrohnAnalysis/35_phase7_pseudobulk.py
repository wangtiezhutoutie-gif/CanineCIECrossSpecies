"""Phase7C: aggregate Human TI raw counts at donor/sample × broad-lineage level."""

from __future__ import annotations

import logging
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
META = ROOT / "00_metadata"
BASE = ROOT / "07_human_crohn" / "GSE266616"
ANNOTATED_SAMPLE_DIR = BASE / "h5ad" / "per_sample_annotated"
PB = BASE / "pseudobulk"
MATRICES = PB / "matrices"
PB_META = PB / "metadata"
LOG = ROOT / "logs" / "phase7_pseudobulk.log"

SAMPLE_MASTER = META / "GSE266616_HUMAN_SAMPLE_MASTER.csv"
PRIMARY_METADATA = PB_META / "HUMAN_PSEUDOBULK_SAMPLE_METADATA.csv"
GENE_ANNOTATION = PB_META / "HUMAN_PSEUDOBULK_GENE_ANNOTATION.csv"
ELIGIBILITY = PB_META / "HUMAN_PRIMARY_LINEAGE_ELIGIBILITY.csv"
PAIRED_PLAN = PB_META / "HUMAN_PAIRED_SECONDARY_PLAN.csv"
PAIRED_METADATA = PB_META / "HUMAN_PAIRED_PSEUDOBULK_SAMPLE_METADATA.csv"
GATE_OUT = PB_META / "PHASE7C_PSEUDOBULK_GATE.csv"

PRIMARY_LINEAGES = ("Epithelial", "T", "Myeloid")
MIN_CELLS = 20
MIN_PAIRED_DONORS = 3
EXPECTED_PRIMARY = {"Healthy": 16, "ActiveCrohn": 11}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def setup_logging() -> logging.Logger:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("phase7_pseudobulk")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)sZ %(levelname)s %(message)s", datefmt="%Y-%m-%dT%H:%M:%S")
    file_handler = logging.FileHandler(LOG, mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def atomic_csv(frame: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix + ".tmp"
    temp = path.with_suffix(suffix)
    frame.to_csv(temp, index=False, compression=compression, encoding="utf-8-sig")
    temp.replace(path)


def audit_samples() -> pd.DataFrame:
    master = pd.read_csv(SAMPLE_MASTER, dtype=str, keep_default_na=False)
    selected = master.loc[master["region"].eq("TI")].copy()
    if len(selected) != 48:
        raise RuntimeError(f"Expected 48 frozen TI samples, observed {len(selected)}")
    for gsm in selected["GSM"]:
        path = ANNOTATED_SAMPLE_DIR / f"{gsm}_annotated.h5ad"
        if not path.is_file():
            raise RuntimeError(f"Missing annotated TI sample: {path}")
    primary = selected.loc[selected["primary_contrast_eligible"].eq("YES_SELECTED")]
    counts = primary.assign(
        primary_group=np.where(primary["healthy_or_crohn"].eq("Healthy"), "Healthy", "ActiveCrohn")
    )["primary_group"].value_counts().to_dict()
    if counts != EXPECTED_PRIMARY:
        raise RuntimeError(f"Primary donor/sample counts changed: {counts}")
    if primary["patient_id"].duplicated().any():
        raise RuntimeError("Primary table contains more than one sample per donor")
    return selected


def aggregate_all(
    selected: pd.DataFrame, logger: logging.Logger
) -> tuple[dict[tuple[str, str], np.ndarray], pd.DataFrame, pd.DataFrame]:
    sums: dict[tuple[str, str], np.ndarray] = {}
    cell_rows: list[dict[str, object]] = []
    reference_var: pd.DataFrame | None = None
    for index, row in enumerate(selected.itertuples(index=False), start=1):
        path = ANNOTATED_SAMPLE_DIR / f"{row.GSM}_annotated.h5ad"
        adata = ad.read_h5ad(path)
        if "counts" not in adata.layers:
            raise RuntimeError(f"Missing layers['counts'] in {path}")
        counts = adata.layers["counts"]
        if not sparse.issparse(counts):
            raise RuntimeError(f"Non-sparse counts layer in {path}")
        values = counts.data
        if values.size and (values.min() < 0 or not np.allclose(values, np.rint(values))):
            raise RuntimeError(f"Nonnegative integer count invariant failed in {path}")
        if reference_var is None:
            reference_var = adata.var.copy()
        elif not adata.var.equals(reference_var):
            raise RuntimeError(f"Gene identity/order changed in {row.GSM}")
        for lineage in PRIMARY_LINEAGES:
            mask = adata.obs["human_level1"].astype(str).eq(lineage).to_numpy()
            n_cells = int(mask.sum())
            summed = np.asarray(counts[mask].sum(axis=0)).ravel().astype(np.int64)
            sums[(str(row.GSM), lineage)] = summed
            cell_rows.append(
                {
                    "GSM": row.GSM,
                    "patient_id": row.patient_id,
                    "disease_group": row.healthy_or_crohn,
                    "state": row.inflammation_state,
                    "lineage": lineage,
                    "n_cells": n_cells,
                    "library_size": int(summed.sum()),
                    "detected_genes": int((summed > 0).sum()),
                    "eligible_ge20_cells": "YES" if n_cells >= MIN_CELLS else "NO",
                }
            )
        logger.info("AGGREGATE %d/%d %s cells=%d", index, len(selected), row.GSM, adata.n_obs)
        del adata
    if reference_var is None:
        raise RuntimeError("No annotated samples were aggregated")
    gene_annotation = reference_var[["gene_id", "gene_symbol", "feature_type"]].copy().reset_index(drop=True)
    gene_annotation.insert(0, "unique_gene_id", gene_annotation["gene_id"].astype(str))
    if gene_annotation["unique_gene_id"].duplicated().any():
        raise RuntimeError("Human pseudobulk feature IDs are not unique")
    return sums, pd.DataFrame(cell_rows), gene_annotation


def matrix_frame(
    gene_annotation: pd.DataFrame, columns: list[str], vectors: list[np.ndarray]
) -> pd.DataFrame:
    if len(columns) != len(vectors) or len(columns) != len(set(columns)):
        raise RuntimeError("Pseudobulk matrix column identity failure")
    frame = gene_annotation[["unique_gene_id", "gene_symbol"]].copy()
    for name, vector in zip(columns, vectors, strict=True):
        if len(vector) != len(frame):
            raise RuntimeError(f"Pseudobulk vector length mismatch: {name}")
        frame[name] = vector
    return frame


def build_primary(
    selected: pd.DataFrame,
    sums: dict[tuple[str, str], np.ndarray],
    cell_summary: pd.DataFrame,
    gene_annotation: pd.DataFrame,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary = selected.loc[selected["primary_contrast_eligible"].eq("YES_SELECTED")].copy()
    primary["condition"] = np.where(primary["healthy_or_crohn"].eq("Healthy"), "Healthy", "ActiveCrohn")
    primary["condition_order"] = primary["condition"].map({"Healthy": 0, "ActiveCrohn": 1})
    primary = primary.sort_values(["condition_order", "patient_id"])
    metadata_rows: list[dict[str, object]] = []
    eligibility_rows: list[dict[str, object]] = []
    for lineage in PRIMARY_LINEAGES:
        columns = primary["patient_id"].astype(str).tolist()
        vectors = [sums[(str(gsm), lineage)] for gsm in primary["GSM"]]
        frame = matrix_frame(gene_annotation, columns, vectors)
        output = MATRICES / f"human_{lineage}_counts.csv.gz"
        atomic_csv(frame, output, compression="gzip")
        lineage_meta = cell_summary.loc[
            cell_summary["GSM"].isin(primary["GSM"]) & cell_summary["lineage"].eq(lineage)
        ].set_index("GSM")
        for row in primary.itertuples(index=False):
            qc = lineage_meta.loc[row.GSM]
            metadata_rows.append(
                {
                    "analysis": "Primary_ActiveCrohn_vs_Healthy",
                    "lineage": lineage,
                    "patient_id": row.patient_id,
                    "sample_id": row.GSM,
                    "condition": row.condition,
                    "reference_level": "Healthy",
                    "state": row.inflammation_state,
                    "n_cells": int(qc.n_cells),
                    "library_size": int(qc.library_size),
                    "detected_genes": int(qc.detected_genes),
                    "eligible": "YES" if int(qc.n_cells) >= MIN_CELLS else "NO",
                    "platform": row.platform,
                    "biologic_treatment": row.biologic_treatment,
                    "biologic_class": row.biologic_class,
                    "age": row.age,
                    "race": row.race,
                    "design": "~condition",
                    "logFC_direction": "positive=ActiveCrohn higher",
                    "matrix_path": str(output),
                }
            )
        group_counts = pd.DataFrame(metadata_rows).loc[lambda x: x["lineage"].eq(lineage)]
        for condition in ("Healthy", "ActiveCrohn"):
            subset = group_counts.loc[group_counts["condition"].eq(condition)]
            eligibility_rows.append(
                {
                    "lineage": lineage,
                    "condition": condition,
                    "expected_donors": EXPECTED_PRIMARY[condition],
                    "donors_total": subset["patient_id"].nunique(),
                    "donors_ge20_cells": subset.loc[subset["n_cells"] >= MIN_CELLS, "patient_id"].nunique(),
                    "min_cells": int(subset["n_cells"].min()),
                    "min_library_size": int(subset["library_size"].min()),
                    "eligible": "YES" if len(subset) == EXPECTED_PRIMARY[condition] and (subset["n_cells"] >= MIN_CELLS).all() else "NO",
                }
            )
        logger.info("PRIMARY %s donors=%d matrix=%s", lineage, len(primary), output)
    return pd.DataFrame(metadata_rows), pd.DataFrame(eligibility_rows)


def paired_donors(selected: pd.DataFrame, target_state: str) -> list[str]:
    donors = []
    for donor, rows in selected.loc[selected["healthy_or_crohn"].eq("Crohn")].groupby("patient_id"):
        active = rows.loc[rows["inflammation_state"].eq("Active_inflamed") & rows["region"].eq("TI")]
        target = rows.loc[rows["inflammation_state"].eq(target_state) & rows["region"].eq("TI")]
        if len(active) > 1 or len(target) > 1:
            raise RuntimeError(f"Repeated paired TI state for donor {donor}: active={len(active)} target={len(target)}")
        if len(active) == 1 and len(target) == 1:
            donors.append(str(donor))
    return sorted(donors)


def build_paired(
    selected: pd.DataFrame,
    sums: dict[tuple[str, str], np.ndarray],
    cell_summary: pd.DataFrame,
    gene_annotation: pd.DataFrame,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata_rows: list[dict[str, object]] = []
    plan_rows: list[dict[str, object]] = []
    state_specs = (("Adjacent", "active_vs_adjacent"), ("Non_inflamed", "active_vs_noninflamed"))
    for target_state, slug in state_specs:
        donors = paired_donors(selected, target_state)
        executable = len(donors) >= MIN_PAIRED_DONORS
        plan_rows.append(
            {
                "comparison": slug,
                "target_state": target_state,
                "paired_donor_n": len(donors),
                "minimum_required": MIN_PAIRED_DONORS,
                "executable": "YES" if executable else "NO_DESCRIPTIVE_ONLY",
                "design": "~patient+state" if executable else "NA",
                "interpretation": "within-Crohn active-inflammation-associated; cross-sectional paired biopsy states, not progression",
            }
        )
        if not executable:
            logger.info("PAIRED SKIP %s donors=%d", slug, len(donors))
            continue
        comparison_samples = []
        crohn = selected.loc[selected["healthy_or_crohn"].eq("Crohn")]
        for donor in donors:
            rows = crohn.loc[crohn["patient_id"].eq(donor)]
            active = rows.loc[rows["inflammation_state"].eq("Active_inflamed")].iloc[0]
            target = rows.loc[rows["inflammation_state"].eq(target_state)].iloc[0]
            comparison_samples.extend([(active, "Active"), (target, target_state)])
        for lineage in PRIMARY_LINEAGES:
            columns = [str(row.GSM) for row, _ in comparison_samples]
            vectors = [sums[(str(row.GSM), lineage)] for row, _ in comparison_samples]
            frame = matrix_frame(gene_annotation, columns, vectors)
            output = MATRICES / f"human_{lineage}_{slug}_counts.csv.gz"
            atomic_csv(frame, output, compression="gzip")
            qc_lookup = cell_summary.loc[cell_summary["lineage"].eq(lineage)].set_index("GSM")
            for row, state_label in comparison_samples:
                qc = qc_lookup.loc[row.GSM]
                metadata_rows.append(
                    {
                        "comparison": slug,
                        "lineage": lineage,
                        "patient_id": row.patient_id,
                        "sample_id": row.GSM,
                        "state": state_label,
                        "reference_state": target_state,
                        "n_cells": int(qc.n_cells),
                        "library_size": int(qc.library_size),
                        "detected_genes": int(qc.detected_genes),
                        "eligible_ge20_cells": "YES" if int(qc.n_cells) >= MIN_CELLS else "NO",
                        "biologic_class": row.biologic_class,
                        "platform": row.platform,
                        "design": "~patient_id+state",
                        "logFC_direction": "positive=Active higher",
                        "matrix_path": str(output),
                    }
                )
            logger.info("PAIRED %s %s donors=%d matrix=%s", slug, lineage, len(donors), output)
    return pd.DataFrame(metadata_rows), pd.DataFrame(plan_rows)


def write_gate(primary_metadata: pd.DataFrame, eligibility: pd.DataFrame, paired_plan: pd.DataFrame) -> None:
    criteria = [
        ("Three primary lineage matrices created", all((MATRICES / f"human_{lineage}_counts.csv.gz").is_file() for lineage in PRIMARY_LINEAGES), "Epithelial/T/Myeloid"),
        ("Primary one sample per donor", primary_metadata.groupby("lineage")["patient_id"].nunique().eq(27).all(), "27 donors per lineage"),
        ("All primary lineage-donor pseudobulks have >=20 cells", primary_metadata["eligible"].eq("YES").all(), f"min={primary_metadata['n_cells'].min()}"),
        ("Expected Healthy/Active donors retained", eligibility["eligible"].eq("YES").all(), eligibility.to_json(orient="records")),
        ("Raw-count library sizes positive", bool((primary_metadata["library_size"] > 0).all()), f"min={primary_metadata['library_size'].min()}"),
        ("Paired secondary feasibility frozen before DE", len(paired_plan) == 2, paired_plan.to_json(orient="records")),
    ]
    gate_pass = all(status for _, status, _ in criteria)
    gate = pd.DataFrame(
        [
            {
                "criterion": name,
                "status": "PASS" if status else "FAIL",
                "evidence": evidence,
                "gate": "GREEN_PSEUDOBULK_GATE" if gate_pass else "YELLOW_PSEUDOBULK_GATE_NOT_RESOLVED",
            }
            for name, status, evidence in criteria
        ]
    )
    atomic_csv(gate, GATE_OUT)
    if not gate_pass:
        raise RuntimeError("Phase7C pseudobulk gate failed")


def main() -> int:
    logger = setup_logging()
    MATRICES.mkdir(parents=True, exist_ok=True)
    PB_META.mkdir(parents=True, exist_ok=True)
    selected = audit_samples()
    sums, cell_summary, gene_annotation = aggregate_all(selected, logger)
    atomic_csv(gene_annotation, GENE_ANNOTATION)
    primary_metadata, eligibility = build_primary(selected, sums, cell_summary, gene_annotation, logger)
    paired_metadata, paired_plan = build_paired(selected, sums, cell_summary, gene_annotation, logger)
    atomic_csv(primary_metadata, PRIMARY_METADATA)
    atomic_csv(eligibility, ELIGIBILITY)
    atomic_csv(paired_plan, PAIRED_PLAN)
    atomic_csv(paired_metadata, PAIRED_METADATA)
    write_gate(primary_metadata, eligibility, paired_plan)
    logger.info("PSEUDOBULK GATE PASS primary_rows=%d paired_rows=%d generated=%s", len(primary_metadata), len(paired_metadata), utc_now())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
