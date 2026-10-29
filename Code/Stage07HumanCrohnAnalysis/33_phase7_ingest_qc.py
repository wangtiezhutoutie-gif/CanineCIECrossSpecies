"""Phase7B: ingest the frozen TI 10x archives and apply the published Human QC thresholds."""

from __future__ import annotations

import csv
import gzip
import json
import logging
import shutil
import tarfile
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
META = ROOT / "00_metadata"
BASE = ROOT / "07_human_crohn" / "GSE266616"
H5AD = BASE / "h5ad"
QC = BASE / "qc"
WORK = ROOT / "work" / "phase7_human_ingest"
LOG = ROOT / "logs" / "phase7_qc.log"

SAMPLE_MASTER = META / "GSE266616_HUMAN_SAMPLE_MASTER.csv"
DOWNLOAD_MANIFEST = META / "GSE266616_DOWNLOAD_MANIFEST.csv"
DATA12_JSON = BASE / "metadata" / "extracted" / "Supplementary_Data12.json"
ZENODO_CELLS = BASE / "metadata" / "source" / "zenodo" / "scRNA_cell_metadata.csv"

RAW_OUT = H5AD / "human_GSE266616_TI_raw.h5ad"
QC_OUT = H5AD / "human_GSE266616_TI_qc.h5ad"
RAW_SAMPLE_DIR = H5AD / "per_sample_raw"
QC_SAMPLE_DIR = H5AD / "per_sample_qc"
GENE_MASTER = META / "GSE266616_HUMAN_GENE_MASTER.csv"
PER_SAMPLE_QC = QC / "HUMAN_QC_PER_SAMPLE.csv"
PER_DONOR_QC = QC / "HUMAN_QC_PER_DONOR.csv"
ZENODO_MATCH = QC / "HUMAN_ZENODO_CELL_ID_MATCH.csv"
ROUTE_DECISION = QC / "HUMAN_QC_ROUTE_DECISION.csv"
QC_GATE = QC / "PHASE7B_QC_GATE.csv"
QC_REPORT = QC / "PHASE7B_QC_REPORT.md"

QC_MIN_GENES = 200
QC_MAX_GENES = 6000
QC_MIN_UMI_EXCLUSIVE = 1000
QC_MAX_MT_PCT_EXCLUSIVE = 50.0
EXPECTED_SAMPLES = 48
EXPECTED_BYTES = 2_302_117_661
REQUIRED_SUFFIXES = ("features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def setup_logging() -> logging.Logger:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("phase7_qc")
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


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temp, index=False, encoding="utf-8-sig")
    temp.replace(path)


def audit_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_master = pd.read_csv(SAMPLE_MASTER, dtype=str, keep_default_na=False)
    selected = sample_master.loc[sample_master["region"].eq("TI")].copy()
    if len(selected) != EXPECTED_SAMPLES:
        raise RuntimeError(f"Expected {EXPECTED_SAMPLES} frozen TI samples, observed {len(selected)}")
    if not selected["tissue"].eq("terminal_ileum").all():
        raise RuntimeError("Frozen TI table contains a non-terminal-ileum tissue")
    if int(selected["processed_file_size"].astype(np.int64).sum()) != EXPECTED_BYTES:
        raise RuntimeError("Frozen TI archive byte total changed")

    manifest = pd.read_csv(DOWNLOAD_MANIFEST, dtype=str, keep_default_na=False)
    if len(manifest) != EXPECTED_SAMPLES or manifest["GSM"].duplicated().any():
        raise RuntimeError("Formal download manifest is not one row per 48 TI GSM")
    required_gsms = set(selected["GSM"])
    if set(manifest["GSM"]) != required_gsms:
        raise RuntimeError("Download manifest GSM set differs from frozen TI selection")
    if not manifest["required_10x_triplet_present"].eq("YES").all():
        raise RuntimeError("Download manifest contains an archive without the required 10x triplet")
    if not (manifest["expected_bytes"].astype(np.int64) == manifest["actual_bytes"].astype(np.int64)).all():
        raise RuntimeError("Download manifest contains a byte-size mismatch")
    for row in manifest.itertuples(index=False):
        path = Path(row.local_path)
        if not path.is_file() or path.stat().st_size != int(row.actual_bytes):
            raise RuntimeError(f"Downloaded archive missing or changed: {path}")

    selected = selected.merge(
        manifest[["GSM", "local_path", "sha256", "actual_bytes"]], on="GSM", how="left", validate="one_to_one"
    )
    return sample_master, selected


def copy_required_members(archive_path: Path, destination: Path) -> dict[str, Path]:
    outputs: dict[str, Path] = {}
    with tarfile.open(archive_path, "r:gz") as archive:
        files = [member for member in archive.getmembers() if member.isfile()]
        for member in files:
            pure = PurePosixPath(member.name)
            if pure.is_absolute() or ".." in pure.parts:
                raise RuntimeError(f"Unsafe tar member in {archive_path.name}: {member.name}")
        for suffix in REQUIRED_SUFFIXES:
            matches = [member for member in files if member.name.lower().endswith(suffix)]
            if len(matches) != 1:
                raise RuntimeError(f"Expected one {suffix} in {archive_path.name}, observed {len(matches)}")
            source = archive.extractfile(matches[0])
            if source is None:
                raise RuntimeError(f"Unable to read {matches[0].name} from {archive_path.name}")
            output = destination / suffix
            with source, output.open("wb") as handle:
                shutil.copyfileobj(source, handle, length=8 * 1024 * 1024)
            outputs[suffix] = output
    return outputs


def read_10x_archive(archive_path: Path) -> tuple[sparse.csr_matrix, pd.DataFrame, list[str]]:
    WORK.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="gsm_", dir=WORK) as temp_dir:
        files = copy_required_members(archive_path, Path(temp_dir))
        with gzip.open(files["features.tsv.gz"], "rt", encoding="utf-8") as handle:
            features = pd.read_csv(
                handle,
                sep="\t",
                header=None,
                names=["gene_id", "gene_symbol", "feature_type"],
                dtype=str,
                keep_default_na=False,
            )
        with gzip.open(files["barcodes.tsv.gz"], "rt", encoding="utf-8") as handle:
            barcodes = [line.strip() for line in handle if line.strip()]
        matrix = mmread(files["matrix.mtx.gz"]).tocsr()
    if matrix.shape != (len(features), len(barcodes)):
        raise RuntimeError(
            f"10x dimension mismatch for {archive_path.name}: matrix={matrix.shape}, features={len(features)}, barcodes={len(barcodes)}"
        )
    if matrix.data.size and (matrix.data.min() < 0 or not np.allclose(matrix.data, np.rint(matrix.data))):
        raise RuntimeError(f"Non-integer or negative counts in {archive_path.name}")
    if matrix.data.size and matrix.data.max() > np.iinfo(np.int32).max:
        raise RuntimeError(f"Count exceeds int32 range in {archive_path.name}")
    matrix.data = matrix.data.astype(np.int32, copy=False)
    if features["gene_id"].duplicated().any():
        raise RuntimeError(f"Duplicate Ensembl feature IDs in {archive_path.name}")
    if len(barcodes) != len(set(barcodes)):
        raise RuntimeError(f"Duplicate cell barcodes in {archive_path.name}")
    return matrix, features, barcodes


def qc_vectors(matrix: sparse.csr_matrix, features: pd.DataFrame) -> tuple[np.ndarray, ...]:
    total_counts = np.asarray(matrix.sum(axis=0)).ravel().astype(np.int64)
    n_genes = np.asarray(matrix.getnnz(axis=0)).ravel().astype(np.int32)
    mt_mask = features["gene_symbol"].str.upper().str.startswith("MT-").to_numpy()
    mt_counts = np.asarray(matrix[mt_mask].sum(axis=0)).ravel().astype(np.int64)
    mt_pct = np.divide(
        mt_counts,
        total_counts,
        out=np.zeros(len(total_counts), dtype=np.float64),
        where=total_counts > 0,
    ) * 100.0
    passed = (
        (n_genes >= QC_MIN_GENES)
        & (n_genes <= QC_MAX_GENES)
        & (total_counts > QC_MIN_UMI_EXCLUSIVE)
        & (mt_pct < QC_MAX_MT_PCT_EXCLUSIVE)
    )
    return total_counts, n_genes, mt_counts, mt_pct, passed


def make_obs(
    row: pd.Series,
    barcodes: list[str],
    total_counts: np.ndarray,
    n_genes: np.ndarray,
    mt_counts: np.ndarray,
    mt_pct: np.ndarray,
    passed: np.ndarray,
) -> pd.DataFrame:
    index = pd.Index([f"{row.GSM}_{barcode}" for barcode in barcodes], name="cell_id")
    obs = pd.DataFrame(index=index)
    obs["barcode"] = barcodes
    for source, target in (
        ("GSM", "GSM"),
        ("sample_title", "sample_title"),
        ("patient_id", "donor"),
        ("healthy_or_crohn", "disease_group"),
        ("inflammation_state", "state"),
        ("region", "region"),
        ("tissue", "tissue"),
        ("biologic_treatment", "biologic_treatment"),
        ("biologic_class", "biologic_class"),
        ("age", "age"),
        ("race", "race"),
        ("platform", "platform"),
        ("primary_contrast_eligible", "primary_contrast_eligible"),
        ("secondary_paired_eligible", "secondary_paired_eligible"),
    ):
        obs[target] = str(row[source])
    obs["total_counts"] = total_counts
    obs["n_genes_by_counts"] = n_genes
    obs["mt_counts"] = mt_counts
    obs["pct_counts_mt"] = mt_pct.astype(np.float32)
    obs["qc_pass_author_thresholds"] = passed
    obs["qc_route"] = "ROUTE_B_PUBLISHED_THRESHOLDS_FULL_COHORT"
    return obs


def var_frame(features: pd.DataFrame) -> pd.DataFrame:
    frame = features.copy()
    frame.index = pd.Index(frame["gene_id"].astype(str), name="gene_id_index")
    return frame


def quantile(values: np.ndarray, probability: float) -> float:
    return float(np.quantile(values, probability)) if len(values) else float("nan")


def process_samples(selected: pd.DataFrame, logger: logging.Logger) -> tuple[list[Path], list[Path], pd.DataFrame, dict[str, set[str]], dict[str, set[str]]]:
    RAW_SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
    QC_SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
    H5AD.mkdir(parents=True, exist_ok=True)
    QC.mkdir(parents=True, exist_ok=True)
    reference_features: pd.DataFrame | None = None
    raw_files: list[Path] = []
    clean_files: list[Path] = []
    summaries: list[dict[str, object]] = []
    raw_barcodes: dict[str, set[str]] = {}
    clean_barcodes: dict[str, set[str]] = {}

    for index, row in enumerate(selected.itertuples(index=False), start=1):
        row_series = pd.Series(row._asdict())
        archive_path = Path(row.local_path)
        logger.info("INGEST %d/%d %s %s", index, len(selected), row.GSM, archive_path.name)
        matrix, features, barcodes = read_10x_archive(archive_path)
        if reference_features is None:
            reference_features = features.copy()
            gene_master = reference_features.copy()
            gene_master.insert(0, "unique_gene_id", gene_master["gene_id"])
            atomic_csv(gene_master, GENE_MASTER)
        elif not features.equals(reference_features):
            raise RuntimeError(f"Gene feature table differs from first TI sample: {row.GSM}")

        total_counts, n_genes, mt_counts, mt_pct, passed = qc_vectors(matrix, features)
        obs = make_obs(row_series, barcodes, total_counts, n_genes, mt_counts, mt_pct, passed)
        var = var_frame(features)
        raw_counts = matrix.transpose().tocsr()
        raw = ad.AnnData(X=raw_counts, obs=obs, var=var)
        raw.uns["qc_route"] = "GEO Cell Ranger filtered matrix before Phase7 Route B thresholds"
        raw_path = RAW_SAMPLE_DIR / f"{row.GSM}_raw.h5ad"
        raw.write_h5ad(raw_path, compression="gzip", compression_opts=4)
        raw_files.append(raw_path)

        retained_counts = matrix[:, passed].transpose().tocsr()
        clean_obs = obs.loc[passed].copy()
        clean = ad.AnnData(X=retained_counts, obs=clean_obs, var=var.copy())
        clean.uns["qc_route"] = "Route B: 200-6000 genes inclusive; >1000 UMI; mt%<50; no added doublet filter"
        clean_path = QC_SAMPLE_DIR / f"{row.GSM}_qc.h5ad"
        clean.write_h5ad(clean_path, compression="gzip", compression_opts=4)
        clean_files.append(clean_path)

        raw_barcodes[str(row.GSM)] = set(barcodes)
        clean_barcodes[str(row.GSM)] = set(np.asarray(barcodes, dtype=object)[passed].tolist())
        summaries.append(
            {
                "GSM": row.GSM,
                "sample_title": row.sample_title,
                "patient_id": row.patient_id,
                "disease_group": row.healthy_or_crohn,
                "state": row.inflammation_state,
                "primary_contrast_eligible": row.primary_contrast_eligible,
                "raw_cells": len(barcodes),
                "retained_cells": int(passed.sum()),
                "retained_fraction": float(passed.mean()),
                "failed_cells": int((~passed).sum()),
                "failed_n_genes_low": int((n_genes < QC_MIN_GENES).sum()),
                "failed_n_genes_high": int((n_genes > QC_MAX_GENES).sum()),
                "failed_umi": int((total_counts <= QC_MIN_UMI_EXCLUSIVE).sum()),
                "failed_mt": int((mt_pct >= QC_MAX_MT_PCT_EXCLUSIVE).sum()),
                "median_umi_raw": quantile(total_counts, 0.5),
                "median_genes_raw": quantile(n_genes, 0.5),
                "median_mt_pct_raw": quantile(mt_pct, 0.5),
                "min_retained_umi": int(total_counts[passed].min()) if passed.any() else 0,
                "max_retained_genes": int(n_genes[passed].max()) if passed.any() else 0,
                "max_retained_mt_pct": float(mt_pct[passed].max()) if passed.any() else float("nan"),
                "qc_route": "ROUTE_B_PUBLISHED_THRESHOLDS_FULL_COHORT",
            }
        )
        del raw, clean, raw_counts, retained_counts, matrix
    if reference_features is None:
        raise RuntimeError("No TI samples were processed")
    return raw_files, clean_files, pd.DataFrame(summaries), raw_barcodes, clean_barcodes


def concat_objects(raw_files: list[Path], clean_files: list[Path], logger: logging.Logger) -> None:
    for output in (RAW_OUT, QC_OUT):
        if output.exists():
            raise RuntimeError(f"Refusing to overwrite existing combined Human object without audit: {output}")
    logger.info("CONCAT raw per-sample h5ad files -> %s", RAW_OUT)
    ad.experimental.concat_on_disk(
        [str(path) for path in raw_files],
        str(RAW_OUT),
        axis=0,
        join="inner",
        merge="same",
        max_loaded_elems=50_000_000,
    )
    logger.info("CONCAT QC-retained per-sample h5ad files -> %s", QC_OUT)
    ad.experimental.concat_on_disk(
        [str(path) for path in clean_files],
        str(QC_OUT),
        axis=0,
        join="inner",
        merge="same",
        max_loaded_elems=50_000_000,
    )


def audit_zenodo(
    raw_barcodes: dict[str, set[str]], clean_barcodes: dict[str, set[str]]
) -> pd.DataFrame:
    data12 = json.loads(DATA12_JSON.read_text(encoding="utf-8"))["rows"]
    mappings: dict[str, str] = {}
    for row in data12:
        if row["scRNA GSE Series"] == "GSE266616" and row["Biopsy Location"] == "TI":
            mappings[str(row["scRNA GSM ID"])] = f"{row['Donor']}-TI"
    zenodo = pd.read_csv(ZENODO_CELLS, dtype=str, keep_default_na=False, low_memory=False)
    rows: list[dict[str, object]] = []
    for gsm in sorted(raw_barcodes):
        library = mappings.get(gsm)
        if library is None:
            rows.append(
                {
                    "GSM": gsm,
                    "Seq_Library_ID": "NA",
                    "published_cells": 0,
                    "raw_cells": len(raw_barcodes[gsm]),
                    "published_in_raw": 0,
                    "routeB_cells": len(clean_barcodes[gsm]),
                    "published_in_routeB": 0,
                    "routeB_extra_vs_published": len(clean_barcodes[gsm]),
                    "published_missing_after_routeB": 0,
                    "coverage_status": "NO_PUBLISHED_METADATA_COVERAGE",
                }
            )
            continue
        subset = zenodo.loc[zenodo["Seq_Library_ID"].eq(library)]
        published = set(subset["barcode"].astype(str))
        raw = raw_barcodes[gsm]
        clean = clean_barcodes[gsm]
        rows.append(
            {
                "GSM": gsm,
                "Seq_Library_ID": library,
                "published_cells": len(published),
                "raw_cells": len(raw),
                "published_in_raw": len(published & raw),
                "routeB_cells": len(clean),
                "published_in_routeB": len(published & clean),
                "routeB_extra_vs_published": len(clean - published),
                "published_missing_after_routeB": len(published - clean),
                "coverage_status": "HEALTHY_REFERENCE_EXACT_BARCODE_SUBSET" if published <= raw else "PUBLISHED_BARCODE_MISMATCH",
            }
        )
    frame = pd.DataFrame(rows)
    atomic_csv(frame, ZENODO_MATCH)
    return frame


def write_outputs(summary: pd.DataFrame, zenodo_match: pd.DataFrame, logger: logging.Logger) -> None:
    atomic_csv(summary, PER_SAMPLE_QC)
    donor = (
        summary.groupby(["patient_id", "disease_group"], as_index=False)
        .agg(
            n_TI_samples=("GSM", "nunique"),
            raw_cells=("raw_cells", "sum"),
            retained_cells=("retained_cells", "sum"),
            primary_samples=("primary_contrast_eligible", lambda x: int((x == "YES_SELECTED").sum())),
        )
        .sort_values(["disease_group", "patient_id"])
    )
    donor["retained_fraction"] = donor["retained_cells"] / donor["raw_cells"]
    atomic_csv(donor, PER_DONOR_QC)

    route = pd.DataFrame(
        [
            {
                "decision_item": "GEO per-GSM matrix stage",
                "decision": "CELL_RANGER_FILTERED_NOT_AUTHOR_FINAL_QC",
                "evidence": "GEO matrices include cells outside 200-6000 genes, >1000 UMI, and mt<50 thresholds",
            },
            {
                "decision_item": "Zenodo final-cell coverage",
                "decision": "SIX_HEALTHY_DONORS_ONLY",
                "evidence": "12 GSE266616 libraries in Data12/Zenodo (six paired TI/AC Healthy donors); zero Crohn libraries",
            },
            {
                "decision_item": "Formal full-cohort QC route",
                "decision": "ROUTE_B_PUBLISHED_THRESHOLDS_FULL_COHORT",
                "evidence": "Consistent across all Healthy and Crohn samples: 200-6000 genes inclusive; >1000 UMI; mt<50%; no added doublet filter",
            },
            {
                "decision_item": "Zenodo use",
                "decision": "REFERENCE_AND_PROVENANCE_AUDIT_ONLY",
                "evidence": "Exact Healthy barcode subset and published annotation reference; not represented as full-cohort annotation",
            },
        ]
    )
    atomic_csv(route, ROUTE_DECISION)

    reference = zenodo_match.loc[zenodo_match["Seq_Library_ID"].ne("NA")]
    criteria = [
        ("48/48 TI archives ingested", len(summary) == EXPECTED_SAMPLES, f"n={len(summary)}"),
        ("All samples retain cells", bool((summary["retained_cells"] > 0).all()), f"min={summary['retained_cells'].min()}"),
        ("Human raw object created", RAW_OUT.is_file(), str(RAW_OUT)),
        ("Human QC object created", QC_OUT.is_file(), str(QC_OUT)),
        ("Gene identity stable", GENE_MASTER.is_file(), "One identical 10x feature table across 48 TI samples"),
        ("Published Healthy barcodes subset of GEO", bool((reference["coverage_status"] == "HEALTHY_REFERENCE_EXACT_BARCODE_SUBSET").all()), f"n={len(reference)} libraries"),
        ("Uniform Route B applied", summary["qc_route"].eq("ROUTE_B_PUBLISHED_THRESHOLDS_FULL_COHORT").all(), "48/48 samples"),
    ]
    gate_pass = all(status for _, status, _ in criteria)
    gate = pd.DataFrame(
        [
            {
                "criterion": name,
                "status": "PASS" if status else "FAIL",
                "evidence": evidence,
                "gate": "GREEN_QC_GATE" if gate_pass else "YELLOW_QC_GATE_NOT_RESOLVED",
            }
            for name, status, evidence in criteria
        ]
    )
    atomic_csv(gate, QC_GATE)
    if not gate_pass:
        raise RuntimeError("Phase7B QC gate failed")

    primary = summary.loc[summary["primary_contrast_eligible"].eq("YES_SELECTED")]
    report = f"""# Phase7B Human TI QC Report

## Gate

**GREEN — QC GATE PASSED**

- TI archives ingested: {len(summary)}.
- Raw cells: {int(summary['raw_cells'].sum()):,}.
- Route-B retained cells: {int(summary['retained_cells'].sum()):,} ({summary['retained_cells'].sum() / summary['raw_cells'].sum():.1%}).
- Primary raw/retained cells: {int(primary['raw_cells'].sum()):,} / {int(primary['retained_cells'].sum()):,}.
- Minimum retained cells in any TI sample: {int(summary['retained_cells'].min()):,}.

## QC route

GEO per-GSM archives are Cell Ranger filtered matrices, not the author-final retained-cell object. The full cohort therefore uses Route B uniformly: 200–6000 detected genes (inclusive), >1000 UMI, and mitochondrial percentage <50%. No additional doublet pipeline was introduced.

Supplementary Data12 and Zenodo provide exact final-cell barcodes and published annotation for six Healthy GSE266616 donors only. These barcodes are retained as provenance/reference evidence; they are not represented as full-cohort published annotation and are not used to apply a different Healthy-only filter.

## Objects

- Raw counts before Phase7 Route B: `{RAW_OUT}`
- Route-B counts before annotation: `{QC_OUT}`

Generated UTC: {utc_now()}
"""
    QC_REPORT.write_text(report, encoding="utf-8", newline="\n")
    logger.info(
        "QC GATE PASS raw_cells=%d retained_cells=%d retained_fraction=%.4f",
        int(summary["raw_cells"].sum()),
        int(summary["retained_cells"].sum()),
        float(summary["retained_cells"].sum() / summary["raw_cells"].sum()),
    )


def main() -> int:
    logger = setup_logging()
    logger.info("START Phase7B ingest/QC anndata=%s scipy_sparse=%s", ad.__version__, sparse.__name__)
    _, selected = audit_inputs()
    raw_files, clean_files, summary, raw_barcodes, clean_barcodes = process_samples(selected, logger)
    concat_objects(raw_files, clean_files, logger)
    zenodo_match = audit_zenodo(raw_barcodes, clean_barcodes)
    write_outputs(summary, zenodo_match, logger)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
