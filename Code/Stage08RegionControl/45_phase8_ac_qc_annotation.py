from __future__ import annotations

import gzip
import logging
import shutil
import tarfile
import tempfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

import anndata as ad
import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P8 = ROOT / "08_region_control"
META = P8 / "metadata"
BASE = P8 / "human_AC_TI"
QC_DIR = BASE / "per_sample_qc"
ANNOTATED_DIR = BASE / "per_sample_annotated"
WORK = ROOT / "work" / "phase8_ac_ingest"
FIGURES = P8 / "figures"
LOGS = (ROOT / "logs" / "phase8_region_analysis.log", P8 / "logs" / "phase8_region_analysis.log")

AUDIT = META / "HUMAN_AC_REGION_AUDIT.csv"
DOWNLOAD = META / "HUMAN_AC_DOWNLOAD_MANIFEST.csv"
MODEL_PATH = ROOT / "07_human_crohn" / "GSE266616" / "annotation" / "human_level1_sgd_classifier.joblib"
PHASE7_GENE_MASTER = ROOT / "00_metadata" / "GSE266616_HUMAN_GENE_MASTER.csv"

QC_SUMMARY = META / "HUMAN_AC_QC_PER_SAMPLE.csv"
CELL_COUNTS = META / "HUMAN_AC_LEVEL1_CELL_COUNTS.csv"
DONOR_COVERAGE = META / "HUMAN_AC_LINEAGE_DONOR_COVERAGE.csv"
GATE = META / "PHASE8A_AC_QC_ANNOTATION_GATE.csv"

QC_MIN_GENES = 200
QC_MAX_GENES = 6000
QC_MIN_UMI_EXCLUSIVE = 1000
QC_MAX_MT_PCT_EXCLUSIVE = 50.0
MIN_PROBABILITY = 0.60
RELAXED_PROBABILITY_WITH_MARKER = 0.45
MIN_CELLS = 20
PRIMARY_LINEAGES = ("Epithelial", "T", "Myeloid")
MODEL_CLASSES = ("Epithelial", "T", "Myeloid", "B", "Plasma", "Mast", "Stromal", "Endothelial")
REQUIRED_SUFFIXES = ("features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz")

MARKERS: dict[str, tuple[str, ...]] = {
    "Epithelial": ("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT20", "VIL1"),
    "T": ("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "LCK", "CD247"),
    "Myeloid": ("LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "LYZ", "SPI1"),
    "B": ("MS4A1", "CD79A", "CD79B", "CD37", "CD74", "HLA-DRA", "CD22"),
    "Plasma": ("MZB1", "JCHAIN", "SDC1", "CD38", "DERL3", "IGHG1", "IGHA1"),
    "Mast": ("TPSAB1", "TPSB2", "KIT", "CPA3", "MS4A2", "HDC"),
    "Stromal": ("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "COL6A2"),
    "Endothelial": ("PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2", "PLVAP"),
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def setup_logging() -> logging.Logger:
    logger = logging.getLogger("phase8_ac_qc_annotation")
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


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temp, index=False, encoding="utf-8-sig")
    temp.replace(path)


def copy_required_members(archive_path: Path, destination: Path) -> dict[str, Path]:
    outputs: dict[str, Path] = {}
    with tarfile.open(archive_path, "r:gz") as archive:
        files = [member for member in archive.getmembers() if member.isfile()]
        for member in files:
            pure = PurePosixPath(member.name)
            if pure.is_absolute() or ".." in pure.parts:
                raise RuntimeError(f"Unsafe tar member: {member.name}")
        for suffix in REQUIRED_SUFFIXES:
            matches = [member for member in files if member.name.lower().endswith(suffix)]
            if len(matches) != 1:
                raise RuntimeError(f"Expected one {suffix} in {archive_path.name}, observed {len(matches)}")
            source = archive.extractfile(matches[0])
            if source is None:
                raise RuntimeError(f"Could not extract {matches[0].name}")
            output = destination / suffix
            with source, output.open("wb") as handle:
                shutil.copyfileobj(source, handle, length=8 * 1024 * 1024)
            outputs[suffix] = output
    return outputs


def read_10x(archive_path: Path) -> tuple[sparse.csr_matrix, pd.DataFrame, list[str]]:
    WORK.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ac_", dir=WORK) as temp_dir:
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
        raise RuntimeError(f"10x dimension mismatch: {archive_path.name} {matrix.shape}")
    if features["gene_id"].duplicated().any() or len(barcodes) != len(set(barcodes)):
        raise RuntimeError(f"Feature/barcode identity failure: {archive_path.name}")
    if matrix.data.size and (matrix.data.min() < 0 or not np.allclose(matrix.data, np.rint(matrix.data))):
        raise RuntimeError(f"Counts not nonnegative integers: {archive_path.name}")
    matrix.data = matrix.data.astype(np.int32, copy=False)
    return matrix, features, barcodes


def qc_vectors(matrix: sparse.csr_matrix, features: pd.DataFrame) -> tuple[np.ndarray, ...]:
    totals = np.asarray(matrix.sum(axis=0)).ravel().astype(np.int64)
    genes = np.asarray(matrix.getnnz(axis=0)).ravel().astype(np.int32)
    mt_mask = features["gene_symbol"].str.upper().str.startswith("MT-").to_numpy()
    mt_counts = np.asarray(matrix[mt_mask].sum(axis=0)).ravel().astype(np.int64)
    mt_pct = np.divide(mt_counts, totals, out=np.zeros(len(totals), dtype=float), where=totals > 0) * 100
    passed = (
        (genes >= QC_MIN_GENES)
        & (genes <= QC_MAX_GENES)
        & (totals > QC_MIN_UMI_EXCLUSIVE)
        & (mt_pct < QC_MAX_MT_PCT_EXCLUSIVE)
    )
    return totals, genes, mt_counts, mt_pct, passed


def normalize_log(counts: sparse.spmatrix) -> sparse.csr_matrix:
    matrix = counts.tocsr().astype(np.float32, copy=True)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    scale = np.divide(10_000.0, totals, out=np.zeros_like(totals, dtype=np.float32), where=totals > 0)
    matrix = sparse.diags(scale, format="csr") @ matrix
    matrix.data = np.log1p(matrix.data).astype(np.float32, copy=False)
    matrix.eliminate_zeros()
    return matrix.tocsr()


def marker_scores(normalized: sparse.csr_matrix, var: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    symbols = var["gene_symbol"].astype(str).str.upper().to_numpy()
    columns = []
    for label in MODEL_CLASSES:
        indices = np.flatnonzero(np.isin(symbols, np.asarray(MARKERS[label], dtype=object)))
        if len(indices) < 3:
            raise RuntimeError(f"Insufficient marker genes for {label}")
        columns.append(np.asarray(normalized[:, indices].mean(axis=1)).ravel())
    scores = np.column_stack(columns).astype(np.float32, copy=False)
    winners = np.asarray(MODEL_CLASSES, dtype=object)[scores.argmax(axis=1)]
    return winners.astype(str), scores.max(axis=1)


def create_figures(qc: pd.DataFrame, counts: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    ordered = qc.sort_values(["disease", "patient_id"])
    colors = np.where(ordered["disease"].eq("Healthy"), "#2878B5", "#C74440")
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(ordered["patient_id"], ordered["retained_fraction"], color=colors)
    ax.set_ylabel("Retained fraction")
    ax.set_xlabel("Donor")
    ax.set_ylim(0, 1)
    ax.tick_params(axis="x", rotation=70)
    ax.set_title("Phase8A AC Route B QC retention")
    fig.tight_layout()
    fig.savefig(FIGURES / "PHASE8_AC_QC_RETENTION.png", dpi=160)
    plt.close(fig)

    composition = counts.pivot_table(index="patient_id", columns="human_level1", values="n_cells", aggfunc="sum", fill_value=0)
    composition = composition.div(composition.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(13, 6))
    composition.plot(kind="bar", stacked=True, ax=ax, width=0.85, colormap="tab20")
    ax.set_ylabel("Fraction of retained cells")
    ax.set_xlabel("Donor")
    ax.set_title("Phase8A AC broad-lineage composition")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "PHASE8_AC_LEVEL1_COMPOSITION.png", dpi=160)
    plt.close(fig)


def main() -> int:
    logger = setup_logging()
    QC_DIR.mkdir(parents=True, exist_ok=True)
    ANNOTATED_DIR.mkdir(parents=True, exist_ok=True)
    audit = pd.read_csv(AUDIT, dtype=str, keep_default_na=False)
    audit = audit.loc[audit["eligible"].str.startswith("YES_")].copy()
    manifest = pd.read_csv(DOWNLOAD, dtype=str, keep_default_na=False)
    if len(audit) != 22 or len(manifest) != 22 or set(audit["GSM"]) != set(manifest["GSM"]):
        raise RuntimeError("AC QC input cardinality/identity failure")
    selected = audit.merge(manifest[["GSM", "local_path", "local_sha256"]], on="GSM", validate="one_to_one")
    bundle = joblib.load(MODEL_PATH)
    model = bundle["model"]
    use_features = np.asarray(bundle["feature_mask"], dtype=bool)
    model_ids = np.asarray(bundle["feature_ids"], dtype=str)
    phase7_genes = pd.read_csv(PHASE7_GENE_MASTER, dtype=str, keep_default_na=False)
    if not np.array_equal(phase7_genes["gene_id"].to_numpy(dtype=str), model_ids):
        raise RuntimeError("Frozen Phase7 gene master differs from annotation model")

    qc_rows: list[dict[str, object]] = []
    count_rows: list[pd.DataFrame] = []
    for index, row in enumerate(selected.itertuples(index=False), start=1):
        logger.info("INGEST_ANNOTATE %d/22 %s %s", index, row.GSM, row.donor_id)
        matrix, features, barcodes = read_10x(Path(row.local_path))
        if not np.array_equal(features["gene_id"].to_numpy(dtype=str), model_ids):
            raise RuntimeError(f"AC gene identity/order differs from frozen model: {row.GSM}")
        totals, genes, mt_counts, mt_pct, passed = qc_vectors(matrix, features)
        if not passed.any():
            raise RuntimeError(f"No retained AC cells: {row.GSM}")
        retained = matrix[:, passed].transpose().tocsr()
        retained_barcodes = np.asarray(barcodes, dtype=object)[passed].astype(str)
        obs = pd.DataFrame(index=pd.Index([f"{row.GSM}_{barcode}" for barcode in retained_barcodes], name="cell_id"))
        obs["barcode"] = retained_barcodes
        obs["GSM"] = row.GSM
        obs["donor"] = row.donor_id
        obs["disease_group"] = row.disease
        obs["state"] = row.inflammation_state
        obs["region"] = "AC"
        obs["biologic_class"] = row.treatment
        obs["total_counts"] = totals[passed]
        obs["n_genes_by_counts"] = genes[passed]
        obs["mt_counts"] = mt_counts[passed]
        obs["pct_counts_mt"] = mt_pct[passed].astype(np.float32)
        obs["qc_pass_author_thresholds"] = True
        obs["qc_route"] = "ROUTE_B_PHASE7_FROZEN_THRESHOLDS"
        var = features.copy()
        var.index = pd.Index(var["gene_id"].astype(str), name="gene_id_index")
        clean = ad.AnnData(X=retained.copy(), obs=obs.copy(), var=var.copy())
        clean.uns["qc_route"] = "Route B: 200-6000 genes inclusive; >1000 UMI; mt%<50; no added doublet filter"
        clean.write_h5ad(QC_DIR / f"{row.GSM}_qc.h5ad", compression="gzip", compression_opts=4)

        normalized = normalize_log(retained)
        probabilities = model.predict_proba(normalized[:, use_features])
        best_index = probabilities.argmax(axis=1)
        predicted = model.classes_[best_index].astype(str)
        best_probability = probabilities[np.arange(len(obs)), best_index]
        marker_winner, marker_max = marker_scores(normalized, var)
        concordant = predicted == marker_winner
        accepted = (best_probability >= MIN_PROBABILITY) | (concordant & (best_probability >= RELAXED_PROBABILITY_WITH_MARKER))
        level1 = np.where(accepted, predicted, "Uncertain").astype(str)
        confidence = np.where(best_probability >= 0.80, "HIGH", np.where(accepted, "MEDIUM", "LOW"))
        annotated = ad.AnnData(X=normalized, obs=obs.copy(), var=var.copy())
        annotated.layers["counts"] = retained.copy()
        annotated.obs["predicted_level1_raw"] = predicted
        annotated.obs["human_level1"] = level1
        annotated.obs["human_level2"] = "NA_BROAD_ONLY"
        annotated.obs["annotation_probability"] = best_probability.astype(np.float32)
        annotated.obs["annotation_confidence"] = confidence.astype(str)
        annotated.obs["marker_winner"] = marker_winner
        annotated.obs["marker_score_max"] = marker_max
        annotated.obs["marker_support"] = np.where(concordant, "CONCORDANT", "DISCORDANT")
        annotated.obs["annotation_source"] = np.where(
            accepted,
            "FROZEN_PHASE7_HUMAN_NATIVE_BROAD_LINEAGE_TRANSFER",
            "LOW_CONFIDENCE_UNASSIGNED",
        )
        annotated.uns["annotation_method"] = "Frozen Phase7 Human-native sparse SGD broad-lineage transfer; no retraining"
        annotated.uns["annotation_thresholds"] = {
            "primary_probability": MIN_PROBABILITY,
            "relaxed_probability_with_marker_concordance": RELAXED_PROBABILITY_WITH_MARKER,
        }
        annotated.write_h5ad(ANNOTATED_DIR / f"{row.GSM}_annotated.h5ad", compression="gzip", compression_opts=4)

        counts = (
            annotated.obs.groupby(["GSM", "donor", "disease_group", "state", "human_level1"], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
            .rename(columns={"donor": "patient_id", "disease_group": "disease"})
        )
        count_rows.append(counts)
        qc_rows.append({
            "GSM": row.GSM,
            "patient_id": row.donor_id,
            "disease": row.disease,
            "state": row.inflammation_state,
            "treatment": row.treatment,
            "region": "AC",
            "raw_cells": len(barcodes),
            "retained_cells": int(passed.sum()),
            "failed_cells": int((~passed).sum()),
            "retained_fraction": float(passed.mean()),
            "failed_n_genes_low": int((genes < QC_MIN_GENES).sum()),
            "failed_n_genes_high": int((genes > QC_MAX_GENES).sum()),
            "failed_umi": int((totals <= QC_MIN_UMI_EXCLUSIVE).sum()),
            "failed_mt": int((mt_pct >= QC_MAX_MT_PCT_EXCLUSIVE).sum()),
            "median_umi_raw": float(np.median(totals)),
            "median_genes_raw": float(np.median(genes)),
            "median_mt_pct_raw": float(np.median(mt_pct)),
            "uncertain_cells": int((level1 == "Uncertain").sum()),
            "uncertain_fraction": float((level1 == "Uncertain").mean()),
            "qc_route": "ROUTE_B_PHASE7_FROZEN_THRESHOLDS",
            "annotation_model": "FROZEN_PHASE7_HUMAN_NATIVE_SGD",
        })
        del clean, annotated, matrix, retained, normalized

    qc = pd.DataFrame(qc_rows).sort_values(["disease", "patient_id"])
    counts = pd.concat(count_rows, ignore_index=True).sort_values(["disease", "patient_id", "human_level1"])
    atomic_csv(qc, QC_SUMMARY)
    atomic_csv(counts, CELL_COUNTS)
    coverage_rows: list[dict[str, object]] = []
    for lineage in PRIMARY_LINEAGES:
        subset = counts.loc[counts["human_level1"].eq(lineage)]
        for disease, expected in (("Healthy", 16), ("Crohn", 6)):
            values = subset.loc[subset["disease"].eq(disease)].set_index("patient_id")["n_cells"]
            all_donors = set(qc.loc[qc["disease"].eq(disease), "patient_id"])
            eligible = set(values.loc[values >= MIN_CELLS].index)
            coverage_rows.append({
                "lineage": lineage,
                "disease": disease,
                "expected_donors": expected,
                "donors_with_cells": len(set(values.index)),
                "donors_ge20_cells": len(eligible),
                "min_cells_among_present": int(values.min()) if len(values) else 0,
                "excluded_lt20_or_absent": ";".join(sorted(all_donors - eligible)),
                "analysis_capable": "YES" if len(eligible) >= 3 else "NO",
            })
    coverage = pd.DataFrame(coverage_rows)
    atomic_csv(coverage, DONOR_COVERAGE)
    create_figures(qc, counts)

    checks = [
        ("22 AC sample objects completed", len(qc) == 22 and len(list(ANNOTATED_DIR.glob("*_annotated.h5ad"))) == 22, f"samples={len(qc)}"),
        ("All AC samples retain cells", qc["retained_cells"].gt(0).all(), f"raw={qc['raw_cells'].sum()};retained={qc['retained_cells'].sum()}"),
        ("Frozen gene identity preserved", True, f"genes={len(model_ids)} exact order"),
        ("Frozen Phase7 annotation model used", qc["annotation_model"].eq("FROZEN_PHASE7_HUMAN_NATIVE_SGD").all(), f"uncertain={qc['uncertain_cells'].sum()}/{qc['retained_cells'].sum()}"),
        ("Three Primary lineages remain donor-level analyzable", coverage["analysis_capable"].eq("YES").all(), ";".join(f"{r.lineage}_{r.disease}={r.donors_ge20_cells}" for r in coverage.itertuples())),
    ]
    gate = pd.DataFrame([{"criterion": name, "status": "PASS" if okay else "FAIL", "evidence": evidence} for name, okay, evidence in checks])
    gate["gate"] = "GREEN_AC_QC_ANNOTATION_COMPLETE" if gate["status"].eq("PASS").all() else "YELLOW_AC_LINEAGE_COVERAGE_LIMITED"
    atomic_csv(gate, GATE)
    if not gate["status"].eq("PASS").all():
        raise RuntimeError("AC QC/annotation gate failed")
    logger.info(
        "QC_ANNOTATION_GATE %s raw=%d retained=%d uncertain=%d",
        gate["gate"].iloc[0],
        int(qc["raw_cells"].sum()),
        int(qc["retained_cells"].sum()),
        int(qc["uncertain_cells"].sum()),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
