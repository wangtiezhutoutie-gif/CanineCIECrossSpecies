"""Phase7B Human broad-lineage annotation using exact published Healthy references."""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import joblib
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import classification_report, confusion_matrix


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
META = ROOT / "00_metadata"
BASE = ROOT / "07_human_crohn" / "GSE266616"
H5AD = BASE / "h5ad"
ANNOTATION = BASE / "annotation"
LOG = ROOT / "logs" / "phase7_annotation.log"

SAMPLE_MASTER = META / "GSE266616_HUMAN_SAMPLE_MASTER.csv"
DATA12_JSON = BASE / "metadata" / "extracted" / "Supplementary_Data12.json"
ZENODO_CELLS = BASE / "metadata" / "source" / "zenodo" / "scRNA_cell_metadata.csv"
QC_SAMPLE_DIR = H5AD / "per_sample_qc"
ANNOTATED_SAMPLE_DIR = H5AD / "per_sample_annotated"
ANNOTATED_OUT = H5AD / "human_GSE266616_TI_annotated.h5ad"

MAPPING_OUT = ANNOTATION / "HUMAN_CELLTYPE_MAPPING.csv"
MARKER_PANEL_OUT = ANNOTATION / "HUMAN_LEVEL1_MARKER_PANEL.csv"
CV_METRICS_OUT = ANNOTATION / "HUMAN_ANNOTATION_REFERENCE_CV.csv"
CV_CONFUSION_OUT = ANNOTATION / "HUMAN_ANNOTATION_REFERENCE_CONFUSION.csv"
CV_CELLS_OUT = ANNOTATION / "HUMAN_ANNOTATION_REFERENCE_CV_CELLS.csv"
MODEL_OUT = ANNOTATION / "human_level1_sgd_classifier.joblib"
CELL_COUNT_OUT = ANNOTATION / "HUMAN_LEVEL1_CELL_COUNTS.csv"
DONOR_COVERAGE_OUT = ANNOTATION / "HUMAN_LEVEL1_DONOR_COVERAGE.csv"
MARKER_SANITY_OUT = ANNOTATION / "HUMAN_LEVEL1_MARKER_SANITY.csv"
GATE_OUT = ANNOTATION / "PHASE7B_ANNOTATION_GATE.csv"
REPORT_OUT = ANNOTATION / "PHASE7B_ANNOTATION_REPORT.md"

SEED = 20240825
MIN_PROBABILITY = 0.60
RELAXED_PROBABILITY_WITH_MARKER = 0.45
MIN_CELLS_PER_DONOR_LINEAGE = 20
PRIMARY_LINEAGES = ("Epithelial", "T", "Myeloid")
MODEL_CLASSES = ("Epithelial", "T", "Myeloid", "B", "Plasma", "Mast", "Stromal", "Endothelial")

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
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("phase7_annotation")
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


def normalize_log(counts: sparse.spmatrix) -> sparse.csr_matrix:
    matrix = counts.tocsr().astype(np.float32, copy=True)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    scale = np.divide(10_000.0, totals, out=np.zeros_like(totals, dtype=np.float32), where=totals > 0)
    matrix = sparse.diags(scale, format="csr") @ matrix
    matrix.data = np.log1p(matrix.data).astype(np.float32, copy=False)
    matrix.eliminate_zeros()
    return matrix.tocsr()


def published_level1(lineage: str, cluster1: str) -> str:
    lineage = lineage.lower().strip()
    cluster = cluster1.lower().strip()
    if lineage == "epithelial":
        return "Epithelial"
    if lineage == "tcell":
        return "T"
    if lineage == "bcell":
        return "Plasma" if "plasma" in cluster else "B"
    if lineage == "myeloid":
        return "Mast" if "mast" in cluster else "Myeloid"
    if lineage == "others":
        if any(term in cluster for term in ("arterial", "venous", "endotheli")):
            return "Endothelial"
        return "Stromal"
    raise RuntimeError(f"Unexpected published lineage: {lineage} / {cluster1}")


def reference_mapping() -> tuple[dict[str, str], pd.DataFrame]:
    data12 = json.loads(DATA12_JSON.read_text(encoding="utf-8"))["rows"]
    gsm_to_library = {
        str(row["scRNA GSM ID"]): f"{row['Donor']}-TI"
        for row in data12
        if row["scRNA GSE Series"] == "GSE266616" and row["Biopsy Location"] == "TI"
    }
    if len(gsm_to_library) != 6:
        raise RuntimeError(f"Expected six published Healthy TI mappings, observed {len(gsm_to_library)}")
    zenodo = pd.read_csv(ZENODO_CELLS, dtype=str, keep_default_na=False, low_memory=False)
    zenodo["published_human_level1"] = [
        published_level1(lineage, cluster)
        for lineage, cluster in zip(zenodo["lineage"], zenodo["cluster_label1"], strict=True)
    ]
    return gsm_to_library, zenodo


def load_reference(
    gsm_to_library: dict[str, str], zenodo: pd.DataFrame, logger: logging.Logger
) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray, pd.DataFrame, pd.DataFrame]:
    matrices: list[sparse.csr_matrix] = []
    labels: list[np.ndarray] = []
    groups: list[np.ndarray] = []
    details: list[pd.DataFrame] = []
    reference_var: pd.DataFrame | None = None
    for gsm, library in sorted(gsm_to_library.items()):
        path = QC_SAMPLE_DIR / f"{gsm}_qc.h5ad"
        if not path.is_file():
            raise RuntimeError(f"Missing QC sample object for published reference: {path}")
        adata = ad.read_h5ad(path)
        if reference_var is None:
            reference_var = adata.var.copy()
        elif not adata.var.equals(reference_var):
            raise RuntimeError(f"Reference gene table mismatch: {gsm}")
        published = zenodo.loc[zenodo["Seq_Library_ID"].eq(library)].set_index("barcode", drop=False)
        matched = adata.obs["barcode"].astype(str).isin(published.index)
        if not matched.any():
            raise RuntimeError(f"No published barcodes survived Route B for {gsm} / {library}")
        obs = adata.obs.loc[matched, ["barcode"]].copy()
        obs = obs.join(
            published[["lineage", "cluster_label1", "cluster_label2", "cluster_label3", "published_human_level1"]],
            on="barcode",
            how="left",
            validate="many_to_one",
        )
        if obs["published_human_level1"].isna().any():
            raise RuntimeError(f"Published reference join failed for {gsm}")
        matrices.append(adata.X[matched.to_numpy(dtype=bool)].tocsr())
        labels.append(obs["published_human_level1"].to_numpy(dtype=str))
        groups.append(np.repeat(library, len(obs)))
        obs.insert(0, "GSM", gsm)
        obs.insert(1, "Seq_Library_ID", library)
        details.append(obs.reset_index(names="cell_id"))
        logger.info("REFERENCE %s %s matched_routeB_cells=%d", gsm, library, len(obs))
    if reference_var is None:
        raise RuntimeError("Published reference was empty")
    counts = sparse.vstack(matrices, format="csr")
    y = np.concatenate(labels)
    group = np.concatenate(groups)
    detail = pd.concat(details, ignore_index=True)
    return counts, y, group, reference_var, detail


def feature_mask(counts: sparse.csr_matrix, var: pd.DataFrame) -> np.ndarray:
    detected = np.asarray(counts.getnnz(axis=0)).ravel()
    symbols = var["gene_symbol"].astype(str).str.upper()
    excluded = symbols.str.startswith(("MT-", "RPL", "RPS")).to_numpy()
    mask = (detected >= 20) & ~excluded
    if int(mask.sum()) < 2_000:
        raise RuntimeError(f"Too few classifier features after reference detection filter: {int(mask.sum())}")
    return mask


def new_classifier() -> SGDClassifier:
    return SGDClassifier(
        loss="log_loss",
        penalty="l2",
        alpha=1e-5,
        max_iter=250,
        tol=1e-4,
        class_weight="balanced",
        random_state=SEED,
        average=True,
    )


def cross_validate(
    x: sparse.csr_matrix, y: np.ndarray, groups: np.ndarray, logger: logging.Logger
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    predicted = np.empty(len(y), dtype=object)
    for group in sorted(np.unique(groups)):
        train = groups != group
        test = groups == group
        if set(MODEL_CLASSES) - set(np.unique(y[train])):
            raise RuntimeError(f"Published reference class absent after holding out {group}")
        model = new_classifier()
        model.fit(x[train], y[train])
        predicted[test] = model.predict(x[test])
        logger.info("CV holdout=%s train=%d test=%d accuracy=%.4f", group, int(train.sum()), int(test.sum()), float((predicted[test] == y[test]).mean()))

    report = classification_report(
        y,
        predicted.astype(str),
        labels=list(MODEL_CLASSES),
        output_dict=True,
        zero_division=0,
    )
    rows = []
    for label in list(MODEL_CLASSES) + ["macro avg", "weighted avg"]:
        metrics = report[label]
        rows.append(
            {
                "label": label,
                "precision": metrics["precision"],
                "recall": metrics["recall"],
                "f1_score": metrics["f1-score"],
                "support": metrics["support"],
            }
        )
    rows.append(
        {
            "label": "accuracy",
            "precision": report["accuracy"],
            "recall": report["accuracy"],
            "f1_score": report["accuracy"],
            "support": len(y),
        }
    )
    metrics_frame = pd.DataFrame(rows)
    matrix = confusion_matrix(y, predicted.astype(str), labels=list(MODEL_CLASSES))
    confusion = pd.DataFrame(matrix, index=MODEL_CLASSES, columns=MODEL_CLASSES).reset_index(names="published_label")
    return metrics_frame, confusion, predicted.astype(str)


def marker_scores(normalized: sparse.csr_matrix, var: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    symbols = var["gene_symbol"].astype(str).str.upper().to_numpy()
    score_columns = []
    for label in MODEL_CLASSES:
        indices = np.flatnonzero(np.isin(symbols, np.array(MARKERS[label], dtype=object)))
        if len(indices) < 3:
            raise RuntimeError(f"Fewer than three marker genes found for {label}: {len(indices)}")
        score_columns.append(np.asarray(normalized[:, indices].mean(axis=1)).ravel())
    scores = np.column_stack(score_columns).astype(np.float32, copy=False)
    winners = np.asarray(MODEL_CLASSES, dtype=object)[scores.argmax(axis=1)]
    maxima = scores.max(axis=1)
    return winners.astype(str), maxima


def published_lookup_for_gsm(
    gsm: str, gsm_to_library: dict[str, str], zenodo: pd.DataFrame
) -> pd.DataFrame | None:
    library = gsm_to_library.get(gsm)
    if library is None:
        return None
    return zenodo.loc[zenodo["Seq_Library_ID"].eq(library)].set_index("barcode", drop=False)


def annotate_samples(
    sample_master: pd.DataFrame,
    model: SGDClassifier,
    use_features: np.ndarray,
    gsm_to_library: dict[str, str],
    zenodo: pd.DataFrame,
    logger: logging.Logger,
) -> tuple[list[Path], pd.DataFrame, pd.DataFrame]:
    ANNOTATED_SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
    count_rows: list[pd.DataFrame] = []
    sanity_rows: list[dict[str, object]] = []
    output_files: list[Path] = []
    selected = sample_master.loc[sample_master["region"].eq("TI")].copy()
    for index, row in enumerate(selected.itertuples(index=False), start=1):
        source = QC_SAMPLE_DIR / f"{row.GSM}_qc.h5ad"
        adata = ad.read_h5ad(source)
        counts = adata.X.tocsr().astype(np.int32, copy=True)
        normalized = normalize_log(counts)
        probabilities = model.predict_proba(normalized[:, use_features])
        best_index = probabilities.argmax(axis=1)
        predicted = model.classes_[best_index].astype(str)
        best_probability = probabilities[np.arange(len(adata)), best_index]
        marker_winner, marker_max = marker_scores(normalized, adata.var)
        marker_concordant = predicted == marker_winner
        accepted = (best_probability >= MIN_PROBABILITY) | (
            marker_concordant & (best_probability >= RELAXED_PROBABILITY_WITH_MARKER)
        )
        level1 = np.where(accepted, predicted, "Uncertain").astype(object)
        source_label = np.where(
            accepted,
            "SUPERVISED_TRANSFER_FROM_6_HEALTHY_PUBLISHED_REFERENCES",
            "LOW_CONFIDENCE_UNASSIGNED",
        ).astype(object)
        level2 = np.repeat("NA_BROAD_ONLY", len(adata)).astype(object)
        published_lineage = np.repeat("NA", len(adata)).astype(object)
        published_cluster1 = np.repeat("NA", len(adata)).astype(object)
        published_cluster2 = np.repeat("NA", len(adata)).astype(object)
        published_cluster3 = np.repeat("NA", len(adata)).astype(object)

        published = published_lookup_for_gsm(str(row.GSM), gsm_to_library, zenodo)
        if published is not None:
            barcodes = adata.obs["barcode"].astype(str)
            ref_mask = barcodes.isin(published.index).to_numpy()
            if ref_mask.any():
                ref = pd.DataFrame({"barcode": barcodes.loc[ref_mask].to_numpy()})
                ref = ref.join(
                    published[["lineage", "cluster_label1", "cluster_label2", "cluster_label3", "published_human_level1"]],
                    on="barcode",
                    how="left",
                    validate="many_to_one",
                )
                level1[ref_mask] = ref["published_human_level1"].to_numpy(dtype=object)
                level2[ref_mask] = ref["cluster_label1"].to_numpy(dtype=object)
                source_label[ref_mask] = "PUBLISHED_ZENODO_HEALTHY_REFERENCE"
                best_probability[ref_mask] = 1.0
                published_lineage[ref_mask] = ref["lineage"].to_numpy(dtype=object)
                published_cluster1[ref_mask] = ref["cluster_label1"].to_numpy(dtype=object)
                published_cluster2[ref_mask] = ref["cluster_label2"].to_numpy(dtype=object)
                published_cluster3[ref_mask] = ref["cluster_label3"].to_numpy(dtype=object)

        confidence = np.where(
            source_label == "PUBLISHED_ZENODO_HEALTHY_REFERENCE",
            "PUBLISHED",
            np.where(best_probability >= 0.80, "HIGH", np.where(accepted, "MEDIUM", "LOW")),
        )
        adata.layers["counts"] = counts
        adata.X = normalized
        adata.obs["predicted_level1_raw"] = predicted
        adata.obs["human_level1"] = level1.astype(str)
        adata.obs["human_level2"] = level2.astype(str)
        adata.obs["annotation_probability"] = best_probability.astype(np.float32)
        adata.obs["annotation_confidence"] = confidence.astype(str)
        adata.obs["marker_winner"] = marker_winner
        adata.obs["marker_score_max"] = marker_max
        adata.obs["marker_support"] = np.where(marker_concordant, "CONCORDANT", "DISCORDANT")
        adata.obs["annotation_source"] = source_label.astype(str)
        adata.obs["published_lineage"] = published_lineage.astype(str)
        adata.obs["published_cluster_label1"] = published_cluster1.astype(str)
        adata.obs["published_cluster_label2"] = published_cluster2.astype(str)
        adata.obs["published_cluster_label3"] = published_cluster3.astype(str)
        adata.uns["annotation_method"] = "Published Healthy reference labels plus donor-LOO-audited sparse SGD broad-lineage transfer"
        adata.uns["annotation_thresholds"] = {
            "primary_probability": MIN_PROBABILITY,
            "relaxed_probability_with_marker_concordance": RELAXED_PROBABILITY_WITH_MARKER,
        }
        output = ANNOTATED_SAMPLE_DIR / f"{row.GSM}_annotated.h5ad"
        adata.write_h5ad(output, compression="gzip", compression_opts=4)
        output_files.append(output)

        counts_frame = (
            adata.obs.groupby(
                ["GSM", "donor", "disease_group", "state", "primary_contrast_eligible", "human_level1"],
                observed=True,
            )
            .size()
            .rename("n_cells")
            .reset_index()
        )
        count_rows.append(counts_frame)
        for label in PRIMARY_LINEAGES:
            mask = adata.obs["human_level1"].eq(label).to_numpy()
            for marker in MARKERS[label]:
                gene_mask = adata.var["gene_symbol"].astype(str).str.upper().eq(marker).to_numpy()
                if gene_mask.any() and mask.any():
                    value = float(np.asarray(normalized[mask][:, gene_mask].mean()).ravel()[0])
                else:
                    value = float("nan")
                sanity_rows.append({"GSM": row.GSM, "human_level1": label, "marker": marker, "mean_log_normalized_expression": value})
        logger.info(
            "ANNOTATE %d/%d %s cells=%d uncertain=%d published=%d",
            index,
            len(selected),
            row.GSM,
            len(adata),
            int((adata.obs["human_level1"] == "Uncertain").sum()),
            int((adata.obs["annotation_source"] == "PUBLISHED_ZENODO_HEALTHY_REFERENCE").sum()),
        )
        del adata, counts, normalized, probabilities
    return output_files, pd.concat(count_rows, ignore_index=True), pd.DataFrame(sanity_rows)


def mapping_table(zenodo: pd.DataFrame) -> pd.DataFrame:
    published = zenodo[["lineage", "cluster_label1", "cluster_label2", "cluster_label3", "published_human_level1"]].drop_duplicates()
    rows = []
    for row in published.itertuples(index=False):
        rows.append(
            {
                "published_cell_type": row.lineage,
                "published_subtype": row.cluster_label1,
                "human_level1": row.published_human_level1,
                "human_level2": row.cluster_label1,
                "confidence": "PUBLISHED_REFERENCE",
                "marker_support": "Author-published annotation; exact barcode mapping for covered Healthy cells",
                "source": "Zenodo scRNA_cell_metadata.csv / Supplementary Data12",
                "notes": f"cluster_label2={row.cluster_label2}; cluster_label3={row.cluster_label3}",
            }
        )
    for label in MODEL_CLASSES:
        rows.append(
            {
                "published_cell_type": "NA_UNCOVERED_CELL",
                "published_subtype": "NA_BROAD_ONLY",
                "human_level1": label,
                "human_level2": "NA_BROAD_ONLY",
                "confidence": "MODEL_PROBABILITY_AND_MARKER_AUDITED",
                "marker_support": ";".join(MARKERS[label]),
                "source": "Sparse SGD transfer trained on six published Healthy TI references",
                "notes": "Broad lineage only; no forced Dog-Human subtype ontology",
            }
        )
    return pd.DataFrame(rows)


def finalize(
    output_files: list[Path],
    counts: pd.DataFrame,
    sanity: pd.DataFrame,
    cv_metrics: pd.DataFrame,
    logger: logging.Logger,
) -> None:
    if ANNOTATED_OUT.exists():
        raise RuntimeError(f"Refusing to overwrite existing annotated Human object without audit: {ANNOTATED_OUT}")
    logger.info("CONCAT annotated per-sample objects -> %s", ANNOTATED_OUT)
    ad.experimental.concat_on_disk(
        [str(path) for path in output_files],
        str(ANNOTATED_OUT),
        axis=0,
        join="inner",
        merge="same",
        max_loaded_elems=50_000_000,
    )
    atomic_csv(counts, CELL_COUNT_OUT)
    atomic_csv(sanity, MARKER_SANITY_OUT)

    primary = counts.loc[counts["primary_contrast_eligible"].eq("YES_SELECTED")].copy()
    coverage_rows = []
    for lineage in PRIMARY_LINEAGES:
        subset = primary.loc[primary["human_level1"].eq(lineage)]
        by_donor = subset.groupby(["donor", "disease_group"], as_index=False)["n_cells"].sum()
        for group in ("Healthy", "Crohn"):
            values = by_donor.loc[by_donor["disease_group"].eq(group), "n_cells"]
            coverage_rows.append(
                {
                    "human_level1": lineage,
                    "primary_group": "ActiveCrohn" if group == "Crohn" else "Healthy",
                    "donors_with_any_cells": int((values > 0).sum()),
                    "donors_ge20_cells": int((values >= MIN_CELLS_PER_DONOR_LINEAGE).sum()),
                    "min_cells_among_covered_donors": int(values.min()) if len(values) else 0,
                    "total_cells": int(values.sum()),
                }
            )
    coverage = pd.DataFrame(coverage_rows)
    atomic_csv(coverage, DONOR_COVERAGE_OUT)

    cv_lookup = cv_metrics.set_index("label")
    uncertain_fraction = float(counts.loc[counts["human_level1"].eq("Uncertain"), "n_cells"].sum() / counts["n_cells"].sum())
    criteria = [
        ("Annotated object created", ANNOTATED_OUT.is_file(), str(ANNOTATED_OUT)),
        ("Counts layer preserved", True, "Verified during per-sample construction; final backed verification follows"),
        ("Reference donor-LOO Epithelial F1 >=0.85", float(cv_lookup.loc["Epithelial", "f1_score"]) >= 0.85, f"F1={float(cv_lookup.loc['Epithelial', 'f1_score']):.4f}"),
        ("Reference donor-LOO T F1 >=0.85", float(cv_lookup.loc["T", "f1_score"]) >= 0.85, f"F1={float(cv_lookup.loc['T', 'f1_score']):.4f}"),
        ("Reference donor-LOO Myeloid F1 >=0.85", float(cv_lookup.loc["Myeloid", "f1_score"]) >= 0.85, f"F1={float(cv_lookup.loc['Myeloid', 'f1_score']):.4f}"),
        ("Primary Healthy/Active donors have >=20 cells in each primary lineage", bool(((coverage.loc[coverage['primary_group'].eq('Healthy'), 'donors_ge20_cells'] == 16).all()) and ((coverage.loc[coverage['primary_group'].eq('ActiveCrohn'), 'donors_ge20_cells'] == 11).all())), coverage.to_json(orient="records")),
        ("Uncertain fraction <20%", uncertain_fraction < 0.20, f"fraction={uncertain_fraction:.4f}"),
    ]
    combined = ad.read_h5ad(ANNOTATED_OUT, backed="r")
    counts_layer_present = "counts" in combined.layers
    obs_count = combined.n_obs
    var_count = combined.n_vars
    combined.file.close()
    criteria[1] = ("Counts layer preserved", counts_layer_present, f"layers['counts']={counts_layer_present}")
    gate_pass = all(status for _, status, _ in criteria)
    gate = pd.DataFrame(
        [
            {
                "criterion": name,
                "status": "PASS" if status else "FAIL",
                "evidence": evidence,
                "gate": "GREEN_ANNOTATION_GATE" if gate_pass else "YELLOW_ANNOTATION_GATE_NOT_RESOLVED",
            }
            for name, status, evidence in criteria
        ]
    )
    atomic_csv(gate, GATE_OUT)
    if not gate_pass:
        raise RuntimeError("Phase7B annotation gate failed")

    primary_counts = counts.loc[counts["human_level1"].isin(PRIMARY_LINEAGES)].groupby("human_level1")["n_cells"].sum()
    report = f"""# Phase7B Human Broad-Lineage Annotation Report

## Gate

**GREEN — ANNOTATION GATE PASSED**

- Final annotated object: {obs_count:,} cells × {var_count:,} Human-native genes.
- Epithelial cells: {int(primary_counts.get('Epithelial', 0)):,}.
- T cells: {int(primary_counts.get('T', 0)):,}.
- Myeloid cells: {int(primary_counts.get('Myeloid', 0)):,}.
- Uncertain broad-lineage fraction: {uncertain_fraction:.2%}.

## Provenance

Exact Zenodo published labels were retained for covered Healthy reference cells. A sparse multinomial classifier was trained in Human-native gene space on those six Healthy TI references, audited by donor leave-one-out prediction, and transferred only at broad-lineage level. Low-confidence cells were left `Uncertain`. Dog subtype names were not imposed.

## Primary reference CV F1

- Epithelial: {float(cv_lookup.loc['Epithelial', 'f1_score']):.4f}
- T: {float(cv_lookup.loc['T', 'f1_score']):.4f}
- Myeloid: {float(cv_lookup.loc['Myeloid', 'f1_score']):.4f}

Generated UTC: {utc_now()}
"""
    REPORT_OUT.write_text(report, encoding="utf-8", newline="\n")
    logger.info("ANNOTATION GATE PASS cells=%d genes=%d uncertain_fraction=%.4f", obs_count, var_count, uncertain_fraction)


def main() -> int:
    logger = setup_logging()
    ANNOTATION.mkdir(parents=True, exist_ok=True)
    sample_master = pd.read_csv(SAMPLE_MASTER, dtype=str, keep_default_na=False)
    selected = sample_master.loc[sample_master["region"].eq("TI")]
    expected_files = [QC_SAMPLE_DIR / f"{gsm}_qc.h5ad" for gsm in selected["GSM"]]
    missing = [str(path) for path in expected_files if not path.is_file()]
    if missing:
        raise RuntimeError(f"Missing {len(missing)} Phase7 QC per-sample objects")

    gsm_to_library, zenodo = reference_mapping()
    counts, labels, groups, var, reference_detail = load_reference(gsm_to_library, zenodo, logger)
    normalized = normalize_log(counts)
    use_features = feature_mask(counts, var)
    x = normalized[:, use_features]
    cv_metrics, cv_confusion, cv_predicted = cross_validate(x, labels, groups, logger)
    reference_detail["cv_predicted_level1"] = cv_predicted
    atomic_csv(cv_metrics, CV_METRICS_OUT)
    atomic_csv(cv_confusion, CV_CONFUSION_OUT)
    atomic_csv(reference_detail, CV_CELLS_OUT)

    model = new_classifier()
    model.fit(x, labels)
    joblib.dump(
        {
            "model": model,
            "feature_mask": use_features,
            "feature_ids": var.index.astype(str).to_numpy(),
            "feature_symbols": var["gene_symbol"].astype(str).to_numpy(),
            "model_classes": np.asarray(MODEL_CLASSES),
            "training_reference_cells": len(labels),
            "training_reference_libraries": sorted(np.unique(groups).tolist()),
            "normalization": "library-size 10,000 followed by log1p",
            "seed": SEED,
        },
        MODEL_OUT,
        compress=3,
    )
    marker_rows = [
        {"human_level1": lineage, "gene_symbol": symbol, "role": "broad-lineage marker audit"}
        for lineage, symbols in MARKERS.items()
        for symbol in symbols
    ]
    atomic_csv(pd.DataFrame(marker_rows), MARKER_PANEL_OUT)
    atomic_csv(mapping_table(zenodo), MAPPING_OUT)
    output_files, cell_counts, sanity = annotate_samples(
        sample_master, model, use_features, gsm_to_library, zenodo, logger
    )
    finalize(output_files, cell_counts, sanity, cv_metrics, logger)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
