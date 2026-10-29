"""Run preregistered donor-aware QC for GSE254005 without formal cell annotation."""

from __future__ import annotations

import csv
import logging
import os
import platform
import sys
import time
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
import scanpy as sc
import scipy.sparse as sp
import scrublet as scr
import seaborn as sns

CODE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_ROOT))
PROJECT_ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"] )

from config.paths import (
    H5AD_CLEAN_DIR,
    H5AD_RAW_DIR,
    LOG_DIR,
    QC_DIR,
    QC_FIGURES_DIR,
    ensure_phase1_directories,
)

SCRIPT_VERSION = "1.0.0"
RANDOM_SEED = 20240824
EXPECTED_DOUBLET_RATE = 0.06
SIM_DOUBLET_RATIO = 2.0
RAW_H5AD = H5AD_RAW_DIR / "dog_GSE254005_raw.h5ad"
CLEAN_H5AD = H5AD_CLEAN_DIR / "dog_GSE254005_clean.h5ad"
LOG_FILE = LOG_DIR / "phase1_qc.log"

PRE_QC_CSV = QC_DIR / "dog_qc_pre_filter_per_donor.csv"
FILTER_FLOW_CSV = QC_DIR / "dog_qc_filter_flow.csv"
SCRUBLET_CSV = QC_DIR / "dog_scrublet_summary.csv"
MITO_CSV = QC_DIR / "dog_mitochondrial_gene_list.csv"
LINEAGE_RULE_CSV = QC_DIR / "dog_coarse_lineage_qc_rule.csv"
FINAL_DONOR_CSV = QC_DIR / "dog_qc_final_per_donor.csv"

DONOR_ORDER = [
    "GSM8031698_h_1",
    "GSM8031699_h_2",
    "GSM8031700_h_3",
    "GSM8031701_cie_2",
    "GSM8031702_cie_3",
    "GSM8031703_cie_4",
    "GSM8031704_cie_6",
]

EXPECTED_MITO_SYMBOLS = {
    "MT-ND1",
    "MT-ND2",
    "MT-CO1",
    "MT-CO2",
    "MT-ATP8",
    "MT-ATP6",
    "MT-CO3",
    "MT-ND3",
    "ND4L",
    "MT-ND4",
    "MT-ND5",
    "MT-ND6",
    "MT-CYB",
}

EPITHELIAL_MARKERS = ["EPCAM", "KRT18", "KRT20", "VIL1"]
IMMUNE_MARKERS = ["PTPRC"]
STROMAL_VASCULAR_MARKERS = [
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "PECAM1",
    "VWF",
    "KDR",
    "EMCN",
    "RGS5",
    "ACTA2",
]


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)sZ %(levelname)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
        handlers=[logging.FileHandler(LOG_FILE, encoding="utf-8"), logging.StreamHandler()],
        force=True,
    )
    logging.Formatter.converter = time.gmtime


def peak_ram_mb() -> float:
    memory = psutil.Process().memory_info()
    peak_bytes = getattr(memory, "peak_wset", memory.rss)
    return peak_bytes / (1024**2)


def atomic_dataframe(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False, encoding="utf-8")
    os.replace(temporary, path)


def condition_for_donor(obs: pd.DataFrame, donor: str) -> str:
    values = obs.loc[obs["donor_id"].astype(str) == donor, "condition"].astype(str).unique()
    if len(values) != 1:
        raise RuntimeError(f"Donor {donor} has non-unique condition values: {values}")
    return values[0]


def audit_mitochondrial_genes(adata: ad.AnnData) -> np.ndarray:
    symbols = adata.var["gene_symbol"].astype(str)
    mito_mask = symbols.str.startswith("MT-").to_numpy() | symbols.eq("ND4L").to_numpy()
    actual = set(symbols[mito_mask])
    if actual != EXPECTED_MITO_SYMBOLS or int(mito_mask.sum()) != 13:
        raise RuntimeError(
            f"Mitochondrial feature audit failed: expected 13 canonical protein-coding genes, found {len(actual)}: {sorted(actual)}"
        )
    mito = adata.var.loc[mito_mask, ["feature_id", "gene_symbol", "feature_type"]].copy()
    mito.insert(0, "var_name", mito.index.astype(str))
    mito["identification_rule"] = "gene_symbol starts with 'MT-' OR gene_symbol equals 'ND4L'"
    mito["gene_reference"] = "CanFam3.1; Ensembl release 104"
    atomic_dataframe(mito.reset_index(drop=True), MITO_CSV)
    logging.info("MITO_AUDIT rule=^MT-|^ND4L$ matched_genes=%s symbols=%s", len(mito), sorted(actual))
    return mito_mask


def compute_qc_metrics(adata: ad.AnnData, mito_mask: np.ndarray) -> None:
    matrix = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    total = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)
    genes = matrix.getnnz(axis=1).astype(np.int32)
    mito_counts = np.asarray(matrix[:, mito_mask].sum(axis=1)).ravel().astype(np.float64)
    mito_pct = np.divide(mito_counts * 100.0, total, out=np.zeros_like(mito_counts), where=total > 0)
    adata.obs["total_counts"] = total
    adata.obs["n_genes_by_counts"] = genes
    adata.obs["pct_counts_mt"] = mito_pct
    adata.obs["nCount_RNA"] = total
    adata.obs["nFeature_RNA"] = genes
    adata.obs["percent_mt"] = mito_pct


def donor_summary(obs: pd.DataFrame, mask: np.ndarray | None = None) -> pd.DataFrame:
    selected = obs if mask is None else obs.loc[mask]
    rows = []
    for donor in DONOR_ORDER:
        group = selected.loc[selected["donor_id"].astype(str) == donor]
        rows.append(
            {
                "donor": donor,
                "condition": condition_for_donor(obs, donor),
                "raw_cells": int((obs["donor_id"].astype(str) == donor).sum()),
                "selected_cells": len(group),
                "median_counts": float(group["total_counts"].median()),
                "mean_counts": float(group["total_counts"].mean()),
                "median_genes": float(group["n_genes_by_counts"].median()),
                "mean_genes": float(group["n_genes_by_counts"].mean()),
                "median_mt_pct": float(group["pct_counts_mt"].median()),
                "mean_mt_pct": float(group["pct_counts_mt"].mean()),
                "p95_mt_pct": float(group["pct_counts_mt"].quantile(0.95)),
            }
        )
    return pd.DataFrame(rows)


def make_pre_qc_figures(obs: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    palette = dict(zip(DONOR_ORDER, sns.color_palette("tab10", len(DONOR_ORDER))))
    condition_palette = {"Healthy": "#2A9D8F", "CIE": "#E76F51"}

    counts = obs.groupby("donor_id", observed=True).size().reindex(DONOR_ORDER)
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.bar(DONOR_ORDER, counts.values, color=[palette[x] for x in DONOR_ORDER])
    ax.set_ylabel("Raw cells")
    ax.set_title("GSE254005 raw cell count per donor")
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(QC_FIGURES_DIR / "01_cell_counts_per_donor.png", dpi=180)
    plt.close(fig)

    for number, column, ylabel, title in [
        (2, "total_counts", "Total counts", "Raw nCount distribution per donor"),
        (3, "n_genes_by_counts", "Detected features", "Raw nFeature distribution per donor"),
        (4, "pct_counts_mt", "Mitochondrial counts (%)", "Raw mitochondrial percentage per donor"),
    ]:
        fig, ax = plt.subplots(figsize=(12, 5))
        sns.boxenplot(data=obs, x="donor_id", y=column, order=DONOR_ORDER, palette=palette, ax=ax)
        ax.set_title(title)
        ax.set_xlabel("Donor")
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", rotation=35)
        fig.tight_layout()
        fig.savefig(QC_FIGURES_DIR / f"0{number}_{column}_by_donor.png", dpi=180)
        plt.close(fig)

    plot_obs = obs.sample(min(len(obs), 50000), random_state=RANDOM_SEED)
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(
        data=plot_obs,
        x="total_counts",
        y="n_genes_by_counts",
        hue="donor_id",
        hue_order=DONOR_ORDER,
        palette=palette,
        s=6,
        alpha=0.3,
        linewidth=0,
        ax=ax,
    )
    ax.set_xscale("log")
    ax.set_title("Raw nCount vs nFeature (donor-resolved)")
    ax.legend(markerscale=2, fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(QC_FIGURES_DIR / "05_counts_vs_features_by_donor.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(
        data=plot_obs,
        x="total_counts",
        y="pct_counts_mt",
        hue="donor_id",
        hue_order=DONOR_ORDER,
        palette=palette,
        s=6,
        alpha=0.3,
        linewidth=0,
        ax=ax,
    )
    ax.set_xscale("log")
    ax.set_title("Raw mitochondrial percentage vs nCount")
    ax.legend(markerscale=2, fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(QC_FIGURES_DIR / "06_mt_pct_vs_counts_by_donor.png", dpi=180)
    plt.close(fig)

    donor_medians = (
        obs.groupby(["donor_id", "condition"], observed=True)[
            ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
        ]
        .median()
        .reset_index()
    )
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    for ax, column, ylabel in zip(
        axes,
        ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
        ["Median nCount", "Median nFeature", "Median mt%"],
    ):
        sns.stripplot(
            data=donor_medians,
            x="condition",
            y=column,
            hue="donor_id",
            palette=palette,
            order=["Healthy", "CIE"],
            jitter=0.12,
            size=9,
            ax=ax,
        )
        ax.set_ylabel(ylabel)
        ax.get_legend().remove()
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=palette[d], label=d, markersize=7)
        for d in DONOR_ORDER
    ]
    fig.legend(handles=handles, loc="center right", fontsize=7)
    fig.suptitle("Condition grouping with donor-level QC medians retained")
    fig.tight_layout(rect=[0, 0, 0.88, 0.94])
    fig.savefig(QC_FIGURES_DIR / "07_condition_with_donor_medians.png", dpi=180)
    plt.close(fig)


def run_scrublet(adata: ad.AnnData, broad_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    scores = np.full(adata.n_obs, np.nan, dtype=np.float64)
    predicted = np.zeros(adata.n_obs, dtype=bool)
    summaries = []
    plot_data = []
    donors = adata.obs["donor_id"].astype(str).to_numpy()
    conditions = adata.obs["condition"].astype(str).to_numpy()
    for donor in DONOR_ORDER:
        indices = np.flatnonzero(broad_mask & (donors == donor))
        if len(indices) < 500:
            raise RuntimeError(f"Too few cells for donor-specific Scrublet after broad QC: {donor} n={len(indices)}")
        matrix = adata.X[indices].tocsr()
        scrub = scr.Scrublet(
            matrix,
            expected_doublet_rate=EXPECTED_DOUBLET_RATE,
            sim_doublet_ratio=SIM_DOUBLET_RATIO,
            random_state=RANDOM_SEED,
        )
        donor_scores, donor_predicted = scrub.scrub_doublets(
            synthetic_doublet_umi_subsampling=1.0,
            use_approx_neighbors=False,
            distance_metric="euclidean",
            min_counts=3,
            min_cells=3,
            min_gene_variability_pctl=85,
            log_transform=False,
            mean_center=True,
            normalize_variance=True,
            n_prin_comps=30,
            svd_solver="arpack",
            verbose=False,
        )
        if not np.isfinite(scrub.threshold_):
            raise RuntimeError(f"Scrublet did not produce a finite threshold for {donor}")
        scores[indices] = donor_scores
        predicted[indices] = donor_predicted
        condition = np.unique(conditions[indices])
        if len(condition) != 1:
            raise RuntimeError(f"Non-unique condition in Scrublet donor {donor}: {condition}")
        summaries.append(
            {
                "donor": donor,
                "condition": condition[0],
                "input_cells_after_broad_qc": len(indices),
                "expected_doublet_rate": EXPECTED_DOUBLET_RATE,
                "sim_doublet_ratio": SIM_DOUBLET_RATIO,
                "random_seed": RANDOM_SEED,
                "use_approx_neighbors": False,
                "neighbor_backend": "scikit-learn exact nearest neighbors; annoy not used",
                "n_neighbors": int(scrub.n_neighbors),
                "n_prin_comps": 30,
                "min_counts": 3,
                "min_cells": 3,
                "min_gene_variability_pctl": 85,
                "score_threshold": float(scrub.threshold_),
                "predicted_doublets": int(donor_predicted.sum()),
                "predicted_doublet_pct": float(donor_predicted.mean() * 100.0),
                "median_observed_score": float(np.median(donor_scores)),
                "p95_observed_score": float(np.quantile(donor_scores, 0.95)),
            }
        )
        plot_data.append((donor, donor_scores, scrub.doublet_scores_sim_, float(scrub.threshold_)))
        logging.info(
            "SCRUBLET donor=%s input=%s threshold=%.6f predicted=%s rate_pct=%.3f",
            donor,
            len(indices),
            scrub.threshold_,
            int(donor_predicted.sum()),
            donor_predicted.mean() * 100,
        )

    fig, axes = plt.subplots(4, 2, figsize=(11, 13))
    for ax, (donor, observed, simulated, threshold) in zip(axes.ravel(), plot_data):
        ax.hist(observed, bins=50, density=True, alpha=0.65, label="observed")
        ax.hist(simulated, bins=50, density=True, alpha=0.45, label="simulated")
        ax.axvline(threshold, color="black", linestyle="--", linewidth=1, label=f"threshold={threshold:.3f}")
        ax.set_title(donor)
        ax.set_xlabel("Scrublet score")
        ax.set_ylabel("Density")
        ax.legend(fontsize=6)
    axes.ravel()[-1].axis("off")
    fig.suptitle("Per-donor Scrublet score distributions")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(QC_FIGURES_DIR / "08_scrublet_score_distributions.png", dpi=180)
    plt.close(fig)
    return scores, predicted, pd.DataFrame(summaries)


def marker_detected(adata: ad.AnnData, marker: str) -> np.ndarray:
    indices = np.flatnonzero(adata.var["gene_symbol"].astype(str).to_numpy() == marker)
    if len(indices) == 0:
        raise RuntimeError(f"Required coarse-QC marker is absent: {marker}")
    values = adata.X[:, indices]
    if sp.issparse(values):
        return np.asarray((values > 0).sum(axis=1)).ravel() > 0
    return np.asarray(values > 0).sum(axis=1) > 0


def coarse_lineage_for_qc(adata: ad.AnnData, eligible_mask: np.ndarray) -> np.ndarray:
    all_markers = EPITHELIAL_MARKERS + IMMUNE_MARKERS + STROMAL_VASCULAR_MARKERS
    detected = {marker: marker_detected(adata, marker) for marker in all_markers}
    epithelial_count = sum(detected[marker].astype(np.int8) for marker in EPITHELIAL_MARKERS)
    stromal_count = sum(detected[marker].astype(np.int8) for marker in STROMAL_VASCULAR_MARKERS)
    clear_epithelial = eligible_mask & (epithelial_count >= 2)
    clear_non_epithelial = eligible_mask & ~clear_epithelial & (
        detected["PTPRC"] | (stromal_count >= 2)
    )
    labels = np.full(adata.n_obs, "not_evaluated", dtype=object)
    labels[eligible_mask] = "uncertain"
    labels[clear_epithelial] = "epithelial"
    labels[clear_non_epithelial] = "non_epithelial"

    rule_rows = []
    for marker in EPITHELIAL_MARKERS:
        rule_rows.append(
            {
                "qc_class": "epithelial",
                "marker": marker,
                "marker_present": True,
                "rule": "at least 2 of EPCAM/KRT18/KRT20/VIL1 detected",
                "purpose": "coarse lineage split only; not Phase2 annotation",
            }
        )
    for marker in IMMUNE_MARKERS:
        rule_rows.append(
            {
                "qc_class": "non_epithelial",
                "marker": marker,
                "marker_present": True,
                "rule": "PTPRC detected, unless already clear epithelial",
                "purpose": "coarse lineage split only; not Phase2 annotation",
            }
        )
    for marker in STROMAL_VASCULAR_MARKERS:
        rule_rows.append(
            {
                "qc_class": "non_epithelial",
                "marker": marker,
                "marker_present": True,
                "rule": "at least 2 stromal/vascular markers detected, unless already clear epithelial",
                "purpose": "coarse lineage split only; not Phase2 annotation",
            }
        )
    rule_rows.append(
        {
            "qc_class": "uncertain",
            "marker": "NA",
            "marker_present": "NA",
            "rule": "eligible cells satisfying neither conservative rule remain uncertain and are not given 12.5% mt filter",
            "purpose": "avoid forced classification and epithelial cell loss",
        }
    )
    atomic_dataframe(pd.DataFrame(rule_rows), LINEAGE_RULE_CSV)

    eligible_labels = pd.Series(labels[eligible_mask])
    counts = eligible_labels.value_counts()
    logging.info("COARSE_QC_LINEAGE counts=%s", counts.to_dict())
    fig, ax = plt.subplots(figsize=(7, 4.5))
    order = ["epithelial", "non_epithelial", "uncertain"]
    values = [int(counts.get(label, 0)) for label in order]
    ax.bar(order, values, color=["#F4A261", "#457B9D", "#BDBDBD"])
    ax.set_ylabel("Cells after broad QC and doublet removal")
    ax.set_title("Conservative lineage split used only for secondary mt QC")
    fig.tight_layout()
    fig.savefig(QC_FIGURES_DIR / "09_coarse_qc_lineage_counts.png", dpi=180)
    plt.close(fig)
    return labels


def filter_flow(
    obs: pd.DataFrame,
    count_mask: np.ndarray,
    feature_mask: np.ndarray,
    broad_mt_mask: np.ndarray,
    predicted_doublet: np.ndarray,
    secondary_mask: np.ndarray,
) -> pd.DataFrame:
    rows = []
    donor_values = obs["donor_id"].astype(str).to_numpy()
    for donor in DONOR_ORDER:
        donor_mask = donor_values == donor
        condition = condition_for_donor(obs, donor)
        current = donor_mask.copy()
        rows.append(
            {
                "step": "raw",
                "donor": donor,
                "condition": condition,
                "cells_before": int(current.sum()),
                "cells_removed": 0,
                "cells_after": int(current.sum()),
                "percent_removed": 0.0,
                "filter_reason": "starting processed matrix cells",
            }
        )
        steps = [
            ("count filter", count_mask, "100 < nCount_RNA < 25000"),
            ("feature filter", feature_mask, "200 < nFeature_RNA < 3000"),
            ("broad mt filter", broad_mt_mask, "percent.mt < 80%"),
            ("doublet", ~predicted_doublet, "per-donor Scrublet predicted_doublet=False"),
            (
                "secondary non-epithelial mt filter",
                secondary_mask,
                "percent.mt < 12.5% only for conservatively classified non-epithelial cells",
            ),
        ]
        for step, criterion, reason in steps:
            before = int(current.sum())
            current = current & criterion
            after = int(current.sum())
            removed = before - after
            rows.append(
                {
                    "step": step,
                    "donor": donor,
                    "condition": condition,
                    "cells_before": before,
                    "cells_removed": removed,
                    "cells_after": after,
                    "percent_removed": removed / before * 100.0 if before else 0.0,
                    "filter_reason": reason,
                }
            )
    return pd.DataFrame(rows)


def main() -> int:
    ensure_phase1_directories()
    configure_logging()
    os.environ.setdefault("LOKY_MAX_CPU_COUNT", str(max(1, os.cpu_count() or 1)))
    started = time.perf_counter()
    logging.info(
        "START script_version=%s python=%s random_seed=%s scanpy=%s anndata=%s scrublet=%s",
        SCRIPT_VERSION,
        platform.python_version(),
        RANDOM_SEED,
        version("scanpy"),
        version("anndata"),
        version("scrublet"),
    )
    logging.info("INPUT raw_h5ad=%s", RAW_H5AD)
    try:
        if CLEAN_H5AD.exists():
            raise FileExistsError(f"Clean h5ad will not be overwritten: {CLEAN_H5AD}")
        adata = ad.read_h5ad(RAW_H5AD)
        if adata.shape != (49005, 20257):
            raise RuntimeError(f"Unexpected raw h5ad shape: {adata.shape}")
        if set(adata.obs["donor_id"].astype(str)) != set(DONOR_ORDER):
            raise RuntimeError("Raw h5ad donor set does not match the seven audited donors")

        mito_mask = audit_mitochondrial_genes(adata)
        compute_qc_metrics(adata, mito_mask)
        pre = donor_summary(adata.obs).drop(columns="selected_cells")
        atomic_dataframe(pre, PRE_QC_CSV)
        make_pre_qc_figures(adata.obs)

        count_mask = ((adata.obs["total_counts"] > 100) & (adata.obs["total_counts"] < 25000)).to_numpy()
        feature_mask = (
            (adata.obs["n_genes_by_counts"] > 200) & (adata.obs["n_genes_by_counts"] < 3000)
        ).to_numpy()
        broad_mt_mask = (adata.obs["pct_counts_mt"] < 80.0).to_numpy()
        broad_mask = count_mask & feature_mask & broad_mt_mask
        logging.info(
            "BROAD_QC raw=%s after=%s removed=%s",
            adata.n_obs,
            int(broad_mask.sum()),
            int((~broad_mask).sum()),
        )

        scrublet_scores, predicted_doublet, scrublet_summary = run_scrublet(adata, broad_mask)
        atomic_dataframe(scrublet_summary, SCRUBLET_CSV)
        adata.obs["scrublet_score"] = scrublet_scores
        adata.obs["predicted_doublet"] = predicted_doublet
        adata.obs["scrublet_evaluated"] = broad_mask

        after_doublet = broad_mask & ~predicted_doublet
        lineage_labels = coarse_lineage_for_qc(adata, after_doublet)
        adata.obs["coarse_lineage_qc"] = pd.Categorical(
            lineage_labels,
            categories=["epithelial", "non_epithelial", "uncertain", "not_evaluated"],
        )
        secondary_mask = (lineage_labels != "non_epithelial") | (
            adata.obs["pct_counts_mt"].to_numpy() < 12.5
        )
        final_mask = after_doublet & secondary_mask
        flow = filter_flow(
            adata.obs,
            count_mask,
            feature_mask,
            broad_mt_mask,
            predicted_doublet,
            secondary_mask,
        )
        atomic_dataframe(flow, FILTER_FLOW_CSV)

        final = donor_summary(adata.obs, final_mask)
        final["removed_cells"] = final["raw_cells"] - final["selected_cells"]
        final["percent_removed"] = final["removed_cells"] / final["raw_cells"] * 100.0
        final = final.rename(columns={"selected_cells": "post_qc_cells"})
        atomic_dataframe(final, FINAL_DONOR_CSV)

        fig, ax = plt.subplots(figsize=(12, 5))
        positions = np.arange(len(DONOR_ORDER))
        raw_values = final.set_index("donor").loc[DONOR_ORDER, "raw_cells"].to_numpy()
        final_values = final.set_index("donor").loc[DONOR_ORDER, "post_qc_cells"].to_numpy()
        ax.bar(positions - 0.2, raw_values, width=0.4, label="raw")
        ax.bar(positions + 0.2, final_values, width=0.4, label="post-QC")
        ax.set_xticks(positions, DONOR_ORDER, rotation=35, ha="right")
        ax.set_ylabel("Cells")
        ax.set_title("Raw and final cells per donor")
        ax.legend()
        fig.tight_layout()
        fig.savefig(QC_FIGURES_DIR / "10_raw_vs_final_cells_per_donor.png", dpi=180)
        plt.close(fig)

        clean = adata[final_mask].copy()
        clean.layers["counts"] = clean.X.copy()
        clean.uns["matrix_state"] = "X=log1p(normalize_total target_sum=10000); layers['counts']=raw processed counts"
        clean.uns["normalization"] = {
            "method": "scanpy.pp.normalize_total followed by scanpy.pp.log1p",
            "target_sum": 10000.0,
            "performed_after_qc": True,
            "scaled_or_regressed": False,
        }
        clean.uns["qc_thresholds"] = {
            "count_min_strict": 100,
            "count_max_strict": 25000,
            "feature_min_strict": 200,
            "feature_max_strict": 3000,
            "broad_mt_pct_max_strict": 80.0,
            "non_epithelial_mt_pct_max_strict": 12.5,
        }
        clean.uns["scrublet"] = {
            "per_donor": True,
            "expected_doublet_rate": EXPECTED_DOUBLET_RATE,
            "sim_doublet_ratio": SIM_DOUBLET_RATIO,
            "random_seed": RANDOM_SEED,
            "use_approx_neighbors": False,
            "neighbor_backend": "scikit-learn exact nearest neighbors; annoy not installed or used",
        }
        clean.uns["coarse_lineage_warning"] = (
            "epithelial/non_epithelial/uncertain is a conservative Phase1 QC split only, not formal cell annotation"
        )
        sc.pp.normalize_total(clean, target_sum=10000.0)
        sc.pp.log1p(clean)

        temporary_h5ad = CLEAN_H5AD.with_suffix(".h5ad.tmp")
        clean.write_h5ad(temporary_h5ad, compression="gzip")
        check = ad.read_h5ad(temporary_h5ad, backed="r")
        required_obs = {
            "sample_id",
            "donor_id",
            "condition",
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
            "scrublet_score",
            "predicted_doublet",
            "coarse_lineage_qc",
        }
        if check.shape != clean.shape or "counts" not in check.layers or not required_obs.issubset(check.obs.columns):
            raise RuntimeError("Clean h5ad reload validation failed")
        check.file.close()
        os.replace(temporary_h5ad, CLEAN_H5AD)

        logging.info(
            "OUTPUT clean_h5ad=%s raw_cells=%s final_cells=%s healthy_final=%s cie_final=%s",
            CLEAN_H5AD,
            adata.n_obs,
            clean.n_obs,
            int((clean.obs["condition"].astype(str) == "Healthy").sum()),
            int((clean.obs["condition"].astype(str) == "CIE").sum()),
        )
        logging.info(
            "END status=success elapsed_seconds=%.3f peak_ram_mb=%.1f",
            time.perf_counter() - started,
            peak_ram_mb(),
        )
        return 0
    except Exception:
        logging.exception(
            "END status=error elapsed_seconds=%.3f peak_ram_mb=%.1f",
            time.perf_counter() - started,
            peak_ram_mb(),
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
