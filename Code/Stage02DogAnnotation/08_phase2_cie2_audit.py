"""Audit GSM8031701_cie_2 and rebuild a no-cie_2 sensitivity embedding.

This script is descriptive only.  It performs no CIE-vs-Healthy test, DEG,
pseudobulk, enrichment, integration, or donor regression.  The annotated main
object is read-only; cie_2 is excluded only from the separate sensitivity run.
"""

from __future__ import annotations

import gc
import logging
import os
import platform
import sys
import time
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

CODE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_ROOT))
PROJECT_ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"] )

from config.paths import LOG_DIR  # noqa: E402
from config.phase2_config import (  # noqa: E402
    ANNOTATED_H5AD,
    ANNOTATION_FIGURES_DIR,
    ANNOTATIONS_DIR,
    N_HVG,
    N_NEIGHBORS,
    N_PCA_COMPUTE,
    N_PCS_NEIGHBORS,
    PHASE2_SEED,
    PHASE2_WORK_DIR,
    UMAP_MIN_DIST,
    ensure_phase2_directories,
)
from Stage02DogAnnotation.phase2_utils import (  # noqa: E402
    EPITHELIAL_LEVEL2_MARKERS,
    LEVEL1_MARKERS,
    MYELOID_LEVEL2_MARKERS,
    T_LEVEL2_MARKERS,
    classify_groups,
    configure_logger,
    group_signature_scores,
    marker_var_names,
    peak_ram_mb,
    validate_phase1_clean,
)

SCRIPT_VERSION = "1.0.0"
CIE2 = "GSM8031701_cie_2"
SENSITIVITY_CLUSTER_KEY = "leiden_no_cie2_0_4"


def atomic_csv(frame: pd.DataFrame, path: Path, **kwargs: object) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temp, index=False, **kwargs)
    os.replace(temp, path)


def expression_summary(adata: ad.AnnData, genes: list[str]) -> tuple[np.ndarray, np.ndarray, int]:
    variables = marker_var_names(adata, genes)
    if not variables:
        zeros = np.zeros(adata.n_obs, dtype=np.float32)
        return zeros, zeros, 0
    matrix = adata[:, variables].X
    if sparse.issparse(matrix):
        mean = np.asarray(matrix.mean(axis=1)).ravel()
        detected = np.asarray((matrix > 0).mean(axis=1)).ravel()
    else:
        mean = np.asarray(matrix).mean(axis=1)
        detected = (np.asarray(matrix) > 0).mean(axis=1)
    return mean.astype(np.float32), detected.astype(np.float32), len(variables)


def add_level1_coherence(adata: ad.AnnData) -> None:
    adata.obs["level1_marker_mean"] = np.nan
    adata.obs["level1_marker_detection"] = np.nan
    adata.obs["level1_marker_n_available"] = 0
    for level1, markers in LEVEL1_MARKERS.items():
        mask = adata.obs["level1"].astype(str).to_numpy() == level1
        if not mask.any():
            continue
        mean, detected, n_available = expression_summary(adata[mask], markers)
        adata.obs.loc[mask, "level1_marker_mean"] = mean
        adata.obs.loc[mask, "level1_marker_detection"] = detected
        adata.obs.loc[mask, "level1_marker_n_available"] = n_available


def quantile_summary(values: pd.Series) -> dict[str, float]:
    clean = pd.to_numeric(values, errors="coerce").dropna()
    return {
        "q25": float(clean.quantile(0.25)),
        "median": float(clean.median()),
        "q75": float(clean.quantile(0.75)),
    }


def donor_level1_qc(adata: ad.AnnData) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    obs = adata.obs.copy()
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["level1"] = obs["level1"].astype(str)
    for (donor, level1), group in obs.groupby(["donor_id", "level1"], observed=True):
        mt = quantile_summary(group["pct_counts_mt"])
        features = quantile_summary(group["n_genes_by_counts"])
        counts = quantile_summary(group["total_counts"])
        rows.append(
            {
                "donor": donor,
                "condition": str(group["condition"].iloc[0]),
                "level1": level1,
                "n_cells": len(group),
                "fraction_of_donor": len(group) / (obs["donor_id"] == donor).sum(),
                "mt_q25": mt["q25"],
                "mt_median": mt["median"],
                "mt_q75": mt["q75"],
                "nFeature_q25": features["q25"],
                "nFeature_median": features["median"],
                "nFeature_q75": features["q75"],
                "nCount_q25": counts["q25"],
                "nCount_median": counts["median"],
                "nCount_q75": counts["q75"],
                "marker_mean_median": float(group["level1_marker_mean"].median()),
                "marker_detection_median": float(group["level1_marker_detection"].median()),
                "scrublet_score_median": float(group["scrublet_score"].median()),
                "predicted_doublet_fraction": float(group["predicted_doublet"].astype(bool).mean()),
                "pct_mt_ge_70_fraction": float((group["pct_counts_mt"] >= 70).mean()),
                "pct_mt_ge_80_fraction": float((group["pct_counts_mt"] >= 80).mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(["level1", "donor"]).reset_index(drop=True)


def quality_thresholds_and_flags(adata: ad.AnnData) -> pd.DataFrame:
    donors = adata.obs["donor_id"].astype(str)
    levels = adata.obs["level1"].astype(str)
    flag = np.zeros(adata.n_obs, dtype=bool)
    rows: list[dict[str, object]] = []
    for level1 in sorted(levels.unique()):
        reference = (levels == level1) & (donors != CIE2)
        target = levels == level1
        if not reference.any():
            continue
        mt_q90 = float(adata.obs.loc[reference, "pct_counts_mt"].quantile(0.90))
        feature_q10 = float(adata.obs.loc[reference, "n_genes_by_counts"].quantile(0.10))
        mt_threshold = max(70.0, mt_q90)
        flag[target] = (
            (adata.obs.loc[target, "pct_counts_mt"].to_numpy() >= mt_threshold)
            & (adata.obs.loc[target, "n_genes_by_counts"].to_numpy() <= feature_q10)
        )
        rows.append(
            {
                "level1": level1,
                "other_donors_mt_q90": mt_q90,
                "operational_mt_threshold": mt_threshold,
                "other_donors_nFeature_q10": feature_q10,
                "definition": "pct_mt >= max(70, other-donor lineage q90) AND nFeature <= other-donor lineage q10",
            }
        )
    adata.obs["operational_low_quality_like"] = flag
    return pd.DataFrame(rows)


def cluster_audit(adata: ad.AnnData) -> pd.DataFrame:
    obs = adata.obs.copy()
    obs["annotation_cluster"] = obs["annotation_cluster"].astype(str)
    obs["donor_id"] = obs["donor_id"].astype(str)
    rows: list[dict[str, object]] = []
    for cluster, group in obs.groupby("annotation_cluster", observed=True):
        donor_counts = group["donor_id"].value_counts()
        cie = group[group["donor_id"] == CIE2]
        non = group[group["donor_id"] != CIE2]
        cie_fraction = len(cie) / len(group)
        rows.append(
            {
                "cluster": cluster,
                "level1": str(group["level1"].mode().iloc[0]),
                "level2": str(group["level2"].mode().iloc[0]),
                "annotation_confidence": str(group["annotation_confidence"].mode().iloc[0]),
                "n_cells": len(group),
                "n_donors": group["donor_id"].nunique(),
                "dominant_donor": str(donor_counts.index[0]),
                "dominant_donor_fraction": float(donor_counts.iloc[0] / len(group)),
                "cie2_n_cells": len(cie),
                "cie2_fraction_of_cluster": cie_fraction,
                "cie2_dominant_ge_0_80": bool(cie_fraction >= 0.80),
                "cie2_mt_median": float(cie["pct_counts_mt"].median()) if len(cie) else np.nan,
                "other_mt_median": float(non["pct_counts_mt"].median()) if len(non) else np.nan,
                "cie2_nFeature_median": float(cie["n_genes_by_counts"].median()) if len(cie) else np.nan,
                "other_nFeature_median": float(non["n_genes_by_counts"].median()) if len(non) else np.nan,
                "cie2_nCount_median": float(cie["total_counts"].median()) if len(cie) else np.nan,
                "cie2_marker_detection_median": float(cie["level1_marker_detection"].median()) if len(cie) else np.nan,
                "other_marker_detection_median": float(non["level1_marker_detection"].median()) if len(non) else np.nan,
                "cie2_stress_score_median": float(cie["epithelial_stress_score"].median()) if len(cie) else np.nan,
                "cie2_doublet_fraction": float(cie["predicted_doublet"].astype(bool).mean()) if len(cie) else np.nan,
                "cie2_low_quality_like_fraction": float(cie["operational_low_quality_like"].mean()) if len(cie) else np.nan,
                "biological_plausibility": (
                    "Low/mixed annotation"
                    if ("Uncertain" in str(group["level2"].mode().iloc[0]) or str(group["annotation_confidence"].mode().iloc[0]) in {"Low", "Uncertain"})
                    else "Supported by retained lineage/state annotation"
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("cie2_fraction_of_cluster", ascending=False).reset_index(drop=True)


def choose_hvgs(adata: ad.AnnData, logger: logging.Logger) -> tuple[pd.Index, str]:
    last_error: Exception | None = None
    method = ""
    for span in (0.3, 0.5, 0.7, 1.0):
        try:
            sc.pp.highly_variable_genes(
                adata,
                flavor="seurat_v3",
                n_top_genes=min(N_HVG + 400, adata.n_vars),
                layer="counts",
                batch_key="donor_id",
                span=span,
                subset=False,
                check_values=True,
            )
            method = f"seurat_v3_donor_span_{span:.1f}"
            break
        except ValueError as exc:
            last_error = exc
            logger.warning("HVG_LOESS_RETRY span=%.1f error=%s", span, exc)
    if not method:
        raise RuntimeError("All fixed donor-aware HVG spans failed") from last_error
    symbols = adata.var["gene_symbol"].astype(str)
    mt = symbols.str.startswith("MT-") | symbols.eq("ND4L")
    ranked = adata.var.loc[~mt & adata.var["highly_variable_rank"].notna()].sort_values(
        ["highly_variable_rank", "highly_variable_nbatches"], ascending=[True, False]
    )
    return ranked.index[: min(N_HVG, len(ranked))], method


def top_cluster_markers(adata: ad.AnnData, key: str, n_top: int = 15) -> pd.DataFrame:
    groups = adata.obs[key].astype(str).to_numpy()
    matrix = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
    total = np.asarray(matrix.sum(axis=0)).ravel()
    symbols = adata.var["gene_symbol"].astype(str).to_numpy()
    keep = np.array([not (x.startswith("MT-") or x.startswith("RPL") or x.startswith("RPS")) for x in symbols])
    rows: list[dict[str, object]] = []
    for group in sorted(set(groups), key=int):
        mask = groups == group
        inside_sum = np.asarray(matrix[mask].sum(axis=0)).ravel()
        difference = inside_sum / mask.sum() - (total - inside_sum) / (~mask).sum()
        candidates = np.flatnonzero(keep)
        indices = candidates[np.argsort(difference[candidates])[-n_top:][::-1]]
        rows.append({"sensitivity_cluster": group, "top_markers": ";".join(symbols[indices])})
    return pd.DataFrame(rows)


def sensitivity_embedding(adata: ad.AnnData, logger: logging.Logger) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mask = adata.obs["donor_id"].astype(str).to_numpy() != CIE2
    sens = adata[mask].copy()
    logger.info("SENSITIVITY_START cells=%d donors=%d excluded=%d", sens.n_obs, sens.obs["donor_id"].nunique(), int((~mask).sum()))
    chosen, hvg_method = choose_hvgs(sens, logger)
    work = sens[:, chosen].copy()
    work.X = work.layers["lognorm"].copy()
    for layer in list(work.layers.keys()):
        del work.layers[layer]
    sc.pp.scale(work, max_value=10)
    sc.tl.pca(work, n_comps=N_PCA_COMPUTE, svd_solver="arpack", random_state=PHASE2_SEED)
    sens.obsm["X_pca_no_cie2"] = np.asarray(work.obsm["X_pca"], dtype=np.float32)
    variance = np.asarray(work.uns["pca"]["variance_ratio"], dtype=float)
    del work
    gc.collect()
    sc.pp.neighbors(
        sens,
        n_neighbors=N_NEIGHBORS,
        n_pcs=N_PCS_NEIGHBORS,
        use_rep="X_pca_no_cie2",
        metric="cosine",
        random_state=PHASE2_SEED,
        key_added="neighbors_no_cie2",
    )
    sc.tl.umap(sens, neighbors_key="neighbors_no_cie2", min_dist=UMAP_MIN_DIST, random_state=PHASE2_SEED)
    sc.tl.leiden(
        sens,
        resolution=0.4,
        random_state=PHASE2_SEED,
        key_added=SENSITIVITY_CLUSTER_KEY,
        neighbors_key="neighbors_no_cie2",
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    score, raw, detection = group_signature_scores(sens, SENSITIVITY_CLUSTER_KEY, LEVEL1_MARKERS)
    predicted = classify_groups(score, raw, detection)
    cross = pd.crosstab(sens.obs[SENSITIVITY_CLUSTER_KEY].astype(str), sens.obs["level1"].astype(str))
    majority = cross.idxmax(axis=1)
    purity = cross.max(axis=1) / cross.sum(axis=1)
    marker_preview = top_cluster_markers(sens, SENSITIVITY_CLUSTER_KEY)
    cluster_table = predicted.reset_index().rename(columns={"cluster": "sensitivity_cluster"})
    cluster_table["reference_majority_level1"] = cluster_table["sensitivity_cluster"].map(majority)
    cluster_table["reference_level1_purity"] = cluster_table["sensitivity_cluster"].map(purity)
    cluster_table["n_cells"] = cluster_table["sensitivity_cluster"].map(cross.sum(axis=1))
    cluster_table["n_donors"] = cluster_table["sensitivity_cluster"].map(
        sens.obs.groupby(SENSITIVITY_CLUSTER_KEY, observed=True)["donor_id"].nunique().astype(int)
    )
    cluster_table = cluster_table.merge(marker_preview, on="sensitivity_cluster", how="left")
    state_sets = {
        "IL7R-high T": ["IL7R-high T"],
        "GZMA-high T": ["GZMA-high T"],
        "Myeloid major states": ["Neutrophil", "Eosinophil", "Macrophage", "Monocyte", "cDC1", "IL22RA2 DC", "Cycling macrophage"],
        "Epithelial major states": ["Enterocyte 1", "Enterocyte 1 (IEG-high)", "Enterocyte 2", "Enterocyte 3", "BEST4+ enterocyte", "Goblet", "Tuft", "IFN enterocyte", "Enteroendocrine", "LTF/KRT14 epithelial"],
    }
    state_rows: list[dict[str, object]] = []
    for state, labels in state_sets.items():
        state_mask = sens.obs["level2"].astype(str).isin(labels)
        cluster_counts = sens.obs.loc[state_mask, SENSITIVITY_CLUSTER_KEY].astype(str).value_counts()
        state_rows.append(
            {
                "state_group": state,
                "retained_cells": int(state_mask.sum()),
                "retained_donors": int(sens.obs.loc[state_mask, "donor_id"].nunique()),
                "sensitivity_clusters": ";".join(cluster_counts.index.tolist()),
                "largest_cluster_fraction": float(cluster_counts.iloc[0] / state_mask.sum()) if state_mask.any() else np.nan,
                "status": "Retained" if state_mask.any() and sens.obs.loc[state_mask, "donor_id"].nunique() >= 2 else "Weak/not retained",
                "interpretation": "Full annotation projected onto independently rebuilt no-cie2 graph; not a disease test.",
            }
        )
    state_table = pd.DataFrame(state_rows)
    ari = adjusted_rand_score(sens.obs["leiden_0_4"].astype(str), sens.obs[SENSITIVITY_CLUSTER_KEY].astype(str))
    nmi = normalized_mutual_info_score(sens.obs["leiden_0_4"].astype(str), sens.obs[SENSITIVITY_CLUSTER_KEY].astype(str))
    weighted_purity = float(np.average(purity, weights=cross.sum(axis=1)))
    metrics = pd.DataFrame(
        [
            {
                "excluded_donor": CIE2,
                "retained_cells": sens.n_obs,
                "retained_donors": sens.obs["donor_id"].nunique(),
                "n_hvg": len(chosen),
                "hvg_method": hvg_method,
                "n_pcs_computed": N_PCA_COMPUTE,
                "n_pcs_neighbors": N_PCS_NEIGHBORS,
                "n_neighbors": N_NEIGHBORS,
                "umap_min_dist": UMAP_MIN_DIST,
                "leiden_resolution": 0.4,
                "random_seed": PHASE2_SEED,
                "sensitivity_clusters": sens.obs[SENSITIVITY_CLUSTER_KEY].nunique(),
                "full_vs_no_cie2_ARI_common_cells": ari,
                "full_vs_no_cie2_NMI_common_cells": nmi,
                "weighted_reference_level1_purity": weighted_purity,
                "cumulative_pca_variance": float(variance.sum()),
                "forced_integration": False,
                "disease_contrast": False,
            }
        ]
    )
    logger.info("SENSITIVITY_END clusters=%d ARI=%.4f NMI=%.4f level1_purity=%.4f", sens.obs[SENSITIVITY_CLUSTER_KEY].nunique(), ari, nmi, weighted_purity)
    return sens, cluster_table, state_table, metrics


def figures(adata: ad.AnnData, sens: ad.AnnData, cluster_table: pd.DataFrame) -> None:
    obs = adata.obs.copy()
    obs["cie2_status"] = np.where(obs["donor_id"].astype(str) == CIE2, "cie_2", "other 6 donors")
    order = [x for x in ["Epithelial", "T", "Myeloid", "B", "Plasma", "Mast", "Stromal", "Uncertain"] if x in set(obs["level1"].astype(str))]
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    data = [[obs.loc[(obs["level1"].astype(str) == level) & (obs["cie2_status"] == status), "pct_counts_mt"].to_numpy() for status in ["cie_2", "other 6 donors"]] for level in order]
    positions = np.arange(len(order))
    for index, status in enumerate(["cie_2", "other 6 donors"]):
        values = [item[index] for item in data]
        axes[0].boxplot(values, positions=positions + (index - 0.5) * 0.28, widths=0.24, showfliers=False, patch_artist=True, boxprops={"facecolor": ["#C4473A", "#4C78A8"][index], "alpha": 0.7})
    axes[0].set_xticks(positions, order, rotation=35, ha="right")
    axes[0].set_ylabel("pct_counts_mt")
    axes[0].set_title("cie_2 vs other donors within Level1")
    axes[0].legend(handles=[plt.Line2D([0], [0], color="#C4473A", lw=8, label="cie_2"), plt.Line2D([0], [0], color="#4C78A8", lw=8, label="other 6 donors")])
    cie = obs[obs["cie2_status"] == "cie_2"]
    axes[1].scatter(cie["n_genes_by_counts"], cie["pct_counts_mt"], s=8, alpha=0.35, c=np.where(cie["operational_low_quality_like"], "#C4473A", "#59A14F"))
    axes[1].set_xlabel("n_genes_by_counts")
    axes[1].set_ylabel("pct_counts_mt")
    axes[1].set_title("cie_2 complexity and operational low-quality-like flag")
    fig.tight_layout()
    fig.savefig(ANNOTATION_FIGURES_DIR / "14_cie2_mt_and_complexity_audit.png", dpi=180)
    plt.close(fig)

    top = cluster_table.sort_values("cie2_fraction_of_cluster", ascending=False).head(20)
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = np.where(top["cie2_dominant_ge_0_80"], "#C4473A", "#4C78A8")
    ax.barh(top["cluster"].astype(str), top["cie2_fraction_of_cluster"], color=colors)
    ax.axvline(0.8, ls="--", color="black", lw=1)
    ax.set_xlim(0, 1)
    ax.invert_yaxis()
    ax.set_xlabel("fraction of annotation cluster from cie_2")
    ax.set_title("cie_2 contribution to annotation clusters")
    fig.tight_layout()
    fig.savefig(ANNOTATION_FIGURES_DIR / "15_cie2_cluster_contribution.png", dpi=180)
    plt.close(fig)

    sens.obs["reference_level1"] = sens.obs["level1"].astype(str).astype("category")
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    sc.pl.umap(sens, color="reference_level1", ax=axes[0], show=False, title="No-cie_2 graph: retained Level1 labels", frameon=False)
    sc.pl.umap(sens, color=SENSITIVITY_CLUSTER_KEY, ax=axes[1], show=False, title="No-cie_2 de novo Leiden 0.4", frameon=False, legend_loc="on data")
    fig.tight_layout()
    fig.savefig(ANNOTATION_FIGURES_DIR / "16_no_cie2_sensitivity_umap.png", dpi=180)
    plt.close(fig)


def main() -> None:
    start = time.time()
    ensure_phase2_directories()
    logger = configure_logger(LOG_DIR / "phase2_cie2_audit.log", "phase2_cie2_audit")
    logger.info("START script_version=%s python=%s scanpy=%s seed=%d input=%s", SCRIPT_VERSION, platform.python_version(), sc.__version__, PHASE2_SEED, ANNOTATED_H5AD)
    try:
        adata = ad.read_h5ad(ANNOTATED_H5AD)
        validate_phase1_clean(adata)
        if CIE2 not in set(adata.obs["donor_id"].astype(str)):
            raise RuntimeError(f"Required donor missing: {CIE2}")
        add_level1_coherence(adata)
        thresholds = quality_thresholds_and_flags(adata)
        qc_table = donor_level1_qc(adata)
        cluster_table = cluster_audit(adata)
        sens, sens_clusters, state_table, metrics = sensitivity_embedding(adata, logger)
        overview = pd.DataFrame(
            [
                {"metric": "cie2_cells", "value": int((adata.obs["donor_id"].astype(str) == CIE2).sum()), "interpretation": "Main object retained; sensitivity only excludes a copy."},
                {"metric": "cie2_overall_mt_median", "value": float(adata.obs.loc[adata.obs["donor_id"].astype(str) == CIE2, "pct_counts_mt"].median()), "interpretation": "Post-Phase1-QC median."},
                {"metric": "cie2_epithelial_fraction", "value": float((adata.obs.loc[adata.obs["donor_id"].astype(str) == CIE2, "level1"].astype(str) == "Epithelial").mean()), "interpretation": "Fraction of cie2 cells annotated Epithelial."},
                {"metric": "cie2_low_quality_like_fraction", "value": float(adata.obs.loc[adata.obs["donor_id"].astype(str) == CIE2, "operational_low_quality_like"].mean()), "interpretation": "Joint high-mt and lineage-relative low-complexity flag; descriptive only."},
                {"metric": "cie2_dominant_annotation_clusters_ge_0_80", "value": int(cluster_table["cie2_dominant_ge_0_80"].sum()), "interpretation": "Requires biological/technical review; not auto-removed."},
                {"metric": "no_cie2_sensitivity_ARI", "value": float(metrics["full_vs_no_cie2_ARI_common_cells"].iloc[0]), "interpretation": "Global cluster agreement on common cells; label permutation-adjusted."},
                {"metric": "no_cie2_sensitivity_NMI", "value": float(metrics["full_vs_no_cie2_NMI_common_cells"].iloc[0]), "interpretation": "Global cluster information agreement on common cells."},
                {"metric": "no_cie2_weighted_level1_purity", "value": float(metrics["weighted_reference_level1_purity"].iloc[0]), "interpretation": "Retained formal Level1 labels within no-cie2 clusters."},
            ]
        )
        atomic_csv(qc_table, ANNOTATIONS_DIR / "cie2_level1_qc_comparison.csv")
        atomic_csv(thresholds, ANNOTATIONS_DIR / "cie2_low_quality_thresholds.csv")
        atomic_csv(cluster_table, ANNOTATIONS_DIR / "cie2_cluster_audit.csv")
        atomic_csv(sens_clusters, ANNOTATIONS_DIR / "cie2_exclusion_cluster_annotation.csv")
        atomic_csv(state_table, ANNOTATIONS_DIR / "cie2_exclusion_state_stability.csv")
        atomic_csv(metrics, ANNOTATIONS_DIR / "cie2_exclusion_embedding_metrics.csv")
        atomic_csv(overview, ANNOTATIONS_DIR / "cie2_sensitivity_summary.csv")
        cell_table = sens.obs[["donor_id", "condition", "level1", "level2", "annotation_cluster", SENSITIVITY_CLUSTER_KEY]].copy()
        cell_table.insert(0, "cell_id", cell_table.index.astype(str))
        cell_table["sensitivity_umap_1"] = sens.obsm["X_umap"][:, 0]
        cell_table["sensitivity_umap_2"] = sens.obsm["X_umap"][:, 1]
        atomic_csv(cell_table, PHASE2_WORK_DIR / "cie2_exclusion_sensitivity_cells.csv.gz", compression="gzip")
        figures(adata, sens, cluster_table)
        logger.info("OUTPUT summary=%s dominant_clusters=%d low_quality_like_fraction=%.4f", ANNOTATIONS_DIR / "cie2_sensitivity_summary.csv", int(cluster_table["cie2_dominant_ge_0_80"].sum()), float(overview.loc[overview["metric"] == "cie2_low_quality_like_fraction", "value"].iloc[0]))
        logger.info("END status=success elapsed_seconds=%.3f peak_ram_mb=%.1f", time.time() - start, peak_ram_mb())
    except Exception:
        logger.exception("END status=error elapsed_seconds=%.3f peak_ram_mb=%.1f", time.time() - start, peak_ram_mb())
        raise


if __name__ == "__main__":
    main()
