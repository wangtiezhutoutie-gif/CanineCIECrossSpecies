"""Independently annotate the canine atlas using paper/canine markers.

Use ``--preview`` to inspect deterministic lineage clusters and marker evidence without
writing any formal Phase 2 annotation deliverables.
"""

from __future__ import annotations

import argparse
import gc
import logging
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
import scanpy as sc
from scipy import sparse

CODE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_ROOT))
PROJECT_ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"] )

from config.paths import LOG_DIR  # noqa: E402
from config.phase2_config import (  # noqa: E402
    ANNOTATED_H5AD,
    ANNOTATION_FIGURES_DIR,
    ANNOTATIONS_DIR,
    EMBEDDED_H5AD,
    PHASE2_SEED,
    PHASE2_WORK_DIR,
    ensure_phase2_directories,
)
from Stage02DogAnnotation.phase2_utils import (  # noqa: E402
    EPITHELIAL_LEVEL2_MARKERS,
    EPITHELIAL_METABOLIC_MARKERS,
    EPITHELIAL_STRESS_MARKERS,
    LEVEL1_MARKERS,
    MYELOID_LEVEL2_MARKERS,
    T_LEVEL2_MARKERS,
    atomic_write_h5ad,
    classify_groups,
    configure_logger,
    group_signature_scores,
    marker_var_names,
    peak_ram_mb,
    validate_phase1_clean,
)


SCRIPT_VERSION = "1.1.0"
GLOBAL_CLUSTER_KEY = "leiden_0_4"

# Frozen after inspection of independent, unintegrated global markers. No Human labels used.
GLOBAL_LEVEL1 = {
    "0": "B",
    "1": "Myeloid",
    "2": "T",
    "3": "T",
    "4": "T",
    "5": "T",
    "6": "T",
    "7": "T",
    "8": "Plasma",
    "9": "Myeloid",
    "10": "Mast",
    "11": "Epithelial",
    "12": "Epithelial",
    "13": "Epithelial",
    "14": "Epithelial",
    "15": "Stromal",
    "16": "T",
    "17": "Epithelial",
    "18": "Uncertain",
    "19": "Epithelial",
    "20": "Myeloid",
    "21": "Myeloid",
    "22": "Epithelial",
    "23": "Epithelial",
}

GLOBAL_CONFIDENCE = {
    **{cluster: "High" for cluster in GLOBAL_LEVEL1},
    "1": "Medium",
    "3": "Medium",
    "6": "Medium",
    "18": "Uncertain",
    "21": "Medium",
    "22": "Medium",
}

# Frozen after the preview marker audit. Values are (label, confidence, rationale).
# Mixed clusters remain in their frozen Level-1 lineage but are not forced into a
# paper subtype. This distinction is retained in the annotation audit table.
T_LEVEL2_OVERRIDES: dict[str, tuple[str, str, str]] = {
    "2": ("GZMA-high T", "High", "CCL5/CCL4/GZMA/GZMB/CD3D/CD8A support cytotoxic T; no leading IFN genes"),
    "3": ("GZMA-high T", "Medium", "CRTAM/CD160/CTSW/CD7/CD96 support cytotoxic T rather than IL7R-high state"),
    "4": ("IL7R-high T", "Medium", "IL7R/LTB and MHC-II-high lymphoid profile; insufficient ILC2-specific support"),
    "11": ("Uncertain epithelial-mixed", "Low", "FABP1/RBP2/APOA1/APOA4/SI dominate despite frozen T Level-1 origin"),
}
EPITHELIAL_LEVEL2_OVERRIDES: dict[str, tuple[str, str, str]] = {
    "4": ("Enterocyte 1 (IEG-high)", "Medium", "FOS/FOSB/DUSP1/EGR1 with FAM13A supports stressed immediate-early enterocyte"),
    "9": ("Uncertain lymphoid-mixed", "Low", "CD79A/CD3E/CCL5 dominate and epithelial subtype evidence is weak"),
    "11": ("Uncertain lymphoid-mixed", "Low", "CD79A-associated profile with weak enterocyte subtype evidence"),
    "17": ("LTF/KRT14 epithelial", "Medium", "KRT18/KRT14/LTF/PIGR/DMBT1 epithelial state without coherent IFN signature"),
}
MYELOID_LEVEL2_OVERRIDES: dict[str, tuple[str, str, str]] = {
    "3": ("Uncertain epithelial-mixed", "Low", "FABP1/RBP2/APOA1/APOA4/SI dominate"),
    "6": ("Cycling macrophage", "Medium", "STMN1/TOP2A/TYMS cycling program with C1QA/C1QB/C1QC/CTSS"),
    "8": ("Macrophage", "High", "C1QA/C1QB/C1QC/CTSS/FCGRT support macrophage identity"),
    "11": ("Uncertain T-mixed", "Low", "CD3D/CD3E/CD7/CD2/CCL5/GZMB dominate"),
    "12": ("Uncertain lymphoid-mixed", "Low", "GATA3/IL7R/IL17RB/CD3D lymphoid profile lacks coherent myeloid support"),
}

LINEAGE_PARAMETERS = {
    "Epithelial": {"prefix": "EPI", "n_hvg": 1_800, "n_pcs": 30, "n_neighbors": 25, "resolution": 0.8},
    "T": {"prefix": "T", "n_hvg": 1_800, "n_pcs": 30, "n_neighbors": 30, "resolution": 0.6},
    "Myeloid": {"prefix": "MY", "n_hvg": 1_200, "n_pcs": 25, "n_neighbors": 25, "resolution": 0.5},
}

LINEAGE_MARKERS = {
    "Epithelial": EPITHELIAL_LEVEL2_MARKERS,
    "T": T_LEVEL2_MARKERS,
    "Myeloid": MYELOID_LEVEL2_MARKERS,
}

LINEAGE_OVERRIDES = {
    "Epithelial": EPITHELIAL_LEVEL2_OVERRIDES,
    "T": T_LEVEL2_OVERRIDES,
    "Myeloid": MYELOID_LEVEL2_OVERRIDES,
}


def top_marker_differences(adata: ad.AnnData, groupby: str, n_top: int = 20) -> pd.DataFrame:
    groups = adata.obs[groupby].astype(str).to_numpy()
    X = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
    symbols = adata.var["gene_symbol"].astype(str).to_numpy()
    keep = np.array(
        [not (gene.startswith("MT-") or gene.startswith("RPL") or gene.startswith("RPS")) for gene in symbols]
    )
    total = np.asarray(X.sum(axis=0)).ravel()
    records: list[dict[str, object]] = []
    for group in sorted(set(groups), key=lambda item: int(item)):
        mask = groups == group
        inside_sum = np.asarray(X[mask].sum(axis=0)).ravel()
        inside = inside_sum / mask.sum()
        outside = (total - inside_sum) / (adata.n_obs - mask.sum())
        difference = inside - outside
        candidates = np.flatnonzero(keep)
        indices = candidates[np.argsort(difference[candidates])[-n_top:][::-1]]
        records.append(
            {
                "cluster": group,
                "n_cells": int(mask.sum()),
                "top_marker_symbols": ";".join(symbols[indices]),
                "top_mean_difference": ";".join(f"{difference[index]:.5f}" for index in indices),
            }
        )
    return pd.DataFrame(records)


def run_lineage(
    adata: ad.AnnData,
    lineage: str,
    logger: logging.Logger,
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    params = LINEAGE_PARAMETERS[lineage]
    subset = adata[adata.obs["level1"].astype(str) == lineage].copy()
    logger.info(
        "LINEAGE_START lineage=%s cells=%d donors=%d n_hvg=%d n_pcs=%d neighbors=%d resolution=%.2f",
        lineage,
        subset.n_obs,
        subset.obs["donor_id"].nunique(),
        params["n_hvg"],
        params["n_pcs"],
        params["n_neighbors"],
        params["resolution"],
    )
    # scikit-misc LOESS can be nearly singular for a lineage with many genes
    # having almost identical low means.  Keep a fixed, logged span escalation;
    # this changes only the smoothing bandwidth, not the donor-aware ranking.
    hvg_method = None
    last_hvg_error: Exception | None = None
    for span in (0.3, 0.5, 0.7, 1.0):
        try:
            sc.pp.highly_variable_genes(
                subset,
                flavor="seurat_v3",
                n_top_genes=min(params["n_hvg"] + 300, subset.n_vars),
                layer="counts",
                batch_key="donor_id",
                subset=False,
                check_values=True,
                span=span,
            )
            hvg_method = f"seurat_v3_donor_span_{span:.1f}"
            break
        except ValueError as exc:
            last_hvg_error = exc
            logger.warning("HVG_LOESS_RETRY lineage=%s span=%.1f error=%s", lineage, span, exc)
    if hvg_method is None:
        raise RuntimeError(f"All fixed donor-aware HVG LOESS spans failed for {lineage}") from last_hvg_error
    mt = subset.var["gene_symbol"].astype(str).str.startswith("MT-") | subset.var["gene_symbol"].astype(str).eq("ND4L")
    ranked = subset.var.loc[~mt & subset.var["highly_variable_rank"].notna()].sort_values(
        ["highly_variable_rank", "highly_variable_nbatches"], ascending=[True, False]
    )
    chosen = ranked.index[: min(params["n_hvg"], len(ranked))]
    subset.var[f"highly_variable_{params['prefix'].lower()}"] = subset.var_names.isin(chosen)
    work = subset[:, chosen].copy()
    work.X = work.layers["lognorm"].copy()
    for layer in list(work.layers.keys()):
        del work.layers[layer]
    sc.pp.scale(work, zero_center=True, max_value=10)
    sc.tl.pca(work, n_comps=params["n_pcs"], svd_solver="arpack", random_state=PHASE2_SEED)
    subset.obsm[f"X_pca_{params['prefix'].lower()}"] = np.asarray(work.obsm["X_pca"], dtype=np.float32)
    variance_ratio = np.asarray(work.uns["pca"]["variance_ratio"], dtype=float)
    del work
    gc.collect()
    neighbor_key = f"{params['prefix'].lower()}_neighbors"
    sc.pp.neighbors(
        subset,
        n_neighbors=params["n_neighbors"],
        n_pcs=params["n_pcs"],
        use_rep=f"X_pca_{params['prefix'].lower()}",
        metric="cosine",
        random_state=PHASE2_SEED,
        key_added=neighbor_key,
    )
    sc.tl.umap(
        subset,
        neighbors_key=neighbor_key,
        min_dist=0.30,
        random_state=PHASE2_SEED,
    )
    cluster_key = f"{params['prefix'].lower()}_cluster"
    sc.tl.leiden(
        subset,
        resolution=params["resolution"],
        random_state=PHASE2_SEED,
        key_added=cluster_key,
        neighbors_key=neighbor_key,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    score, raw, detection = group_signature_scores(subset, cluster_key, LINEAGE_MARKERS[lineage])
    classification = classify_groups(score, raw, detection)
    overrides = LINEAGE_OVERRIDES[lineage]
    classification["auto_label"] = classification["predicted_label"]
    classification["override_applied"] = False
    classification["override_rationale"] = ""
    for cluster, (label, confidence, rationale) in overrides.items():
        if cluster not in classification.index:
            raise RuntimeError(f"Override cluster {lineage}:{cluster} not present")
        classification.loc[cluster, "predicted_label"] = label
        classification.loc[cluster, "annotation_confidence"] = confidence
        classification.loc[cluster, "override_applied"] = True
        classification.loc[cluster, "override_rationale"] = rationale
    classification["lineage"] = lineage
    classification["lineage_cluster"] = [f"{params['prefix']}_{cluster}" for cluster in classification.index]
    classification["n_cells"] = subset.obs[cluster_key].astype(str).value_counts().reindex(classification.index).to_numpy()
    marker_preview = top_marker_differences(subset, cluster_key)
    marker_preview.insert(0, "lineage", lineage)
    marker_preview["lineage_cluster"] = params["prefix"] + "_" + marker_preview["cluster"].astype(str)
    audit = pd.crosstab(subset.obs[cluster_key].astype(str), subset.obs["donor_id"].astype(str))
    dominance = audit.max(axis=1) / audit.sum(axis=1)
    metrics = {
        "lineage": lineage,
        "n_cells": subset.n_obs,
        "n_donors": subset.obs["donor_id"].nunique(),
        "n_hvg": len(chosen),
        "hvg_method": hvg_method,
        "n_pcs": params["n_pcs"],
        "n_neighbors": params["n_neighbors"],
        "resolution": params["resolution"],
        "n_clusters": subset.obs[cluster_key].nunique(),
        "cumulative_pca_variance": float(variance_ratio.sum()),
        "donor_dominant_clusters_gt_0_8": int((dominance > 0.8).sum()),
        "max_dominant_donor_fraction": float(dominance.max()),
    }
    logger.info(
        "LINEAGE_END lineage=%s clusters=%d donor_dominant_gt80=%d auto_labels=%s",
        lineage,
        metrics["n_clusters"],
        metrics["donor_dominant_clusters_gt_0_8"],
        classification["predicted_label"].value_counts().to_dict(),
    )
    return subset, classification.reset_index(), marker_preview, metrics


def gene_expression(adata: ad.AnnData, symbol: str) -> np.ndarray:
    indices = np.flatnonzero(adata.var["gene_symbol"].astype(str).to_numpy() == symbol)
    if not len(indices):
        return np.zeros(adata.n_obs, dtype=np.float32)
    values = adata.X[:, indices]
    if sparse.issparse(values):
        return np.asarray(values.sum(axis=1)).ravel().astype(np.float32)
    return np.asarray(values).sum(axis=1).astype(np.float32)


def module_mean(adata: ad.AnnData, genes: list[str]) -> np.ndarray:
    variables = marker_var_names(adata, genes)
    if not variables:
        return np.zeros(adata.n_obs, dtype=np.float32)
    values = adata[:, variables].X
    result = np.asarray(values.mean(axis=1)).ravel()
    return result.astype(np.float32)


def extract_rank_markers(adata: ad.AnnData) -> pd.DataFrame:
    key = "rank_genes_annotation_cluster"
    sc.tl.rank_genes_groups(
        adata,
        groupby="annotation_cluster",
        method="wilcoxon",
        use_raw=False,
        n_genes=min(500, adata.n_vars),
        pts=True,
        tie_correct=True,
        key_added=key,
    )
    result = sc.get.rank_genes_groups_df(adata, group=None, key=key)
    rename = {
        "group": "cluster",
        "names": "unique_var_name",
        "logfoldchanges": "logFC",
        "scores": "score",
        "pvals_adj": "adjusted_p",
        "pct_nz_group": "pct_in",
        "pct_nz_reference": "pct_out",
    }
    result = result.rename(columns=rename)
    if "pct_in" not in result or "pct_out" not in result:
        pts = adata.uns[key].get("pts")
        pts_rest = adata.uns[key].get("pts_rest")
        if pts is None or pts_rest is None:
            raise RuntimeError("Scanpy did not expose marker detection proportions")
        result["pct_in"] = [pts.loc[gene, group] for gene, group in zip(result["unique_var_name"], result["cluster"])]
        result["pct_out"] = [pts_rest.loc[gene, group] for gene, group in zip(result["unique_var_name"], result["cluster"])]
    symbol_lookup = adata.var["gene_symbol"].astype(str).to_dict()
    result["gene"] = result["unique_var_name"].map(symbol_lookup)
    result["marker_rank"] = result.groupby("cluster", observed=True).cumcount() + 1
    columns = [
        "cluster", "gene", "unique_var_name", "logFC", "score", "pct_in", "pct_out",
        "adjusted_p", "marker_rank",
    ]
    return result[columns]


def plot_categorical(
    coordinates: np.ndarray,
    values: pd.Series,
    destination: Path,
    title: str,
    point_size: float = 0.9,
) -> None:
    categories = sorted(values.astype(str).unique())
    cmap = plt.get_cmap("tab20", max(1, len(categories)))
    fig, ax = plt.subplots(figsize=(8.5, 7))
    for index, category in enumerate(categories):
        mask = values.astype(str).to_numpy() == category
        ax.scatter(coordinates[mask, 0], coordinates[mask, 1], s=point_size, alpha=0.65, linewidths=0, color=cmap(index), label=category, rasterized=True)
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(markerscale=5, fontsize=7, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(destination, dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_positive_controls(adata: ad.AnnData) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    def add_row(name: str, lineage: str, mask: np.ndarray, supporting: str, status: str, notes: str) -> None:
        selected = adata.obs.loc[mask]
        clusters = selected["annotation_cluster"].astype(str).value_counts().head(5).index.tolist()
        healthy_donors = selected.loc[selected["condition"].astype(str) == "Healthy", "donor_id"].astype(str).nunique()
        cie_donors = selected.loc[selected["condition"].astype(str) == "CIE", "donor_id"].astype(str).nunique()
        rows.append(
            {
                "positive_control": name,
                "expected_lineage": lineage,
                "detected": bool(mask.sum() > 0),
                "cluster_or_state": ";".join(clusters),
                "supporting_markers": supporting,
                "n_cells": int(mask.sum()),
                "n_donors": int(selected["donor_id"].astype(str).nunique()),
                "healthy_donors_detected": int(healthy_donors),
                "cie_donors_detected": int(cie_donors),
                "reproduction_status": status,
                "notes": notes,
            }
        )

    level1 = adata.obs["level1"].astype(str).to_numpy()
    level2 = adata.obs["level2"].astype(str).to_numpy()
    add_row(
        "IL7R-high T", "T", (level1 == "T") & (level2 == "IL7R-high T"),
        "IL7R;TCF7;CXCR4;SELL;CCR7;LTB", "Reproduced",
        "Identity/state presence only; no disease abundance test.",
    )
    add_row(
        "GZMA-high T", "T", (level1 == "T") & (level2 == "GZMA-high T"),
        "GZMA;GZMB;ITGB7;ITGAE;CCL5;CD7;CD2", "Reproduced",
        "Identity/state presence only; no disease abundance test.",
    )
    sod2 = gene_expression(adata, "SOD2") > 0
    il1a = gene_expression(adata, "IL1A") > 0
    add_row(
        "SOD2-associated myeloid state", "Myeloid", (level1 == "Myeloid") & sod2,
        "SOD2;S100A12;S100A8;NFKBIA;SELL", "Reproduced",
        "Marker-positive myeloid cells described without CIE-vs-Healthy testing.",
    )
    add_row(
        "IL1A-associated myeloid state", "Myeloid", (level1 == "Myeloid") & il1a,
        "IL1A;SOD2;TNFAIP3;S100A12;S100A8", "Reproduced" if ((level1 == "Myeloid") & il1a).sum() >= 10 else "Partially reproduced",
        "Marker-positive myeloid cells described without CIE-vs-Healthy testing.",
    )
    epithelial = level1 == "Epithelial"
    stress = adata.obs["epithelial_stress_score"].to_numpy()
    metabolic = adata.obs["epithelial_metabolic_score"].to_numpy()
    stress_threshold = float(np.quantile(stress[epithelial], 0.90))
    metabolic_threshold = float(np.quantile(metabolic[epithelial], 0.90))
    add_row(
        "epithelial stress-related state", "Epithelial", epithelial & (stress >= stress_threshold),
        ";".join(EPITHELIAL_STRESS_MARKERS), "Partially reproduced",
        "Top-decile paper-compatible stress program; no pathway enrichment or disease DEG performed.",
    )
    add_row(
        "epithelial metabolic state", "Epithelial", epithelial & (metabolic >= metabolic_threshold),
        ";".join(EPITHELIAL_METABOLIC_MARKERS), "Partially reproduced",
        "Top-decile enterocyte metabolic program; no pathway enrichment or disease DEG performed.",
    )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preview", action="store_true")
    args = parser.parse_args()
    ensure_phase2_directories()
    logger = configure_logger(LOG_DIR / "phase2_annotation.log", "phase2_annotation")
    process = psutil.Process(os.getpid())
    started = datetime.now(timezone.utc)
    np.random.seed(PHASE2_SEED)
    if not args.preview:
        formal = [
            ANNOTATED_H5AD,
            ANNOTATIONS_DIR / "dog_cell_annotation.csv",
            ANNOTATIONS_DIR / "cluster_markers.csv",
        ]
        existing = [path for path in formal if path.exists()]
        if existing:
            raise FileExistsError(f"Refusing to overwrite formal annotation outputs: {existing}")
    logger.info(
        "START script_version=%s python=%s scanpy=%s seed=%d preview=%s input=%s",
        SCRIPT_VERSION,
        platform.python_version(),
        sc.__version__,
        PHASE2_SEED,
        args.preview,
        EMBEDDED_H5AD,
    )
    try:
        adata = ad.read_h5ad(EMBEDDED_H5AD)
        validate_phase1_clean(adata)
        if GLOBAL_CLUSTER_KEY not in adata.obs:
            raise RuntimeError(f"Missing selected global cluster key: {GLOBAL_CLUSTER_KEY}")
        observed_clusters = set(adata.obs[GLOBAL_CLUSTER_KEY].astype(str))
        if observed_clusters != set(GLOBAL_LEVEL1):
            raise RuntimeError("Frozen global cluster map does not match embedded object")
        adata.obs["global_cluster"] = adata.obs[GLOBAL_CLUSTER_KEY].astype(str)
        adata.obs["level1"] = adata.obs["global_cluster"].map(GLOBAL_LEVEL1).astype(str)
        adata.obs["annotation_confidence"] = adata.obs["global_cluster"].map(GLOBAL_CONFIDENCE).astype(str)
        logger.info("LEVEL1 counts=%s", adata.obs["level1"].value_counts().to_dict())

        lineage_subsets: dict[str, ad.AnnData] = {}
        classifications: list[pd.DataFrame] = []
        previews: list[pd.DataFrame] = []
        metrics: list[dict[str, object]] = []
        for lineage in ("Epithelial", "T", "Myeloid"):
            subset, classification, marker_preview, metric = run_lineage(adata, lineage, logger)
            lineage_subsets[lineage] = subset
            classifications.append(classification)
            previews.append(marker_preview)
            metrics.append(metric)
        classification_table = pd.concat(classifications, ignore_index=True)
        marker_preview_table = pd.concat(previews, ignore_index=True)
        metrics_table = pd.DataFrame(metrics)
        if args.preview:
            classification_table.to_csv(PHASE2_WORK_DIR / "lineage_preview_assignments.csv", index=False)
            marker_preview_table.to_csv(PHASE2_WORK_DIR / "lineage_preview_cluster_markers.csv", index=False)
            metrics_table.to_csv(PHASE2_WORK_DIR / "lineage_preview_metrics.csv", index=False)
            logger.info(
                "OUTPUT preview_assignments=%s preview_markers=%s",
                PHASE2_WORK_DIR / "lineage_preview_assignments.csv",
                PHASE2_WORK_DIR / "lineage_preview_cluster_markers.csv",
            )
            logger.info(
                "END status=success_preview elapsed_seconds=%.3f peak_ram_mb=%.1f",
                (datetime.now(timezone.utc) - started).total_seconds(),
                peak_ram_mb(process),
            )
            return

        adata.obs["lineage_cluster"] = "NA"
        adata.obs["level2"] = adata.obs["level1"].astype(str)
        for lineage, subset in lineage_subsets.items():
            params = LINEAGE_PARAMETERS[lineage]
            prefix = params["prefix"]
            cluster_key = f"{prefix.lower()}_cluster"
            table = classification_table.loc[classification_table["lineage"] == lineage].set_index("cluster")
            clusters = subset.obs[cluster_key].astype(str)
            lineage_cluster = prefix + "_" + clusters
            level2 = clusters.map(table["predicted_label"]).astype(str)
            confidence = clusters.map(table["annotation_confidence"]).astype(str)
            adata.obs.loc[subset.obs_names, "lineage_cluster"] = lineage_cluster.to_numpy()
            adata.obs.loc[subset.obs_names, "level2"] = level2.to_numpy()
            adata.obs.loc[subset.obs_names, "annotation_confidence"] = confidence.to_numpy()
            full_coordinates = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)
            positions = adata.obs_names.get_indexer(subset.obs_names)
            full_coordinates[positions] = np.asarray(subset.obsm["X_umap"], dtype=np.float32)
            adata.obsm[f"X_umap_{prefix.lower()}"] = full_coordinates
            adata.obs[f"{prefix.lower()}_cluster"] = "NA"
            adata.obs.loc[subset.obs_names, f"{prefix.lower()}_cluster"] = clusters.to_numpy()

        non_primary = ~adata.obs["level1"].isin(["Epithelial", "T", "Myeloid"])
        adata.obs["annotation_cluster"] = adata.obs["lineage_cluster"].astype(str)
        adata.obs.loc[non_primary, "annotation_cluster"] = (
            "G" + adata.obs.loc[non_primary, "global_cluster"].astype(str) + "_" + adata.obs.loc[non_primary, "level1"].astype(str)
        )
        adata.obs["epithelial_stress_score"] = module_mean(adata, EPITHELIAL_STRESS_MARKERS)
        adata.obs["epithelial_metabolic_score"] = module_mean(adata, EPITHELIAL_METABOLIC_MARKERS)
        adata.obs["SOD2_detected"] = gene_expression(adata, "SOD2") > 0
        adata.obs["IL1A_detected"] = gene_expression(adata, "IL1A") > 0

        logger.info("MARKERS_START method=wilcoxon groupby=annotation_cluster n_genes=500 disease_contrast=False")
        markers = extract_rank_markers(adata)
        markers.to_csv(ANNOTATIONS_DIR / "cluster_markers.csv", index=False)
        logger.info("MARKERS_END rows=%d clusters=%d", len(markers), markers["cluster"].nunique())

        donor_composition = (
            adata.obs.groupby(["annotation_cluster", "donor_id", "condition"], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        donor_composition["fraction_of_cluster"] = donor_composition["n_cells"] / donor_composition.groupby(
            "annotation_cluster", observed=True
        )["n_cells"].transform("sum")
        donor_composition["fraction_of_donor"] = donor_composition["n_cells"] / donor_composition.groupby(
            "donor_id", observed=True
        )["n_cells"].transform("sum")
        donor_composition.rename(columns={"annotation_cluster": "cluster", "donor_id": "donor"}).to_csv(
            ANNOTATIONS_DIR / "cluster_donor_composition.csv", index=False
        )

        annotation_rows: list[dict[str, object]] = []
        top_markers = markers.loc[~markers["gene"].astype(str).str.match(r"^(MT-|RPL|RPS)")]
        for cluster, group in adata.obs.groupby("annotation_cluster", observed=True):
            donors = group["donor_id"].astype(str).value_counts()
            cluster_markers = top_markers.loc[top_markers["cluster"].astype(str) == str(cluster)].head(12)
            level1 = group["level1"].astype(str).mode().iloc[0]
            level2 = group["level2"].astype(str).mode().iloc[0]
            confidence = group["annotation_confidence"].astype(str).mode().iloc[0]
            paper_match = "Reproduced" if level1 not in {"Uncertain", "Endothelial"} else "Not reproduced"
            paper_state = level2 if paper_match == "Reproduced" else "NA"
            notes: list[str] = []
            if donors.iloc[0] / donors.sum() > 0.8:
                notes.append("single-donor dominated; descriptive only")
            if level1 == "Uncertain":
                notes.append("mixed plasma/epithelial markers; retained as uncertain")
            annotation_rows.append(
                {
                    "cluster_id": cluster,
                    "level1": level1,
                    "level2": level2,
                    "annotation_confidence": confidence,
                    "positive_markers": ";".join(cluster_markers["gene"].astype(str)),
                    "negative_or_conflicting_markers": "See marker score audit and cluster_markers.csv",
                    "original_paper_match": paper_match,
                    "paper_cell_state": paper_state,
                    "n_cells": len(group),
                    "n_donors": int(group["donor_id"].astype(str).nunique()),
                    "dominant_donor": donors.index[0],
                    "dominant_donor_fraction": float(donors.iloc[0] / donors.sum()),
                    "notes": "; ".join(notes) or "Marker-coherent independent annotation",
                }
            )
        annotation_table = pd.DataFrame(annotation_rows).sort_values(["level1", "level2", "cluster_id"])
        annotation_table.to_csv(ANNOTATIONS_DIR / "dog_cell_annotation.csv", index=False)

        proportions = (
            adata.obs.groupby(["donor_id", "condition", "level1", "level2"], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        proportions["fraction_of_donor"] = proportions["n_cells"] / proportions.groupby(
            "donor_id", observed=True
        )["n_cells"].transform("sum")
        proportions.rename(columns={"donor_id": "donor"}).to_csv(
            ANNOTATIONS_DIR / "celltype_proportion_per_donor.csv", index=False
        )
        positive_controls = build_positive_controls(adata)
        positive_controls.to_csv(ANNOTATIONS_DIR / "positive_control_summary.csv", index=False)
        classification_table.to_csv(ANNOTATIONS_DIR / "lineage_marker_assignment_audit.csv", index=False)
        metrics_table.to_csv(ANNOTATIONS_DIR / "lineage_embedding_parameters.csv", index=False)

        adata.uns["phase2_annotation"] = {
            "script_version": SCRIPT_VERSION,
            "random_seed": PHASE2_SEED,
            "primary_global_cluster_key": GLOBAL_CLUSTER_KEY,
            "primary_annotation_source": "GSE254005 paper and author v0.1.0 canine markers plus current cluster markers",
            "human_reference_used_as_primary": False,
            "disease_deg_performed": False,
            "formal_level1": sorted(adata.obs["level1"].astype(str).unique()),
            "lineage_specific_reclustering": ["Epithelial", "T", "Myeloid"],
            "annotation_rule": "multiple coherent positive markers; uncertain/mixed cells retained",
        }

        plot_categorical(adata.obsm["X_umap"], adata.obs["level1"], ANNOTATION_FIGURES_DIR / "07_umap_by_level1.png", "Primary UMAP by formal Level1")
        plot_categorical(adata.obsm["X_umap"], adata.obs["annotation_cluster"], ANNOTATION_FIGURES_DIR / "08_umap_by_annotation_cluster.png", "Primary UMAP by annotation cluster", 0.7)
        sc.pl.umap(adata, color=["IL7R", "GZMA", "SOD2", "IL1A"], ncols=2, show=False, use_raw=False, cmap="magma")
        plt.savefig(ANNOTATION_FIGURES_DIR / "09_positive_control_feature_umaps.png", dpi=220, bbox_inches="tight")
        plt.close("all")
        dot_markers = {
            label: [gene for gene in genes if gene in set(adata.var_names)]
            for label, genes in LEVEL1_MARKERS.items()
        }
        dot_markers = {label: genes[:5] for label, genes in dot_markers.items() if genes}
        dot = sc.pl.dotplot(adata, var_names=dot_markers, groupby="level1", use_raw=False, show=False, return_fig=True)
        dot.savefig(ANNOTATION_FIGURES_DIR / "10_level1_marker_dotplot.png", dpi=220, bbox_inches="tight")
        plt.close("all")
        for index, lineage in enumerate(("Epithelial", "T", "Myeloid"), start=11):
            subset = lineage_subsets[lineage]
            plot_categorical(
                subset.obsm["X_umap"],
                adata.obs.loc[subset.obs_names, "level2"],
                ANNOTATION_FIGURES_DIR / f"{index:02d}_{lineage.lower()}_lineage_umap.png",
                f"{lineage} lineage-specific UMAP by Level2",
                1.2,
            )

        atomic_write_h5ad(adata, ANNOTATED_H5AD)
        check = ad.read_h5ad(ANNOTATED_H5AD, backed="r")
        try:
            required = {"level1", "level2", "annotation_cluster", "annotation_confidence"}
            if check.shape != adata.shape or not required.issubset(check.obs.columns) or "counts" not in check.layers or "lognorm" not in check.layers:
                raise RuntimeError("Annotated h5ad read-back validation failed")
        finally:
            check.file.close()
        logger.info(
            "OUTPUT annotated_h5ad=%s cells=%d level1=%s level2_n=%d annotation_clusters=%d",
            ANNOTATED_H5AD,
            adata.n_obs,
            adata.obs["level1"].value_counts().to_dict(),
            adata.obs["level2"].nunique(),
            adata.obs["annotation_cluster"].nunique(),
        )
        logger.info(
            "END status=success elapsed_seconds=%.3f peak_ram_mb=%.1f",
            (datetime.now(timezone.utc) - started).total_seconds(),
            peak_ram_mb(process),
        )
    except Exception:
        logger.exception(
            "END status=error elapsed_seconds=%.3f peak_ram_mb=%.1f",
            (datetime.now(timezone.utc) - started).total_seconds(),
            peak_ram_mb(process),
        )
        raise


if __name__ == "__main__":
    main()
