"""Build the primary unintegrated canine atlas embedding and multi-resolution clusters."""

from __future__ import annotations

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

CODE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_ROOT))
PROJECT_ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"] )

from config.paths import LOG_DIR  # noqa: E402
from config.phase2_config import (  # noqa: E402
    ANNOTATION_FIGURES_DIR,
    ANNOTATIONS_DIR,
    CLEAN_H5AD,
    EMBEDDED_H5AD,
    LEIDEN_RESOLUTIONS,
    N_HVG,
    N_NEIGHBORS,
    N_PCA_COMPUTE,
    N_PCS_NEIGHBORS,
    PHASE2_SEED,
    PRIMARY_LEIDEN_KEY,
    UMAP_MIN_DIST,
    ensure_phase2_directories,
)
from Stage02DogAnnotation.phase2_utils import (  # noqa: E402
    LEVEL1_MARKERS,
    atomic_write_h5ad,
    classify_groups,
    configure_logger,
    group_signature_scores,
    peak_ram_mb,
    validate_phase1_clean,
)


SCRIPT_VERSION = "1.0.0"


def categorical_umap(
    coordinates: np.ndarray,
    values: pd.Series,
    destination: Path,
    title: str,
    point_size: float = 1.0,
) -> None:
    categories = sorted(values.astype(str).unique())
    cmap = plt.get_cmap("tab20", max(1, len(categories)))
    fig, ax = plt.subplots(figsize=(8, 7))
    for index, category in enumerate(categories):
        mask = values.astype(str).to_numpy() == category
        ax.scatter(
            coordinates[mask, 0],
            coordinates[mask, 1],
            s=point_size,
            alpha=0.65,
            linewidths=0,
            color=cmap(index),
            label=category,
            rasterized=True,
        )
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(markerscale=5, fontsize=7, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(destination, dpi=220, bbox_inches="tight")
    plt.close(fig)


def continuous_umap(
    coordinates: np.ndarray,
    values: pd.Series,
    destination: Path,
    title: str,
) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    scatter = ax.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        c=values.to_numpy(dtype=float),
        s=1.0,
        alpha=0.65,
        linewidths=0,
        cmap="viridis",
        rasterized=True,
    )
    fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    fig.tight_layout()
    fig.savefig(destination, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ensure_phase2_directories()
    if EMBEDDED_H5AD.exists():
        raise FileExistsError(f"Refusing to overwrite formal output: {EMBEDDED_H5AD}")
    logger = configure_logger(LOG_DIR / "phase2_embedding.log", "phase2_embedding")
    process = psutil.Process(os.getpid())
    started = datetime.now(timezone.utc)
    np.random.seed(PHASE2_SEED)
    logger.info(
        "START script_version=%s python=%s scanpy=%s seed=%d input=%s",
        SCRIPT_VERSION,
        platform.python_version(),
        sc.__version__,
        PHASE2_SEED,
        CLEAN_H5AD,
    )
    try:
        adata = ad.read_h5ad(CLEAN_H5AD)
        validate_phase1_clean(adata)
        if "lognorm" in adata.layers:
            raise RuntimeError("Unexpected pre-existing lognorm layer in Phase1 object")
        adata.layers["lognorm"] = adata.X.copy()

        logger.info(
            "HVG_START flavor=seurat_v3 batch_key=donor_id layer=counts requested_non_mt=%d",
            N_HVG,
        )
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat_v3",
            n_top_genes=3_000,
            layer="counts",
            batch_key="donor_id",
            subset=False,
            check_values=True,
            span=0.3,
        )
        mt = adata.var["gene_symbol"].astype(str).str.startswith("MT-") | adata.var["gene_symbol"].astype(str).eq("ND4L")
        ranked = adata.var.loc[~mt & adata.var["highly_variable_rank"].notna()].sort_values(
            ["highly_variable_rank", "highly_variable_nbatches"], ascending=[True, False]
        )
        if len(ranked) < N_HVG:
            raise RuntimeError(f"Only {len(ranked)} non-mito ranked HVGs available")
        selected = ranked.index[:N_HVG]
        adata.var["highly_variable_phase2"] = adata.var_names.isin(selected)
        adata.var["mitochondrial_phase2"] = mt.to_numpy()
        logger.info(
            "HVG_END selected=%d mt_excluded=True donor_batches=%d median_nbatches=%.1f",
            int(adata.var["highly_variable_phase2"].sum()),
            adata.obs["donor_id"].nunique(),
            float(adata.var.loc[selected, "highly_variable_nbatches"].median()),
        )

        work = adata[:, adata.var["highly_variable_phase2"]].copy()
        work.X = work.layers["lognorm"].copy()
        for layer in list(work.layers.keys()):
            del work.layers[layer]
        logger.info("PCA_START scale=HVG_only max_value=10 regressors=none")
        sc.pp.scale(work, zero_center=True, max_value=10)
        sc.tl.pca(
            work,
            n_comps=N_PCA_COMPUTE,
            svd_solver="arpack",
            random_state=PHASE2_SEED,
        )
        variance = np.asarray(work.uns["pca"]["variance"])
        variance_ratio = np.asarray(work.uns["pca"]["variance_ratio"])
        pca_table = pd.DataFrame(
            {
                "pc": np.arange(1, len(variance) + 1),
                "variance": variance,
                "variance_ratio": variance_ratio,
                "cumulative_variance_ratio": np.cumsum(variance_ratio),
                "used_for_neighbors": np.arange(1, len(variance) + 1) <= N_PCS_NEIGHBORS,
            }
        )
        pca_table.to_csv(ANNOTATIONS_DIR / "pca_variance.csv", index=False)
        adata.obsm["X_pca"] = np.asarray(work.obsm["X_pca"], dtype=np.float32)
        loadings = np.zeros((adata.n_vars, N_PCA_COMPUTE), dtype=np.float32)
        loadings[adata.var["highly_variable_phase2"].to_numpy(), :] = np.asarray(
            work.varm["PCs"], dtype=np.float32
        )
        adata.varm["PCs"] = loadings
        adata.uns["pca"] = {
            "variance": variance.astype(np.float32),
            "variance_ratio": variance_ratio.astype(np.float32),
            "params": {
                "zero_center": True,
                "use_highly_variable": True,
                "mask_var": "highly_variable_phase2",
            },
        }
        del work, loadings
        gc.collect()
        logger.info(
            "PCA_END computed=%d used=%d cumulative_variance_used=%.6f",
            N_PCA_COMPUTE,
            N_PCS_NEIGHBORS,
            float(pca_table.loc[pca_table["used_for_neighbors"], "variance_ratio"].sum()),
        )

        sc.pp.neighbors(
            adata,
            n_neighbors=N_NEIGHBORS,
            n_pcs=N_PCS_NEIGHBORS,
            use_rep="X_pca",
            metric="cosine",
            random_state=PHASE2_SEED,
            key_added="phase2_neighbors",
        )
        sc.tl.umap(
            adata,
            neighbors_key="phase2_neighbors",
            min_dist=UMAP_MIN_DIST,
            random_state=PHASE2_SEED,
        )
        resolution_rows: list[dict[str, object]] = []
        marker_audits: list[pd.DataFrame] = []
        for resolution in LEIDEN_RESOLUTIONS:
            key = f"leiden_{str(resolution).replace('.', '_')}"
            sc.tl.leiden(
                adata,
                resolution=resolution,
                random_state=PHASE2_SEED,
                key_added=key,
                neighbors_key="phase2_neighbors",
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
            counts = adata.obs[key].astype(str).value_counts()
            donor_table = pd.crosstab(adata.obs[key].astype(str), adata.obs["donor_id"].astype(str))
            dominance = donor_table.max(axis=1) / donor_table.sum(axis=1)
            score, raw, detection = group_signature_scores(adata, key, LEVEL1_MARKERS)
            classification = classify_groups(score, raw, detection).reset_index()
            classification.insert(0, "leiden_key", key)
            classification.insert(1, "resolution", resolution)
            marker_audits.append(classification)
            resolution_rows.append(
                {
                    "leiden_key": key,
                    "resolution": resolution,
                    "n_clusters": len(counts),
                    "min_cluster_cells": int(counts.min()),
                    "median_cluster_cells": float(counts.median()),
                    "max_cluster_cells": int(counts.max()),
                    "clusters_lt_20_cells": int((counts < 20).sum()),
                    "clusters_dominant_donor_gt_0_8": int((dominance > 0.8).sum()),
                    "max_dominant_donor_fraction": float(dominance.max()),
                    "low_or_uncertain_marker_clusters": int(
                        classification["annotation_confidence"].isin(["Low", "Uncertain"]).sum()
                    ),
                    "selected_primary": key == PRIMARY_LEIDEN_KEY,
                }
            )
            logger.info(
                "LEIDEN key=%s resolution=%.2f clusters=%d min_cells=%d donor_dominant_gt80=%d marker_low_or_uncertain=%d",
                key,
                resolution,
                len(counts),
                int(counts.min()),
                int((dominance > 0.8).sum()),
                int(classification["annotation_confidence"].isin(["Low", "Uncertain"]).sum()),
            )
        pd.DataFrame(resolution_rows).to_csv(
            ANNOTATIONS_DIR / "leiden_resolution_summary.csv", index=False
        )
        pd.concat(marker_audits, ignore_index=True).to_csv(
            ANNOTATIONS_DIR / "leiden_resolution_marker_audit.csv", index=False
        )

        parameters = pd.DataFrame(
            [
                ("random_seed", PHASE2_SEED, "all stochastic steps"),
                ("hvg_flavor", "seurat_v3", "raw counts; batch_key=donor_id"),
                ("hvg_non_mito", N_HVG, "mitochondrial genes retained in object but excluded from HVG"),
                ("scale", "HVG only; max_value=10", "not saved over counts/lognorm"),
                ("pca_components_computed", N_PCA_COMPUTE, "ARPACK"),
                ("pca_components_neighbors", N_PCS_NEIGHBORS, "fixed before clustering"),
                ("neighbors", N_NEIGHBORS, "cosine on unintegrated PCA"),
                ("umap_min_dist", UMAP_MIN_DIST, "primary embedding"),
                ("leiden_resolutions", ";".join(map(str, LEIDEN_RESOLUTIONS)), "low/mid/high audit"),
                ("primary_leiden_key", PRIMARY_LEIDEN_KEY, "mid-resolution candidate; marker/donor audited"),
                ("donor_batch_correction", "none", "condition/donor not regressed or integrated"),
                ("mt_regression", "none", "cie_2 mt signal intentionally preserved for audit"),
            ],
            columns=["parameter", "value", "notes"],
        )
        parameters.to_csv(ANNOTATIONS_DIR / "embedding_parameters.csv", index=False)
        adata.uns["phase2_embedding"] = {
            "script_version": SCRIPT_VERSION,
            "random_seed": PHASE2_SEED,
            "hvg_flavor": "seurat_v3",
            "hvg_batch_key": "donor_id",
            "n_hvg_non_mito": N_HVG,
            "n_pca_computed": N_PCA_COMPUTE,
            "n_pcs_neighbors": N_PCS_NEIGHBORS,
            "n_neighbors": N_NEIGHBORS,
            "neighbor_metric": "cosine",
            "umap_min_dist": UMAP_MIN_DIST,
            "leiden_resolutions": list(LEIDEN_RESOLUTIONS),
            "primary_leiden_key": PRIMARY_LEIDEN_KEY,
            "donor_batch_correction": "none",
            "regressors": "none",
            "mt_genes_excluded_from_hvg_only": True,
        }

        coordinates = adata.obsm["X_umap"]
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(pca_table["pc"], pca_table["variance_ratio"], marker="o", markersize=2)
        ax.axvline(N_PCS_NEIGHBORS, color="firebrick", linestyle="--", label=f"Used PCs = {N_PCS_NEIGHBORS}")
        ax.set(xlabel="PC", ylabel="Explained variance ratio", title="Phase2 PCA variance")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(ANNOTATION_FIGURES_DIR / "01_pca_variance.png", dpi=220)
        plt.close(fig)

        fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
        for ax, resolution in zip(axes, LEIDEN_RESOLUTIONS):
            key = f"leiden_{str(resolution).replace('.', '_')}"
            values = adata.obs[key].astype(str)
            categories = sorted(values.unique(), key=lambda item: int(item))
            cmap = plt.get_cmap("tab20", max(1, len(categories)))
            for index, category in enumerate(categories):
                mask = values.to_numpy() == category
                ax.scatter(coordinates[mask, 0], coordinates[mask, 1], s=0.7, alpha=0.55, linewidths=0, color=cmap(index), rasterized=True)
            ax.set_title(f"Leiden resolution {resolution} (n={len(categories)})")
            ax.set_xticks([])
            ax.set_yticks([])
        fig.tight_layout()
        fig.savefig(ANNOTATION_FIGURES_DIR / "02_umap_multi_resolution.png", dpi=220)
        plt.close(fig)
        categorical_umap(coordinates, adata.obs["donor_id"], ANNOTATION_FIGURES_DIR / "03_umap_by_donor.png", "Primary UMAP by donor")
        categorical_umap(coordinates, adata.obs["condition"], ANNOTATION_FIGURES_DIR / "04_umap_by_condition.png", "Primary UMAP by condition (descriptive only)")
        continuous_umap(coordinates, adata.obs["pct_counts_mt"], ANNOTATION_FIGURES_DIR / "05_umap_by_mt_pct.png", "Primary UMAP by mitochondrial percentage")
        continuous_umap(coordinates, adata.obs["n_genes_by_counts"], ANNOTATION_FIGURES_DIR / "06_umap_by_nfeature.png", "Primary UMAP by detected genes")

        atomic_write_h5ad(adata, EMBEDDED_H5AD)
        check = ad.read_h5ad(EMBEDDED_H5AD, backed="r")
        try:
            if check.shape != adata.shape or "lognorm" not in check.layers or "X_umap" not in check.obsm:
                raise RuntimeError("Embedded h5ad read-back validation failed")
        finally:
            check.file.close()
        logger.info(
            "OUTPUT embedded_h5ad=%s shape=%s primary_key=%s primary_clusters=%d",
            EMBEDDED_H5AD,
            adata.shape,
            PRIMARY_LEIDEN_KEY,
            adata.obs[PRIMARY_LEIDEN_KEY].nunique(),
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
