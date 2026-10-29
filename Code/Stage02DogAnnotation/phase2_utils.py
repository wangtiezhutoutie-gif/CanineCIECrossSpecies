"""Shared, canine-only utilities for Phase 2 scripts."""

from __future__ import annotations

import hashlib
import logging
import os
import time
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
import psutil
from scipy import sparse


EXPECTED_SHAPE = (38_186, 20_257)
EXPECTED_CONDITIONS = {"Healthy": 15_589, "CIE": 22_597}
EXPECTED_DONORS = {
    "GSM8031698_h_1",
    "GSM8031699_h_2",
    "GSM8031700_h_3",
    "GSM8031701_cie_2",
    "GSM8031702_cie_3",
    "GSM8031703_cie_4",
    "GSM8031704_cie_6",
}


LEVEL1_MARKERS: dict[str, list[str]] = {
    "Epithelial": ["EPCAM", "KRT18", "KRT20", "VIL1", "ANPEP", "SI", "APOA4", "APOC3", "FABP1"],
    "T": ["CD3D", "CD3E", "LCK", "TRAT1", "CD2", "CD7", "IL7R", "GZMA"],
    "B": ["MS4A1", "CD79A", "CD37", "CD22", "PAX5", "TNFRSF13C", "FCRLA"],
    "Plasma": ["JCHAIN", "MZB1", "SDC1", "TXNDC5", "DERL3", "PRDX4", "IGHM"],
    "Myeloid": ["AIF1", "FCER1G", "CTSS", "CSF1R", "C1QA", "C1QC", "S100A8", "S100A12", "FLT3"],
    "Mast": ["KIT", "CPA3", "MS4A2", "IGF1", "F2RL3", "SVIL"],
    "Stromal": ["COL1A1", "COL1A2", "COL3A1", "DCN", "THY1", "ACTA2", "TPM2", "IGFBP7"],
    "Endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "ESAM", "ENG", "PLVAP"],
}

T_LEVEL2_MARKERS: dict[str, list[str]] = {
    "IL7R-high T": ["IL7R", "TCF7", "CXCR4", "SELL", "CCR7", "LTB"],
    "GZMA-high T": ["GZMA", "GZMB", "ITGB7", "ITGAE", "CCL5", "CD7", "CD2"],
    "Memory T": ["CCL4", "GZMK", "FOS", "SH2D1A", "KLF4"],
    "Treg": ["TNFRSF4", "TNFRSF18", "CTLA4", "ICOS", "NFKBIA"],
    "T IFN": ["ISG15", "MX2", "OAS1", "IFI44", "IFI44L"],
    "ILC2": ["IL13", "IL17RB", "GATA3", "PLAC8A", "ARL4C"],
    "Cycling T": ["MKI67", "TOP2A", "STMN1", "CENPF", "TUBA1B"],
}

EPITHELIAL_LEVEL2_MARKERS: dict[str, list[str]] = {
    "Enterocyte 1": ["UGT1A6", "SI", "MYO1A", "FAM13A", "SLC22A4", "SLC43A2"],
    "Enterocyte 2": ["APOC3", "APOA4", "APOA1", "APOB", "TM4SF20", "CLCA4", "SLC40A1"],
    "Enterocyte 3": ["FABP1", "MGST3", "COX4I1", "RBP2", "GSTP1", "SLC25A5", "SLC25A6"],
    "BEST4+ enterocyte": ["GUCA2A", "GUCA2B", "CFTR", "STOM", "CA7", "BEST4"],
    "Goblet": ["ZG16", "SPINK4", "CLCA1", "BCAS1", "LYPD8", "PNLIP", "AGR2"],
    "Tuft": ["ANXA4", "CA2", "TRPM5", "IRAG2", "RYR1", "IL17RB", "ATP2A3"],
    "IFN enterocyte": ["ISG15", "RSAD2", "IFIT1", "DDX60", "RNF213", "OAS1"],
    "Enteroendocrine": ["CHGB", "ADGRG4", "SCG2", "PCSK1N", "TPH1"],
}

MYELOID_LEVEL2_MARKERS: dict[str, list[str]] = {
    "Neutrophil": ["S100A12", "S100A8", "SERPINA1", "SOD2", "NFKBIA", "SELL", "ISG20"],
    "Eosinophil": ["MMP1", "MS4A2", "CAT", "PGF", "CHI3L1", "ADAMDEC1", "FABP3"],
    "Monocyte": ["MT1E", "MT2A", "HMOX1", "CTSS", "DLA-DRA"],
    "Macrophage": ["CCL3", "MAFB", "C1QA", "C1QC", "STAB1", "CTSS", "C1QB", "MSR1", "CD163"],
    "cDC1": ["BATF3", "ECRG4", "CLEC1B", "NAPSA", "DNASE1L3", "FLT3"],
    "IL22RA2 DC": ["IL22RA2", "CD52", "ECRG4", "MAL", "S100A5", "FBP1", "CA2", "LSP1", "FLT3"],
}

EPITHELIAL_STRESS_MARKERS = [
    "HSPA1A", "HSPA1B", "HSPA2", "HSPA5", "HSP90AA1", "ATF3", "JUN", "FOS", "DDIT3", "XBP1"
]
EPITHELIAL_METABOLIC_MARKERS = [
    "APOA4", "APOA1", "APOC3", "FABP1", "FABP2", "RBP2", "ALDOB", "SI", "MGAM", "SLC25A5"
]


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def configure_logger(path: Path, name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter(
        fmt="%(asctime)sZ %(levelname)s %(message)s", datefmt="%Y-%m-%dT%H:%M:%S"
    )
    formatter.converter = time.gmtime
    file_handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def peak_ram_mb(process: psutil.Process | None = None) -> float:
    process = process or psutil.Process(os.getpid())
    memory = process.memory_info()
    return getattr(memory, "peak_wset", memory.rss) / 1024**2


def validate_phase1_clean(adata: ad.AnnData) -> None:
    if adata.shape != EXPECTED_SHAPE:
        raise RuntimeError(f"Phase1 clean shape mismatch: {adata.shape}")
    conditions = adata.obs["condition"].astype(str).value_counts().to_dict()
    if conditions != EXPECTED_CONDITIONS:
        raise RuntimeError(f"Phase1 condition counts mismatch: {conditions}")
    donors = set(adata.obs["donor_id"].astype(str))
    if donors != EXPECTED_DONORS:
        raise RuntimeError(f"Phase1 donor set mismatch: {donors}")
    required_obs = {
        "sample_id", "donor_id", "condition", "total_counts", "n_genes_by_counts",
        "pct_counts_mt", "scrublet_score", "predicted_doublet",
    }
    required_var = {"feature_id", "gene_symbol", "feature_type"}
    if not required_obs.issubset(adata.obs.columns):
        raise RuntimeError(f"Missing obs fields: {required_obs - set(adata.obs.columns)}")
    if not required_var.issubset(adata.var.columns):
        raise RuntimeError(f"Missing var fields: {required_var - set(adata.var.columns)}")
    if "counts" not in adata.layers or not sparse.issparse(adata.layers["counts"]):
        raise RuntimeError("Missing sparse counts layer")
    counts_data = adata.layers["counts"].data
    if counts_data.size and (
        counts_data.min() < 0 or not np.allclose(counts_data, np.rint(counts_data))
    ):
        raise RuntimeError("Counts layer is not non-negative integer-valued")
    if not sparse.issparse(adata.X):
        raise RuntimeError("Expected sparse log-normalized X")
    normalization = adata.uns.get("normalization", {})
    if normalization.get("method") != "scanpy.pp.normalize_total followed by scanpy.pp.log1p":
        raise RuntimeError(f"Unexpected X normalization record: {normalization}")


def marker_var_names(adata: ad.AnnData, genes: Iterable[str]) -> list[str]:
    requested = set(genes)
    return adata.var_names[adata.var["gene_symbol"].astype(str).isin(requested)].tolist()


def marker_availability_rows(
    adata: ad.AnnData,
    source: str,
    marker_sets: dict[str, list[str]],
) -> list[dict[str, object]]:
    symbols = set(adata.var["gene_symbol"].astype(str))
    rows: list[dict[str, object]] = []
    for label, genes in marker_sets.items():
        for gene in genes:
            rows.append(
                {
                    "marker_source": source,
                    "marker_set": label,
                    "gene": gene,
                    "present_in_canFam3_1_features": gene in symbols,
                    "n_matching_features": int((adata.var["gene_symbol"].astype(str) == gene).sum()),
                }
            )
    return rows


def group_signature_scores(
    adata: ad.AnnData,
    groupby: str,
    marker_sets: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return discriminative z-score, mean expression, and detection per group/signature."""
    groups = pd.Index(sorted(adata.obs[groupby].astype(str).unique()), name=groupby)
    codes = pd.Categorical(adata.obs[groupby].astype(str), categories=groups).codes
    union = sorted({gene for genes in marker_sets.values() for gene in genes})
    gene_vars: dict[str, list[int]] = {}
    symbols = adata.var["gene_symbol"].astype(str).to_numpy()
    for gene in union:
        indices = np.flatnonzero(symbols == gene).tolist()
        if indices:
            gene_vars[gene] = indices
    gene_order = sorted(gene_vars)
    means = np.zeros((len(groups), len(gene_order)), dtype=np.float64)
    detection = np.zeros_like(means)
    sizes = np.bincount(codes, minlength=len(groups)).astype(np.float64)
    for gene_index, gene in enumerate(gene_order):
        indices = gene_vars[gene]
        vector = adata.X[:, indices]
        if sparse.issparse(vector):
            values = np.asarray(vector.sum(axis=1)).ravel()
        else:
            values = np.asarray(vector).sum(axis=1)
        means[:, gene_index] = np.bincount(codes, weights=values, minlength=len(groups)) / sizes
        detection[:, gene_index] = np.bincount(
            codes, weights=(values > 0).astype(float), minlength=len(groups)
        ) / sizes
    mean_df = pd.DataFrame(means, index=groups, columns=gene_order)
    det_df = pd.DataFrame(detection, index=groups, columns=gene_order)
    std = mean_df.std(axis=0, ddof=0).replace(0, np.nan)
    z = mean_df.subtract(mean_df.mean(axis=0), axis=1).divide(std, axis=1).fillna(0)
    score = pd.DataFrame(index=groups)
    raw = pd.DataFrame(index=groups)
    detect = pd.DataFrame(index=groups)
    for label, genes in marker_sets.items():
        present = [gene for gene in genes if gene in mean_df.columns]
        if not present:
            score[label] = -np.inf
            raw[label] = 0.0
            detect[label] = 0.0
        else:
            score[label] = z[present].mean(axis=1)
            raw[label] = mean_df[present].mean(axis=1)
            detect[label] = det_df[present].mean(axis=1)
    return score, raw, detect


def classify_groups(
    score: pd.DataFrame,
    raw: pd.DataFrame,
    detection: pd.DataFrame,
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for group in score.index:
        ranked = score.loc[group].sort_values(ascending=False)
        best, second = ranked.index[:2]
        margin = float(ranked.iloc[0] - ranked.iloc[1])
        best_detection = float(detection.loc[group, best])
        label = best
        confidence = "High" if margin >= 0.60 and best_detection >= 0.15 else "Medium"
        if margin < 0.15 or best_detection < 0.05:
            confidence = "Low"
        if margin < 0.03 and best_detection < 0.10:
            label = "Uncertain"
            confidence = "Uncertain"
        records.append(
            {
                "cluster": str(group),
                "predicted_label": label,
                "annotation_confidence": confidence,
                "top_score": float(ranked.iloc[0]),
                "second_label": second,
                "second_score": float(ranked.iloc[1]),
                "score_margin": margin,
                "mean_marker_expression": float(raw.loc[group, best]),
                "mean_marker_detection": best_detection,
            }
        )
    return pd.DataFrame(records).set_index("cluster")


def atomic_write_h5ad(adata: ad.AnnData, destination: Path) -> None:
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite formal output: {destination}")
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    if temporary.exists():
        temporary.unlink()
    adata.write_h5ad(temporary, compression="gzip")
    check = ad.read_h5ad(temporary, backed="r")
    try:
        if check.shape != adata.shape:
            raise RuntimeError("Temporary h5ad shape mismatch after write")
    finally:
        check.file.close()
    os.replace(temporary, destination)

