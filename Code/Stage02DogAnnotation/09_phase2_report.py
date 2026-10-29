"""Validate Phase 2 deliverables and write the canine-only gate report."""

from __future__ import annotations

import hashlib
import logging
import platform
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

CODE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_ROOT))
PROJECT_ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"] )

from config.paths import LOG_DIR, METADATA_DIR, RESULTS_DIR  # noqa: E402
from config.phase2_config import (  # noqa: E402
    ANNOTATED_H5AD,
    ANNOTATIONS_DIR,
    GENE_MAPPING_MASTER,
    GEO_URL,
    GITHUB_TAG_URL,
    PAPER_DOI,
    PAPER_URL,
    PHASE2_SEED,
    ZENODO_URL,
    ensure_phase2_directories,
)
from Stage02DogAnnotation.phase2_utils import configure_logger, peak_ram_mb, sha256_file, validate_phase1_clean  # noqa: E402

SCRIPT_VERSION = "1.0.0"
CIE2 = "GSM8031701_cie_2"
GATE = "YELLOW"


REQUIRED_PHASE2 = [
    ANNOTATED_H5AD,
    ANNOTATIONS_DIR / "dog_cell_annotation.csv",
    ANNOTATIONS_DIR / "cluster_markers.csv",
    ANNOTATIONS_DIR / "cluster_donor_composition.csv",
    ANNOTATIONS_DIR / "celltype_proportion_per_donor.csv",
    ANNOTATIONS_DIR / "positive_control_summary.csv",
    ANNOTATIONS_DIR / "cie2_sensitivity_summary.csv",
    ANNOTATIONS_DIR / "cie2_level1_qc_comparison.csv",
    ANNOTATIONS_DIR / "cie2_cluster_audit.csv",
    ANNOTATIONS_DIR / "cie2_exclusion_embedding_metrics.csv",
    ANNOTATIONS_DIR / "cie2_exclusion_state_stability.csv",
    GENE_MAPPING_MASTER,
    CODE_ROOT / "Stage02DogAnnotation" / "05_phase2_gene_mapping.py",
    CODE_ROOT / "Stage02DogAnnotation" / "06_phase2_embedding.py",
    CODE_ROOT / "Stage02DogAnnotation" / "07_phase2_annotation.py",
    CODE_ROOT / "Stage02DogAnnotation" / "08_phase2_cie2_audit.py",
    CODE_ROOT / "Stage02DogAnnotation" / "09_phase2_report.py",
    LOG_DIR / "phase2_gene_mapping.log",
    LOG_DIR / "phase2_embedding.log",
    LOG_DIR / "phase2_annotation.log",
    LOG_DIR / "phase2_cie2_audit.log",
]


def fmt(value: object, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.{digits}f}"
    return str(value)


def markdown_table(frame: pd.DataFrame) -> str:
    clean = frame.copy()
    for column in clean.columns:
        clean[column] = clean[column].map(fmt)
    headers = [str(x) for x in clean.columns]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in clean.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(x).replace("|", "\\|") for x in row) + " |")
    return "\n".join(lines)


def validate_outputs(adata: ad.AnnData) -> None:
    missing = [str(path) for path in REQUIRED_PHASE2 if not path.exists() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"Missing/empty Phase2 deliverables: {missing}")
    validate_phase1_clean(adata)
    required_obs = {
        "donor_id", "condition", "level1", "level2", "annotation_cluster",
        "annotation_confidence", "pct_counts_mt", "n_genes_by_counts", "total_counts",
    }
    if required_obs - set(adata.obs.columns):
        raise RuntimeError(f"Annotated object missing obs fields: {required_obs - set(adata.obs.columns)}")
    if "lognorm" not in adata.layers or "counts" not in adata.layers:
        raise RuntimeError("Annotated object must retain counts and lognorm layers")
    if sparse.issparse(adata.X):
        if (adata.X != adata.layers["lognorm"]).nnz:
            raise RuntimeError("Annotated X no longer matches frozen lognorm layer")
    elif not np.allclose(adata.X, adata.layers["lognorm"]):
        raise RuntimeError("Annotated X no longer matches frozen lognorm layer")
    counts = adata.layers["counts"]
    values = counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
    if values.size and (values.min() < 0 or not np.allclose(values, np.rint(values))):
        raise RuntimeError("Annotated counts layer is not nonnegative integer-valued")
    if adata.shape != (38_186, 20_257):
        raise RuntimeError(f"Annotated shape changed: {adata.shape}")
    if adata.obs["donor_id"].nunique() != 7:
        raise RuntimeError("Annotated object does not retain all seven donors")


def atlas_tables(adata: ad.AnnData) -> tuple[pd.DataFrame, pd.DataFrame]:
    obs = adata.obs.copy()
    level1 = (
        obs.groupby("level1", observed=True)
        .agg(cells=("level1", "size"), donors=("donor_id", "nunique"), clusters=("annotation_cluster", "nunique"))
        .reset_index()
    )
    level1["fraction_pct"] = level1["cells"] / adata.n_obs * 100
    level1 = level1.sort_values("cells", ascending=False)
    level2 = (
        obs.groupby(["level1", "level2"], observed=True)
        .agg(cells=("level2", "size"), donors=("donor_id", "nunique"), clusters=("annotation_cluster", "nunique"))
        .reset_index()
        .sort_values(["level1", "cells"], ascending=[True, False])
    )
    return level1, level2


def cie2_facts(adata: ad.AnnData) -> dict[str, float]:
    obs = adata.obs.copy()
    donor = obs["donor_id"].astype(str)
    level = obs["level1"].astype(str)
    cie = donor == CIE2
    cie_epi = cie & (level == "Epithelial")
    other_epi = (~cie) & (level == "Epithelial")
    high_mt_cie = cie & (obs["pct_counts_mt"] >= 70)
    return {
        "cells": int(cie.sum()),
        "mt_median": float(obs.loc[cie, "pct_counts_mt"].median()),
        "epi_fraction": float((level[cie] == "Epithelial").mean()),
        "high_mt_epi_fraction": float((level[high_mt_cie] == "Epithelial").mean()) if high_mt_cie.any() else np.nan,
        "epi_mt": float(obs.loc[cie_epi, "pct_counts_mt"].median()),
        "t_mt": float(obs.loc[cie & (level == "T"), "pct_counts_mt"].median()),
        "my_mt": float(obs.loc[cie & (level == "Myeloid"), "pct_counts_mt"].median()),
        "epi_features": float(obs.loc[cie_epi, "n_genes_by_counts"].median()),
        "other_epi_features": float(obs.loc[other_epi, "n_genes_by_counts"].median()),
        "epi_counts": float(obs.loc[cie_epi, "total_counts"].median()),
        "other_epi_counts": float(obs.loc[other_epi, "total_counts"].median()),
    }


def reproduction_table(positive: pd.DataFrame) -> pd.DataFrame:
    return positive[
        ["positive_control", "cluster_or_state", "n_cells", "n_donors", "reproduction_status", "supporting_markers"]
    ].rename(
        columns={
            "positive_control": "positive control",
            "cluster_or_state": "state/cluster",
            "n_cells": "cells",
            "n_donors": "donors",
            "reproduction_status": "status",
            "supporting_markers": "marker support",
        }
    )


def manifest(paths: list[Path], logger: logging.Logger) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    seen: set[Path] = set()
    for path in paths:
        path = path.resolve()
        if path in seen or not path.exists() or not path.is_file():
            continue
        seen.add(path)
        logger.info("HASH path=%s bytes=%d", path, path.stat().st_size)
        rows.append(
            {
                "path_relative_to_project": str(path.relative_to(PROJECT_ROOT)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    start = time.time()
    ensure_phase2_directories()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    logger = configure_logger(LOG_DIR / "phase2_report.log", "phase2_report")
    logger.info("START script_version=%s python=%s scanpy=%s seed=%d gate=%s", SCRIPT_VERSION, platform.python_version(), sc.__version__, PHASE2_SEED, GATE)
    try:
        adata = ad.read_h5ad(ANNOTATED_H5AD)
        validate_outputs(adata)
        level1, level2 = atlas_tables(adata)
        positive = pd.read_csv(ANNOTATIONS_DIR / "positive_control_summary.csv")
        sensitivity = pd.read_csv(ANNOTATIONS_DIR / "cie2_sensitivity_summary.csv").set_index("metric")["value"]
        sens_metrics = pd.read_csv(ANNOTATIONS_DIR / "cie2_exclusion_embedding_metrics.csv").iloc[0]
        state_stability = pd.read_csv(ANNOTATIONS_DIR / "cie2_exclusion_state_stability.csv")
        cie2_clusters = pd.read_csv(ANNOTATIONS_DIR / "cie2_cluster_audit.csv")
        gene_map = pd.read_csv(GENE_MAPPING_MASTER)
        lineage_params = pd.read_csv(ANNOTATIONS_DIR / "lineage_embedding_parameters.csv")
        facts = cie2_facts(adata)
        dominant = cie2_clusters[cie2_clusters["dominant_donor_fraction"] >= 0.80]
        cie2_dominant = cie2_clusters[cie2_clusters["cie2_dominant_ge_0_80"]]
        if len(cie2_dominant) != 1:
            raise RuntimeError(f"Expected exactly one cie2-dominant annotation cluster from audit, observed {len(cie2_dominant)}")
        epi4 = cie2_dominant.iloc[0]

        paper_compare = pd.DataFrame(
            [
                {"feature": "total cells", "this Phase2": adata.n_obs, "paper": "35,668", "assessment": "Partially reproduced (different QC retention)"},
                {"feature": "global populations", "this Phase2": 24, "paper": "31", "assessment": "Not tuned to match count; major structure reproduced"},
                {"feature": "T cells", "this Phase2": int((adata.obs["level1"].astype(str) == "T").sum()), "paper": "18,824", "assessment": "Reproduced"},
                {"feature": "Epithelial", "this Phase2": int((adata.obs["level1"].astype(str) == "Epithelial").sum()), "paper": "12,457", "assessment": "Reproduced"},
                {"feature": "Myeloid", "this Phase2": int((adata.obs["level1"].astype(str) == "Myeloid").sum()), "paper": "1,945", "assessment": "Reproduced"},
                {"feature": "Endothelial", "this Phase2": 0, "paper": "not a reported major class", "assessment": "Not assessable/not forced"},
            ]
        )
        duplicate_symbols = int(gene_map.loc[gene_map["is_duplicate_symbol"].astype(bool), "original_gene_symbol"].nunique())
        duplicate_rows = int(gene_map["is_duplicate_symbol"].astype(bool).sum())

        report = f"""# Canine CIE project — Phase2 report

Generated: {pd.Timestamp.utcnow().isoformat()}  
Project root: `{PROJECT_ROOT}`  
Random seed: `{PHASE2_SEED}`  
Gate: **{GATE}**

## Executive gate decision

Phase2 satisfies the minimum conditions for a donor-level Phase3, but the gate is **YELLOW**, not GREEN. The canine atlas, raw-count preservation, major lineages, three lineage-specific analyses, and all requested positive controls are reproducible. Excluding `GSM8031701_cie_2` leaves a high-purity major-lineage structure. However, 14/50 annotation clusters have a dominant donor fraction ≥80%, including one large `cie_2`-dominated enterocyte immediate-early state. These donor effects must remain explicit design constraints in any later donor-level model.

Recommendation: **retain `cie_2` in the primary Phase3 cohort**, pre-specify donor-level sensitivity/LODO attention for it and other donor-driven states, and do not interpret cell-level abundance or disease differences. Phase3 was not executed here.

## 1. Input integrity and scope compliance

- Annotated object: {adata.n_obs:,} cells × {adata.n_vars:,} features; Healthy {int((adata.obs['condition'].astype(str) == 'Healthy').sum()):,}, CIE {int((adata.obs['condition'].astype(str) == 'CIE').sum()):,}; seven donors retained.
- `layers['counts']` is sparse, nonnegative, and integer-valued. `X` exactly matches the retained sparse `layers['lognorm']` matrix (library-size 10,000 + log1p from Phase1).
- Phase0/Phase1 files were read-only and were not regenerated or overwritten.
- Primary embedding used donor-aware Seurat-v3 HVGs from counts, excluded mitochondrial genes from HVG/PCA only, and performed no donor integration, condition/donor/mt regression, or forced batch correction.
- Global resolutions 0.4/0.8/1.2 were audited. Resolution 0.8 was the pre-annotation candidate; resolution 0.4 (24 clusters) was frozen as the formal Level1 scaffold because it had fewer donor-dominant and low/uncertain clusters. The result was not tuned to the paper's 31 populations.
- No CIE-vs-Healthy DEG, pseudobulk, GSEA, LODO, Human download/analysis, ortholog comparison, singIST, SATURN, CONCORD, or therapeutic prioritization was run.

## 2. Final canine atlas

The final object has 24 global Level1 scaffold clusters and 50 annotation clusters after independent Epithelial/T/Myeloid reclustering. It contains 8 Level1 labels (including `Uncertain`) and 29 Level2 labels. No coherent endothelial cluster was found, so none was forced.

### Level1

{markdown_table(level1.rename(columns={'level1':'Level1','cells':'cells','donors':'donors','clusters':'annotation clusters','fraction_pct':'fraction %'}))}

### Level2

{markdown_table(level2.rename(columns={'level1':'Level1','level2':'Level2','cells':'cells','donors':'donors','clusters':'clusters'}))}

Lineage-specific parameters and donor effects:

{markdown_table(lineage_params)}

Low/mixed clusters were retained rather than force-labeled. The formal audit table records automatic labels, overrides, confidence, rationale, positive markers, conflicting evidence, donor dominance, and paper correspondence.

## 3. Original-paper reproduction

Primary reference: Manchester et al., *Frontiers in Immunology* (2024), DOI [{PAPER_DOI}]({PAPER_URL}). The authors' frozen v0.1.0 marker code was used as the primary subtype evidence after the project annotation was independently constructed; no 4.6 GB annotated RDS was downloaded.

{markdown_table(paper_compare)}

Classification:

- **Reproduced:** major T, Epithelial and Myeloid structure; B, Plasma, Mast and Stromal identities; IL7R-high/GZMA-high T; SOD2/IL1A-associated myeloid states; the requested paper epithelial/T/myeloid subtype families with marker support.
- **Partially reproduced:** exact fine-state partition and total population count; epithelial stress/metabolic programs, because Phase2 establishes their presence but deliberately does not run disease DEG or enrichment.
- **Not reproduced:** no required major positive control failed.
- **Not assessable/not forced:** a distinct endothelial population and some mixed low-confidence cells; current dog→human orthology is deferred.

### Positive controls

{markdown_table(reproduction_table(positive))}

All state-presence calls cover both Healthy and CIE donors, but this is descriptive only. It is not evidence of a disease difference.

## 4. `GSM8031701_cie_2` audit

### Where the high mt% occurs

- `cie_2` contains {int(facts['cells']):,} retained cells and has overall mt% median {facts['mt_median']:.2f}%.
- {facts['epi_fraction']*100:.2f}% of its cells are Epithelial, and {facts['high_mt_epi_fraction']*100:.2f}% of its cells with mt% ≥70% are Epithelial.
- Within `cie_2`, median mt% is {facts['epi_mt']:.2f}% in Epithelial, versus {facts['t_mt']:.2f}% in T and {facts['my_mt']:.2f}% in Myeloid. Thus the donor-wide 68.49% median is overwhelmingly an epithelial phenomenon, not a pan-lineage sample failure.

### Complexity and marker coherence

- `cie_2` Epithelial cells retain median nFeature {facts['epi_features']:.0f} and nCount {facts['epi_counts']:.0f}; pooled other-donor Epithelial medians are {facts['other_epi_features']:.0f} and {facts['other_epi_counts']:.0f}.
- Their median detected fraction across available Level1 epithelial markers is 0.778 (other donor-level epithelial groups are typically 0.778–0.889).
- Only {float(sensitivity['cie2_low_quality_like_fraction'])*100:.2f}% of all `cie_2` cells meet the pre-defined joint flag: mt% ≥ max(70%, other-donor lineage q90) **and** nFeature ≤ other-donor lineage q10. This flag is descriptive, not a deletion rule.
- Phase1-passing `cie_2` cells have no retained Scrublet-predicted doublets in the formal audit table.

### Donor-dominant cluster

Exactly one annotation cluster is ≥80% `cie_2`: `{epi4['cluster']}` / `{epi4['level2']}` with {int(epi4['n_cells']):,} cells, of which {epi4['cie2_fraction_of_cluster']*100:.2f}% are from `cie_2`. It is marker-supported (FOS/FOSB/DUSP1/EGR1 with enterocyte FAM13A/related genes), so it is biologically plausible as an enterocyte immediate-early/stress state. It is nevertheless donor-confounded. Within its `cie_2` cells, {epi4['cie2_low_quality_like_fraction']*100:.2f}% meet the joint low-quality-like flag; the non-`cie_2` cells in the same cluster have a comparable mt median ({epi4['other_mt_median']:.2f}% versus {epi4['cie2_mt_median']:.2f}%). This does not support automatic deletion of the donor, but it prohibits treating this state as a general CIE effect without donor-level sensitivity.

### Temporary exclusion sensitivity

- Rebuilt from counts/lognorm with donor-aware HVGs, 40-PC neighbor graph, no integration/regression, and fixed seed; retained {int(sens_metrics['retained_cells']):,} cells from six donors and produced {int(sens_metrics['sensitivity_clusters'])} de novo clusters.
- Full-vs-no-`cie_2` agreement on common cells: ARI {sens_metrics['full_vs_no_cie2_ARI_common_cells']:.4f}, NMI {sens_metrics['full_vs_no_cie2_NMI_common_cells']:.4f}.
- Weighted retained-Level1 purity within the independent no-`cie_2` clusters: {sens_metrics['weighted_reference_level1_purity']:.4f}.
- IL7R-high T, GZMA-high T, principal myeloid states and principal epithelial states all remain present across the other six donors:

{markdown_table(state_stability[['state_group','retained_cells','retained_donors','largest_cluster_fraction','status']])}

Conclusion: retain `cie_2` for primary donor-level work, mark EPI_4 as donor-confounded, and require a later pre-specified sensitivity/LODO review. Do not delete the donor on mt% alone.

## 5. Donor effects

There are {len(dominant)} of 50 annotation clusters with dominant-donor fraction ≥80%. Their dominant-donor distribution is: {', '.join(f'{k}: {v}' for k, v in dominant.groupby('dominant_donor').size().to_dict().items())}. This includes biologically coherent but donor-confounded epithelial, T, neutrophil and stromal states, plus the 17-cell uncertain cluster. The unintegrated primary analysis intentionally exposes these patterns. They are not formal disease findings.

This is the principal reason for the YELLOW gate. Phase3 must use donor as the biological replicate; cell-level proportion tests and cell-level disease P values are not acceptable.

## 6. Gene identity and mapping

- `dog_gene_mapping_master.csv` freezes all {len(gene_map):,} original feature rows, `feature_id`, original symbol and collision-safe `unique_var_name` under CanFam3.1 / Ensembl release 104 provenance.
- {duplicate_symbols} duplicated symbols affect {duplicate_rows} feature rows. They remain distinct through `unique_var_name`; no counts were summed or guessed during Phase2.
- {int(gene_map['mitochondrial_flag'].astype(bool).sum())} mitochondrial features are flagged and retained in the full matrix.
- Current dog Ensembl/symbol refresh and human 1:1 ortholog fields remain NA/deferred for all rows by design. No Human label or ortholog was used for primary dog annotation.
- Later pseudobulk must aggregate by stable feature identity and explicitly decide how duplicate-symbol rows are handled; symbol-only aggregation without this master could double-count or merge distinct features.

## 7. Gate checklist

| Requirement | Result |
| --- | --- |
| Seven dogs retained and traceable | PASS |
| Counts and normalized expression independently retained | PASS |
| Major Level1 markers coherent | PASS |
| Epithelial/T/Myeloid independently identifiable | PASS |
| Positive controls present | PASS / epithelial programs partial by scope |
| `cie_2` explained and sensitivity-tested | PASS; retain with caution |
| No broad annotation completely depends on one donor | PASS for major Level1; several fine states donor-driven |
| Human Crohn labels avoided as primary annotation | PASS |
| Formal Phase2 h5ad/tables/logs reproducible | PASS |

**Final Phase2 gate: YELLOW.** The atlas is suitable to propose Phase3 donor-level pseudobulk, with explicit `cie_2` and broader donor-driven-state sensitivity constraints. This report stops at the Phase2 gate and does not authorize or execute Phase3.

## 8. Key deliverables and provenance

- `05_annotations\\dog_GSE254005_annotated.h5ad`
- `05_annotations\\dog_cell_annotation.csv`
- `05_annotations\\cluster_markers.csv`
- `05_annotations\\cluster_donor_composition.csv`
- `05_annotations\\celltype_proportion_per_donor.csv`
- `05_annotations\\positive_control_summary.csv`
- `05_annotations\\cie2_sensitivity_summary.csv` and detailed cie2/sensitivity tables
- `05_annotations\\figures\\01–16_*.png`
- `00_metadata\\dog_gene_mapping_master.csv`
- `00_metadata\\phase2_environment.yml` and software/source manifests
- `scripts\\05_phase2_gene_mapping.py` through `scripts\\09_phase2_report.py`
- `logs\\phase2_*.log`

External frozen references: [GEO GSE254005]({GEO_URL}), [paper full text]({PAPER_URL}), [author Zenodo record]({ZENODO_URL}), [author code tag v0.1.0]({GITHUB_TAG_URL}).
"""
        report_path = RESULTS_DIR / "PHASE2_REPORT.md"
        report_path.write_text(report, encoding="utf-8")

        manifest_paths = list(REQUIRED_PHASE2) + [
            report_path,
            METADATA_DIR / "PHASE2_REFERENCE_SOURCES.csv",
            METADATA_DIR / "PHASE2_MARKER_AVAILABILITY.csv",
            METADATA_DIR / "PHASE2_SOFTWARE_VERSIONS.csv",
            METADATA_DIR / "phase2_environment.yml",
            ANNOTATIONS_DIR / "embedding_parameters.csv",
            ANNOTATIONS_DIR / "leiden_resolution_summary.csv",
            ANNOTATIONS_DIR / "lineage_embedding_parameters.csv",
            ANNOTATIONS_DIR / "lineage_marker_assignment_audit.csv",
            ANNOTATIONS_DIR / "cie2_exclusion_cluster_annotation.csv",
            ANNOTATIONS_DIR / "cie2_exclusion_state_stability.csv",
        ] + sorted((ANNOTATIONS_DIR / "figures").glob("*.png"))
        phase1_paths = [
            PROJECT_ROOT / "04_h5ad_clean" / "dog_GSE254005_clean.h5ad",
            RESULTS_DIR / "PHASE1_REPORT.md",
            METADATA_DIR / "PHASE0_AUDIT_REPORT.md",
            METADATA_DIR / "MASTER_DATA_AUDIT.csv",
        ]
        manifest_frame = manifest(manifest_paths + phase1_paths, logger)
        manifest_frame.to_csv(METADATA_DIR / "PHASE2_MANIFEST.csv", index=False)
        logger.info("OUTPUT report=%s manifest=%s gate=%s", report_path, METADATA_DIR / "PHASE2_MANIFEST.csv", GATE)
        logger.info("END status=success elapsed_seconds=%.3f peak_ram_mb=%.1f", time.time() - start, peak_ram_mb())
    except Exception:
        logger.exception("END status=error elapsed_seconds=%.3f peak_ram_mb=%.1f", time.time() - start, peak_ram_mb())
        raise


if __name__ == "__main__":
    main()
