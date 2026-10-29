"""Freeze original canine feature identities and Phase 2 primary-source evidence."""

from __future__ import annotations

import importlib.metadata
import logging
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd
import psutil

CODE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_ROOT))
PROJECT_ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"] )

from config.paths import LOG_DIR, METADATA_DIR  # noqa: E402
from config.phase2_config import (  # noqa: E402
    AUTHOR_SIGNATURES,
    CLEAN_H5AD,
    GENE_MAPPING_MASTER,
    GEO_URL,
    GITHUB_TAG_URL,
    MARKER_AVAILABILITY,
    PAPER_DOI,
    PAPER_URL,
    PHASE2_ENV_DIR,
    REFERENCE_SOURCES,
    ZENODO_URL,
    ensure_phase2_directories,
)
from Stage02DogAnnotation.phase2_utils import (  # noqa: E402
    EPITHELIAL_LEVEL2_MARKERS,
    LEVEL1_MARKERS,
    MYELOID_LEVEL2_MARKERS,
    T_LEVEL2_MARKERS,
    configure_logger,
    marker_availability_rows,
    peak_ram_mb,
    sha256_file,
    validate_phase1_clean,
)


SCRIPT_VERSION = "1.0.0"


def version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "NOT INSTALLED"


def main() -> None:
    ensure_phase2_directories()
    outputs = [GENE_MAPPING_MASTER, REFERENCE_SOURCES, MARKER_AVAILABILITY]
    existing = [path for path in outputs if path.exists()]
    if existing:
        raise FileExistsError(f"Refusing to overwrite Phase2 outputs: {existing}")
    logger = configure_logger(LOG_DIR / "phase2_gene_mapping.log", "phase2_gene_mapping")
    process = psutil.Process(os.getpid())
    started = datetime.now(timezone.utc)
    logger.info(
        "START script_version=%s python=%s seed=NA input=%s",
        SCRIPT_VERSION,
        platform.python_version(),
        CLEAN_H5AD,
    )
    try:
        adata = ad.read_h5ad(CLEAN_H5AD)
        validate_phase1_clean(adata)
        symbols = adata.var["gene_symbol"].astype(str)
        duplicate_n = symbols.map(symbols.value_counts()).astype(int)
        mitochondrial = symbols.str.startswith("MT-") | symbols.eq("ND4L")
        mapping = pd.DataFrame(
            {
                "feature_id": adata.var["feature_id"].astype(str).to_numpy(),
                "original_gene_symbol": symbols.to_numpy(),
                "unique_var_name": adata.var_names.astype(str),
                "reference_version": "CanFam3.1",
                "ensembl_release": "104",
                "is_duplicate_symbol": duplicate_n.gt(1).to_numpy(),
                "duplicate_symbol_n": duplicate_n.to_numpy(),
                "mitochondrial_flag": mitochondrial.to_numpy(),
                "current_dog_ensembl_id": pd.NA,
                "current_dog_symbol": pd.NA,
                "human_1to1_ortholog": pd.NA,
                "mapping_status": "original_feature_frozen; current/ortholog mapping deferred",
                "mapping_source": "GSE254005 processed features; raw and clean h5ad var metadata",
                "mapping_version": "Phase2_v1_CanFam3.1_Ensembl104",
                "notes": "No Human ortholog guessed; duplicate symbols retained as separate feature IDs",
            }
        )
        mapping.to_csv(GENE_MAPPING_MASTER, index=False)

        author = pd.read_csv(AUTHOR_SIGNATURES)
        author_sets = {
            str(label).strip(): group["gene"].astype(str).tolist()
            for label, group in author.groupby("celltype.l2", observed=True)
        }
        availability: list[dict[str, object]] = []
        for source, marker_sets in (
            ("Phase2 predefined Level1 canine/paper markers", LEVEL1_MARKERS),
            ("Manchester et al. 2024 / author v0.1.0 signatures", author_sets),
            ("Paper-derived T Level2 markers", T_LEVEL2_MARKERS),
            ("Paper-derived Epithelial Level2 markers", EPITHELIAL_LEVEL2_MARKERS),
            ("Paper-derived Myeloid Level2 markers", MYELOID_LEVEL2_MARKERS),
        ):
            availability.extend(marker_availability_rows(adata, source, marker_sets))
        pd.DataFrame(availability).to_csv(MARKER_AVAILABILITY, index=False)

        code_zip = METADATA_DIR / "phase2_reference_code" / "canine_duodenal_atlas-v0.1.0.zip"
        sources = pd.DataFrame(
            [
                {
                    "source_type": "primary_article",
                    "title": "Single cell transcriptomic analysis of the canine duodenum in chronic inflammatory enteropathy and health",
                    "url": PAPER_URL,
                    "version_or_accessed": "Frontiers in Immunology 15 (2024); DOI " + PAPER_DOI,
                    "local_path": pd.NA,
                    "sha256": pd.NA,
                    "use_in_phase2": "primary annotation and positive-control evidence",
                },
                {
                    "source_type": "official_dataset",
                    "title": "NCBI GEO GSE254005",
                    "url": GEO_URL,
                    "version_or_accessed": "accessed 2026-08-24",
                    "local_path": pd.NA,
                    "sha256": pd.NA,
                    "use_in_phase2": "sample and reference metadata",
                },
                {
                    "source_type": "author_code",
                    "title": "dyammons/canine_duodenal_atlas",
                    "url": GITHUB_TAG_URL,
                    "version_or_accessed": "v0.1.0",
                    "local_path": str(code_zip),
                    "sha256": sha256_file(code_zip),
                    "use_in_phase2": "author canine marker definitions; no RDS labels copied",
                },
                {
                    "source_type": "author_archive_metadata",
                    "title": "Zenodo record 11434076",
                    "url": ZENODO_URL,
                    "version_or_accessed": "v1.0.0",
                    "local_path": pd.NA,
                    "sha256": pd.NA,
                    "use_in_phase2": "resource provenance only; 4.6 GB annotated RDS not downloaded",
                },
            ]
        )
        sources.to_csv(REFERENCE_SOURCES, index=False)

        software = pd.DataFrame(
            [
                ("python", platform.python_version()),
                ("scanpy", version("scanpy")),
                ("anndata", version("anndata")),
                ("numpy", version("numpy")),
                ("pandas", version("pandas")),
                ("scipy", version("scipy")),
                ("h5py", version("h5py")),
                ("igraph", version("igraph")),
                ("leidenalg", version("leidenalg")),
                ("scikit-misc", version("scikit-misc")),
                ("scikit-learn", version("scikit-learn")),
                ("umap-learn", version("umap-learn")),
            ],
            columns=["software", "version"],
        )
        software.to_csv(METADATA_DIR / "PHASE2_SOFTWARE_VERSIONS.csv", index=False)
        freeze = subprocess.run(
            [sys.executable, "-m", "pip", "freeze", "--all"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.splitlines()
        environment_lines = [
            "name: canine_cie_phase2",
            "channels:",
            "  - conda-forge",
            "dependencies:",
            f"  - python={platform.python_version()}",
            "  - pip",
            "  - pip:",
            *[f"      - {item}" for item in freeze if item.strip()],
        ]
        (METADATA_DIR / "phase2_environment.yml").write_text(
            "\n".join(environment_lines) + "\n", encoding="utf-8"
        )
        duplicate_symbols = int(symbols[symbols.duplicated(keep=False)].nunique())
        logger.info(
            "OUTPUT mapping=%s rows=%d duplicate_symbols=%d mt_features=%d marker_audit_rows=%d",
            GENE_MAPPING_MASTER,
            len(mapping),
            duplicate_symbols,
            int(mitochondrial.sum()),
            len(availability),
        )
        logger.info(
            "REFERENCE author_code_sha256=%s annotated_RDS_downloaded=False human_ortholog_mapping=deferred",
            sha256_file(code_zip),
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
