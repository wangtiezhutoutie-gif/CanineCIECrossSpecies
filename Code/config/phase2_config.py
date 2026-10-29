"""Canonical configuration for the authorized canine-only Phase 2 analysis."""

from __future__ import annotations

from config.paths import LOG_DIR, METADATA_DIR, PROJECT_ROOT, RESULTS_DIR


PHASE2_SEED = 20240824
ANNOTATIONS_DIR = PROJECT_ROOT / "05_annotations"
ANNOTATION_FIGURES_DIR = ANNOTATIONS_DIR / "figures"
PHASE2_WORK_DIR = PROJECT_ROOT / "work" / "phase2"
PHASE2_ENV_DIR = PROJECT_ROOT / "envs" / "phase2"

CLEAN_H5AD = PROJECT_ROOT / "04_h5ad_clean" / "dog_GSE254005_clean.h5ad"
EMBEDDED_H5AD = ANNOTATIONS_DIR / "dog_GSE254005_embedded.h5ad"
ANNOTATED_H5AD = ANNOTATIONS_DIR / "dog_GSE254005_annotated.h5ad"

GENE_MAPPING_MASTER = METADATA_DIR / "dog_gene_mapping_master.csv"
REFERENCE_SOURCES = METADATA_DIR / "PHASE2_REFERENCE_SOURCES.csv"
MARKER_AVAILABILITY = METADATA_DIR / "PHASE2_MARKER_AVAILABILITY.csv"
AUTHOR_CODE_DIR = (
    METADATA_DIR
    / "phase2_reference_code"
    / "canine_duodenal_atlas-0.1.0"
)
AUTHOR_SIGNATURES = AUTHOR_CODE_DIR / "scripts" / "input" / "genesig_long.csv"

N_HVG = 2_500
N_PCA_COMPUTE = 50
N_PCS_NEIGHBORS = 40
N_NEIGHBORS = 30
UMAP_MIN_DIST = 0.30
LEIDEN_RESOLUTIONS = (0.4, 0.8, 1.2)
PRIMARY_LEIDEN_KEY = "leiden_0_8"

PAPER_DOI = "10.3389/fimmu.2024.1397590"
PAPER_URL = (
    "https://www.frontiersin.org/journals/immunology/articles/"
    "10.3389/fimmu.2024.1397590/full"
)
GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254005"
ZENODO_URL = "https://zenodo.org/records/11434076"
GITHUB_TAG_URL = "https://github.com/dyammons/canine_duodenal_atlas/tree/v0.1.0"


def ensure_phase2_directories() -> None:
    for path in (ANNOTATIONS_DIR, ANNOTATION_FIGURES_DIR, PHASE2_WORK_DIR):
        path.mkdir(parents=True, exist_ok=True)
