"""Canonical configuration for authorized canine-only Phase 3 pseudobulk."""

from __future__ import annotations

from config.paths import LOG_DIR, METADATA_DIR, PROJECT_ROOT, RESULTS_DIR


PHASE3_SEED = 20240824
ANNOTATED_H5AD = PROJECT_ROOT / "05_annotations" / "dog_GSE254005_annotated.h5ad"
GENE_MAPPING_MASTER = METADATA_DIR / "dog_gene_mapping_master.csv"
PHASE2_MANIFEST = METADATA_DIR / "PHASE2_MANIFEST.csv"
PHASE2_REPORT = RESULTS_DIR / "PHASE2_REPORT.md"

PSEUDOBULK_DIR = PROJECT_ROOT / "06_pseudobulk"
PB_MATRICES_DIR = PSEUDOBULK_DIR / "matrices"
PB_METADATA_DIR = PSEUDOBULK_DIR / "metadata"
PB_DE_DIR = PSEUDOBULK_DIR / "de_results"
PB_SENSITIVITY_DIR = PSEUDOBULK_DIR / "sensitivity"
PB_FIGURES_DIR = PSEUDOBULK_DIR / "figures"
PB_RANKED_DIR = PSEUDOBULK_DIR / "ranked_lists"
PB_WORK_DIR = PROJECT_ROOT / "work" / "phase3"

R_PORTABLE_ROOT = PROJECT_ROOT / "envs" / "phase3_r_portable" / "app"
R_USER_LIBRARY = PROJECT_ROOT / "envs" / "phase3_r_portable" / "library_user"
RSCRIPT = R_PORTABLE_ROOT / "bin" / "x64" / "Rscript.exe"

PRIMARY_LINEAGES = ("Epithelial", "T", "Myeloid")
SECONDARY_LINEAGES = ("B", "Plasma", "Mast", "Stromal")
CONDITION_LEVELS = ("Healthy", "CIE")
CIE2_DONOR = "GSM8031701_cie_2"
MIN_STATE_CELLS_PER_DONOR = 20
MEANINGFUL_EFFECT_THRESHOLD = 0.25
TOP_EFFECT_N = 100


def ensure_phase3_directories() -> None:
    for path in (
        PSEUDOBULK_DIR,
        PB_MATRICES_DIR,
        PB_METADATA_DIR,
        PB_DE_DIR,
        PB_SENSITIVITY_DIR,
        PB_FIGURES_DIR,
        PB_RANKED_DIR,
        PB_WORK_DIR,
        R_USER_LIBRARY,
    ):
        path.mkdir(parents=True, exist_ok=True)
