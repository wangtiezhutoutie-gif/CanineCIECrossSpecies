"""Canonical project paths. All Phase 1 scripts import paths from here."""

import os
from pathlib import Path

PROJECT_ROOT = Path(os.environ["CANINE_CIE_PROJECT_ROOT"])
METADATA_DIR = PROJECT_ROOT / "00_metadata"
DOG_RAW_DIR = PROJECT_ROOT / "01_raw_processed" / "dog_GSE254005"
H5AD_RAW_DIR = PROJECT_ROOT / "02_h5ad_raw"
QC_DIR = PROJECT_ROOT / "03_qc"
QC_FIGURES_DIR = QC_DIR / "figures"
H5AD_CLEAN_DIR = PROJECT_ROOT / "04_h5ad_clean"
RESULTS_DIR = PROJECT_ROOT / "08_results"
LOG_DIR = PROJECT_ROOT / "logs"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
CONFIG_DIR = PROJECT_ROOT / "config"


def ensure_phase1_directories() -> None:
    """Create only the directories authorized for Phase 1."""
    for path in (
        METADATA_DIR,
        DOG_RAW_DIR,
        H5AD_RAW_DIR,
        QC_FIGURES_DIR,
        H5AD_CLEAN_DIR,
        RESULTS_DIR,
        LOG_DIR,
        SCRIPTS_DIR,
        CONFIG_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)
