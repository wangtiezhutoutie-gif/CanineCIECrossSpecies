"""Audit, extract, and ingest the seven GSE254005 processed 10x matrices."""

from __future__ import annotations

import csv
import gzip
import hashlib
import logging
import os
import platform
import re
import shutil
import sys
import tarfile
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import psutil
import scipy.sparse as sp
from scipy.io import mmread

CODE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_ROOT))
PROJECT_ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"] )

from config.paths import DOG_RAW_DIR, H5AD_RAW_DIR, LOG_DIR, METADATA_DIR, ensure_phase1_directories

SCRIPT_VERSION = "1.1.0"
ARCHIVE = DOG_RAW_DIR / "GSE254005_RAW.tar"
EXTRACTED_DIR = DOG_RAW_DIR / "extracted"
RAW_H5AD = H5AD_RAW_DIR / "dog_GSE254005_raw.h5ad"
SAMPLE_MANIFEST = METADATA_DIR / "GSE254005_SAMPLE_MANIFEST.csv"
INGEST_SUMMARY = METADATA_DIR / "GSE254005_INGEST_SUMMARY.csv"
LOG_FILE = LOG_DIR / "phase1_ingest.log"

EXPECTED = {
    "h_1": ("GSM8031698", "Healthy"),
    "h_2": ("GSM8031699", "Healthy"),
    "h_3": ("GSM8031700", "Healthy"),
    "cie_2": ("GSM8031701", "CIE"),
    "cie_3": ("GSM8031702", "CIE"),
    "cie_4": ("GSM8031703", "CIE"),
    "cie_6": ("GSM8031704", "CIE"),
}
FILE_KINDS = ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
MEMBER_RE = re.compile(
    r"^(?P<gsm>GSM\d+)_(?P<sample>h_\d+|cie_\d+)_(?P<kind>barcodes\.tsv\.gz|features\.tsv\.gz|matrix\.mtx\.gz)$"
)


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


def atomic_csv(rows: list[dict[str, object]], path: Path, fields: list[str]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def audit_archive() -> tuple[dict[str, dict[str, tarfile.TarInfo]], list[tarfile.TarInfo]]:
    if not ARCHIVE.exists():
        raise FileNotFoundError(ARCHIVE)
    grouped: dict[str, dict[str, tarfile.TarInfo]] = {}
    with tarfile.open(ARCHIVE, "r:") as archive:
        members = archive.getmembers()
    if len(members) != 21:
        raise RuntimeError(f"Expected exactly 21 TAR members, found {len(members)}")
    for member in members:
        if not member.isfile():
            raise RuntimeError(f"Unexpected non-file TAR member: {member.name}")
        if member.size <= 0:
            raise RuntimeError(f"Zero-byte TAR member: {member.name}")
        if Path(member.name).name != member.name or Path(member.name).is_absolute():
            raise RuntimeError(f"Unsafe or nested TAR member: {member.name}")
        match = MEMBER_RE.fullmatch(member.name)
        if not match:
            raise RuntimeError(f"Unexpected TAR member name: {member.name}")
        gsm, sample, kind = match.group("gsm", "sample", "kind")
        if sample not in EXPECTED:
            raise RuntimeError(f"Unexpected sample prefix: {sample}")
        if gsm != EXPECTED[sample][0]:
            raise RuntimeError(f"GEO/sample mismatch: {gsm} vs {sample}")
        if kind in grouped.setdefault(sample, {}):
            raise RuntimeError(f"Duplicate member for {sample}/{kind}")
        grouped[sample][kind] = member
    if set(grouped) != set(EXPECTED):
        raise RuntimeError(f"Sample set mismatch: {sorted(grouped)}")
    for sample, files in grouped.items():
        if set(files) != set(FILE_KINDS):
            raise RuntimeError(f"Incomplete three-file set for {sample}: {sorted(files)}")
    logging.info("ARCHIVE_AUDIT donors=%s members=%s extras=0 zero_bytes=0", len(grouped), len(members))
    return grouped, members


def extract_verified_members(members: list[tarfile.TarInfo]) -> None:
    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)
    expected_names = {member.name for member in members}
    unexpected = [path.name for path in EXTRACTED_DIR.iterdir() if path.is_file() and path.name not in expected_names]
    if unexpected:
        raise RuntimeError(f"Unexpected files already present in extraction directory: {unexpected}")
    with tarfile.open(ARCHIVE, "r:") as archive:
        for member in members:
            target = EXTRACTED_DIR / member.name
            member_digest = sha256_tar_member(archive, member)
            if target.exists():
                if target.stat().st_size != member.size:
                    raise RuntimeError(f"Existing extracted file has wrong size: {target}")
                target_digest = sha256_path(target)
                if target_digest == member_digest:
                    logging.info("EXTRACT_SKIP existing_verified_sha256 path=%s", target)
                    continue
                quarantine = target.with_name(f"{target.name}.corrupt.{target_digest[:12]}")
                os.replace(target, quarantine)
                logging.warning(
                    "EXTRACT_QUARANTINE path=%s replacement_source_sha256=%s",
                    quarantine,
                    member_digest,
                )
            source = archive.extractfile(member)
            if source is None:
                raise RuntimeError(f"Could not open TAR member: {member.name}")
            temporary = target.with_suffix(target.suffix + ".tmp")
            with source, temporary.open("xb") as output:
                shutil.copyfileobj(source, output, length=8 * 1024 * 1024)
                output.flush()
                os.fsync(output.fileno())
            if temporary.stat().st_size != member.size:
                raise RuntimeError(f"Extracted size mismatch for {member.name}")
            os.replace(temporary, target)
            if sha256_path(target) != member_digest:
                raise RuntimeError(f"Extracted SHA-256 mismatch for {member.name}")
    logging.info("EXTRACT_COMPLETE directory=%s files=%s", EXTRACTED_DIR, len(members))


def count_gzip_lines(path: Path) -> int:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_tar_member(archive: tarfile.TarFile, member: tarfile.TarInfo) -> str:
    source = archive.extractfile(member)
    if source is None:
        raise RuntimeError(f"Could not open TAR member for hashing: {member.name}")
    digest = hashlib.sha256()
    with source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_sample(sample: str, gsm: str, condition: str, files: dict[str, tarfile.TarInfo]):
    barcode_path = EXTRACTED_DIR / files["barcodes.tsv.gz"].name
    feature_path = EXTRACTED_DIR / files["features.tsv.gz"].name
    matrix_path = EXTRACTED_DIR / files["matrix.mtx.gz"].name

    barcodes = pd.read_csv(barcode_path, sep="\t", header=None, compression="gzip", dtype=str)
    features = pd.read_csv(feature_path, sep="\t", header=None, compression="gzip", dtype=str)
    if features.shape[1] < 2:
        raise RuntimeError(f"Feature table for {sample} has fewer than two columns")
    matrix_feature_by_cell = mmread(matrix_path)
    matrix = sp.csr_matrix(matrix_feature_by_cell.T, dtype=np.float32)
    if matrix.shape != (len(barcodes), len(features)):
        raise RuntimeError(
            f"Shape mismatch for {sample}: matrix={matrix.shape}, barcodes={len(barcodes)}, features={len(features)}"
        )
    if np.any(matrix.data < 0) or not np.allclose(matrix.data, np.round(matrix.data)):
        raise RuntimeError(f"Matrix for {sample} is not non-negative integer-valued counts")

    full_sample_id = f"{gsm}_{sample}"
    obs = pd.DataFrame(index=pd.Index([f"{full_sample_id}:{value}" for value in barcodes.iloc[:, 0]], name="cell_id"))
    obs["sample_id"] = full_sample_id
    obs["donor_id"] = full_sample_id
    obs["geo_sample_id"] = gsm
    obs["sample_prefix"] = sample
    obs["condition"] = condition
    obs["dataset"] = "GSE254005"
    obs["species"] = "Canis lupus familiaris"
    obs["tissue"] = "proximal duodenum"

    var = pd.DataFrame(index=pd.Index(features.iloc[:, 1].astype(str), name="gene_name"))
    var["feature_id"] = features.iloc[:, 0].to_numpy()
    var["gene_symbol"] = features.iloc[:, 1].to_numpy()
    var["feature_type"] = features.iloc[:, 2].to_numpy() if features.shape[1] >= 3 else "Gene Expression"
    sample_adata = ad.AnnData(X=matrix, obs=obs, var=var)
    sample_adata.var_names_make_unique()

    summary = {
        "sample_prefix": sample,
        "sample_id": full_sample_id,
        "geo_sample_id": gsm,
        "condition": condition,
        "cells": sample_adata.n_obs,
        "genes_features": sample_adata.n_vars,
        "matrix_rows_cells": sample_adata.shape[0],
        "matrix_columns_features": sample_adata.shape[1],
        "nnz": sample_adata.X.nnz,
        "duplicate_gene_symbols": int(features.iloc[:, 1].duplicated().sum()),
    }
    manifest = {
        "sample_prefix": sample,
        "geo_sample_id": gsm,
        "condition": condition,
        "species": "Canis lupus familiaris",
        "tissue": "proximal duodenum",
        "matrix_file": files["matrix.mtx.gz"].name,
        "features_file": files["features.tsv.gz"].name,
        "barcodes_file": files["barcodes.tsv.gz"].name,
        "n_barcodes_raw": count_gzip_lines(barcode_path),
        "n_features_raw": count_gzip_lines(feature_path),
        "gene_reference": "CanFam3.1; Ensembl release 104",
        "notes": "donor_id=sample; paper dog ID, breed, age and medication remain NA/unmapped",
    }
    logging.info(
        "READ_SAMPLE sample=%s condition=%s shape=%s nnz=%s duplicate_symbols=%s",
        full_sample_id,
        condition,
        sample_adata.shape,
        sample_adata.X.nnz,
        summary["duplicate_gene_symbols"],
    )
    return sample_adata, summary, manifest


def main() -> int:
    ensure_phase1_directories()
    configure_logging()
    started = time.perf_counter()
    logging.info("START script_version=%s python=%s", SCRIPT_VERSION, platform.python_version())
    logging.info("INPUT archive=%s", ARCHIVE)
    try:
        grouped, members = audit_archive()
        extract_verified_members(members)
        if RAW_H5AD.exists():
            raise FileExistsError(f"Raw h5ad will not be overwritten: {RAW_H5AD}")

        sample_objects = []
        summaries = []
        manifests = []
        reference_feature_ids: np.ndarray | None = None
        reference_var_names: np.ndarray | None = None
        reference_var: pd.DataFrame | None = None
        for sample, (gsm, condition) in EXPECTED.items():
            sample_adata, summary, manifest = read_sample(sample, gsm, condition, grouped[sample])
            current_ids = sample_adata.var["feature_id"].to_numpy()
            current_names = sample_adata.var_names.to_numpy()
            if reference_feature_ids is None:
                reference_feature_ids = current_ids.copy()
                reference_var_names = current_names.copy()
                reference_var = sample_adata.var.copy()
            elif not np.array_equal(current_ids, reference_feature_ids) or not np.array_equal(
                current_names, reference_var_names
            ):
                raise RuntimeError(f"Feature order differs in {sample}; outer join was intentionally not attempted")
            sample_objects.append(sample_adata)
            summaries.append(summary)
            manifests.append(manifest)

        assert reference_var is not None
        combined_x = sp.vstack([obj.X for obj in sample_objects], format="csr")
        combined_obs = pd.concat([obj.obs for obj in sample_objects], axis=0)
        combined = ad.AnnData(X=combined_x, obs=combined_obs, var=reference_var)
        combined.uns["dataset"] = "GSE254005"
        combined.uns["matrix_state"] = "pre-QC processed filtered counts; not normalized"
        combined.uns["gene_reference"] = "CanFam3.1; Ensembl release 104"
        combined.uns["script_version_ingest"] = SCRIPT_VERSION
        combined.uns["condition_definition"] = "Healthy=h_*; CIE=cie_* from GEO sample names"
        combined.uns["metadata_warning"] = (
            "Do not map GEO cie_2/cie_3/cie_4/cie_6 to publication CIE_1-CIE_4 without author verification"
        )

        temporary_h5ad = RAW_H5AD.with_suffix(".h5ad.tmp")
        combined.write_h5ad(temporary_h5ad, compression="gzip")
        check = ad.read_h5ad(temporary_h5ad, backed="r")
        if check.shape != combined.shape:
            raise RuntimeError(f"Raw h5ad reload shape mismatch: {check.shape} vs {combined.shape}")
        check.file.close()
        os.replace(temporary_h5ad, RAW_H5AD)

        atomic_csv(
            manifests,
            SAMPLE_MANIFEST,
            [
                "sample_prefix",
                "geo_sample_id",
                "condition",
                "species",
                "tissue",
                "matrix_file",
                "features_file",
                "barcodes_file",
                "n_barcodes_raw",
                "n_features_raw",
                "gene_reference",
                "notes",
            ],
        )
        atomic_csv(
            summaries,
            INGEST_SUMMARY,
            [
                "sample_prefix",
                "sample_id",
                "geo_sample_id",
                "condition",
                "cells",
                "genes_features",
                "matrix_rows_cells",
                "matrix_columns_features",
                "nnz",
                "duplicate_gene_symbols",
            ],
        )
        logging.info("OUTPUT raw_h5ad=%s shape=%s", RAW_H5AD, combined.shape)
        logging.info("OUTPUT sample_manifest=%s ingest_summary=%s", SAMPLE_MANIFEST, INGEST_SUMMARY)
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
