#!/usr/bin/env python3
"""Download only Phase7A lightweight official metadata.

Explicitly excluded: GSE266616_RAW.tar, SRA/FASTQ, Zenodo 6.1 GB count matrix,
and every per-GSM expression archive.  The latter are permitted only after the
Phase7A metadata gate passes.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import shutil
import ssl
import tarfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
BASE = ROOT / "07_human_crohn" / "GSE266616" / "metadata"
GEO_DIR = BASE / "source" / "geo"
PAPER_DIR = BASE / "source" / "paper"
ZENODO_DIR = BASE / "source" / "zenodo"
EXTRACT_DIR = BASE / "extracted"
LOG = ROOT / "logs" / "phase7_metadata.log"
MANIFEST = ROOT / "00_metadata" / "GSE266616_METADATA_DOWNLOAD_MANIFEST.csv"

GEO_SOFT = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266616/soft/GSE266616_family.soft.gz"
GEO_MINIML = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266616/miniml/GSE266616_family.xml.tgz"
EPMC_FULLTEXT = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11751449/fullTextXML"
SUPP_DATA12 = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41467-025-56212-2/MediaObjects/"
    "41467_2025_56212_MOESM14_ESM.xlsx"
)
ZENODO_API = "https://zenodo.org/api/records/14165540"
ZENODO_ALLOWED = {
    "scRNA_cell_metadata.csv",
    "scRNA_gene_counts_col_cellIDs.txt",
    "scRNA_gene_counts_row_genes.txt",
}
FORBIDDEN_SUBSTRINGS = ("GSE266616_RAW.tar", ".fastq", ".sra", "scRNA_gene_counts.mtx")


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def log(message: str) -> None:
    line = f"[{now_utc()}] {message}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(line + "\n")


def digest(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def assert_allowed(url: str, filename: str) -> None:
    joined = f"{url} {filename}".lower()
    for forbidden in FORBIDDEN_SUBSTRINGS:
        if forbidden.lower() in joined:
            raise RuntimeError(f"Forbidden Phase7A download requested: {url}")


def download(url: str, destination: Path) -> str:
    assert_allowed(url, destination.name)
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    if partial.exists():
        partial.unlink()
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Canine-CIE-Phase7A-metadata/1.0", "Accept": "*/*"},
    )
    with urllib.request.urlopen(request, context=ssl.create_default_context(), timeout=300) as response:
        with partial.open("wb") as output:
            shutil.copyfileobj(response, output)
    os.replace(partial, destination)
    timestamp = now_utc()
    log(f"Downloaded {url} -> {destination} ({destination.stat().st_size} bytes)")
    return timestamp


def ensure_download(url: str, destination: Path, expected_size: int | None = None) -> str:
    if destination.is_file() and destination.stat().st_size > 0:
        if expected_size is None or destination.stat().st_size == expected_size:
            log(f"Using existing verified-size candidate: {destination}")
            return "EXISTING_FILE_TIME_UNKNOWN"
    return download(url, destination)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    for directory in (GEO_DIR, PAPER_DIR, ZENODO_DIR, EXTRACT_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("", encoding="utf-8")
    log("Phase7A lightweight metadata download started")
    rows: list[dict[str, object]] = []

    static_jobs = [
        ("GEO", GEO_SOFT, GEO_DIR / "GSE266616_family.soft.gz", None, "NA_NOT_PROVIDED"),
        ("GEO", GEO_MINIML, GEO_DIR / "GSE266616_family.xml.tgz", None, "NA_NOT_PROVIDED"),
        ("Europe_PMC", EPMC_FULLTEXT, PAPER_DIR / "PMC11751449_fulltext.xml", None, "NA_NOT_PROVIDED"),
        (
            "Springer_Nature",
            SUPP_DATA12,
            PAPER_DIR / "41467_2025_56212_MOESM14_ESM.xlsx",
            9476,
            "etag_md5:37f1a6421ac1a5bf16b1690448bccba0",
        ),
        ("Zenodo", ZENODO_API, ZENODO_DIR / "zenodo_record_14165540.json", None, "NA_NOT_PROVIDED"),
    ]
    for source, url, destination, expected_size, checksum in static_jobs:
        timestamp = ensure_download(url, destination, expected_size)
        if expected_size is not None and destination.stat().st_size != expected_size:
            raise RuntimeError(
                f"Size mismatch for {destination}: expected {expected_size}, observed {destination.stat().st_size}"
            )
        if checksum.startswith("etag_md5:"):
            expected_md5 = checksum.split(":", 1)[1]
            if digest(destination, "md5") != expected_md5:
                raise RuntimeError(f"Nature Supplementary Data 12 MD5 mismatch: {destination}")
        rows.append(
            {
                "source": source,
                "url": url,
                "filename": destination.name,
                "expected_size": expected_size if expected_size is not None else "NA_NOT_PROVIDED",
                "actual_size": destination.stat().st_size,
                "repository_checksum": checksum,
                "local_sha256": digest(destination),
                "download_time_utc": timestamp,
                "status": "VERIFIED",
                "local_path": str(destination),
                "notes": "Phase7A lightweight metadata",
            }
        )

    zenodo = json.loads((ZENODO_DIR / "zenodo_record_14165540.json").read_text(encoding="utf-8"))
    available_names = {str(item["key"]) for item in zenodo.get("files", [])}
    missing = sorted(ZENODO_ALLOWED - available_names)
    if missing:
        raise RuntimeError(f"Required Zenodo lightweight metadata missing: {missing}")
    for item in zenodo["files"]:
        name = str(item["key"])
        if name not in ZENODO_ALLOWED:
            continue
        url = str(item.get("links", {}).get("content") or item.get("links", {}).get("self"))
        destination = ZENODO_DIR / name
        expected_size = int(item["size"])
        timestamp = ensure_download(url, destination, expected_size)
        checksum = str(item.get("checksum") or "NA_NOT_PROVIDED")
        if checksum.startswith("md5:"):
            expected_md5 = checksum.split(":", 1)[1]
            if digest(destination, "md5").lower() != expected_md5.lower():
                raise RuntimeError(f"Zenodo MD5 mismatch for {destination}")
        rows.append(
            {
                "source": "Zenodo",
                "url": url,
                "filename": name,
                "expected_size": expected_size,
                "actual_size": destination.stat().st_size,
                "repository_checksum": checksum,
                "local_sha256": digest(destination),
                "download_time_utc": timestamp,
                "status": "VERIFIED",
                "local_path": str(destination),
                "notes": "Allowed lightweight metadata; 6.1 GB matrix explicitly excluded",
            }
        )

    soft_gz = GEO_DIR / "GSE266616_family.soft.gz"
    soft_txt = EXTRACT_DIR / "GSE266616_family.soft"
    with gzip.open(soft_gz, "rb") as source, soft_txt.open("wb") as output:
        shutil.copyfileobj(source, output)
    rows.append(
        {
            "source": "derived_from_GEO",
            "url": GEO_SOFT,
            "filename": soft_txt.name,
            "expected_size": "NA",
            "actual_size": soft_txt.stat().st_size,
            "repository_checksum": "NA_DERIVED",
            "local_sha256": digest(soft_txt),
            "download_time_utc": now_utc(),
            "status": "EXTRACTED_VERIFIED",
            "local_path": str(soft_txt),
            "notes": "Decompressed SOFT; source gzip retained unchanged",
        }
    )

    miniml_tgz = GEO_DIR / "GSE266616_family.xml.tgz"
    miniml_dir = EXTRACT_DIR / "miniml"
    miniml_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(miniml_tgz, "r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            target = (miniml_dir / member.name).resolve()
            if miniml_dir.resolve() not in target.parents:
                raise RuntimeError(f"Unsafe MINiML archive member: {member.name}")
            source = archive.extractfile(member)
            if source is None:
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with source, target.open("wb") as output:
                shutil.copyfileobj(source, output)
            rows.append(
                {
                    "source": "derived_from_GEO",
                    "url": GEO_MINIML,
                    "filename": str(target.relative_to(miniml_dir)).replace("\\", "/"),
                    "expected_size": member.size,
                    "actual_size": target.stat().st_size,
                    "repository_checksum": "NA_DERIVED",
                    "local_sha256": digest(target),
                    "download_time_utc": now_utc(),
                    "status": "EXTRACTED_VERIFIED",
                    "local_path": str(target),
                    "notes": "Extracted MINiML member; source tgz retained unchanged",
                }
            )

    forbidden_found = []
    for path in BASE.rglob("*"):
        if not path.is_file():
            continue
        lower = path.name.lower()
        if any(term.lower() in lower for term in FORBIDDEN_SUBSTRINGS):
            forbidden_found.append(str(path))
    if forbidden_found:
        raise RuntimeError(f"Forbidden payload present in Phase7A metadata tree: {forbidden_found}")

    fields = [
        "source",
        "url",
        "filename",
        "expected_size",
        "actual_size",
        "repository_checksum",
        "local_sha256",
        "download_time_utc",
        "status",
        "local_path",
        "notes",
    ]
    write_csv(MANIFEST, fields, rows)
    log(f"Phase7A metadata download complete: {len(rows)} manifest rows")
    log("Confirmed exclusions: no GEO RAW tar, no per-GSM matrices, no SRA/FASTQ, no Zenodo count matrix")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        log(f"FATAL: {type(exc).__name__}: {exc}")
        raise
