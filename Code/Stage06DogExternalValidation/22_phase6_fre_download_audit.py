#!/usr/bin/env python3
"""Phase6 FRE data-source audit and minimal processed-data download.

This script intentionally downloads only the published Zenodo artifact and the
official Europe PMC supplementary archive.  It never downloads FASTQ files.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import ssl
import sys
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
RAW_ROOT = ROOT / "01_raw_processed" / "dog_FRE_2026"
ZENODO_DIR = RAW_ROOT / "downloaded" / "zenodo"
EPMC_DIR = RAW_ROOT / "downloaded" / "europepmc"
EXTRACT_DIR = RAW_ROOT / "extracted" / "PMC13293934_SupplementaryFiles"
REFERENCE_DIR = RAW_ROOT / "reference"
METADATA_DIR = ROOT / "00_metadata"
LOG_PATH = ROOT / "logs" / "phase6_fre_download_audit.log"

ZENODO_RECORD_API = "https://zenodo.org/api/records/20330471"
FRONTIERS_XML = (
    "https://www.frontiersin.org/journals/veterinary-science/articles/"
    "10.3389/fvets.2026.1829399/xml/nlm"
)
EPMC_SEARCH_API = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
    "query=EXT_ID:42368340&format=json"
)
EPMC_FULLTEXT_XML = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13293934/fullTextXML"
)
EPMC_SUPPLEMENT_ZIP = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13293934/"
    "supplementaryFiles"
)
NCBI_OA_API = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC13293934"
DOI = "10.5281/zenodo.20330471"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def log(message: str) -> None:
    timestamped = f"[{utc_now()}] {message}"
    print(timestamped, flush=True)
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LOG_PATH.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(timestamped + "\n")


def digest(path: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    if temporary.exists():
        temporary.unlink()
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Canine-CIE-Phase6-audit/1.0",
            "Accept": "*/*",
        },
    )
    context = ssl.create_default_context()
    with urllib.request.urlopen(request, context=context, timeout=180) as response:
        with temporary.open("wb") as output:
            shutil.copyfileobj(response, output)
    os.replace(temporary, destination)
    log(f"Downloaded {url} -> {destination} ({destination.stat().st_size} bytes)")


def fetch_reference(url: str, destination: Path) -> None:
    if destination.exists() and destination.stat().st_size > 0:
        log(f"Reference already present: {destination}")
        return
    download(url, destination)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    for directory in (ZENODO_DIR, EPMC_DIR, EXTRACT_DIR, REFERENCE_DIR, METADATA_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    LOG_PATH.write_text("", encoding="utf-8")
    log("Phase6 FRE data audit started; FASTQ download is prohibited by design.")

    reference_jobs = [
        (ZENODO_RECORD_API, REFERENCE_DIR / "zenodo_record_20330471.json"),
        (FRONTIERS_XML, REFERENCE_DIR / "frontiers_2026_1829399.xml"),
        (EPMC_SEARCH_API, REFERENCE_DIR / "europepmc_search_PMID42368340.json"),
        (EPMC_FULLTEXT_XML, REFERENCE_DIR / "europepmc_PMC13293934_fulltext.xml"),
        (NCBI_OA_API, REFERENCE_DIR / "ncbi_oa_PMC13293934.xml"),
    ]
    for url, destination in reference_jobs:
        fetch_reference(url, destination)

    zenodo_metadata = json.loads(
        (REFERENCE_DIR / "zenodo_record_20330471.json").read_text(encoding="utf-8")
    )
    files = zenodo_metadata.get("files", [])
    if len(files) != 1:
        raise RuntimeError(f"Expected one current Zenodo file; observed {len(files)}")
    zfile = files[0]
    zname = str(zfile["key"])
    zurl = str(zfile.get("links", {}).get("content") or zfile.get("links", {}).get("self"))
    zpath = ZENODO_DIR / zname
    if not zpath.exists() or zpath.stat().st_size != int(zfile["size"]):
        download(zurl, zpath)
    official_checksum = str(zfile.get("checksum", ""))
    if official_checksum.startswith("md5:"):
        observed_md5 = digest(zpath, "md5")
        expected_md5 = official_checksum.split(":", 1)[1]
        if observed_md5.lower() != expected_md5.lower():
            raise RuntimeError(
                f"Zenodo MD5 mismatch for {zpath}: expected {expected_md5}, got {observed_md5}"
            )
        log(f"Verified Zenodo MD5 for {zname}: {observed_md5}")

    zip_path = EPMC_DIR / "PMC13293934_SupplementaryFiles.zip"
    if not zip_path.exists() or zip_path.stat().st_size == 0:
        download(EPMC_SUPPLEMENT_ZIP, zip_path)
    if not zipfile.is_zipfile(zip_path):
        raise RuntimeError(f"Official Europe PMC response is not a ZIP archive: {zip_path}")

    with zipfile.ZipFile(zip_path) as archive:
        members = [member for member in archive.infolist() if not member.is_dir()]
        for member in members:
            target = (EXTRACT_DIR / member.filename).resolve()
            if EXTRACT_DIR.resolve() not in target.parents:
                raise RuntimeError(f"Unsafe ZIP member path: {member.filename}")
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists() or target.stat().st_size != member.file_size:
                with archive.open(member) as source, target.open("wb") as output:
                    shutil.copyfileobj(source, output)
        log(f"Extracted/verified {len(members)} supplementary archive members")

    audit_rows: list[dict[str, object]] = []
    download_rows: list[dict[str, object]] = []
    timestamp = utc_now()

    audit_rows.append(
        {
            "source": "Zenodo record 20330471",
            "source_url": zurl,
            "file_name": zname,
            "file_type": zpath.suffix.lower().lstrip("."),
            "file_size": zpath.stat().st_size,
            "checksum": official_checksum,
            "local_sha256": digest(zpath, "sha256"),
            "likely_content": "TPM values for selected bile acid-associated genes",
            "is_raw_count": "NO",
            "is_tpm": "YES_TARGETED_ONLY",
            "is_de_table": "NO",
            "is_sample_metadata": "NO",
            "is_fastq": "NO",
            "is_pathway_results": "NO",
            "required_for_phase6": "YES_CONTEXT_ONLY",
            "download_status": "VERIFIED",
            "local_path": str(zpath),
            "notes": "Only file exposed by the live Zenodo record; not a whole-transcriptome matrix.",
        }
    )
    download_rows.append(
        {
            "url": zurl,
            "doi": DOI,
            "filename": zname,
            "bytes": zpath.stat().st_size,
            "repository_checksum": official_checksum,
            "local_sha256": digest(zpath, "sha256"),
            "download_utc_time": timestamp,
            "status": "VERIFIED",
            "local_path": str(zpath),
        }
    )

    audit_rows.append(
        {
            "source": "Europe PMC supplementary archive",
            "source_url": EPMC_SUPPLEMENT_ZIP,
            "file_name": zip_path.name,
            "file_type": "zip",
            "file_size": zip_path.stat().st_size,
            "checksum": "NA_NOT_PROVIDED",
            "local_sha256": digest(zip_path, "sha256"),
            "likely_content": "Official archive of ten supplementary PDF/DOCX files",
            "is_raw_count": "NO",
            "is_tpm": "YES_TARGETED_ONLY",
            "is_de_table": "NO_FULL_DE_TABLE_IDENTIFIED",
            "is_sample_metadata": "YES_SUPPLEMENTARY_TABLE_1",
            "is_fastq": "NO",
            "is_pathway_results": "YES_SUPPLEMENTARY_TABLE_4",
            "required_for_phase6": "YES_DATA_GATE_AUDIT",
            "download_status": "VERIFIED",
            "local_path": str(zip_path),
            "notes": "Archive contents are itemized below; content classification is verified during ingest.",
        }
    )
    download_rows.append(
        {
            "url": EPMC_SUPPLEMENT_ZIP,
            "doi": "10.3389/fvets.2026.1829399",
            "filename": zip_path.name,
            "bytes": zip_path.stat().st_size,
            "repository_checksum": "NA_NOT_PROVIDED",
            "local_sha256": digest(zip_path, "sha256"),
            "download_utc_time": timestamp,
            "status": "VERIFIED",
            "local_path": str(zip_path),
        }
    )

    for member_path in sorted(EXTRACT_DIR.rglob("*")):
        if not member_path.is_file():
            continue
        audit_rows.append(
            {
                "source": "Europe PMC supplementary archive member",
                "source_url": EPMC_SUPPLEMENT_ZIP,
                "file_name": str(member_path.relative_to(EXTRACT_DIR)).replace("\\", "/"),
                "file_type": member_path.suffix.lower().lstrip("."),
                "file_size": member_path.stat().st_size,
                "checksum": "NA_NOT_PROVIDED",
                "local_sha256": digest(member_path, "sha256"),
                "likely_content": "PENDING_DOCUMENT_INGEST",
                "is_raw_count": "PENDING_DOCUMENT_INGEST",
                "is_tpm": "PENDING_DOCUMENT_INGEST",
                "is_de_table": "PENDING_DOCUMENT_INGEST",
                "is_sample_metadata": "PENDING_DOCUMENT_INGEST",
                "is_fastq": "NO",
                "is_pathway_results": "PENDING_DOCUMENT_INGEST",
                "required_for_phase6": "PENDING_DOCUMENT_INGEST",
                "download_status": "EXTRACTED_VERIFIED",
                "local_path": str(member_path),
                "notes": "Original extracted member retained without modification.",
            }
        )

    audit_fields = [
        "source",
        "source_url",
        "file_name",
        "file_type",
        "file_size",
        "checksum",
        "local_sha256",
        "likely_content",
        "is_raw_count",
        "is_tpm",
        "is_de_table",
        "is_sample_metadata",
        "is_fastq",
        "is_pathway_results",
        "required_for_phase6",
        "download_status",
        "local_path",
        "notes",
    ]
    download_fields = [
        "url",
        "doi",
        "filename",
        "bytes",
        "repository_checksum",
        "local_sha256",
        "download_utc_time",
        "status",
        "local_path",
    ]
    write_csv(METADATA_DIR / "PHASE6_FRE_DATA_AUDIT.csv", audit_fields, audit_rows)
    write_csv(METADATA_DIR / "FRE2026_DOWNLOAD_MANIFEST.csv", download_fields, download_rows)
    log(f"Wrote initial audit with {len(audit_rows)} rows")
    log("Phase6 FRE data audit completed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        log(f"FATAL: {type(exc).__name__}: {exc}")
        raise
