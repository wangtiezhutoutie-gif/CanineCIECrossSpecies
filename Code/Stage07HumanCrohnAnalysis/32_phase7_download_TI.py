from __future__ import annotations

import csv
import hashlib
import os
import subprocess
import sys
import tarfile
import time
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
SAMPLE_MASTER = ROOT / "00_metadata" / "GSE266616_HUMAN_SAMPLE_MASTER.csv"
OUT_ROOT = ROOT / "07_human_crohn" / "GSE266616" / "raw_processed"
MANIFEST = ROOT / "00_metadata" / "GSE266616_DOWNLOAD_MANIFEST.csv"
LOG = ROOT / "logs" / "phase7_TI_download.log"

FORBIDDEN = (
    "GSE266616_RAW.tar",
    ".fastq",
    ".sra",
    "scRNA_gene_counts.mtx",
)
EXPECTED_HOST = "ftp.ncbi.nlm.nih.gov"
REQUIRED_SUFFIXES = ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def log(message: str) -> None:
    line = f"[{utc_now()}] {message}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(line + "\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def validate_url(url: str) -> str:
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https" or parsed.hostname != EXPECTED_HOST:
        raise RuntimeError(f"Refusing unexpected processed URL: {url}")
    if "/geo/samples/" not in parsed.path or "/suppl/" not in parsed.path:
        raise RuntimeError(f"Refusing non-per-GSM supplementary URL: {url}")
    filename = urllib.parse.unquote(PurePosixPath(parsed.path).name)
    lower = filename.lower()
    if not lower.endswith(".tar.gz"):
        raise RuntimeError(f"Refusing non-archive processed file: {filename}")
    for token in FORBIDDEN:
        if token.lower() in lower:
            raise RuntimeError(f"Forbidden download target encountered: {filename}")
    return filename


def subdir_for(row: dict[str, str]) -> str:
    if row["region"] != "TI":
        raise RuntimeError(f"Non-TI row reached downloader: {row['GSM']} {row['region']}")
    state = row["inflammation_state"]
    if row["healthy_or_crohn"] == "Healthy" and state == "Healthy":
        return "healthy_TI"
    mapping = {
        "Active_inflamed": "crohn_TI_active",
        "Adjacent": "crohn_TI_adjacent",
        "Non_inflamed": "crohn_TI_noninflamed",
    }
    if row["healthy_or_crohn"] == "Crohn" and state in mapping:
        return mapping[state]
    raise RuntimeError(f"Unexpected TI cohort state: {row['GSM']} {row['healthy_or_crohn']} {state}")


def validate_tar(path: Path) -> tuple[int, str]:
    required = {suffix: False for suffix in REQUIRED_SUFFIXES}
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        if not members:
            raise RuntimeError(f"Empty archive: {path}")
        for member in members:
            posix = PurePosixPath(member.name)
            if posix.is_absolute() or ".." in posix.parts:
                raise RuntimeError(f"Unsafe member path in {path.name}: {member.name}")
            if not member.isfile():
                continue
            lower = member.name.lower()
            for suffix in REQUIRED_SUFFIXES:
                if lower.endswith(suffix):
                    required[suffix] = True
    missing = [name for name, present in required.items() if not present]
    if missing:
        raise RuntimeError(f"Processed archive missing required 10x files: {path.name}: {missing}")
    return len(members), "YES"


def download_once(transfer_url: str, part: Path, expected_size: int) -> None:
    command = [
        "curl.exe",
        "--silent",
        "--show-error",
        "--location",
        "--fail",
        "--http1.1",
        "--retry",
        "10",
        "--retry-all-errors",
        "--retry-delay",
        "2",
        "--continue-at",
        "-",
        "--connect-timeout",
        "30",
        "--max-time",
        "7200",
        "--output",
        str(part),
        transfer_url,
    ]
    completed = subprocess.run(command, capture_output=True, text=True, timeout=7300, check=False)
    if completed.returncode != 0:
        detail = completed.stderr.strip().replace("\n", " | ")[-1000:]
        raise RuntimeError(f"curl exit {completed.returncode}: {detail}")
    downloaded = part.stat().st_size if part.exists() else 0
    if downloaded != expected_size:
        raise RuntimeError(f"Downloaded byte mismatch for {transfer_url}: {downloaded} != {expected_size}")


def process_row(row: dict[str, str]) -> dict[str, object]:
    url = row["processed_download_url"]
    filename = validate_url(url)
    transfer_url = url
    expected_size = int(row["processed_file_size"])
    directory = OUT_ROOT / subdir_for(row)
    directory.mkdir(parents=True, exist_ok=True)
    final = directory / filename
    part = directory / f".{filename}.part"

    status = "REUSED_VERIFIED"
    if not final.exists() or final.stat().st_size != expected_size:
        status = "DOWNLOADED"
        error: Exception | None = None
        if part.exists() and part.stat().st_size == expected_size:
            validate_tar(part)
            os.replace(part, final)
        else:
            if part.exists() and part.stat().st_size > expected_size:
                raise RuntimeError(f"Partial file exceeds expected size and requires manual audit: {part}")
            for attempt in range(1, 4):
                try:
                    download_once(transfer_url, part, expected_size)
                    validate_tar(part)
                    os.replace(part, final)
                    error = None
                    break
                except Exception as exc:  # noqa: BLE001 - retain exact download failure
                    error = exc
                    log(f"Retry {attempt}/3 for {row['GSM']} after {type(exc).__name__}: {exc}")
                    time.sleep(attempt * 2)
        if error is not None:
            raise error

    actual_size = final.stat().st_size
    if actual_size != expected_size:
        raise RuntimeError(f"Final byte mismatch for {final}: {actual_size} != {expected_size}")
    member_count, triplet = validate_tar(final)
    checksum = sha256_file(final)
    verified_time = utc_now()
    return {
        "GSM": row["GSM"],
        "sample_title": row["sample_title"],
        "patient_id": row["patient_id"],
        "state": row["inflammation_state"],
        "region": "terminal ileum",
        "url": url,
        "filename": final.name,
        "expected_size": expected_size,
        "actual_size": actual_size,
        "local_sha256": checksum,
        "download_time": verified_time,
        "archive_valid": triplet,
        "status": status,
        "healthy_or_crohn": row["healthy_or_crohn"],
        "inflammation_state": row["inflammation_state"],
        "primary_contrast_eligible": row["primary_contrast_eligible"],
        "secondary_paired_eligible": row["secondary_paired_eligible"],
        "processed_download_url": url,
        "transfer_url": transfer_url,
        "local_path": str(final),
        "expected_bytes": expected_size,
        "actual_bytes": actual_size,
        "repository_checksum": row.get("repository_checksum", "NA") or "NA",
        "sha256": checksum,
        "archive_member_count": member_count,
        "required_10x_triplet_present": triplet,
        "download_status": status,
        "verified_utc": verified_time,
    }


def write_manifest(rows: list[dict[str, object]]) -> None:
    fields = [
        "GSM",
        "patient_id",
        "state",
        "region",
        "url",
        "filename",
        "expected_size",
        "actual_size",
        "repository_checksum",
        "local_sha256",
        "download_time",
        "archive_valid",
        "status",
        "sample_title",
        "healthy_or_crohn",
        "inflammation_state",
        "primary_contrast_eligible",
        "secondary_paired_eligible",
        "processed_download_url",
        "transfer_url",
        "local_path",
        "expected_bytes",
        "actual_bytes",
        "sha256",
        "archive_member_count",
        "required_10x_triplet_present",
        "download_status",
        "verified_utc",
    ]
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("", encoding="utf-8")
    with SAMPLE_MASTER.open("r", encoding="utf-8-sig", newline="") as handle:
        all_rows = list(csv.DictReader(handle))
    selected = [row for row in all_rows if row["region"] == "TI"]
    if len(selected) != 48:
        raise RuntimeError(f"Expected 48 explicit TI samples, observed {len(selected)}")
    if any(row["tissue"] != "terminal_ileum" for row in selected):
        raise RuntimeError("TI selection contains non-terminal-ileum tissue")
    if sum(int(row["processed_file_size"]) for row in selected) != 2_302_117_661:
        raise RuntimeError("Frozen TI byte total changed from Phase7A")

    log("Starting 48 per-GSM TI processed archive downloads; AC/RAW/SRA/FASTQ/Zenodo matrix excluded")
    completed: dict[str, dict[str, object]] = {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(process_row, row): row for row in selected}
        for future in as_completed(futures):
            row = futures[future]
            result = future.result()
            completed[str(result["GSM"])] = result
            log(
                f"Verified {result['GSM']} {result['inflammation_state']} "
                f"{int(result['actual_bytes']) / (1024 * 1024):.1f} MiB ({len(completed)}/48)"
            )

    ordered = [completed[row["GSM"]] for row in selected]
    write_manifest(ordered)
    if len(ordered) != 48 or any(row["required_10x_triplet_present"] != "YES" for row in ordered):
        raise RuntimeError("TI download manifest completeness failure")
    log(f"Phase7B TI download gate PASS: 48/48 archives, {sum(int(row['actual_bytes']) for row in ordered)} bytes")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # noqa: BLE001 - top-level provenance
        log(f"FATAL: {type(exc).__name__}: {exc}")
        raise
