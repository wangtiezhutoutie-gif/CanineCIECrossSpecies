"""Safely download and freeze the GSE254005 processed count-matrix archive."""

from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import ftplib
import gzip
import hashlib
import logging
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.error
import urllib.request
import zlib
from datetime import datetime, timezone
from pathlib import Path

import psutil

CODE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_ROOT))
PROJECT_ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"] )

from config.paths import DOG_RAW_DIR, LOG_DIR, METADATA_DIR, ensure_phase1_directories

SCRIPT_VERSION = "1.6.1"
DATASET = "GSE254005"
ACCESSION = "GSE254005"
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254005/suppl/GSE254005_RAW.tar"
FTP_HOST = "ftp.ncbi.nlm.nih.gov"
FTP_PATH = "/geo/series/GSE254nnn/GSE254005/suppl/GSE254005_RAW.tar"
FTP_URL = f"ftp://{FTP_HOST}{FTP_PATH}"
FILENAME = "GSE254005_RAW.tar"
TARGET = DOG_RAW_DIR / FILENAME
PART = DOG_RAW_DIR / f"{FILENAME}.part"
MANIFEST = METADATA_DIR / "GSE254005_DOWNLOAD_MANIFEST.csv"
LOG_FILE = LOG_DIR / "phase1_download.log"
SEGMENT_DIR = DOG_RAW_DIR / "GSE254005_RAW.tar.https_segments"
USER_AGENT = "Canine-CIE-Phase1/1.0 (scientific dataset audit)"
CHUNK_BYTES = 8 * 1024 * 1024
SEGMENT_BYTES = 2 * 1024 * 1024
SEGMENT_WORKERS = 8

MANIFEST_FIELDS = [
    "dataset",
    "accession",
    "url",
    "filename",
    "local_path",
    "download_time_utc",
    "http_status",
    "content_length",
    "actual_bytes",
    "sha256",
    "archive_readable",
    "status",
    "notes",
    "script_version",
    "python_version",
]


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)sZ %(levelname)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
        handlers=[logging.FileHandler(LOG_FILE, encoding="utf-8"), logging.StreamHandler()],
        force=True,
    )
    logging.Formatter.converter = time.gmtime


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()


def archive_readable(path: Path) -> tuple[bool, str]:
    try:
        with tarfile.open(path, mode="r:") as archive:
            members = archive.getmembers()
            if not members:
                return False, "TAR contains no members"
            for member in members:
                if not member.isfile() or member.size <= 0:
                    return False, f"Invalid TAR member: {member.name} size={member.size}"
                if not member.name.endswith(".gz"):
                    return False, f"Unexpected non-gzip TAR member: {member.name}"
                source = archive.extractfile(member)
                if source is None:
                    return False, f"Could not open TAR member: {member.name}"
                with source, gzip.GzipFile(fileobj=source, mode="rb") as stream:
                    while stream.read(CHUNK_BYTES):
                        pass
        return True, f"TAR readable; {len(members)} members; all gzip streams passed full decompression/CRC"
    except (tarfile.TarError, OSError, EOFError, gzip.BadGzipFile, zlib.error) as exc:
        return False, f"Archive/member validation failed: {type(exc).__name__}: {exc}"


def get_remote_metadata() -> tuple[int | str, int | str]:
    request = urllib.request.Request(URL, method="HEAD", headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        status = getattr(response, "status", response.getcode())
        length_header = response.headers.get("Content-Length")
        return status, int(length_header) if length_header else "NA"


def download_to_part(row: dict[str, object], max_attempts: int = 3) -> int:
    """Download with bounded retries and safe Range resume into the .part file."""
    last_error: Exception | None = None
    for attempt in range(1, max_attempts + 1):
        offset = PART.stat().st_size if PART.exists() else 0
        headers = {"User-Agent": USER_AGENT, "Connection": "close"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        request = urllib.request.Request(URL, method="GET", headers=headers)
        try:
            logging.info("DOWNLOAD_ATTEMPT attempt=%s offset=%s", attempt, offset)
            with urllib.request.urlopen(request, timeout=120) as response:
                status = getattr(response, "status", response.getcode())
                row["http_status"] = status
                if offset and status != 206:
                    stale = PART.with_name(
                        f"{PART.name}.range_ignored.{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
                    )
                    os.replace(PART, stale)
                    raise RuntimeError(
                        f"Server ignored Range request (status={status}); preserved partial as {stale.name}"
                    )
                if not offset and status != 200:
                    raise RuntimeError(f"Unexpected GET status: {status}")

                mode = "ab" if PART.exists() else "xb"
                actual = offset
                with PART.open(mode) as output:
                    while True:
                        block = response.read(CHUNK_BYTES)
                        if not block:
                            break
                        output.write(block)
                        actual += len(block)
                        if actual and actual % (64 * 1024 * 1024) < CHUNK_BYTES:
                            logging.info("DOWNLOAD_PROGRESS bytes=%s", actual)
                    output.flush()
                    os.fsync(output.fileno())

            expected = row["content_length"]
            actual = PART.stat().st_size
            if expected != "NA" and actual < int(expected):
                raise ConnectionError(f"Short transfer: expected {expected}, currently {actual}")
            if expected != "NA" and actual > int(expected):
                raise RuntimeError(f"Oversized transfer: expected {expected}, currently {actual}")
            return actual
        except (OSError, urllib.error.URLError, ConnectionError, RuntimeError) as exc:
            last_error = exc
            current = PART.stat().st_size if PART.exists() else 0
            logging.warning(
                "DOWNLOAD_RETRY attempt=%s current_bytes=%s error=%s: %s",
                attempt,
                current,
                type(exc).__name__,
                exc,
            )
            if attempt < max_attempts:
                time.sleep(min(2**attempt, 30))
    raise RuntimeError(f"Download failed after {max_attempts} attempts: {last_error}")


def download_to_part_ftp(row: dict[str, object], max_attempts: int = 4) -> int:
    """Fallback to the same NCBI host/path over FTP, with REST resume."""
    last_error: Exception | None = None
    for attempt in range(1, max_attempts + 1):
        offset = PART.stat().st_size if PART.exists() else 0
        try:
            logging.info("FTP_DOWNLOAD_ATTEMPT attempt=%s offset=%s", attempt, offset)
            with ftplib.FTP(FTP_HOST, timeout=120) as ftp:
                ftp.login()
                ftp.voidcmd("TYPE I")
                remote_size = ftp.size(FTP_PATH)
                if remote_size is not None:
                    if row["content_length"] != "NA" and int(row["content_length"]) != int(remote_size):
                        raise RuntimeError(
                            f"HTTPS HEAD size {row['content_length']} differs from FTP size {remote_size}"
                        )
                    row["content_length"] = int(remote_size)

                mode = "ab" if PART.exists() else "xb"
                actual = offset
                with PART.open(mode) as output:
                    def write_block(block: bytes) -> None:
                        nonlocal actual
                        output.write(block)
                        actual += len(block)
                        if actual and actual % (64 * 1024 * 1024) < 64 * 1024:
                            logging.info("FTP_DOWNLOAD_PROGRESS bytes=%s", actual)

                    ftp.retrbinary(f"RETR {FTP_PATH}", write_block, blocksize=64 * 1024, rest=offset or None)
                    output.flush()
                    os.fsync(output.fileno())

            expected = row["content_length"]
            actual = PART.stat().st_size
            if expected != "NA" and actual < int(expected):
                raise ConnectionError(f"Short FTP transfer: expected {expected}, currently {actual}")
            if expected != "NA" and actual > int(expected):
                raise RuntimeError(f"Oversized FTP transfer: expected {expected}, currently {actual}")
            return actual
        except ftplib.all_errors + (RuntimeError,) as exc:
            last_error = exc
            current = PART.stat().st_size if PART.exists() else 0
            logging.warning(
                "FTP_DOWNLOAD_RETRY attempt=%s current_bytes=%s error=%s: %s",
                attempt,
                current,
                type(exc).__name__,
                exc,
            )
            if attempt < max_attempts:
                time.sleep(min(2**attempt, 30))
    raise RuntimeError(f"FTP download failed after {max_attempts} attempts: {last_error}")


def download_to_part_curl_ftp(row: dict[str, object], max_attempts: int = 12) -> int:
    """Resume the same NCBI FTP object using curl's mature Windows transfer backend."""
    expected = row["content_length"]
    last_error = "NA"
    for attempt in range(1, max_attempts + 1):
        offset = PART.stat().st_size if PART.exists() else 0
        resume_arguments = ["--continue-at", "-"] if offset > 0 else []
        command = [
            "curl.exe",
            "--fail",
            "--silent",
            "--show-error",
            "--location",
            "--connect-timeout",
            "30",
            "--max-time",
            "300",
            "--speed-limit",
            "1024",
            "--speed-time",
            "60",
            *resume_arguments,
            "--output",
            str(PART),
            FTP_URL,
        ]
        logging.info("CURL_FTP_ATTEMPT attempt=%s offset=%s url=%s", attempt, offset, FTP_URL)
        result = subprocess.run(command, text=True, capture_output=True, check=False)
        actual = PART.stat().st_size if PART.exists() else 0
        if result.stderr.strip():
            logging.warning("CURL_STDERR %s", result.stderr.strip().replace("\n", " | "))
        logging.info(
            "CURL_FTP_RESULT attempt=%s returncode=%s bytes_before=%s bytes_after=%s",
            attempt,
            result.returncode,
            offset,
            actual,
        )
        if expected != "NA" and actual > int(expected):
            raise RuntimeError(f"curl FTP oversized partial: expected {expected}, actual {actual}")
        if expected != "NA" and actual == int(expected):
            logging.info("CURL_FTP_END bytes=%s", actual)
            return actual
        last_error = f"returncode={result.returncode}; bytes={actual}; stderr={result.stderr.strip()}"
        if attempt < max_attempts:
            time.sleep(min(2**min(attempt, 4), 15))
    raise RuntimeError(f"curl FTP transfer incomplete after {max_attempts} attempts: {last_error}")


def download_https_segment(index: int, start: int, end: int, total: int, max_attempts: int = 8) -> Path:
    """Download one exact HTTPS byte range; TLS protects each accepted segment in transit."""
    SEGMENT_DIR.mkdir(parents=True, exist_ok=True)
    expected_bytes = end - start + 1
    segment = SEGMENT_DIR / f"segment_{index:04d}_{start}_{end}.bin"
    temporary = segment.with_suffix(".bin.part")
    if segment.exists():
        if segment.stat().st_size == expected_bytes:
            return segment
        raise RuntimeError(f"Existing segment has wrong size: {segment}")
    if temporary.exists():
        temporary.unlink()

    last_error: Exception | None = None
    for attempt in range(1, max_attempts + 1):
        request = urllib.request.Request(
            URL,
            method="GET",
            headers={
                "User-Agent": USER_AGENT,
                "Connection": "close",
                "Range": f"bytes={start}-{end}",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=240) as response:
                status = getattr(response, "status", response.getcode())
                content_range = response.headers.get("Content-Range", "")
                expected_range = f"bytes {start}-{end}/{total}"
                if status != 206 or content_range != expected_range:
                    raise RuntimeError(
                        f"Range response mismatch: status={status}, Content-Range={content_range}, expected={expected_range}"
                    )
                written = 0
                with temporary.open("xb") as output:
                    while True:
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        output.write(block)
                        written += len(block)
                    output.flush()
                    os.fsync(output.fileno())
            if written != expected_bytes:
                raise ConnectionError(f"Segment {index} short: expected {expected_bytes}, wrote {written}")
            os.replace(temporary, segment)
            logging.info(
                "HTTPS_SEGMENT_OK index=%s range=%s-%s bytes=%s attempt=%s",
                index,
                start,
                end,
                written,
                attempt,
            )
            return segment
        except (OSError, urllib.error.URLError, ConnectionError, RuntimeError) as exc:
            last_error = exc
            if temporary.exists():
                temporary.unlink()
            logging.warning(
                "HTTPS_SEGMENT_RETRY index=%s attempt=%s error=%s: %s",
                index,
                attempt,
                type(exc).__name__,
                exc,
            )
            if attempt < max_attempts:
                time.sleep(min(2**min(attempt, 4), 15))
    raise RuntimeError(f"HTTPS segment {index} failed after {max_attempts} attempts: {last_error}")


def download_to_part_https_segments(total: int) -> int:
    """Download fixed TLS-protected ranges concurrently, then assemble in byte order."""
    if PART.exists():
        raise RuntimeError(f"Refusing to assemble over an existing partial file: {PART}")
    ranges = []
    for index, start in enumerate(range(0, total, SEGMENT_BYTES)):
        end = min(start + SEGMENT_BYTES - 1, total - 1)
        ranges.append((index, start, end))
    logging.info(
        "HTTPS_SEGMENT_PLAN total_bytes=%s segments=%s segment_bytes=%s workers=%s",
        total,
        len(ranges),
        SEGMENT_BYTES,
        SEGMENT_WORKERS,
    )
    completed: dict[int, Path] = {}
    with ThreadPoolExecutor(max_workers=SEGMENT_WORKERS, thread_name_prefix="geo-range") as executor:
        futures = {
            executor.submit(download_https_segment, index, start, end, total): index
            for index, start, end in ranges
        }
        for future in as_completed(futures):
            index = futures[future]
            completed[index] = future.result()
            logging.info("HTTPS_SEGMENT_PROGRESS completed=%s total=%s", len(completed), len(ranges))
    if len(completed) != len(ranges):
        raise RuntimeError(f"Segment completion mismatch: {len(completed)} vs {len(ranges)}")

    with PART.open("xb") as output:
        for index, start, end in ranges:
            segment = completed[index]
            if segment.stat().st_size != end - start + 1:
                raise RuntimeError(f"Segment changed before assembly: {segment}")
            with segment.open("rb") as source:
                shutil.copyfileobj(source, output, length=8 * 1024 * 1024)
        output.flush()
        os.fsync(output.fileno())
    actual = PART.stat().st_size
    if actual != total:
        raise RuntimeError(f"Assembled HTTPS file size mismatch: expected {total}, actual {actual}")
    logging.info("HTTPS_SEGMENT_ASSEMBLY_COMPLETE bytes=%s", actual)
    return actual


def write_manifest(row: dict[str, object]) -> None:
    temporary = MANIFEST.with_suffix(".csv.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerow(row)
    os.replace(temporary, MANIFEST)


def base_row() -> dict[str, object]:
    return {
        "dataset": DATASET,
        "accession": ACCESSION,
        "url": URL,
        "filename": FILENAME,
        "local_path": str(TARGET),
        "download_time_utc": datetime.now(timezone.utc).isoformat(),
        "http_status": "NA",
        "content_length": "NA",
        "actual_bytes": "NA",
        "sha256": "NA",
        "archive_readable": False,
        "status": "ERROR",
        "notes": "",
        "script_version": SCRIPT_VERSION,
        "python_version": platform.python_version(),
    }


def main() -> int:
    ensure_phase1_directories()
    configure_logging()
    started = time.perf_counter()
    process = psutil.Process(os.getpid())

    def peak_ram_mb() -> float:
        memory = process.memory_info()
        return getattr(memory, "peak_wset", memory.rss) / 1024**2

    row = base_row()
    logging.info("START script_version=%s python=%s", SCRIPT_VERSION, sys.version.replace("\n", " "))
    logging.info("INPUT url=%s target=%s", URL, TARGET)
    try:
        remote_status, remote_length = get_remote_metadata()
        row["http_status"] = remote_status
        row["content_length"] = remote_length
        logging.info("REMOTE http_status=%s content_length=%s", remote_status, remote_length)

        quarantine_notes: list[str] = []
        if TARGET.exists():
            actual = TARGET.stat().st_size
            readable, note = archive_readable(TARGET)
            digest = sha256_file(TARGET)
            size_ok = remote_length == "NA" or actual == remote_length
            row.update(
                {
                    "download_time_utc": datetime.fromtimestamp(
                        TARGET.stat().st_mtime, timezone.utc
                    ).isoformat(),
                    "actual_bytes": actual,
                    "sha256": digest,
                    "archive_readable": readable,
                    "status": "EXISTING_VALID" if readable and size_ok else "EXISTING_INVALID",
                    "notes": f"Existing target was not overwritten. {note}; remote size match={size_ok}",
                }
            )
            write_manifest(row)
            logging.info("EXISTING actual_bytes=%s sha256=%s readable=%s size_ok=%s", actual, digest, readable, size_ok)
            if not readable or not size_ok:
                quarantine = TARGET.with_name(
                    f"{TARGET.name}.corrupt.{digest[:12]}.{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
                )
                os.replace(TARGET, quarantine)
                quarantine_notes.append(f"Prior invalid target preserved as {quarantine.name}: {note}")
                logging.warning("QUARANTINE path=%s reason=%s", quarantine, note)
            else:
                logging.info(
                    "END status=success_existing_valid elapsed_seconds=%.3f peak_ram_mb=%.1f",
                    time.perf_counter() - started,
                    peak_ram_mb(),
                )
                return 0

        backend_note = "parallel fixed-range HTTPS GET with TLS-protected 2-MiB segments"
        if remote_length == "NA":
            raise RuntimeError("Segmented HTTPS download requires a verified remote Content-Length")
        if PART.exists():
            part_digest = sha256_file(PART)
            quarantine = PART.with_name(
                f"{PART.name}.untrusted_ftp.{part_digest[:12]}.{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
            )
            os.replace(PART, quarantine)
            quarantine_notes.append(
                f"Untrusted FTP partial preserved as {quarantine.name}; fixed-range FTP was not byte-stable"
            )
            logging.warning("QUARANTINE path=%s reason=untrusted FTP partial", quarantine)

        actual = download_to_part_https_segments(int(remote_length))
        readable, note = archive_readable(PART)
        if not readable:
            raise RuntimeError(note)

        digest = sha256_file(PART)
        os.replace(PART, TARGET)
        row.update(
            {
                "download_time_utc": datetime.now(timezone.utc).isoformat(),
                "actual_bytes": actual,
                "sha256": digest,
                "archive_readable": True,
                "status": "DOWNLOADED_VALID",
                "notes": (
                    f"{note}; atomically renamed from .part; transfer backend={backend_note}; "
                    f"no repository checksum was asserted; {' | '.join(quarantine_notes) or 'no prior invalid target'}"
                ),
            }
        )
        write_manifest(row)
        logging.info("OUTPUT path=%s bytes=%s sha256=%s", TARGET, actual, row["sha256"])
        logging.info(
            "END status=success elapsed_seconds=%.3f peak_ram_mb=%.1f",
            time.perf_counter() - started,
            peak_ram_mb(),
        )
        return 0
    except Exception as exc:
        row["notes"] = f"{row.get('notes', '')} ERROR {type(exc).__name__}: {exc}".strip()
        if PART.exists():
            row["actual_bytes"] = PART.stat().st_size
        write_manifest(row)
        logging.exception(
            "END status=error elapsed_seconds=%.3f peak_ram_mb=%.1f",
            time.perf_counter() - started,
            peak_ram_mb(),
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
