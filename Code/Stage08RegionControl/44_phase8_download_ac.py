from __future__ import annotations

import csv
import hashlib
import os
import tarfile
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P8 = ROOT / "08_region_control"
AUDIT = P8 / "metadata" / "HUMAN_AC_REGION_AUDIT.csv"
OUT = P8 / "human_AC_TI" / "raw_processed"
MANIFEST = P8 / "metadata" / "HUMAN_AC_DOWNLOAD_MANIFEST.csv"
LOGS = (ROOT / "logs" / "phase8_region_analysis.log", P8 / "logs" / "phase8_region_analysis.log")

EXPECTED_HOST = "ftp.ncbi.nlm.nih.gov"
REQUIRED_SUFFIXES = ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def log(message: str) -> None:
    line = f"[{utc_now()}] {message}"
    print(line, flush=True)
    for path in LOGS:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(line + "\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_url(url: str) -> str:
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https" or parsed.hostname != EXPECTED_HOST:
        raise RuntimeError(f"Unexpected URL host/scheme: {url}")
    if "/geo/samples/" not in parsed.path or "/suppl/" not in parsed.path:
        raise RuntimeError(f"Not a per-GSM GEO supplementary URL: {url}")
    filename = urllib.parse.unquote(PurePosixPath(parsed.path).name)
    if not filename.lower().endswith(".tar.gz"):
        raise RuntimeError(f"Not a processed tar.gz: {filename}")
    forbidden = ("gse266616_raw", ".fastq", ".sra", "scrna_gene_counts.mtx")
    if any(token in filename.lower() for token in forbidden):
        raise RuntimeError(f"Forbidden target: {filename}")
    return filename


def validate_tar(path: Path) -> int:
    found = {suffix: 0 for suffix in REQUIRED_SUFFIXES}
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        if not members:
            raise RuntimeError(f"Empty archive: {path}")
        for member in members:
            pure = PurePosixPath(member.name)
            if pure.is_absolute() or ".." in pure.parts:
                raise RuntimeError(f"Unsafe member path: {member.name}")
            if member.isfile():
                for suffix in REQUIRED_SUFFIXES:
                    if member.name.lower().endswith(suffix):
                        found[suffix] += 1
    if any(value != 1 for value in found.values()):
        raise RuntimeError(f"10x triplet cardinality failure in {path.name}: {found}")
    return len(members)


def urllib_download(url: str, part: Path, expected: int) -> None:
    current = part.stat().st_size if part.exists() else 0
    if current > expected:
        raise RuntimeError(f"Partial file exceeds expected bytes: {part}")
    headers = {"User-Agent": "Canine-CIE-Phase8A/1.0"}
    if current:
        headers["Range"] = f"bytes={current}-"
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=120) as response:
        status = int(getattr(response, "status", response.getcode()))
        if current and status == 206:
            mode = "ab"
        elif current and status == 200:
            mode = "wb"
            current = 0
        elif not current and status == 200:
            mode = "wb"
        else:
            raise RuntimeError(f"Unexpected HTTP status {status} for resume offset {current}")
        with part.open(mode) as handle:
            while True:
                chunk = response.read(4 * 1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
                if handle.tell() > expected:
                    raise RuntimeError(f"Remote response exceeded expected size for {part.name}")
    if part.stat().st_size != expected:
        raise RuntimeError(f"Byte mismatch after download: {part.stat().st_size} != {expected}")


def process(row: dict[str, str]) -> dict[str, object]:
    filename = validate_url(row["processed_download_url"])
    expected = int(row["processed_file_size"])
    group = "healthy_AC" if row["disease"] == "Healthy" else "crohn_AC_active"
    directory = OUT / group
    directory.mkdir(parents=True, exist_ok=True)
    final = directory / filename
    part = directory / f".{filename}.part"
    status = "REUSED_VERIFIED"
    if not final.is_file() or final.stat().st_size != expected:
        status = "DOWNLOADED"
        last: Exception | None = None
        for attempt in range(1, 5):
            try:
                urllib_download(row["processed_download_url"], part, expected)
                validate_tar(part)
                os.replace(part, final)
                last = None
                break
            except Exception as exc:  # noqa: BLE001
                last = exc
                log(f"RETRY {attempt}/4 {row['GSM']} {type(exc).__name__}: {exc}")
                time.sleep(2 * attempt)
        if last is not None:
            raise last
    if final.stat().st_size != expected:
        raise RuntimeError(f"Final size mismatch: {final}")
    members = validate_tar(final)
    checksum = sha256(final)
    return {
        "GSM": row["GSM"],
        "donor_id": row["donor_id"],
        "disease": row["disease"],
        "state": row["inflammation_state"],
        "region": "AC",
        "treatment": row["treatment"],
        "url": row["processed_download_url"],
        "filename": filename,
        "expected_size": expected,
        "actual_size": final.stat().st_size,
        "repository_checksum": "NA_NOT_PROVIDED",
        "local_sha256": checksum,
        "download_time": utc_now(),
        "archive_valid": "YES",
        "archive_member_count": members,
        "required_10x_triplet_present": "YES",
        "local_path": str(final),
        "status": status,
    }


def write_manifest(rows: list[dict[str, object]]) -> None:
    fields = list(rows[0])
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    temp = MANIFEST.with_suffix(MANIFEST.suffix + ".tmp")
    with temp.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp, MANIFEST)


def main() -> int:
    for path in LOGS:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("", encoding="utf-8")
    with AUDIT.open("r", encoding="utf-8-sig", newline="") as handle:
        all_rows = list(csv.DictReader(handle))
    selected = [row for row in all_rows if row["eligible"].startswith("YES_")]
    expected_total = sum(int(row["processed_file_size"]) for row in selected)
    if len(selected) != 22 or expected_total != 944_416_204:
        raise RuntimeError(f"Frozen AC download scope changed: n={len(selected)} bytes={expected_total}")
    if any(row["donor_id"] == "OR95" for row in selected):
        raise RuntimeError("OR95 reached the AC download scope")
    log(f"START 22 per-GSM AC processed archives expected_bytes={expected_total}; RAW/SRA/FASTQ excluded")
    completed: dict[str, dict[str, object]] = {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(process, row): row for row in selected}
        for future in as_completed(futures):
            result = future.result()
            completed[str(result["GSM"])] = result
            log(f"VERIFIED {result['GSM']} {result['disease']} {result['state']} {len(completed)}/22")
    ordered = [completed[row["GSM"]] for row in selected]
    write_manifest(ordered)
    if len(ordered) != 22 or any(row["archive_valid"] != "YES" for row in ordered):
        raise RuntimeError("AC download manifest completeness failure")
    if sum(int(row["actual_size"]) for row in ordered) != expected_total:
        raise RuntimeError("AC final byte total mismatch")
    if list(OUT.rglob("*.part")):
        raise RuntimeError("Partial files remain after completed AC download")
    log(f"DOWNLOAD_GATE PASS files=22 bytes={expected_total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
