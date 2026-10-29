from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import tempfile
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
BASE = ROOT / "07_human_crohn" / "IBDverse"
RAW = BASE / "raw_summary"
META = BASE / "metadata"
LOG = ROOT / "logs" / "phase7_ibdverse.log"

STUDY_API = {
    "S-BSST2944": "https://www.ebi.ac.uk/biostudies/api/v1/studies/S-BSST2944",
    "E-MTAB-16999": "https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-16999",
}
ASSETS = {
    "DGE.zip": (
        "https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/944/S-BSST2944/Files/DGE.zip",
        292_385_193,
        "IBDverse precomputed DGE summary statistics",
    ),
    "individual_metadata.tsv": (
        "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/999/E-MTAB-16999/Files/individual_metadata.tsv",
        22_203,
        "public anonymized individual metadata",
    ),
    "sample_metadata.tsv": (
        "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/999/E-MTAB-16999/Files/sample_metadata.tsv",
        59_724,
        "public sample metadata",
    ),
}
USER_AGENT = "Mozilla/5.0 Phase7-IBDverse-Lightweight-Audit/1.0"


def setup_logger() -> logging.Logger:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("phase7_ibdverse")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)sZ %(levelname)s %(message)s", datefmt="%Y-%m-%dT%H:%M:%S")
    file_handler = logging.FileHandler(LOG, encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_bytes(data: bytes, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    tmp = Path(name)
    try:
        tmp.write_bytes(data)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    tmp = Path(name)
    try:
        frame.to_csv(tmp, index=False, quoting=csv.QUOTE_MINIMAL)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def request(url: str, *, range_start: int | None = None) -> urllib.request.Request:
    headers = {"User-Agent": USER_AGENT, "Accept": "*/*"}
    if range_start is not None and range_start > 0:
        headers["Range"] = f"bytes={range_start}-"
    return urllib.request.Request(url, headers=headers)


def api_snapshots(logger: logging.Logger) -> None:
    META.mkdir(parents=True, exist_ok=True)
    for accession, url in STUDY_API.items():
        with urllib.request.urlopen(request(url), timeout=120) as response:
            data = response.read()
        payload = json.loads(data)
        if payload.get("accno") != accession:
            raise RuntimeError(f"BioStudies accession mismatch for {accession}")
        atomic_bytes(json.dumps(payload, ensure_ascii=False, indent=2).encode("utf-8"), META / f"{accession}_api.json")
        logger.info("API %s bytes=%d", accession, len(data))


def download_asset(name: str, url: str, expected: int, logger: logging.Logger) -> Path:
    if ".h5ad" in name.lower() or ".h5ad" in url.lower():
        raise RuntimeError("Full-cell object is prohibited in the Phase7D lightweight route")
    RAW.mkdir(parents=True, exist_ok=True)
    final = RAW / name
    if final.is_file() and final.stat().st_size == expected:
        logger.info("REUSE %s bytes=%d", name, expected)
        return final
    part = RAW / f"{name}.part"
    start = part.stat().st_size if part.is_file() else 0
    req = request(url, range_start=start if start else None)
    with urllib.request.urlopen(req, timeout=300) as response:
        status = response.status
        if start and status != 206:
            part.unlink(missing_ok=True)
            start = 0
            response.close()
            with urllib.request.urlopen(request(url), timeout=300) as fresh:
                with part.open("wb") as handle:
                    while True:
                        chunk = fresh.read(4 * 1024 * 1024)
                        if not chunk:
                            break
                        handle.write(chunk)
        else:
            mode = "ab" if start else "wb"
            with part.open(mode) as handle:
                while True:
                    chunk = response.read(4 * 1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
                    if part.stat().st_size % (64 * 1024 * 1024) < len(chunk):
                        logger.info("DOWNLOAD %s bytes=%d/%d", name, part.stat().st_size, expected)
    if part.stat().st_size != expected:
        raise RuntimeError(f"Size mismatch for {name}: {part.stat().st_size} != {expected}")
    os.replace(part, final)
    logger.info("COMPLETE %s bytes=%d", name, expected)
    return final


def main() -> int:
    logger = setup_logger()
    logger.info("START Phase7D lightweight IBDverse audit")
    api_snapshots(logger)
    rows: list[dict[str, object]] = []
    paths: dict[str, Path] = {}
    for name, (url, expected, role) in ASSETS.items():
        path = download_asset(name, url, expected, logger)
        paths[name] = path
        rows.append({
            "asset": name, "url": url, "expected_bytes": expected,
            "actual_bytes": path.stat().st_size, "sha256": sha256(path),
            "role": role, "status": "PASS", "downloaded_utc": datetime.now(timezone.utc).isoformat(),
        })
    with zipfile.ZipFile(paths["DGE.zip"]) as archive:
        bad = archive.testzip()
        if bad is not None:
            raise RuntimeError(f"DGE.zip CRC failure: {bad}")
        member_rows = [
            {
                "member": item.filename,
                "compressed_bytes": item.compress_size,
                "uncompressed_bytes": item.file_size,
                "crc32": f"{item.CRC:08x}",
                "is_directory": item.is_dir(),
            }
            for item in archive.infolist()
        ]
    atomic_csv(pd.DataFrame(rows), META / "IBDVERSE_LIGHTWEIGHT_DOWNLOAD_MANIFEST.csv")
    atomic_csv(pd.DataFrame(member_rows), META / "IBDVERSE_DGE_ZIP_INVENTORY.csv")
    atomic_csv(pd.DataFrame([{
        "route": "LIGHTWEIGHT_SUMMARY_LEVEL",
        "DGE_zip_available": "YES",
        "DGE_zip_bytes": paths["DGE.zip"].stat().st_size,
        "public_metadata_available": "YES",
        "full_cohort_h5ad_downloaded": "NO",
        "full_cohort_h5ad_repository_bytes": 25_476_124_776,
        "website_legacy_directory_links": "HTTP_404",
        "biostudies_fire_endpoint": "HTTP_200",
        "status": "READY_FOR_SUMMARY_REPLICATION",
    }]), META / "IBDVERSE_AVAILABILITY_GATE.csv")
    logger.info("END Phase7D download audit members=%d", len(member_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
