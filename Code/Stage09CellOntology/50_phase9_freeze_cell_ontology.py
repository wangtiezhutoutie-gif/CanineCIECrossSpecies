#!/usr/bin/env python3
"""Freeze the official Cell Ontology v2026-06-08 asset for Phase9."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
PHASE9 = ROOT / "09_cell_ontology"
REFERENCE = PHASE9 / "reference"
AUDIT = PHASE9 / "audit"
METADATA = ROOT / "00_metadata"
LOGS = ROOT / "logs"

VERSION = "v2026-06-08"
RELEASE_DATE = "2026-06-08"
RELEASE_API = (
    "https://api.github.com/repos/obophenotype/cell-ontology/"
    "releases/tags/v2026-06-08"
)
RELEASE_PAGE = (
    "https://github.com/obophenotype/cell-ontology/releases/tag/v2026-06-08"
)
ASSET_NAME = "cl-base.obo"
EXPECTED_SIZE = 3_757_082
EXPECTED_SHA256 = "c6274fb65cdf740b983a06b74cb678814249d7f8e580c2ec9de54bf204a2b1c7"


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def fetch_release() -> dict:
    request = urllib.request.Request(
        RELEASE_API,
        headers={
            "Accept": "application/vnd.github+json",
            "User-Agent": "Canine-CIE-Phase9",
        },
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.load(response)


def download_asset(url: str, target: Path) -> None:
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/octet-stream",
            "User-Agent": "Canine-CIE-Phase9",
        },
    )
    temp = target.with_suffix(target.suffix + ".part")
    with urllib.request.urlopen(request, timeout=300) as response, temp.open("wb") as out:
        while chunk := response.read(1024 * 1024):
            out.write(chunk)
    os.replace(temp, target)


def main() -> int:
    for directory in [
        REFERENCE,
        PHASE9 / "dog",
        PHASE9 / "human",
        PHASE9 / "harmonized",
        AUDIT,
        PHASE9 / "figures",
        PHASE9 / "summary",
        LOGS,
    ]:
        directory.mkdir(parents=True, exist_ok=True)

    required = [
        METADATA / "PHASE8_MANIFEST.csv",
        ROOT / "05_annotations" / "dog_GSE254005_annotated.h5ad",
        ROOT / "05_annotations" / "dog_cell_annotation.csv",
        ROOT
        / "07_human_crohn"
        / "GSE266616"
        / "h5ad"
        / "human_GSE266616_TI_annotated.h5ad",
        ROOT
        / "07_human_crohn"
        / "GSE266616"
        / "annotation"
        / "HUMAN_CELLTYPE_MAPPING.csv",
        ROOT / "08_region_control" / "summary" / "REGION_INTERPRETATION_MATRIX.csv",
        ROOT
        / "07_human_crohn"
        / "GSE266616"
        / "gsea"
        / "DOG_FORWARD_IN_HUMAN_OBSERVATION.csv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required Phase9 inputs: " + "; ".join(missing))

    release = fetch_release()
    if release.get("tag_name") != VERSION:
        raise RuntimeError(f"Unexpected release tag: {release.get('tag_name')!r}")
    assets = {item["name"]: item for item in release.get("assets", [])}
    if ASSET_NAME not in assets:
        raise RuntimeError(f"Official release does not contain {ASSET_NAME}")
    asset = assets[ASSET_NAME]
    api_digest = str(asset.get("digest", ""))
    if api_digest != "sha256:" + EXPECTED_SHA256:
        raise RuntimeError(f"Official API digest changed: {api_digest}")
    if int(asset.get("size", -1)) != EXPECTED_SIZE:
        raise RuntimeError(f"Official API size changed: {asset.get('size')}")

    target = REFERENCE / ASSET_NAME
    if not target.exists() or target.stat().st_size != EXPECTED_SIZE:
        # The release-asset API endpoint is more reliable than the public
        # browser redirect in the managed Windows environment and resolves to
        # the identical immutable release blob.
        download_asset(asset["url"], target)
    actual_sha = sha256_file(target)
    if target.stat().st_size != EXPECTED_SIZE or actual_sha != EXPECTED_SHA256:
        raise RuntimeError(
            f"Downloaded ontology verification failed: size={target.stat().st_size}, "
            f"sha256={actual_sha}"
        )

    retrieval_utc = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    audit_rows = [
        {
            "ontology_name": "Cell Ontology",
            "version": VERSION,
            "release_date": RELEASE_DATE,
            "release_published_utc": release.get("published_at", ""),
            "source_url": asset["browser_download_url"],
            "release_page": RELEASE_PAGE,
            "downloaded_filename": ASSET_NAME,
            "size": target.stat().st_size,
            "SHA256": actual_sha,
            "release_api_digest": api_digest,
            "retrieval_utc": retrieval_utc,
            "verification_status": "PASS",
        }
    ]
    fields = list(audit_rows[0])
    write_csv(METADATA / "PHASE9_ONTOLOGY_SOURCE_AUDIT.csv", audit_rows, fields)

    log = LOGS / "phase9_ontology_freeze.log"
    log.write_text(
        "\n".join(
            [
                f"retrieval_utc={retrieval_utc}",
                f"release={VERSION}",
                f"release_page={RELEASE_PAGE}",
                f"asset_url={asset['browser_download_url']}",
                f"asset_size={target.stat().st_size}",
                f"asset_sha256={actual_sha}",
                "verification_status=PASS",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps(audit_rows[0], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
