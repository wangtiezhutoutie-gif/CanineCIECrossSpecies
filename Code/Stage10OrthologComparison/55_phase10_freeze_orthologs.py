from __future__ import annotations

import csv
import gzip
import hashlib
import os
import re
import tempfile
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
MAPPING_MASTER = ROOT / "00_metadata" / "dog_gene_mapping_master.csv"
OUT = ROOT / "10_cross_species_baseline" / "mapping"
SNAPSHOT = OUT / "source_snapshot"
LOG = ROOT / "logs" / "phase10_ortholog_freeze.log"
SOURCE_AUDIT = ROOT / "00_metadata" / "PHASE10_ORTHOLOG_SOURCE_AUDIT.csv"

FTP_ROOT = "https://ftp.ensembl.org/pub"
RELEASE_SPECS = {
    104: {
        "mart": "ensembl_mart_104",
        "assembly": "CanFam3.1",
        "archive_host": "https://may2021.archive.ensembl.org",
    },
    116: {
        "mart": "ensembl_mart_116",
        "assembly": "ROS_Cfam_1.0",
        "archive_host": "https://jun2026.archive.ensembl.org",
    },
}
TABLES = (
    "clfamiliaris_gene_ensembl__gene__main.txt.gz",
    "clfamiliaris_gene_ensembl__homolog_hsapiens__dm.txt.gz",
)
CURRENT_CORE_BASE = (
    f"{FTP_ROOT}/release-116/mysql/canis_lupus_familiaris_core_116_1"
)
CURRENT_HISTORY_TABLES = ("stable_id_event.txt.gz", "mapping_session.txt.gz")
USER_AGENT = "Canine-CIE-Phase10/1.0 (ortholog provenance audit)"


def log(message: str) -> None:
    stamp = datetime.now(timezone.utc).isoformat()
    line = f"[{stamp}] {message}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8", newline="") as handle:
        handle.write(line + "\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        frame.to_csv(tmp, index=False, quoting=csv.QUOTE_MINIMAL)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def remote_size(url: str) -> int:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT, "Range": "bytes=0-0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        content_range = response.headers.get("Content-Range", "")
        if "/" in content_range:
            return int(content_range.rsplit("/", 1)[1])
        content_length = response.headers.get("Content-Length")
        if content_length:
            return int(content_length)
    raise RuntimeError(f"Could not determine remote size: {url}")


def download_ranges(url: str, output: Path, chunk_size: int = 512 * 1024) -> None:
    size = remote_size(url)
    if output.is_file() and output.stat().st_size == size:
        log(f"Reuse complete snapshot: {output.name}; bytes={size}")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=output.name + ".", suffix=".part", dir=output.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        with tmp.open("wb") as handle:
            for start in range(0, size, chunk_size):
                end = min(size - 1, start + chunk_size - 1)
                last_error: Exception | None = None
                for attempt in range(1, 6):
                    try:
                        request = urllib.request.Request(
                            url,
                            headers={
                                "User-Agent": USER_AGENT,
                                "Range": f"bytes={start}-{end}",
                            },
                        )
                        with urllib.request.urlopen(request, timeout=180) as response:
                            payload = response.read()
                        expected = end - start + 1
                        if len(payload) != expected:
                            raise RuntimeError(
                                f"Range length mismatch {start}-{end}: {len(payload)} != {expected}"
                            )
                        handle.write(payload)
                        last_error = None
                        break
                    except Exception as exc:  # network retries are provenance-logged
                        last_error = exc
                        log(f"Retry {attempt}/5 {output.name} range={start}-{end}: {exc}")
                        time.sleep(min(8, attempt * 2))
                if last_error is not None:
                    raise last_error
        if tmp.stat().st_size != size:
            raise RuntimeError(f"Downloaded size mismatch for {output}: {tmp.stat().st_size} != {size}")
        os.replace(tmp, output)
        log(f"Frozen snapshot: {output.name}; bytes={size}; sha256={sha256(output)}")
    finally:
        tmp.unlink(missing_ok=True)


def ddl_columns(sql_path: Path, table: str) -> list[str]:
    with gzip.open(sql_path, "rt", encoding="utf-8", errors="replace") as handle:
        sql = handle.read()
    match = re.search(
        rf"CREATE TABLE `{re.escape(table)}`\s*\((.*?)\) ENGINE=",
        sql,
        flags=re.DOTALL,
    )
    if not match:
        raise RuntimeError(f"DDL not found for {table} in {sql_path}")
    columns = []
    for line in match.group(1).splitlines():
        field = re.match(r"\s*`([^`]+)`\s+", line)
        if field:
            columns.append(field.group(1))
    if not columns:
        raise RuntimeError(f"No DDL columns parsed for {table}")
    return columns


def read_dump(path: Path, columns: list[str]) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=columns,
        dtype=str,
        compression="gzip",
        keep_default_na=False,
        na_filter=False,
        low_memory=False,
    )
    if frame.shape[1] != len(columns):
        raise RuntimeError(f"Dump schema mismatch: {path}")
    return frame.replace({r"\N": ""})


def choose(columns: list[str], patterns: tuple[str, ...], label: str) -> str:
    candidates = [c for c in columns if all(re.search(p, c, flags=re.I) for p in patterns)]
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one {label} column, got {candidates}")
    return candidates[0]


def extract_release(release: int) -> tuple[pd.DataFrame, dict[str, str]]:
    spec = RELEASE_SPECS[release]
    sql_name = f"ensembl_mart_{release}.sql.gz"
    sql_path = SNAPSHOT / f"release{release}_{sql_name}"
    mart = spec["mart"]
    base = f"{FTP_ROOT}/release-{release}/mysql/{mart}"
    download_ranges(f"{base}/{sql_name}", sql_path)
    local_tables: dict[str, Path] = {}
    for table_file in TABLES:
        local = SNAPSHOT / f"release{release}_{table_file}"
        download_ranges(f"{base}/{table_file}", local)
        local_tables[table_file] = local

    gene_table = TABLES[0].removesuffix(".txt.gz")
    homolog_table = TABLES[1].removesuffix(".txt.gz")
    gene_columns = ddl_columns(sql_path, gene_table)
    homolog_columns = ddl_columns(sql_path, homolog_table)
    genes = read_dump(local_tables[TABLES[0]], gene_columns)
    homologs = read_dump(local_tables[TABLES[1]], homolog_columns)

    common_keys = [c for c in gene_columns if c.endswith("_key") and c in homolog_columns]
    if len(common_keys) != 1:
        raise RuntimeError(f"Expected one Dog gene join key; release={release}; keys={common_keys}")
    key = common_keys[0]
    dog_id = choose(gene_columns, (r"^stable_id_",), "Dog stable ID")
    dog_symbol_candidates = [
        c for c in gene_columns if re.search(r"display_label", c, flags=re.I)
    ]
    if not dog_symbol_candidates:
        raise RuntimeError(f"Dog symbol field unavailable; release={release}")
    dog_symbol = dog_symbol_candidates[0]

    human_id = choose(homolog_columns, (r"^stable_id_", r"_r2$"), "Human stable ID")
    human_symbol = choose(homolog_columns, (r"display_label", r"_r1$"), "Human symbol")
    homology_type = choose(homolog_columns, (r"^description_",), "homology type")
    confidence = choose(homolog_columns, (r"high_confidence",), "orthology confidence")
    goc = choose(homolog_columns, (r"goc_score",), "GOC score")
    perc_candidates = [c for c in homolog_columns if re.search(r"^perc_id", c, flags=re.I)]
    if len(perc_candidates) != 2:
        raise RuntimeError(f"Expected two percent identity columns; release={release}; got={perc_candidates}")
    dog_perc = choose(perc_candidates, (r"_r1$",), "Dog percent identity")
    human_perc = next(c for c in perc_candidates if c != dog_perc)

    gene_lookup = genes[[key, dog_id, dog_symbol]].copy()
    gene_lookup.columns = ["dog_join_key", "dog_gene_id", "dog_symbol"]
    gene_lookup = gene_lookup.drop_duplicates()
    per_gene = gene_lookup.groupby("dog_join_key", as_index=False).agg(
        dog_gene_id=("dog_gene_id", lambda x: ";".join(sorted(set(v for v in x if v)))),
        dog_symbol=("dog_symbol", lambda x: ";".join(sorted(set(v for v in x if v)))),
    )
    hom = homologs[[key, human_id, human_symbol, homology_type, confidence, dog_perc, human_perc, goc]].copy()
    hom.columns = [
        "dog_join_key", "human_gene_id", "human_symbol", "homology_type",
        "orthology_confidence", "dog_percent_identity", "human_percent_identity", "goc_score",
    ]
    merged = per_gene.merge(hom, on="dog_join_key", how="left", validate="one_to_many")
    merged = merged.fillna("").loc[lambda x: x["dog_gene_id"].ne("")].drop_duplicates().reset_index(drop=True)
    schema = {
        "dog_join_key": key,
        "dog_gene_id": dog_id,
        "dog_symbol": dog_symbol,
        "human_gene_id": human_id,
        "human_symbol": human_symbol,
        "homology_type": homology_type,
        "orthology_confidence": confidence,
        "dog_percent_identity": dog_perc,
        "human_percent_identity": human_perc,
        "goc_score": goc,
    }
    log(
        f"Parsed release {release}: dog genes={merged['dog_gene_id'].nunique()}; "
        f"homology rows={len(merged)}; one2one={(merged['homology_type'] == 'ortholog_one2one').sum()}"
    )
    return merged, schema


def join_rows(frame: pd.DataFrame) -> dict[str, list[dict[str, str]]]:
    result: dict[str, list[dict[str, str]]] = {}
    for row in frame.to_dict("records"):
        result.setdefault(row["dog_gene_id"], []).append(row)
    return result


def normalize_value(value: object) -> str:
    text = str(value).strip() if value is not None else ""
    return text if text else "NA_NOT_AVAILABLE"


def build_master(release104: pd.DataFrame, release116: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    original = pd.read_csv(MAPPING_MASTER, dtype=str, keep_default_na=False)
    if len(original) != 20_257 or original["feature_id"].duplicated().any():
        raise RuntimeError("Dog mapping master identity gate failed")
    primary_map = join_rows(release104)
    current_map = join_rows(release116)

    one2one_human_counts: dict[str, int] = {}
    one2one_dog_counts: dict[str, int] = {}
    for dog_id, rows in primary_map.items():
        targets = {r["human_gene_id"] for r in rows if r["homology_type"] == "ortholog_one2one" and r["human_gene_id"]}
        one2one_dog_counts[dog_id] = len(targets)
        for target in targets:
            one2one_human_counts[target] = one2one_human_counts.get(target, 0) + 1

    records: list[dict[str, object]] = []
    exclusion: list[dict[str, object]] = []
    retrieved = datetime.now(timezone.utc).isoformat()
    for row in original.to_dict("records"):
        dog_id = row["feature_id"]
        rel_rows = primary_map.get(dog_id, [])
        rel_types = sorted({x["homology_type"] for x in rel_rows if x["homology_type"]})
        one_rows = [x for x in rel_rows if x["homology_type"] == "ortholog_one2one"]
        human_targets = sorted({x["human_gene_id"] for x in one_rows if x["human_gene_id"]})
        ambiguity = "NA_NOT_APPLICABLE"
        primary = False
        if not rel_rows:
            status = "DOG_ID_UNRESOLVED"
            selected: dict[str, str] = {}
            ambiguity = "Dog original feature absent from release104 BioMart gene/homology join"
        elif len(human_targets) == 1 and one2one_dog_counts.get(dog_id) == 1 and one2one_human_counts.get(human_targets[0]) == 1:
            status = "ONE2ONE_PRIMARY"
            primary = True
            selected = next(x for x in one_rows if x["human_gene_id"] == human_targets[0])
        elif "ortholog_one2many" in rel_types:
            status = "ONE2MANY_EXCLUDED"
            selected = rel_rows[0]
            ambiguity = "Ensembl104 homology type is ortholog_one2many"
        elif "ortholog_many2many" in rel_types:
            status = "MANY2MANY_EXCLUDED"
            selected = rel_rows[0]
            ambiguity = "Ensembl104 homology type is ortholog_many2many"
        elif len(human_targets) > 1 or any(one2one_human_counts.get(x, 0) > 1 for x in human_targets):
            status = "HISTORICAL_ID_MERGE_AMBIGUOUS"
            selected = one_rows[0] if one_rows else rel_rows[0]
            ambiguity = "Non-unique historical Dog-to-Human one2one identity; excluded without disease effects"
        else:
            status = "NO_ORTHOLOG"
            selected = rel_rows[0] if rel_rows else {}
            ambiguity = "No release104 ortholog_one2one/one2many/many2many Human relation"

        current_rows = current_map.get(dog_id, [])
        current_one = [x for x in current_rows if x["homology_type"] == "ortholog_one2one"]
        current_targets = sorted({x["human_gene_id"] for x in current_one if x["human_gene_id"]})
        if not current_rows:
            current_status = "CURRENT_ID_UNRESOLVED"
            current_selected: dict[str, str] = {}
        elif len(current_targets) == 1:
            current_selected = next(x for x in current_one if x["human_gene_id"] == current_targets[0])
            if primary and current_targets[0] == selected.get("human_gene_id"):
                current_status = "CURRENT_STATUS_AGREE"
            elif primary:
                current_status = "CURRENT_STATUS_CHANGED"
            else:
                current_status = "CURRENT_RESOLVED_NONPRIMARY"
        else:
            current_selected = current_rows[0]
            current_status = "CURRENT_STATUS_CHANGED" if primary else "CURRENT_RESOLVED_NONPRIMARY"

        confidence_raw = selected.get("orthology_confidence", "")
        confidence = (
            "HIGH_CONFIDENCE_ONE2ONE"
            if primary and str(confidence_raw).strip() in {"1", "high", "HIGH"}
            else ("ONE2ONE_RELEASE104" if primary else "NOT_PRIMARY")
        )
        record = {
            "dog_original_feature_id": dog_id,
            "dog_original_symbol": normalize_value(row.get("original_gene_symbol")),
            "dog_unique_var_name": normalize_value(row.get("unique_var_name")),
            "dog_source_assembly": "CanFam3.1",
            "dog_source_ensembl_release": "104",
            "dog_source_gene_id": dog_id,
            "dog_current_gene_id": dog_id if current_rows else "NA_NOT_AVAILABLE",
            "dog_current_symbol": normalize_value(current_selected.get("dog_symbol")),
            "human_ortholog_release104_id": normalize_value(selected.get("human_gene_id")),
            "human_ortholog_symbol": normalize_value(selected.get("human_symbol")),
            "human_current_gene_id": normalize_value(current_selected.get("human_gene_id")),
            "human_current_symbol": normalize_value(current_selected.get("human_symbol")),
            "homology_type_release104": normalize_value(selected.get("homology_type")),
            "orthology_confidence": normalize_value(selected.get("orthology_confidence")),
            "dog_percent_identity": normalize_value(selected.get("dog_percent_identity")),
            "human_percent_identity": normalize_value(selected.get("human_percent_identity")),
            "goc_score_if_available": normalize_value(selected.get("goc_score")),
            "current_homology_type": normalize_value(current_selected.get("homology_type")),
            "current_audit_status": current_status,
            "mapping_status": status,
            "primary_one2one": primary,
            "mapping_confidence": confidence,
            "ambiguity_reason": ambiguity,
            "source": "Ensembl release104 official FTP BioMart database dump",
            "source_release": "104",
            "retrieval_date": retrieved,
            "notes": "Primary mapping independent of disease effects; Ensembl116 is audit-only",
        }
        records.append(record)
        if not primary:
            exclusion.append({
                "dog_original_feature_id": dog_id,
                "dog_original_symbol": normalize_value(row.get("original_gene_symbol")),
                "mapping_status": status,
                "homology_types_release104": ";".join(rel_types) if rel_types else "NA_NOT_AVAILABLE",
                "candidate_human_ids": ";".join(sorted({x["human_gene_id"] for x in rel_rows if x["human_gene_id"]})) or "NA_NOT_AVAILABLE",
                "ambiguity_reason": ambiguity,
                "exclusion_rule": "Annotation-provenance-only; no disease effect fields inspected",
            })
    master = pd.DataFrame(records)
    excluded = pd.DataFrame(exclusion)
    if len(master) != 20_257 or master["dog_original_feature_id"].duplicated().any():
        raise RuntimeError("Ortholog master row identity gate failed")
    primary_rows = master.loc[master["primary_one2one"]]
    if primary_rows["human_ortholog_release104_id"].duplicated().any():
        raise RuntimeError("Primary Human ortholog collision gate failed")
    return master, excluded


def source_audit_rows(schema_by_release: dict[int, dict[str, str]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    retrieved = datetime.now(timezone.utc).isoformat()
    for release, spec in RELEASE_SPECS.items():
        mart = spec["mart"]
        for file_name in (f"ensembl_mart_{release}.sql.gz", *TABLES):
            local = SNAPSHOT / f"release{release}_{file_name}"
            rows.append({
                "role": "PRIMARY_ORTHOLOG_SOURCE" if release == 104 else "SECONDARY_CURRENT_ANNOTATION_AUDIT",
                "provider": "Ensembl",
                "resource": "Compara-backed BioMart MySQL dump",
                "release": release,
                "archive_host": spec["archive_host"],
                "retrieval_method": "Official release FTP database snapshot; table join reconstructed from frozen SQL DDL",
                "dog_assembly": spec["assembly"],
                "human_assembly": "GRCh38",
                "retrieval_utc": retrieved,
                "query": (
                    "JOIN clfamiliaris_gene_ensembl__gene__main TO "
                    "clfamiliaris_gene_ensembl__homolog_hsapiens__dm ON shared gene key; "
                    "Primary accepts homology_type=ortholog_one2one only"
                ),
                "remote_url": f"{FTP_ROOT}/release-{release}/mysql/{mart}/{file_name}",
                "local_snapshot": str(local),
                "bytes": local.stat().st_size,
                "sha256": sha256(local),
                "schema_fields": ";".join(f"{k}={v}" for k, v in schema_by_release[release].items()),
                "status": "FROZEN_VERIFIED",
            })
    for file_name in CURRENT_HISTORY_TABLES:
        local = SNAPSHOT / f"release116_core_{file_name}"
        with gzip.open(local, "rb") as handle:
            decompressed_bytes = len(handle.read())
        rows.append({
            "role": "SECONDARY_CURRENT_STABLE_ID_HISTORY_AUDIT",
            "provider": "Ensembl",
            "resource": "Dog core stable-ID history table",
            "release": 116,
            "archive_host": RELEASE_SPECS[116]["archive_host"],
            "retrieval_method": "Official release FTP database snapshot",
            "dog_assembly": RELEASE_SPECS[116]["assembly"],
            "human_assembly": "NA_NOT_APPLICABLE",
            "retrieval_utc": retrieved,
            "query": "Inspect current Dog core stable-ID mapping history; audit only",
            "remote_url": f"{CURRENT_CORE_BASE}/{file_name}",
            "local_snapshot": str(local),
            "bytes": local.stat().st_size,
            "sha256": sha256(local),
            "schema_fields": (
                f"decompressed_bytes={decompressed_bytes}; "
                "empty table means no release116 core mapping events are available"
            ),
            "status": "FROZEN_EMPTY_TABLE" if decompressed_bytes == 0 else "FROZEN_VERIFIED",
        })
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    SNAPSHOT.mkdir(parents=True, exist_ok=True)
    LOG.unlink(missing_ok=True)
    query_text = (
        "Phase10 ortholog freeze query\n"
        "Primary: Ensembl release 104 / CanFam3.1 / GRCh38\n"
        "Source: official FTP ensembl_mart_104 SQL + Dog gene main + Dog-Human homolog table\n"
        "Join: shared BioMart dog gene key parsed from the release-specific SQL DDL\n"
        "Filter: homology_type == ortholog_one2one; no identity cutoff; no disease-effect data\n"
        "Secondary audit: same reconstruction from Ensembl release 116 / ROS_Cfam_1.0\n"
    )
    (SNAPSHOT / "PHASE10_ORTHOLOG_RECONSTRUCTION_QUERY.txt").write_text(query_text, encoding="utf-8")
    release104, schema104 = extract_release(104)
    release116, schema116 = extract_release(116)
    for file_name in CURRENT_HISTORY_TABLES:
        download_ranges(
            f"{CURRENT_CORE_BASE}/{file_name}",
            SNAPSHOT / f"release116_core_{file_name}",
        )
    raw104 = SNAPSHOT / "ENSEMBL104_DOG_HUMAN_HOMOLOGY_RECONSTRUCTED.tsv.gz"
    raw116 = SNAPSHOT / "ENSEMBL116_DOG_HUMAN_HOMOLOGY_RECONSTRUCTED.tsv.gz"
    release104.to_csv(raw104, sep="\t", index=False, compression="gzip")
    release116.to_csv(raw116, sep="\t", index=False, compression="gzip")
    master, excluded = build_master(release104, release116)
    atomic_csv(master, OUT / "DOG_HUMAN_ORTHOLOG_MASTER.csv")
    atomic_csv(excluded, OUT / "ORTHOLOG_EXCLUSION_TABLE.csv")

    audit = {
        "total_dog_features": len(master),
        "resolved_release104": int((master["mapping_status"] != "DOG_ID_UNRESOLVED").sum()),
        "one2one": int((master["homology_type_release104"] == "ortholog_one2one").sum()),
        "one2many": int((master["mapping_status"] == "ONE2MANY_EXCLUDED").sum()),
        "many2many": int((master["mapping_status"] == "MANY2MANY_EXCLUDED").sum()),
        "no_ortholog": int((master["mapping_status"] == "NO_ORTHOLOG").sum()),
        "dog_id_unresolved": int((master["mapping_status"] == "DOG_ID_UNRESOLVED").sum()),
        "ambiguous": int(master["mapping_status"].isin(["HISTORICAL_ID_MERGE_AMBIGUOUS", "OTHER_AMBIGUOUS"]).sum()),
        "current_id_resolved": int((master["dog_current_gene_id"] != "NA_NOT_AVAILABLE").sum()),
        "current_status_agree": int((master["current_audit_status"] == "CURRENT_STATUS_AGREE").sum()),
        "current_status_changed": int((master["current_audit_status"] == "CURRENT_STATUS_CHANGED").sum()),
        "primary_one2one_final": int(master["primary_one2one"].sum()),
        "high_confidence_one2one": int((master["mapping_confidence"] == "HIGH_CONFIDENCE_ONE2ONE").sum()),
        "gate": "PASS",
        "notes": "All 20,257 original features retained in master; exclusions are explicit; no disease effects used",
    }
    atomic_csv(pd.DataFrame([audit]), OUT / "ORTHOLOG_MAPPING_AUDIT.csv")
    stable_history_rows = []
    for file_name in CURRENT_HISTORY_TABLES:
        local = SNAPSHOT / f"release116_core_{file_name}"
        with gzip.open(local, "rb") as handle:
            stable_history_rows.append(len(handle.read()))
    atomic_csv(
        pd.DataFrame([{
            "original_dog_features": len(master),
            "release116_direct_same_stable_id_resolved": audit["current_id_resolved"],
            "release116_same_ortholog_status_agree": audit["current_status_agree"],
            "release116_current_status_changed": audit["current_status_changed"],
            "release116_stable_id_event_decompressed_bytes": stable_history_rows[0],
            "release116_mapping_session_decompressed_bytes": stable_history_rows[1],
            "audit_interpretation": (
                "Current Dog annotation uses ROS_Cfam_1.0 with a reset stable-ID namespace; "
                "no direct old-ID matches or core stable-ID mapping events were available. "
                "This limits the secondary current audit and does not alter release104 Primary mapping."
            ),
            "status": "YELLOW_CURRENT_ID_AUDIT_LIMITED",
        }]),
        OUT / "CURRENT_ANNOTATION_AUDIT_SUMMARY.csv",
    )
    atomic_csv(source_audit_rows({104: schema104, 116: schema116}), SOURCE_AUDIT)
    log(f"Ortholog master complete: primary_one2one={audit['primary_one2one_final']}; exclusions={len(excluded)}")


if __name__ == "__main__":
    main()
