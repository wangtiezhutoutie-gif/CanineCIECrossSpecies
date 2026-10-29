#!/usr/bin/env python3
"""Build the Phase7A GSM and donor masters and adjudicate the metadata gate."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import ssl
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
BASE = ROOT / "07_human_crohn" / "GSE266616" / "metadata"
SOFT = BASE / "extracted" / "GSE266616_family.soft"
DATA12_JSON = BASE / "extracted" / "Supplementary_Data12.json"
ZENODO_CELLS = BASE / "source" / "zenodo" / "scRNA_cell_metadata.csv"
ZENODO_CELL_IDS = BASE / "source" / "zenodo" / "scRNA_gene_counts_col_cellIDs.txt"
ZENODO_GENES = BASE / "source" / "zenodo" / "scRNA_gene_counts_row_genes.txt"
META = ROOT / "00_metadata"
REPORT = ROOT / "08_results" / "PHASE7A_METADATA_REPORT.md"
LOG = ROOT / "logs" / "phase7_metadata_audit.log"

SAMPLE_MASTER = META / "GSE266616_HUMAN_SAMPLE_MASTER.csv"
DONOR_MASTER = META / "GSE266616_HUMAN_DONOR_MASTER.csv"
COUNT_AUDIT = META / "GSE266616_79_VS_85_AUDIT.csv"
ZENODO_COVERAGE = META / "GSE266616_ZENODO_METADATA_COVERAGE.csv"
PLATFORM_AUDIT = META / "GSE266616_PRIMARY_PLATFORM_AUDIT.csv"
TREATMENT_PLAN = META / "HUMAN_TREATMENT_SENSITIVITY_PLAN.csv"
GATE = META / "PHASE7A_METADATA_GATE.csv"

USER_AGENT = "Canine-CIE-Phase7A-audit/1.0"


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def log(message: str) -> None:
    line = f"[{now_utc()}] {message}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(line + "\n")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_soft() -> list[dict[str, object]]:
    lines = SOFT.read_text(encoding="utf-8").splitlines()
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    for line in lines:
        if line.startswith("^SAMPLE = "):
            if current is not None:
                samples.append(current)
            current = {"GSM": line.split("=", 1)[1].strip(), "characteristics": {}, "data_processing": []}
            continue
        if current is None or not line.startswith("!Sample_"):
            continue
        key, value = line[1:].split(" = ", 1)
        if key == "Sample_characteristics_ch1":
            if ": " not in value:
                raise RuntimeError(f"Malformed characteristic: {line}")
            characteristic, characteristic_value = value.split(": ", 1)
            current["characteristics"][characteristic] = characteristic_value
        elif key == "Sample_data_processing":
            current["data_processing"].append(value)
        elif key in current:
            existing = current[key]
            if isinstance(existing, list):
                existing.append(value)
            else:
                current[key] = [existing, value]
        else:
            current[key] = value
    if current is not None:
        samples.append(current)
    if len(samples) != 85:
        raise RuntimeError(f"Expected 85 GEO samples, observed {len(samples)}")
    return samples


def canonical_patient(sample: dict[str, object]) -> tuple[str, str]:
    characteristics = sample["characteristics"]
    raw_patient = str(characteristics["patient"])
    title = str(sample["Sample_title"])
    patient_match = re.match(r"^(HA\d+|OR\d+)", raw_patient, flags=re.IGNORECASE)
    title_match = re.match(r"^(HA\d+|OR\d+)", title, flags=re.IGNORECASE)
    if not patient_match or not title_match:
        raise RuntimeError(f"Cannot canonicalize patient for {sample['GSM']}: patient={raw_patient}, title={title}")
    patient = patient_match.group(1).upper()
    title_patient = title_match.group(1).upper()
    if patient != title_patient:
        raise RuntimeError(f"Patient/title mismatch for {sample['GSM']}: {patient} vs {title_patient}")
    note = "DIRECT_GEO_PATIENT_FIELD" if raw_patient.upper() == patient else "CANONICALIZED_SUFFIXED_GEO_PATIENT_FIELD"
    return patient, note


def normalize_treatment(raw: str) -> tuple[str, str]:
    normalized = raw.strip().lower().replace("-", "_")
    mapping = {
        "no": "NO_BIOLOGIC",
        "anti_tnf": "ANTI_TNF",
        "anti_tnfalpha": "ANTI_TNF",
        "anti_il12_il23": "ANTI_IL12_23",
    }
    if normalized not in mapping:
        raise RuntimeError(f"Unexpected biologic treatment category: {raw}")
    return mapping[normalized], "NA_NOT_PROVIDED"


def https_url(url: str) -> str:
    if url.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
        return "https://ftp.ncbi.nlm.nih.gov/" + url.split("ftp.ncbi.nlm.nih.gov/", 1)[1]
    return url


def head_size(url: str) -> tuple[int, str]:
    request = urllib.request.Request(https_url(url), method="HEAD", headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, context=ssl.create_default_context(), timeout=120) as response:
        length = response.headers.get("Content-Length")
        if length is None:
            raise RuntimeError(f"No Content-Length for {url}")
        return int(length), str(response.geturl())


def freeze_samples(samples: list[dict[str, object]]) -> list[dict[str, object]]:
    size_by_gsm: dict[str, tuple[int, str]] = {}
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {
            executor.submit(head_size, str(sample["Sample_supplementary_file_1"])): str(sample["GSM"])
            for sample in samples
        }
        for future in as_completed(futures):
            gsm = futures[future]
            size_by_gsm[gsm] = future.result()
    if len(size_by_gsm) != 85:
        raise RuntimeError(f"Processed-file HEAD audit incomplete: {len(size_by_gsm)}/85")
    log("Verified 85/85 per-GSM processed archive URLs and Content-Length values")

    interim: list[dict[str, object]] = []
    by_patient: dict[str, list[dict[str, object]]] = defaultdict(list)
    for sample in samples:
        characteristics = sample["characteristics"]
        gsm = str(sample["GSM"])
        title = str(sample["Sample_title"])
        patient, patient_note = canonical_patient(sample)
        condition_raw = str(characteristics["condition"])
        tissue_raw = str(characteristics["tissue"])
        sample_region = str(characteristics["sample region"])
        is_healthy = condition_raw.lower() == "control"
        region = "TI" if tissue_raw == "terminal_ileum" else "AC"
        if is_healthy:
            inflammation_state = "Healthy"
            disease = "Healthy"
            healthy_or_crohn = "Healthy"
        else:
            disease = "Crohn"
            healthy_or_crohn = "Crohn"
            inflammation_state = {
                "inf": "Active_inflamed",
                "adjacent": "Adjacent",
                "non": "Non_inflamed",
            }[sample_region]
        treatment_raw = str(characteristics["treatment (biological)"])
        treatment_class, specific_drug = normalize_treatment(treatment_raw)
        url_raw = str(sample["Sample_supplementary_file_1"])
        processed_size, resolved_url = size_by_gsm[gsm]
        accession_number = int(gsm.replace("GSM", ""))
        cohort_record_version = "ORIGINAL_79" if 8252464 <= accession_number <= 8252542 else "POST_SUMMARY_ADDITION_6"
        row = {
            "GSM": gsm,
            "sample_title": title,
            "patient_id": patient,
            "disease": disease,
            "condition": condition_raw,
            "region": region,
            "tissue": tissue_raw,
            "inflammation_state": inflammation_state,
            "healthy_or_crohn": healthy_or_crohn,
            "age": str(characteristics["age"]),
            "race": str(characteristics["race"]),
            "biologic_treatment": treatment_raw,
            "biologic_class": treatment_class,
            "specific_drug_if_available": specific_drug,
            "sample_region_label": sample_region,
            "platform": str(sample["Sample_instrument_model"]),
            "platform_accession": str(sample["Sample_platform_id"]),
            "library_protocol": "10x Genomics 3' RNA v3.1; Cell Ranger v7.0.1; hg38",
            "GEO_processed_filename": Path(urllib.parse.urlparse(url_raw).path).name,
            "processed_file_size": processed_size,
            "processed_download_url": resolved_url,
            "repository_checksum": "NA_NOT_PROVIDED",
            "cohort_record_version": cohort_record_version,
            "patient_id_audit": patient_note,
            "paired_donor": "PENDING",
            "paired_sample_id": "PENDING",
            "primary_contrast_eligible": "PENDING",
            "secondary_paired_eligible": "PENDING",
            "exclusion_reason": "PENDING",
            "metadata_source": "GEO family SOFT and MINiML; Supplementary Data12/Zenodo provide subset cross-check only",
            "notes": "",
        }
        interim.append(row)
        by_patient[patient].append(row)

    primary_candidates = [
        row
        for row in interim
        if row["region"] == "TI"
        and (
            row["healthy_or_crohn"] == "Healthy"
            or (row["healthy_or_crohn"] == "Crohn" and row["inflammation_state"] == "Active_inflamed")
        )
    ]
    candidate_counts = Counter(str(row["patient_id"]) for row in primary_candidates)
    duplicates = {patient: n for patient, n in candidate_counts.items() if n > 1}
    if duplicates:
        raise RuntimeError(f"Primary donor has multiple candidate samples and requires manual audit: {duplicates}")
    primary_gsms = {str(row["GSM"]) for row in primary_candidates}

    for patient, patient_rows in by_patient.items():
        ti_active = [row for row in patient_rows if row["region"] == "TI" and row["inflammation_state"] == "Active_inflamed"]
        ti_secondary = [
            row
            for row in patient_rows
            if row["region"] == "TI" and row["inflammation_state"] in {"Adjacent", "Non_inflamed"}
        ]
        paired = bool(ti_active and ti_secondary)
        for row in patient_rows:
            other_ids = [str(other["GSM"]) for other in patient_rows if other is not row]
            row["paired_donor"] = "YES" if len(patient_rows) > 1 else "NO"
            row["paired_sample_id"] = ";".join(other_ids) if other_ids else "NA"
            row["primary_contrast_eligible"] = "YES_SELECTED" if row["GSM"] in primary_gsms else "NO"
            row["secondary_paired_eligible"] = (
                "YES"
                if paired
                and row["region"] == "TI"
                and row["inflammation_state"] in {"Active_inflamed", "Adjacent", "Non_inflamed"}
                else "NO"
            )
            exclusions = []
            if row["GSM"] not in primary_gsms:
                if row["region"] == "AC":
                    exclusions.append("AC_EXCLUDED_FROM_PRIMARY")
                elif row["inflammation_state"] == "Adjacent":
                    exclusions.append("ADJACENT_SECONDARY_ONLY")
                elif row["inflammation_state"] == "Non_inflamed":
                    exclusions.append("NONINFLAMED_SECONDARY_ONLY")
                else:
                    exclusions.append("NOT_PRIMARY_CONTRAST")
            row["exclusion_reason"] = ";".join(exclusions) if exclusions else "NA"
            notes = []
            if row["cohort_record_version"] == "POST_SUMMARY_ADDITION_6":
                notes.append("One of six AC-only samples added after the original 79-sample submission")
            if row["patient_id_audit"] != "DIRECT_GEO_PATIENT_FIELD":
                notes.append("GEO patient field contained sample suffix; canonical donor prefix agrees with sample title")
            row["notes"] = "; ".join(notes) if notes else "NA"

    fields = [
        "GSM",
        "sample_title",
        "patient_id",
        "disease",
        "condition",
        "region",
        "tissue",
        "inflammation_state",
        "healthy_or_crohn",
        "age",
        "race",
        "biologic_treatment",
        "biologic_class",
        "specific_drug_if_available",
        "sample_region_label",
        "platform",
        "platform_accession",
        "library_protocol",
        "GEO_processed_filename",
        "processed_file_size",
        "processed_download_url",
        "repository_checksum",
        "cohort_record_version",
        "patient_id_audit",
        "paired_donor",
        "paired_sample_id",
        "primary_contrast_eligible",
        "secondary_paired_eligible",
        "exclusion_reason",
        "metadata_source",
        "notes",
    ]
    write_csv(SAMPLE_MASTER, fields, interim)
    return interim


def build_donor_master(samples: list[dict[str, object]]) -> list[dict[str, object]]:
    by_patient: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in samples:
        by_patient[str(row["patient_id"])].append(row)
    donors: list[dict[str, object]] = []
    for patient in sorted(by_patient):
        rows = by_patient[patient]
        group_values = {str(row["healthy_or_crohn"]) for row in rows}
        age_values = {str(row["age"]) for row in rows}
        race_values = {str(row["race"]) for row in rows}
        treatment_raw_values = {str(row["biologic_treatment"]) for row in rows}
        treatment_values = {str(row["biologic_class"]) for row in rows}
        if len(group_values) != 1 or len(age_values) != 1 or len(treatment_values) != 1:
            raise RuntimeError(
                f"Donor metadata conflict for {patient}: group={group_values}, age={age_values}, treatment={treatment_values}"
            )
        primary = [row for row in rows if row["primary_contrast_eligible"] == "YES_SELECTED"]
        if len(primary) > 1:
            raise RuntimeError(f"Multiple primary samples for donor {patient}")
        notes = []
        if len(race_values) > 1:
            notes.append(f"race_metadata_variants={'|'.join(sorted(race_values))}")
        if any(row["cohort_record_version"] == "POST_SUMMARY_ADDITION_6" for row in rows):
            notes.append("post-summary AC-only donor")
        donors.append(
            {
                "patient_id": patient,
                "healthy_or_crohn": next(iter(group_values)),
                "n_total_samples": len(rows),
                "n_TI_samples": sum(row["region"] == "TI" for row in rows),
                "n_AC_samples": sum(row["region"] == "AC" for row in rows),
                "n_active_TI": sum(row["region"] == "TI" and row["inflammation_state"] == "Active_inflamed" for row in rows),
                "n_adjacent_TI": sum(row["region"] == "TI" and row["inflammation_state"] == "Adjacent" for row in rows),
                "n_noninflamed_TI": sum(row["region"] == "TI" and row["inflammation_state"] == "Non_inflamed" for row in rows),
                "age": next(iter(age_values)),
                "race": next(iter(race_values)) if len(race_values) == 1 else ";".join(sorted(race_values)),
                "biologic_treatment": next(iter(treatment_raw_values)) if len(treatment_raw_values) == 1 else ";".join(sorted(treatment_raw_values)),
                "biologic_class": next(iter(treatment_values)),
                "primary_sample_id": str(primary[0]["GSM"]) if primary else "NA",
                "primary_eligible": "YES" if primary else "NO",
                "paired_secondary_eligible": "YES" if any(row["secondary_paired_eligible"] == "YES" for row in rows) else "NO",
                "notes": "; ".join(notes) if notes else "NA",
            }
        )
    donor_groups = Counter(str(row["healthy_or_crohn"]) for row in donors)
    if donor_groups != Counter({"Crohn": 29, "Healthy": 17}):
        raise RuntimeError(
            "Unexpected current sample-derived donor composition: "
            f"{dict(donor_groups)}; expected 17 Healthy + 29 Crohn from the 85 explicit GSM records"
        )
    write_csv(
        DONOR_MASTER,
        [
            "patient_id",
            "healthy_or_crohn",
            "n_total_samples",
            "n_TI_samples",
            "n_AC_samples",
            "n_active_TI",
            "n_adjacent_TI",
            "n_noninflamed_TI",
            "age",
            "race",
            "biologic_treatment",
            "biologic_class",
            "primary_sample_id",
            "primary_eligible",
            "paired_secondary_eligible",
            "notes",
        ],
        donors,
    )
    return donors


def audit_zenodo_coverage(samples: list[dict[str, object]]) -> dict[str, object]:
    data12 = json.loads(DATA12_JSON.read_text(encoding="utf-8"))
    data12_rows = data12["rows"]
    gse_rows = [row for row in data12_rows if row["scRNA GSE Series"] == "GSE266616"]
    gse_gsms = {str(row["scRNA GSM ID"]) for row in gse_rows}
    all_gsms = {str(row["GSM"]) for row in samples}
    if not gse_gsms <= all_gsms or len(gse_rows) != 12:
        raise RuntimeError("Supplementary Data12 GSE266616 mapping did not resolve to 12 GEO samples")

    cell_rows = 0
    libraries: Counter[str] = Counter()
    full_ids: set[str] = set()
    with ZENODO_CELLS.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        unnamed = reader.fieldnames[0]
        for row in reader:
            cell_rows += 1
            libraries[row["Seq_Library_ID"]] += 1
            full_ids.add(row[unnamed])
    with ZENODO_CELL_IDS.open("r", encoding="utf-8") as handle:
        cell_id_rows = sum(1 for line in handle if line.strip())
    with ZENODO_GENES.open("r", encoding="utf-8") as handle:
        gene_rows = sum(1 for line in handle if line.strip())
    if cell_rows != cell_id_rows or cell_rows != len(full_ids):
        raise RuntimeError(
            f"Zenodo cell ID mismatch: metadata={cell_rows}, column IDs={cell_id_rows}, unique metadata IDs={len(full_ids)}"
        )
    gse_library_set = {f"C{i}-{region}" for i in range(1, 7) for region in ("AC", "TI")}
    if not gse_library_set <= set(libraries):
        raise RuntimeError(f"Expected GSE-linked Zenodo libraries missing: {sorted(gse_library_set - set(libraries))}")
    coverage_rows = [
        {
            "audit_item": "Supplementary Data12 total rows",
            "value": len(data12_rows),
            "status": "PASS",
            "interpretation": "Includes 12 GSE266616 healthy samples and 12 GSE242086 pouch/prepouch samples",
        },
        {
            "audit_item": "Supplementary Data12 GSE266616 rows",
            "value": len(gse_rows),
            "status": "PASS_LIMITED_COVERAGE",
            "interpretation": "Six healthy donors with paired TI/AC; no Crohn samples",
        },
        {
            "audit_item": "Zenodo cell metadata rows",
            "value": cell_rows,
            "status": "PASS",
            "interpretation": "Cell metadata and column cell-ID file match one-to-one",
        },
        {
            "audit_item": "Zenodo unique Seq_Library_ID",
            "value": len(libraries),
            "status": "PASS",
            "interpretation": "Contains six GSE266616 healthy donors plus pouch/prepouch study libraries",
        },
        {
            "audit_item": "Zenodo gene-ID rows",
            "value": gene_rows,
            "status": "PASS",
            "interpretation": "Human-native gene symbols for the unavailable 6.1 GB combined matrix",
        },
        {
            "audit_item": "Crohn cell annotation coverage",
            "value": 0,
            "status": "ABSENT",
            "interpretation": "Published Zenodo cell metadata cannot serve as full-cohort Crohn annotation",
        },
    ]
    write_csv(
        ZENODO_COVERAGE,
        ["audit_item", "value", "status", "interpretation"],
        coverage_rows,
    )
    return {
        "data12_total": len(data12_rows),
        "data12_gse": len(gse_rows),
        "zenodo_cells": cell_rows,
        "zenodo_libraries": len(libraries),
        "zenodo_genes": gene_rows,
    }


def write_count_audit(samples: list[dict[str, object]], donors: list[dict[str, object]]) -> None:
    original = [row for row in samples if row["cohort_record_version"] == "ORIGINAL_79"]
    added = [row for row in samples if row["cohort_record_version"] == "POST_SUMMARY_ADDITION_6"]
    rows = [
        {
            "audit_item": "Series overall-design sample count",
            "observed_value": 79,
            "evidence": "GEO Series_overall_design text",
            "interpretation": "Original submission/design count",
            "primary_TI_impact": "NONE",
        },
        {
            "audit_item": "Current GEO GSM count",
            "observed_value": len(samples),
            "evidence": "85 SAMPLE blocks in current SOFT/MINiML",
            "interpretation": "Current live record count",
            "primary_TI_impact": "NONE_AFTER_AUDIT",
        },
        {
            "audit_item": "Original contiguous GSM block",
            "observed_value": len(original),
            "evidence": "GSM8252464-GSM8252542 inclusive; submitted 2024-05-05",
            "interpretation": "Exactly 79 original sample records",
            "primary_TI_impact": "PRIMARY_TI_DERIVED_FROM_CURRENT_METADATA",
        },
        {
            "audit_item": "Post-summary additions",
            "observed_value": len(added),
            "evidence": "GSM8647246-GSM8647251; submitted 2024-11-22",
            "interpretation": "Six AC-only adjacent/inflamed samples from HA98, HA102 and HA104",
            "primary_TI_impact": "NONE_ALL_ASCENDING_COLON",
        },
        {
            "audit_item": "Series overall-design Crohn donor claim",
            "observed_value": 32,
            "evidence": "GEO Series_overall_design text",
            "interpretation": "Not reproduced by the explicit patient fields in the 79 original GSM records",
            "primary_TI_impact": "NONE_PRIMARY_USES_EXPLICIT_SAMPLE_METADATA",
        },
        {
            "audit_item": "Original sample-derived Crohn donors",
            "observed_value": len({str(row["patient_id"]) for row in original if row["healthy_or_crohn"] == "Crohn"}),
            "evidence": "Unique canonical patient IDs across the 79 original GSM records",
            "interpretation": "26 explicit Crohn donors; six fewer than the overall-design claim of 32",
            "primary_TI_impact": "NONE_NO_UNOBSERVED_DONORS_INFERRED",
        },
        {
            "audit_item": "Current unique donors",
            "observed_value": len(donors),
            "evidence": "Canonical patient prefix audit across all 85 GSM",
            "interpretation": "17 Healthy + 29 Crohn; three post-summary Crohn donors are AC-only",
            "primary_TI_impact": "NONE_ADDITIONAL_3_DONORS_AC_ONLY",
        },
    ]
    write_csv(COUNT_AUDIT, ["audit_item", "observed_value", "evidence", "interpretation", "primary_TI_impact"], rows)


def write_platform_and_treatment(samples: list[dict[str, object]], donors: list[dict[str, object]]) -> None:
    primary = [row for row in samples if row["primary_contrast_eligible"] == "YES_SELECTED"]
    platform_rows = []
    for platform in sorted({str(row["platform"]) for row in primary}):
        for condition in ("Healthy", "Crohn"):
            count = sum(row["platform"] == platform and row["healthy_or_crohn"] == condition for row in primary)
            platform_rows.append(
                {
                    "platform": platform,
                    "primary_condition": condition,
                    "donor_n": count,
                    "platform_condition_cell_nonzero": "YES" if count > 0 else "NO",
                    "notes": "Primary model remains ~condition; platform sensitivity only if estimable and non-sparse.",
                }
            )
    write_csv(
        PLATFORM_AUDIT,
        ["platform", "primary_condition", "donor_n", "platform_condition_cell_nonzero", "notes"],
        platform_rows,
    )

    active = [row for row in donors if row["n_active_TI"] == 1]
    treatment_counts = Counter(str(row["biologic_class"]) for row in active)
    plan_rows = []
    for treatment in ("NO_BIOLOGIC", "ANTI_TNF", "ANTI_IL12_23"):
        n = treatment_counts[treatment]
        plan_rows.append(
            {
                "sensitivity": f"{treatment}_active_Crohn_TI_vs_Healthy_TI",
                "active_crohn_donor_n": n,
                "minimum_required": 3,
                "executable": "YES" if n >= 3 else "NO_DESCRIPTIVE_ONLY",
                "design": "~condition" if n >= 3 else "NA",
                "interpretation_limit": "Treatment-confounding sensitivity only; not a causal treatment-effect analysis.",
            }
        )
    for excluded in ("ANTI_TNF", "ANTI_IL12_23"):
        retained_n = len(active) - treatment_counts[excluded]
        plan_rows.append(
            {
                "sensitivity": f"exclude_{excluded}",
                "active_crohn_donor_n": retained_n,
                "minimum_required": 3,
                "executable": "YES" if retained_n >= 3 else "NO_DESCRIPTIVE_ONLY",
                "design": "~condition" if retained_n >= 3 else "NA",
                "interpretation_limit": "Treatment-confounding sensitivity only; not a causal treatment-effect analysis.",
            }
        )
    write_csv(
        TREATMENT_PLAN,
        ["sensitivity", "active_crohn_donor_n", "minimum_required", "executable", "design", "interpretation_limit"],
        plan_rows,
    )


def write_gate_and_report(
    samples: list[dict[str, object]], donors: list[dict[str, object]], coverage: dict[str, object]
) -> None:
    primary = [row for row in samples if row["primary_contrast_eligible"] == "YES_SELECTED"]
    primary_healthy = [row for row in primary if row["healthy_or_crohn"] == "Healthy"]
    primary_active = [row for row in primary if row["healthy_or_crohn"] == "Crohn"]
    active_treatments = Counter(str(row["biologic_class"]) for row in primary_active)
    added = [row for row in samples if row["cohort_record_version"] == "POST_SUMMARY_ADDITION_6"]
    all_urls_valid = all(int(row["processed_file_size"]) > 0 for row in samples)
    unique_primary = len({str(row["patient_id"]) for row in primary}) == len(primary)
    criteria = [
        ("Primary Healthy TI donors uniquely identified", len(primary_healthy) > 0 and len({row['patient_id'] for row in primary_healthy}) == len(primary_healthy), f"n={len(primary_healthy)}"),
        ("Primary active Crohn TI donors uniquely identified", len(primary_active) > 0 and len({row['patient_id'] for row in primary_active}) == len(primary_active), f"n={len(primary_active)}"),
        ("Primary one sample per donor", unique_primary, f"samples={len(primary)} donors={len({row['patient_id'] for row in primary})}"),
        ("Region explicit", all(row["region"] in {"TI", "AC"} for row in samples), "85/85"),
        ("Inflammation state explicit", all(row["inflammation_state"] in {"Healthy", "Active_inflamed", "Adjacent", "Non_inflamed"} for row in samples), "85/85"),
        ("Biologic class traceable", all(row["biologic_class"] in {"NO_BIOLOGIC", "ANTI_TNF", "ANTI_IL12_23"} for row in samples), "85/85"),
        ("Per-GSM processed link and size confirmed", all_urls_valid, "85/85 HEAD checks"),
        ("79-vs-85 discrepancy isolated from Primary TI", len(added) == 6 and all(row["region"] == "AC" for row in added), "6/6 additions are AC"),
        ("Overall-design donor-count discrepancy disclosed", sum(row["healthy_or_crohn"] == "Healthy" for row in donors) == 17 and sum(row["healthy_or_crohn"] == "Crohn" for row in donors) == 29, "sample-derived current count=17 Healthy + 29 Crohn; design text claims 17 + 32"),
    ]
    gate_pass = all(status for _, status, _ in criteria)
    gate_rows = [
        {
            "criterion": criterion,
            "status": "PASS" if status else "FAIL",
            "evidence": evidence,
            "gate": "GREEN_METADATA_GATE" if gate_pass else "YELLOW_HUMAN_METADATA_GATE_NOT_RESOLVED",
        }
        for criterion, status, evidence in criteria
    ]
    write_csv(GATE, ["criterion", "status", "evidence", "gate"], gate_rows)
    if not gate_pass:
        raise RuntimeError("Phase7A metadata gate failed")

    report = f"""# Phase7A Metadata Freeze Report

## Gate

**GREEN — HUMAN METADATA GATE PASSED**

The current GEO record contains 85 GSM records representing {len(donors)} unique donors. Primary terminal-ileum samples were frozen before any expression matrix was downloaded.

## 79 vs 85 resolution

- `GSM8252464` through `GSM8252542` form an exact contiguous block of 79 original samples submitted on 2024-05-05.
- `GSM8647246` through `GSM8647251` are six later samples submitted on 2024-11-22.
- The six additions are paired adjacent/inflamed ascending-colon samples from three Crohn donors: HA98, HA102 and HA104.
- The live Series sample list was updated to 85, while `Series_overall_design` retained the original 79-sample wording.
- All six additions are AC, so they do not change or contaminate the Primary TI contrast.

This is a resolved GEO record-version discrepancy, not an unexplained biopsy split.

Separately, the original 79 explicit GSM records contain 17 Healthy and 26 Crohn donor IDs, whereas `Series_overall_design` claims 17 Healthy and 32 Crohn patients. The six later AC samples add three Crohn donor IDs, yielding the current sample-derived total of 17 Healthy + 29 Crohn. No identities are available for the six claimed-but-unrepresented Crohn donors, so none were inferred. This unresolved design-text/donor-field discrepancy is disclosed in the audit and does not alter the Primary TI samples selected from explicit GSM metadata.

## Donor and sample structure

- Total GSM: {len(samples)}.
- Current unique donors: {len(donors)} ({sum(row['healthy_or_crohn'] == 'Healthy' for row in donors)} Healthy; {sum(row['healthy_or_crohn'] == 'Crohn' for row in donors)} Crohn).
- Donors with at least one TI sample: {sum(int(row['n_TI_samples']) > 0 for row in donors)}.
- Primary Healthy TI: {len(primary_healthy)} independent donors / {len(primary_healthy)} samples.
- Primary active Crohn TI: {len(primary_active)} independent donors / {len(primary_active)} samples.
- Every Primary donor contributes exactly one sample.
- TI paired active/adjacent-or-noninflamed donors: {sum(row['paired_secondary_eligible'] == 'YES' for row in donors)}.

## Active-TI biologic treatment

- No biologic: {active_treatments['NO_BIOLOGIC']} donors.
- Anti-TNF: {active_treatments['ANTI_TNF']} donors.
- Anti-IL12/23: {active_treatments['ANTI_IL12_23']} donors.

Sensitivity feasibility is frozen in `HUMAN_TREATMENT_SENSITIVITY_PLAN.csv`; all analyses are treatment-confounding sensitivities, not causal treatment-effect models.

## Platform

Primary model remains `~condition`. Platform-by-condition counts are frozen in `GSE266616_PRIMARY_PLATFORM_AUDIT.csv`; a platform sensitivity is allowed only if the design is estimable and non-sparse.

## Supplementary Data12 / Zenodo coverage

- Supplementary Data12 contains {coverage['data12_total']} rows, of which {coverage['data12_gse']} map to GSE266616.
- Those 12 rows are paired TI/AC samples from six Healthy donors only.
- Zenodo cell metadata contains {coverage['zenodo_cells']} cells across {coverage['zenodo_libraries']} combined-study libraries and matches the cell-ID list one-to-one.
- Crohn annotation coverage in this Zenodo metadata is zero. It cannot be represented as full GSE266616 published annotation.
- The six Healthy reference donors can be used for annotation/QC provenance checks. The remaining Human cohort requires an explicitly documented independent broad-lineage annotation route.

## Download readiness

- 85/85 per-GSM processed URLs and Content-Length values were confirmed.
- Phase7B is authorized only for TI processed archives selected in `GSE266616_HUMAN_SAMPLE_MASTER.csv`.
- AC, SRA/FASTQ, `GSE266616_RAW.tar`, and the 6.1 GB Zenodo matrix remain excluded.

Generated UTC: {now_utc()}
"""
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(report, encoding="utf-8", newline="\n")
    log(
        f"Phase7A Gate PASS: Primary Healthy TI n={len(primary_healthy)}, Active Crohn TI n={len(primary_active)}"
    )


def main() -> int:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("", encoding="utf-8")
    samples_raw = parse_soft()
    samples = freeze_samples(samples_raw)
    donors = build_donor_master(samples)
    coverage = audit_zenodo_coverage(samples)
    write_count_audit(samples, donors)
    write_platform_and_treatment(samples, donors)
    write_gate_and_report(samples, donors, coverage)
    log("Phase7A metadata audit completed; expression-data download gate may proceed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        log(f"FATAL: {type(exc).__name__}: {exc}")
        raise
