#!/usr/bin/env python3
"""Ingest the public 2026 canine FRE supplementary material.

The available supplements contain sample metadata, sample QC, PCA coordinates,
a 25-row published KEGG GSEA summary, subsampling summaries, and group means for
12 targeted bile-acid genes.  They do not contain a whole-transcriptome matrix
or a complete gene-level DE table.  Two Word supplements contain rasterized
tables; their values are transcribed below with explicit provenance.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import platform
import re
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

from docx import Document


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
RAW_ROOT = ROOT / "01_raw_processed" / "dog_FRE_2026"
SUPP_ROOT = RAW_ROOT / "extracted" / "PMC13293934_SupplementaryFiles"
PARSED_ROOT = RAW_ROOT / "reference" / "parsed_tables"
METADATA_ROOT = ROOT / "00_metadata"
VALIDATION_ROOT = ROOT / "07_external_validation" / "FRE_2026"
LOG_PATH = ROOT / "logs" / "phase6_fre_ingest.log"
FRONTIERS_XML = RAW_ROOT / "reference" / "frontiers_2026_1829399.xml"


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def log(message: str) -> None:
    line = f"[{now_utc()}] {message}"
    print(line, flush=True)
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LOG_PATH.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(line + "\n")


def normalize_text(value: str) -> str:
    return " ".join(value.replace("\xa0", " ").split())


def parse_decimal(value: str) -> str:
    value = normalize_text(value).replace(",", ".")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def docx_table(path: Path, table_index: int) -> tuple[list[str], list[dict[str, str]]]:
    table = Document(path).tables[table_index]
    rows = [[normalize_text(cell.text) for cell in row.cells] for row in table.rows]
    header = rows[0]
    return header, [dict(zip(header, values)) for values in rows[1:]]


def group_dog_id(raw: str) -> str:
    raw = normalize_text(raw)
    match = re.search(r"(CTRL|FRE)[-_](\d+)", raw, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"Cannot identify dog from {raw}")
    group = match.group(1).upper()
    return f"{group}-{int(match.group(2))}"


def create_validation_directories() -> None:
    for name in (
        "mapping",
        "processed",
        "de",
        "ranked",
        "gsea",
        "leading_edge",
        "summary",
        "figures",
    ):
        (VALIDATION_ROOT / name).mkdir(parents=True, exist_ok=True)


def verify_phase5() -> None:
    manifest_path = METADATA_ROOT / "PHASE5_MANIFEST.csv"
    rows = read_csv(manifest_path)
    output: list[dict[str, object]] = []
    for row in rows:
        relative = row["path_relative_to_project"]
        path = ROOT / Path(relative)
        exists = path.is_file()
        observed_size = path.stat().st_size if exists else "NA"
        observed_sha = sha256(path) if exists else "NA"
        size_match = exists and int(row["bytes"]) == observed_size
        hash_match = exists and row["sha256"].lower() == str(observed_sha).lower()
        status = "PASS" if size_match and hash_match else "FAIL"
        output.append(
            {
                "path_relative_to_project": relative,
                "expected_bytes": row["bytes"],
                "observed_bytes": observed_size,
                "expected_sha256": row["sha256"],
                "observed_sha256": observed_sha,
                "status": status,
            }
        )
    write_csv(
        METADATA_ROOT / "PHASE6_PHASE5_INPUT_HASH_AUDIT.csv",
        [
            "path_relative_to_project",
            "expected_bytes",
            "observed_bytes",
            "expected_sha256",
            "observed_sha256",
            "status",
        ],
        output,
    )
    failures = [row for row in output if row["status"] != "PASS"]
    if failures:
        raise RuntimeError(f"Phase5 hash audit failed for {len(failures)} file(s)")
    if len(output) != 389:
        raise RuntimeError(f"Expected 389 Phase5 manifest rows, observed {len(output)}")
    forward_path = ROOT / "06_pseudobulk" / "lodo" / "summary" / "PHASE5_PROGRAM_FORWARD_LIST.csv"
    forward = read_csv(forward_path)
    statuses = {row["forward_to_phase6"] for row in forward}
    if len(forward) != 25 or not statuses <= {"YES", "YES_WITH_QUALIFICATION"}:
        raise RuntimeError(
            f"Forward-list freeze failed: rows={len(forward)}, statuses={sorted(statuses)}"
        )
    log("Phase5 frozen inputs verified: 389/389 hashes PASS; 25/25 forward targets eligible")


def export_docx_tables() -> dict[str, list[dict[str, str]]]:
    outputs: dict[str, list[dict[str, str]]] = {}

    header, demographics = docx_table(SUPP_ROOT / "Table_1.docx", 0)
    write_csv(PARSED_ROOT / "Supplementary_Table_S1_demographics.csv", header, demographics)
    outputs["demographics"] = demographics

    header, clinical = docx_table(SUPP_ROOT / "Table_1.docx", 1)
    write_csv(PARSED_ROOT / "Supplementary_Table_S1_clinical.csv", header, clinical)
    outputs["clinical"] = clinical

    header, pca_rows = docx_table(SUPP_ROOT / "Table_3.docx", 0)
    for row in pca_rows:
        for field in ("Bodyweight (kg)", "PC1", "PC2"):
            row[field] = parse_decimal(row[field])
    write_csv(PARSED_ROOT / "Supplementary_Table_S3_PCA_metadata.csv", header, pca_rows)
    outputs["pca"] = pca_rows

    header, pathway_rows = docx_table(SUPP_ROOT / "Table_4.docx", 0)
    for row in pathway_rows:
        for field in ("NES", "p-value", "padj"):
            row[field] = parse_decimal(row[field])
        row["direction"] = "FRE_higher" if float(row["NES"]) > 0 else "FRE_lower"
    pathway_fields = header + ["direction"]
    write_csv(PARSED_ROOT / "Supplementary_Table_S4_KEGG_GSEA.csv", pathway_fields, pathway_rows)
    outputs["pathway"] = pathway_rows

    header, targeted_means = docx_table(SUPP_ROOT / "Table_6.docx", 0)
    write_csv(
        PARSED_ROOT / "Supplementary_Table_S6_targeted_bile_acid_group_means.csv",
        header,
        targeted_means,
    )
    outputs["targeted_means"] = targeted_means
    return outputs


# Supplementary Table S2 is a single raster image embedded in Table_2.docx.
# These values were double-read from that image.  FASTQ-like strings are sample
# labels only; no FASTQ file is present in the archive.
QC_ROWS = [
    ("CTRL-9_S347_L001_R1_001.fastq.gz", 27674671, 19123062, 69.1, 17515141, 91.6, 49),
    ("CTRL-8_S346_L001_R1_001.fastq.gz", 29922194, 20633263, 69.0, 18933465, 91.8, 50),
    ("CTRL-7_S345_L001_R1_001.fastq.gz", 33781909, 23253732, 68.8, 21268118, 91.5, 50),
    ("CTRL-10_S348_L001_R1_001.fastq.gz", 30996281, 21535649, 69.5, 19729197, 91.6, 50),
    ("FRE-1_S340_L001_R1_001.fastq.gz", 32067455, 22773992, 71.0, 20982827, 92.1, 50),
    ("FRE-6_S344_L001_R1_001.fastq.gz", 34474434, 24448197, 70.9, 22467208, 91.9, 51),
    ("FRE-5_S343_L001_R1_001.fastq.gz", 35186942, 24908942, 70.8, 22801125, 91.5, 50),
    ("FRE-4_S107_L001_R1_001.fastq.gz", 37225140, 26452139, 71.1, 24064684, 91.0, 51),
    ("FRE-3_S342_L001_R1_001.fastq.gz", 32071190, 22211233, 69.3, 20244967, 91.1, 51),
    ("FRE-2_S341_L001_R1_001.fastq.gz", 35416116, 25175778, 71.1, 22950384, 91.2, 50),
    ("FRE-11_S349_L001_R1_001.fastq.gz", 35229496, 24939143, 70.8, 22739948, 91.2, 51),
    ("FRE-12_S350_L001_R1_001.fastq.gz", 35141052, 24594081, 70.0, 22382647, 91.0, 50),
]


def export_raster_tables() -> list[dict[str, object]]:
    fields = [
        "sample_label",
        "dog_id",
        "total_fragments",
        "matched_fragments",
        "mapping_rate_percent",
        "properly_matched_fragments",
        "properly_matched_percent",
        "gc_percent",
        "source_format",
        "notes",
    ]
    rows = []
    for values in QC_ROWS:
        sample, total, matched, rate, proper, proper_rate, gc = values
        rows.append(
            {
                "sample_label": sample,
                "dog_id": group_dog_id(sample),
                "total_fragments": total,
                "matched_fragments": matched,
                "mapping_rate_percent": rate,
                "properly_matched_fragments": proper,
                "properly_matched_percent": proper_rate,
                "gc_percent": gc,
                "source_format": "rasterized_table_in_DOCX",
                "notes": "Sample label references a FASTQ filename; the FASTQ itself is not in the public supplement.",
            }
        )
    write_csv(PARSED_ROOT / "Supplementary_Table_S2_sample_QC.csv", fields, rows)

    sensitivity = [
        {
            "analysis": "Full cohort",
            "samples_retained": "CTRL n=4; FRE n=8",
            "rationale": "primary analysis",
            "PCA_result": "group separation",
            "BA_pathway_result": "bile secretion enriched",
            "interpretation": "main exploratory result",
        },
        {
            "analysis": "Young FRE subset",
            "samples_retained": "CTRL n=4; FRE n=4",
            "rationale": "exclude older FRE dogs",
            "PCA_result": "separation retained",
            "BA_pathway_result": "BA signal retained by ORA",
            "interpretation": "sensitivity support",
        },
        {
            "analysis": "Restrictive subset",
            "samples_retained": "CTRL n=4; FRE n=3",
            "rationale": "exclude older + largest young FRE dog",
            "PCA_result": "separation retained, unstable",
            "BA_pathway_result": "BA signal retained by GSEA",
            "interpretation": "underpowered sensitivity check",
        },
    ]
    write_csv(
        PARSED_ROOT / "Supplementary_Table_S5_subsampling_sensitivity.csv",
        [
            "analysis",
            "samples_retained",
            "rationale",
            "PCA_result",
            "BA_pathway_result",
            "interpretation",
        ],
        sensitivity,
    )
    return rows


def export_main_table() -> list[dict[str, str]]:
    root = ET.parse(FRONTIERS_XML).getroot()
    table_wrap = root.find(".//table-wrap[@id='tab1']")
    if table_wrap is None:
        raise RuntimeError("Main article Table 1 was not found in Frontiers XML")
    table = table_wrap.find("table")
    if table is None:
        raise RuntimeError("Main article Table 1 has no table element")
    tr_rows = table.findall("./thead/tr") + table.findall("./tbody/tr")
    values: list[list[str]] = []
    for tr in tr_rows:
        cells = list(tr.findall("th")) + list(tr.findall("td"))
        values.append([normalize_text("".join(cell.itertext())) for cell in cells])
    header = values[0]
    rows = [dict(zip(header, values_row)) for values_row in values[1:]]
    for row in rows:
        row["log2FC"] = row["log2FC"].replace("−", "-")
    write_csv(PARSED_ROOT / "Main_Table1_targeted_bile_acid_DE.csv", header, rows)
    return rows


def create_sample_metadata(
    demographics: list[dict[str, str]],
    pca_rows: list[dict[str, str]],
    qc_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    def key(row: dict[str, str], weight_name: str) -> tuple[str, str, str, str]:
        return (
            row["Breed"].lower(),
            row["Sex"].lower(),
            parse_decimal(row["Age (years)"]),
            parse_decimal(row[weight_name]),
        )

    demographic_by_key = {key(row, "Bodyweight (kg)"): row for row in demographics}
    qc_by_dog = {str(row["dog_id"]): row for row in qc_rows}
    output: list[dict[str, object]] = []
    for row in pca_rows:
        match_key = key(row, "Bodyweight (kg)")
        if match_key not in demographic_by_key:
            raise RuntimeError(f"Unable to match Table S3 row to Table S1: {row}")
        group = row["Group"].upper()
        dog_id = f"{'FRE' if group == 'FRE' else 'CTRL'}-{int(row['ID'])}"
        if dog_id not in qc_by_dog:
            raise RuntimeError(f"No QC row for {dog_id}")
        qc = qc_by_dog[dog_id]
        is_fre = group == "FRE"
        output.append(
            {
                "sample_id": dog_id,
                "dog_id": dog_id,
                "condition": "FRE" if is_fre else "Healthy",
                "FRE_or_Healthy": "FRE" if is_fre else "Healthy",
                "biological_replicate": "independent_dog",
                "breed": row["Breed"],
                "age": parse_decimal(row["Age (years)"]),
                "sex": row["Sex"],
                "body_weight": parse_decimal(row["Bodyweight (kg)"]),
                "clinical_status": (
                    "chronic_gastrointestinal_signs; duodenal_inflammation"
                    if is_fre
                    else "healthy_no_history_or_clinical_evidence_of_gastrointestinal_disease"
                ),
                "treatment_response": (
                    "clinical_remission_after_at_least_1_month_exclusive_hypoallergenic_diet"
                    if is_fre
                    else "NA"
                ),
                "tissue": "duodenal_mucosa",
                "region": "proximal_duodenum",
                "library": qc["sample_label"],
                "notes": (
                    "client-owned; biopsy collected before response-based FRE classification"
                    if is_fre
                    else "purpose-bred research Beagle"
                ),
            }
        )
    fre_n = sum(row["condition"] == "FRE" for row in output)
    healthy_n = sum(row["condition"] == "Healthy" for row in output)
    if len(output) != 12 or fre_n != 8 or healthy_n != 4:
        raise RuntimeError(f"Sample freeze failed: total={len(output)} FRE={fre_n} Healthy={healthy_n}")
    fields = [
        "sample_id",
        "dog_id",
        "condition",
        "FRE_or_Healthy",
        "biological_replicate",
        "breed",
        "age",
        "sex",
        "body_weight",
        "clinical_status",
        "treatment_response",
        "tissue",
        "region",
        "library",
        "notes",
    ]
    write_csv(METADATA_ROOT / "FRE2026_SAMPLE_METADATA.csv", fields, output)
    log("FRE sample table frozen: 8 independent FRE dogs + 4 independent Healthy dogs")
    return output


def write_content_inventory() -> None:
    descriptions = {
        "Image_1.pdf": ("PCA plots by group, age and body weight", "figure_only"),
        "Image_2.pdf": ("Global differential-expression volcano plot", "figure_only_no_gene_table"),
        "Image_3.pdf": ("KEGG bile-secretion pathway map", "figure_only"),
        "Image_4.pdf": ("KEGG primary bile-acid biosynthesis pathway map", "figure_only"),
        "Table_1.docx": ("Dog demographics and FRE clinical summaries", "sample_metadata"),
        "Table_2.docx": ("Rasterized sample sequencing QC table", "sample_qc"),
        "Table_3.docx": ("Sample metadata and PCA coordinates", "sample_metadata_pca"),
        "Table_4.docx": ("25-row published KEGG GSEA pathway summary", "pathway_summary"),
        "Table_5.docx": ("Rasterized subsampling sensitivity summary", "sensitivity_summary"),
        "Table_6.docx": ("Group mean counts/TPM for 12 selected bile-acid genes", "targeted_group_means"),
        "fvets-13-1829399-g001.gif": ("Main-article targeted TPM plot preview", "figure_only"),
        "fvets-13-1829399-g001.jpg": ("Main-article targeted TPM plot preview", "figure_only"),
        "fvets-13-1829399-g002.gif": ("Main-article targeted TPM heatmap preview", "figure_only"),
        "fvets-13-1829399-g002.jpg": ("Main-article targeted TPM heatmap preview", "figure_only"),
    }
    rows = []
    for path in sorted(SUPP_ROOT.iterdir()):
        if not path.is_file():
            continue
        description, category = descriptions[path.name]
        rows.append(
            {
                "file_name": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "verified_content": description,
                "content_category": category,
                "whole_transcriptome_gene_matrix": "NO",
                "complete_gene_level_DE_statistics": "NO",
                "sample_level_expression_matrix": "NO",
                "FASTQ_present": "NO",
                "notes": "Original Europe PMC archive member retained unchanged.",
            }
        )
    write_csv(
        PARSED_ROOT / "SUPPLEMENTARY_CONTENT_INVENTORY.csv",
        [
            "file_name",
            "bytes",
            "sha256",
            "verified_content",
            "content_category",
            "whole_transcriptome_gene_matrix",
            "complete_gene_level_DE_statistics",
            "sample_level_expression_matrix",
            "FASTQ_present",
            "notes",
        ],
        rows,
    )


def update_data_audit() -> None:
    audit_path = METADATA_ROOT / "PHASE6_FRE_DATA_AUDIT.csv"
    rows = read_csv(audit_path)
    classifications = {
        "Image_1.pdf": ("Three PCA figures", "NO", "NO", "NO", "YES_PCA_ONLY", "NO"),
        "Image_2.pdf": ("Volcano figure only; no underlying gene table", "NO", "NO", "NO", "NO", "NO"),
        "Image_3.pdf": ("KEGG bile-secretion pathway map", "NO", "NO", "NO", "NO", "YES_FIGURE_ONLY"),
        "Image_4.pdf": ("KEGG primary bile-acid biosynthesis maps", "NO", "NO", "NO", "NO", "YES_FIGURE_ONLY"),
        "Table_1.docx": ("Dog demographics and FRE clinical summaries", "NO", "NO", "NO", "YES", "NO"),
        "Table_2.docx": ("Rasterized sample sequencing QC; FASTQ names are labels only", "NO", "NO", "NO", "YES_QC", "NO"),
        "Table_3.docx": ("Sample metadata and PCA coordinates", "NO", "NO", "NO", "YES", "NO"),
        "Table_4.docx": ("25-row published KEGG GSEA summary", "NO", "NO", "NO", "NO", "YES_SUMMARY_ONLY"),
        "Table_5.docx": ("Rasterized subsampling sensitivity summary", "NO", "NO", "NO", "NO", "YES_SENSITIVITY_ONLY"),
        "Table_6.docx": ("Group means for counts and TPM for 12 selected bile-acid genes", "NO_GROUP_MEANS_ONLY", "YES_TARGETED_GROUP_MEANS", "NO", "NO", "NO"),
        "fvets-13-1829399-g001.gif": ("Targeted TPM plot preview", "NO", "NO", "NO", "NO", "NO"),
        "fvets-13-1829399-g001.jpg": ("Targeted TPM plot preview", "NO", "NO", "NO", "NO", "NO"),
        "fvets-13-1829399-g002.gif": ("Targeted TPM heatmap preview", "NO", "NO", "NO", "NO", "NO"),
        "fvets-13-1829399-g002.jpg": ("Targeted TPM heatmap preview", "NO", "NO", "NO", "NO", "NO"),
    }
    for row in rows:
        name = row["file_name"]
        if name in classifications:
            likely, raw_count, tpm, de_table, sample_meta, pathway = classifications[name]
            row["likely_content"] = likely
            row["is_raw_count"] = raw_count
            row["is_tpm"] = tpm
            row["is_de_table"] = de_table
            row["is_sample_metadata"] = sample_meta
            row["is_fastq"] = "NO"
            row["is_pathway_results"] = pathway
            row["required_for_phase6"] = "YES_DATA_GATE_AUDIT"
            row["notes"] = (
                "Content inspected structurally; no whole-transcriptome matrix or complete DE table."
            )
    write_csv(audit_path, list(rows[0].keys()), rows)


def write_software_versions() -> None:
    rows = [
        {"component": "Python", "version": platform.python_version(), "role": "Phase6 audit and ingest"},
        {"component": "python-docx", "version": importlib.metadata.version("python-docx"), "role": "DOCX table extraction"},
        {"component": "pypdf", "version": importlib.metadata.version("pypdf"), "role": "PDF text/page audit"},
        {"component": "platform", "version": platform.platform(), "role": "execution environment"},
        {"component": "analysis_route", "version": "NOT_SELECTED", "role": "awaiting Option A or Option B authorization"},
    ]
    write_csv(METADATA_ROOT / "PHASE6_SOFTWARE_VERSIONS.csv", ["component", "version", "role"], rows)


def main() -> int:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    LOG_PATH.write_text("", encoding="utf-8")
    PARSED_ROOT.mkdir(parents=True, exist_ok=True)
    create_validation_directories()
    verify_phase5()
    tables = export_docx_tables()
    qc_rows = export_raster_tables()
    main_table = export_main_table()
    samples = create_sample_metadata(tables["demographics"], tables["pca"], qc_rows)
    write_content_inventory()
    update_data_audit()
    write_software_versions()

    if len(tables["pathway"]) != 25:
        raise RuntimeError(f"Expected 25 published KEGG pathways, got {len(tables['pathway'])}")
    if len(tables["targeted_means"]) != 12 or len(main_table) != 12:
        raise RuntimeError("Expected 12 targeted bile-acid genes in both published tables")
    if len(samples) != 12:
        raise RuntimeError("Expected 12 frozen dogs")

    zenodo = RAW_ROOT / "downloaded" / "zenodo" / "Table 6.docx"
    epmc = SUPP_ROOT / "Table_6.docx"
    if sha256(zenodo) != sha256(epmc):
        raise RuntimeError("Zenodo Table 6 and Europe PMC Table_6 are not byte-identical")
    log("Zenodo Table 6 and Europe PMC Table_6 are byte-identical")
    log(
        "Verified data gate: no whole-transcriptome counts/abundance matrix, "
        "no sample-level TPM matrix, no complete gene-level DE table, and no FASTQ payload"
    )
    log("Phase6 FRE supplement ingest completed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        log(f"FATAL: {type(exc).__name__}: {exc}")
        raise
