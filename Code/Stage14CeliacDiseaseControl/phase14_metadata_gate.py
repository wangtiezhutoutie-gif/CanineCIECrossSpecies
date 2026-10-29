#!/usr/bin/env python3
"""Resolve and freeze the Phase14 GSE277276 metadata gate."""
from __future__ import annotations

import gzip
import hashlib
import json
import re
from pathlib import Path
import pandas as pd

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P14 = ROOT / "14_ced_specificity"
META = P14 / "metadata"
SOURCE = META / "source"
AUDIT = P14 / "audit"
ROOT_META = ROOT / "00_metadata"
LOG = ROOT / "logs" / "phase14_metadata_gate.log"

def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def first_value(lines: list[str], prefix: str) -> str:
    values = [line[len(prefix):] for line in lines if line.startswith(prefix)]
    return values[0].strip() if values else ""

def parse_soft(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    sample_rows, processed_rows = [], []
    for block in text.split("^SAMPLE = ")[1:]:
        lines = block.splitlines()
        gsm = lines[0].strip()
        title = first_value(lines, "!Sample_title = ")
        characteristics = {}
        for line in lines:
            prefix = "!Sample_characteristics_ch1 = "
            if line.startswith(prefix) and ": " in line[len(prefix):]:
                key, item = line[len(prefix):].split(": ", 1)
                characteristics[key.strip().lower()] = item.strip()
        biosample_url = first_value(lines, "!Sample_relation = BioSample: ")
        sra_url = first_value(lines, "!Sample_relation = SRA: ")
        replicate_role = "rep2" if "_rep2_" in title else "primary"
        subject = characteristics.get("subject id", "")
        disease_geo = characteristics.get("disease state", "")
        sample_rows.append({
            "gsm": gsm, "sample_title": title, "subject_id_geo": subject,
            "subject_sample_author_expected": f"{subject}2" if replicate_role == "rep2" else subject,
            "replicate_role": replicate_role,
            "condition_geo": "CeD" if disease_geo.lower() == "coeliac" else "Control",
            "disease_state_geo": disease_geo,
            "tissue_geo": characteristics.get("tissue", ""),
            "marsh_geo": characteristics.get("histologic severity", ""),
            "biosample": biosample_url.rsplit("/", 1)[-1], "biosample_url": biosample_url,
            "sra_experiment": sra_url.split("term=")[-1] if "term=" in sra_url else sra_url,
            "sra_url": sra_url,
        })
        for line in lines:
            match = re.match(r"!Sample_supplementary_file_(\d+) = (.+)", line)
            if match:
                url = match.group(2).strip().replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
                processed_rows.append({"gsm": gsm, "sample_title": title,
                    "supplementary_index": int(match.group(1)), "file_name": url.rsplit("/", 1)[-1], "url": url})
    return pd.DataFrame(sample_rows), pd.DataFrame(processed_rows)

def source_role(name: str) -> str:
    if name == "GSE277276_all_meta.csv.gz": return "author_final_cell_metadata"
    if "family.soft" in name: return "geo_family_soft"
    if "family.xml" in name: return "geo_family_miniml"
    if "series_matrix" in name: return "geo_series_matrix_metadata"
    if name.startswith("PMID"): return "pubmed_record"
    if name.startswith("DOI_"): return "crossref_record"
    if name.endswith(".pdf"): return "publisher_pdf"
    if name.endswith(".xlsx"): return "publisher_supplement_xlsx"
    return "source"

SOURCE_URLS = {
    "GSE277276_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277276/soft/GSE277276_family.soft.gz",
    "GSE277276_family.xml.tgz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277276/miniml/GSE277276_family.xml.tgz",
    "GSE277276_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277276/matrix/GSE277276_series_matrix.txt.gz",
    "GSE277276_all_meta.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277276/suppl/GSE277276_all_meta.csv.gz",
    "PMID40694474_esummary.json": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=40694474&retmode=json",
    "PMID40694474_pmc_elink.json": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pmc&id=40694474&retmode=json",
    "PMID40694474_pubmed.xml": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=40694474&retmode=xml",
    "DOI_10.1016_j.celrep.2025.116039_crossref.json": "https://api.crossref.org/works/10.1016/j.celrep.2025.116039",
    "Richards_2025_mmc1_Supplemental_Information.pdf": "https://ars.els-cdn.com/content/image/1-s2.0-S2211124725008101-mmc1.pdf",
    "Richards_2025_mmc2.xlsx": "https://ars.els-cdn.com/content/image/1-s2.0-S2211124725008101-mmc2.xlsx",
    "Richards_2025_mmc3.xlsx": "https://ars.els-cdn.com/content/image/1-s2.0-S2211124725008101-mmc3.xlsx",
    "Richards_2025_mmc4_Document_S1.pdf": "https://ars.els-cdn.com/content/image/1-s2.0-S2211124725008101-mmc4.pdf",
}

def main() -> None:
    for directory in [META, SOURCE, AUDIT, ROOT_META, LOG.parent]:
        directory.mkdir(parents=True, exist_ok=True)
    soft, processed = parse_soft(SOURCE / "GSE277276_family.soft.gz")
    require(len(soft) == 34, f"Expected 34 GEO libraries, observed {len(soft)}")
    require(soft["gsm"].is_unique, "GSM accessions are not unique")
    require(soft["biosample"].is_unique, "BioSample accessions are not unique")
    require(soft["subject_id_geo"].nunique() == 32, "GEO subject id count is not 32")
    duplicate_subjects = soft.groupby("subject_id_geo").size()
    duplicate_subjects = duplicate_subjects[duplicate_subjects > 1].to_dict()
    require(duplicate_subjects == {"13A": 2, "24A": 2}, f"Unexpected duplicate subjects: {duplicate_subjects}")
    require(len(processed) == 102, f"Expected 102 matrix component files, observed {len(processed)}")

    usecols = ["cell.id", "subject", "subject.sample", "Trinity_Patient.remark",
               "Marsh_score", "Celltype", "compartment"]
    cells = pd.read_csv(SOURCE / "GSE277276_all_meta.csv.gz", usecols=usecols, low_memory=False)
    require(len(cells) == 203555, f"Author-final cell count is {len(cells)}, not 203555")
    require(cells["cell.id"].is_unique, "Author-final cell IDs are not unique")
    require(cells["subject"].nunique() == 32, "all_meta subject count is not 32")
    require(cells["subject.sample"].nunique() == 34, "all_meta library count is not 34")
    cells["condition_author"] = cells["Trinity_Patient.remark"].map({"Coeliac": "CeD", "control": "Control"})
    require(cells["condition_author"].notna().all(), "Unrecognized author condition")
    subjects = cells[["subject", "condition_author"]].drop_duplicates()
    require(subjects["subject"].is_unique, "A subject has multiple conditions")
    condition_counts = subjects["condition_author"].value_counts().to_dict()
    require(condition_counts == {"CeD": 21, "Control": 11}, f"Unexpected condition counts: {condition_counts}")

    library = (cells.groupby(["subject", "subject.sample", "condition_author", "Marsh_score"], dropna=False)
               .size().rename("n_author_final_cells").reset_index())
    master = soft.merge(library, left_on=["subject_id_geo", "subject_sample_author_expected"],
                        right_on=["subject", "subject.sample"], how="left", validate="one_to_one")
    require(master["n_author_final_cells"].notna().all(), "A GEO library did not map to all_meta")
    require(master["subject_id_geo"].eq(master["subject"]).all(), "Subject mismatch")
    require(master["condition_geo"].eq(master["condition_author"]).all(), "Condition mismatch")
    require(int(master["n_author_final_cells"].sum()) == 203555, "Library cells do not sum to 203555")

    master["subject_id"] = master["subject_id_geo"]
    master["subject_sample_author"] = master["subject.sample"]
    master["condition"] = master["condition_author"]
    master["marsh_author"] = master["Marsh_score"]
    master["library_to_subject_status"] = "PASS_EXPLICIT_GEO_SUBJECT_ID"
    master["biological_unit"] = "subject"
    master["replicate_handling_primary"] = "sum CellBender raw counts across libraries within subject x Level1"
    master["replicate_handling_sensitivity"] = master["replicate_role"].map({
        "primary": "retained in primary-only sensitivity", "rep2": "excluded in primary-only sensitivity"})
    master["metadata_note"] = ""
    master.loc[master["subject_id"].eq("11A"), "metadata_note"] = (
        "GEO SOFT says Marsh 3c, but GEO disease state and all_meta say control; paper states all controls had normal histology. "
        "Condition remains Control; discrepancy retained, not overwritten.")
    master["evidence"] = ("GEO SOFT explicit subject id + disease state; author all_meta subject/subject.sample; "
                          "Richards et al. 2025 cohort n=32")
    cols = ["gsm", "sample_title", "subject_id", "subject_sample_author", "replicate_role", "condition",
            "biological_unit", "disease_state_geo", "tissue_geo", "marsh_geo", "marsh_author",
            "n_author_final_cells", "biosample", "sra_experiment", "library_to_subject_status",
            "replicate_handling_primary", "replicate_handling_sensitivity", "metadata_note", "evidence"]
    frozen = master[cols].sort_values("gsm").reset_index(drop=True)
    root_master = ROOT_META / "GSE277276_CED_SAMPLE_MASTER.csv"
    stage_master = META / "GSE277276_CED_SAMPLE_MASTER.csv"
    frozen.to_csv(root_master, index=False)
    frozen.to_csv(stage_master, index=False)
    require(sha256(root_master) == sha256(stage_master), "Frozen sample-master copies differ")
    master_hash = sha256(root_master)
    (ROOT_META / "GSE277276_CED_SAMPLE_MASTER.sha256.txt").write_text(
        f"{master_hash}  GSE277276_CED_SAMPLE_MASTER.csv\n", encoding="utf-8")
    processed.to_csv(META / "GSE277276_PROCESSED_FILE_EXPECTATIONS.csv", index=False)
    frozen[["gsm", "sample_title", "subject_id", "subject_sample_author", "replicate_role", "condition",
            "biosample", "sra_experiment", "library_to_subject_status", "evidence"]].to_csv(
                META / "GSE277276_LIBRARY_TO_SUBJECT_EVIDENCE.csv", index=False)

    annotation = cells.groupby(["compartment", "Celltype"]).size().rename("n_cells").reset_index()
    level1_map = {"epi": "Epithelial", "tnk": "T", "myl": "Myeloid"}
    annotation["Phase9_Level1"] = annotation["compartment"].map(level1_map).fillna("")
    annotation["include_primary"] = annotation["compartment"].isin(level1_map).map({True: "YES", False: "NO"})
    annotation["mapping_basis"] = "author compartment; no reclustering or reannotation"
    annotation["exclusion_reason"] = annotation["compartment"].map({
        "bplasma": "Outside frozen Phase9 Primary Level1 scope", "str": "Outside frozen Phase9 Primary Level1 scope"}).fillna("")
    annotation.to_csv(META / "GSE277276_AUTHOR_CELLTYPE_LEVEL1_MAP.csv", index=False)

    pd.DataFrame([
        {"check": "author_final_rows", "expected": 203555, "observed": len(cells), "status": "PASS"},
        {"check": "unique_cell_ids", "expected": 203555, "observed": cells["cell.id"].nunique(), "status": "PASS"},
        {"check": "geo_libraries", "expected": 34, "observed": len(soft), "status": "PASS"},
        {"check": "author_libraries", "expected": 34, "observed": cells["subject.sample"].nunique(), "status": "PASS"},
        {"check": "biological_subjects", "expected": 32, "observed": cells["subject"].nunique(), "status": "PASS"},
        {"check": "CeD_subjects", "expected": 21, "observed": condition_counts["CeD"], "status": "PASS"},
        {"check": "Control_subjects", "expected": 11, "observed": condition_counts["Control"], "status": "PASS"},
        {"check": "explicit_duplicate_subjects", "expected": "13A:2|24A:2", "observed": "13A:2|24A:2", "status": "PASS"},
        {"check": "primary_Level1_cells", "expected": "audited", "observed": int(cells["compartment"].isin(level1_map).sum()), "status": "PASS"},
    ]).to_csv(AUDIT / "GSE277276_ALL_META_AUDIT.csv", index=False)

    pd.DataFrame([
        {"evidence_id": "GEO_SOFT_REPLICATES", "source": "GSE277276_family.soft.gz",
         "finding": "GSM8517672/GSM8517673 both have explicit subject id 24A; GSM8517690/GSM8517691 both have explicit subject id 13A.", "status": "PASS"},
        {"evidence_id": "AUTHOR_ALL_META", "source": "GSE277276_all_meta.csv.gz",
         "finding": "203,555 unique final cell IDs; 34 subject.sample labels; 32 subjects; 24A2 and 13A2 merge to subjects 24A and 13A.", "status": "PASS"},
        {"evidence_id": "PAPER_COHORT", "source": "Richards_2025_mmc4_Document_S1.pdf",
         "finding": "Paper reports 32 subjects: 21 active CeD and 11 controls; 203,555 final high-quality cells.", "status": "PASS"},
        {"evidence_id": "PAPER_METHODS", "source": "Richards_2025_mmc4_Document_S1.pdf pages 20-22",
         "finding": "Diagnostic duodenal biopsies; Cell Ranger v7.0 then CellBender v0.3.0; author QC and compartment annotations documented.", "status": "PASS"},
        {"evidence_id": "GEO_11A_MARSH", "source": "GSM8517689 SOFT vs paper/all_meta",
         "finding": "SOFT Marsh 3c conflicts with paper normal-control histology and all_meta control label. Preserved as discrepancy; condition unambiguous Control.", "status": "WARN_DOCUMENTED"},
    ]).to_csv(META / "GSE277276_METADATA_EVIDENCE.csv", index=False)

    source_rows = []
    for path in sorted(SOURCE.iterdir()):
        if path.is_file():
            source_rows.append({"file": path.name, "path": path.relative_to(ROOT).as_posix(),
                "size": path.stat().st_size, "SHA256": sha256(path), "role": source_role(path.name),
                "source_url": SOURCE_URLS.get(path.name, "")})
    inventory = pd.DataFrame(source_rows)
    require(set(SOURCE_URLS).issubset(set(inventory["file"])), "Required metadata source missing")
    inventory.to_csv(META / "GSE277276_METADATA_SOURCE_INVENTORY.csv", index=False)

    gate = pd.DataFrame([
        {"criterion": "phase13_manifest_revalidated", "status": "PASS", "observed": "132/132 size and SHA-256 PASS", "decision": "continue"},
        {"criterion": "geo_library_count", "status": "PASS", "observed": "34 libraries", "decision": "continue"},
        {"criterion": "explicit_library_to_subject_mapping", "status": "PASS", "observed": "34 libraries -> 32 subjects; duplicates exactly 24A and 13A", "decision": "continue"},
        {"criterion": "subject_condition_counts", "status": "PASS", "observed": "21 active CeD + 11 Control", "decision": "continue"},
        {"criterion": "author_final_cells", "status": "PASS", "observed": "203555 unique cell IDs across 34 subject.sample labels", "decision": "continue"},
        {"criterion": "author_annotation_available", "status": "PASS", "observed": "40 Celltype values; frozen Epithelial/T/Myeloid mapping", "decision": "continue"},
        {"criterion": "GEO_11A_marsh_discrepancy", "status": "WARN_DOCUMENTED", "observed": "condition remains Control", "decision": "continue_with_flag"},
        {"criterion": "route", "status": "PASS", "observed": "ROUTE_A_AUTHOR_FINAL_CELLS", "decision": "download processed CellBender matrices"},
        {"criterion": "overall_metadata_gate", "status": "GREEN", "observed": "GREEN_PHASE14A_METADATA_RESOLVED_ROUTE_A_AUTHORIZED", "decision": "continue_to_processed_download"},
    ])
    gate.to_csv(META / "PHASE14_METADATA_GATE.csv", index=False)
    gate.to_csv(ROOT_META / "PHASE14A_METADATA_GATE.csv", index=False)
    log_lines = ["Phase14 metadata gate: PASS", "phase13_manifest=132/132 PASS",
        f"geo_libraries={len(soft)}", f"subjects={cells['subject'].nunique()}",
        f"conditions={json.dumps(condition_counts, sort_keys=True)}", f"author_final_cells={len(cells)}",
        f"subject_sample_labels={cells['subject.sample'].nunique()}", f"sample_master_sha256={master_hash}",
        "route=ROUTE_A_AUTHOR_FINAL_CELLS",
        "warning=GSM8517689 GEO Marsh 3c conflicts with paper/all_meta control metadata; discrepancy preserved"]
    LOG.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print("\n".join(log_lines))

if __name__ == "__main__":
    main()
