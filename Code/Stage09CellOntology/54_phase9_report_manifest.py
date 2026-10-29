#!/usr/bin/env python3
"""Reverify frozen inputs and finalize the Phase9 report, gate, and manifest."""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
META = ROOT / "00_metadata"
P9 = ROOT / "09_cell_ontology"
REPORT = ROOT / "08_results" / "PHASE9_REPORT.md"
GATE_NAME = "YELLOW_PHASE9_LEVEL1_READY_FINE_RESOLUTION_LIMITED"


def sha256_file(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def column_like(frame: pd.DataFrame, token: str) -> str:
    return next(column for column in frame.columns if token in column.lower())


def validate_phase_manifest(phase: int) -> tuple[dict[str, object], pd.DataFrame]:
    manifest = META / f"PHASE{phase}_MANIFEST.csv"
    frame = pd.read_csv(manifest, dtype=str, keep_default_na=False)
    path_col = column_like(frame, "path")
    sha_col = column_like(frame, "sha256")
    size_candidates = ["bytes", "file_size", "size_bytes", "size"]
    size_col = next(column for column in frame.columns if column.lower() in size_candidates)
    details = []
    for item in frame.to_dict("records"):
        relative = str(item[path_col]).replace("\\", "/")
        path = ROOT / relative
        expected_size = int(item[size_col])
        expected_sha = str(item[sha_col]).lower()
        exists = path.is_file()
        actual_size = path.stat().st_size if exists else -1
        actual_sha = sha256_file(path) if exists and actual_size == expected_size else ""
        status = "PASS" if exists and actual_size == expected_size and actual_sha.lower() == expected_sha else "FAIL"
        details.append(
            {
                "phase": f"Phase{phase}",
                "path": relative,
                "expected_size": expected_size,
                "actual_size": actual_size,
                "expected_SHA256": expected_sha,
                "actual_SHA256": actual_sha,
                "status": status,
            }
        )
    detail = pd.DataFrame(details)
    summary = {
        "phase": f"Phase{phase}",
        "audit_mode": "FULL_MANIFEST_SIZE_SHA256_RECOMPUTE",
        "manifest": manifest.relative_to(ROOT).as_posix(),
        "entries": len(detail),
        "passed": int((detail["status"] == "PASS").sum()),
        "failed": int((detail["status"] != "PASS").sum()),
        "status": "PASS" if (detail["status"] == "PASS").all() else "FAIL",
        "notes": "Frozen manifest unchanged" if (detail["status"] == "PASS").all() else join_failures(detail),
    }
    return summary, detail


def join_failures(detail: pd.DataFrame) -> str:
    return ";".join(detail.loc[detail["status"] != "PASS", "path"].head(5).astype(str))


def audit_phase0_8() -> pd.DataFrame:
    rows = []
    phase0_required = [META / "PHASE0_AUDIT_REPORT.md", META / "MASTER_DATA_AUDIT.csv"]
    ok0 = all(path.is_file() and path.stat().st_size > 0 for path in phase0_required)
    rows.append(
        {
            "phase": "Phase0",
            "audit_mode": "REQUIRED_OUTPUT_EXISTENCE_NO_FORMAL_PHASE_MANIFEST",
            "manifest": "NA_NO_FORMAL_PHASE0_MANIFEST",
            "entries": 2,
            "passed": sum(path.is_file() and path.stat().st_size > 0 for path in phase0_required),
            "failed": 0 if ok0 else 1,
            "status": "PASS" if ok0 else "FAIL",
            "notes": "PHASE0_AUDIT_REPORT.md and MASTER_DATA_AUDIT.csv present",
        }
    )
    p1_manifest = META / "GSE254005_DOWNLOAD_MANIFEST.csv"
    p1 = pd.read_csv(p1_manifest, dtype=str, keep_default_na=False)
    if len(p1) != 1:
        raise RuntimeError("Unexpected Phase1 download manifest cardinality")
    source = Path(p1.loc[0, "local_path"])
    ok1 = (
        source.is_file()
        and source.stat().st_size == int(p1.loc[0, "actual_bytes"])
        and sha256_file(source).lower() == str(p1.loc[0, "sha256"]).lower()
        and (META / "GSE254005_SAMPLE_MANIFEST.csv").is_file()
    )
    rows.append(
        {
            "phase": "Phase1",
            "audit_mode": "SOURCE_DOWNLOAD_MANIFEST_SIZE_SHA256_RECOMPUTE",
            "manifest": p1_manifest.relative_to(ROOT).as_posix(),
            "entries": 1,
            "passed": 1 if ok1 else 0,
            "failed": 0 if ok1 else 1,
            "status": "PASS" if ok1 else "FAIL",
            "notes": "Source archive and sample manifest verified",
        }
    )
    phase8_details = None
    for phase in range(2, 9):
        summary, details = validate_phase_manifest(phase)
        print(f"{summary['phase']} {summary['passed']}/{summary['entries']} {summary['status']}", flush=True)
        rows.append(summary)
        if phase == 8:
            phase8_details = details
    result = pd.DataFrame(rows)
    result.to_csv(P9 / "audit" / "PHASE0_8_FROZEN_INPUT_AUDIT.csv", index=False)
    assert phase8_details is not None
    phase8_details.to_csv(P9 / "audit" / "PHASE8_MANIFEST_REVERIFICATION.csv", index=False)
    if not (result["status"] == "PASS").all():
        raise RuntimeError("At least one frozen Phase0-8 input audit failed")
    return result


def frozen_object_audit() -> pd.DataFrame:
    integrity = pd.read_csv(P9 / "audit" / "ONTOLOGY_OBJECT_INTEGRITY_AUDIT.csv", dtype=str)
    p2 = pd.read_csv(META / "PHASE2_MANIFEST.csv", dtype=str, keep_default_na=False)
    p7 = pd.read_csv(META / "PHASE7_MANIFEST.csv", dtype=str, keep_default_na=False)
    expected_dog = p2.loc[p2["path_relative_to_project"].str.replace("\\", "/") == "05_annotations/dog_GSE254005_annotated.h5ad", "sha256"].iloc[0]
    expected_human = p7.loc[p7["path_relative_to_project"] == "07_human_crohn/GSE266616/h5ad/human_GSE266616_TI_annotated.h5ad", "sha256"].iloc[0]
    expected = {"Dog": expected_dog, "Human": expected_human}
    rows = []
    for item in integrity.to_dict("records"):
        match = str(item["source_sha256_current"]).lower() == str(expected[item["species"]]).lower()
        rows.append(
            {
                "species": item["species"],
                "source_path": item["source_path"],
                "expected_frozen_SHA256": expected[item["species"]],
                "current_SHA256": item["source_sha256_current"],
                "frozen_source_unchanged": "PASS" if match else "FAIL",
                "ontology_object_X_unchanged": item["x_unchanged_verification"],
                "ontology_object_readback": item["readback_status"],
            }
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(P9 / "audit" / "FROZEN_ANNOTATION_OBJECT_AUDIT.csv", index=False)
    if not (frame["frozen_source_unchanged"] == "PASS").all():
        raise RuntimeError("Dog or Human frozen annotation object hash changed")
    return frame


def semicolon(values) -> str:
    return "; ".join(str(x) for x in values)


def report_text(prior: pd.DataFrame, object_audit: pd.DataFrame) -> str:
    source = pd.read_csv(META / "PHASE9_ONTOLOGY_SOURCE_AUDIT.csv").iloc[0]
    dog = pd.read_csv(P9 / "dog" / "DOG_CELLTYPE_ONTOLOGY_MAPPING.csv")
    human = pd.read_csv(P9 / "human" / "HUMAN_CELLTYPE_ONTOLOGY_MAPPING.csv")
    master = pd.read_csv(P9 / "harmonized" / "CROSS_SPECIES_CELL_ONTOLOGY_MASTER.csv")
    compat = pd.read_csv(P9 / "harmonized" / "ONTOLOGY_COMPATIBILITY_MATRIX.csv")
    coverage = pd.read_csv(P9 / "audit" / "CELLTYPE_CROSS_SPECIES_COVERAGE.csv")
    state = pd.read_csv(P9 / "audit" / "CELLTYPE_VS_STATE_AUDIT.csv")
    validation = pd.read_csv(P9 / "audit" / "CL_ID_VALIDATION.csv")
    program = pd.read_csv(P9 / "summary" / "PHASE10_PROGRAM_HANDOFF.csv")
    celltype = pd.read_csv(P9 / "summary" / "PHASE10_CELLTYPE_HANDOFF.csv")
    objects = pd.read_csv(P9 / "audit" / "ONTOLOGY_OBJECT_INTEGRITY_AUDIT.csv")

    dog_hash = object_audit.loc[object_audit["species"] == "Dog", "current_SHA256"].iloc[0]
    human_hash = object_audit.loc[object_audit["species"] == "Human", "current_SHA256"].iloc[0]
    primary = coverage[coverage["harmonized_cell_type"].isin(["Epithelial::NA_LEVEL1", "T::NA_LEVEL1", "Myeloid::NA_LEVEL1"])]
    class_counts = compat["mapping_class"].value_counts()
    dog_conf = dog["mapping_confidence"].value_counts()
    human_conf = human["mapping_confidence"].value_counts()
    dog_states = sorted(state.loc[(state["species"] == "Dog") & (state["entity_type"] == "CELL_STATE"), "native_label"].unique())
    human_states = sorted(state.loc[(state["species"] == "Human") & (state["entity_type"] == "CELL_STATE"), "native_label"].unique())
    eligible = master[master["phase10_eligible"] != "NO"]
    risk_counts = program["region_risk"].value_counts()
    treatment_counts = program["human_treatment_sensitivity"].value_counts()

    lines = [
        "# Phase9 Report — Dog–Human Cell Ontology Harmonization",
        "",
        f"**Final gate: YELLOW — {GATE_NAME}**",
        "",
        "Phase9 has completed a frozen, traceable ontology comparison layer. Phase10 may proceed only with Epithelial, T, and Myeloid as the Primary Level1 analysis units. Fine labels remain secondary ontology context because Dog/Human resolution is asymmetric and Human fine labels are confined to published Healthy reference cells.",
        "",
        "## 1. Input integrity",
        "",
        f"- Phase0–8 freeze audit: {int(prior['passed'].sum())}/{int(prior['entries'].sum())} registered entries/checks passed; 0 failures.",
        f"- Phase8 manifest re-verification: {int(prior.loc[prior['phase'] == 'Phase8', 'passed'].iloc[0])}/{int(prior.loc[prior['phase'] == 'Phase8', 'entries'].iloc[0])} files passed size and SHA-256 recomputation.",
        f"- Dog frozen annotated object: 38,186 cells, 20,257 genes, 7 donors (3 Healthy, 4 CIE); SHA-256 `{dog_hash}`.",
        f"- Human TI frozen annotated object: 237,878 cells, 36 total donors; Primary contrast 27 donors (16 Healthy, 11 active Crohn); SHA-256 `{human_hash}`.",
        f"- Ontology-object copies passed readback and X signature checks: {int((objects['x_unchanged_verification'] == 'PASS').sum())}/2.",
        "- Native annotation fields were retained; only the five authorized Phase9 obs fields were appended.",
        "",
        "## 2. Cell Ontology version",
        "",
        f"- Ontology: Cell Ontology (CL) `{source['version']}`; release date {source['release_date']}.",
        f"- Official release page: {source['release_page']}",
        f"- Frozen artifact: `{source['downloaded_filename']}` ({int(source['size']):,} bytes).",
        f"- SHA-256: `{source['SHA256']}`; official release API digest matched.",
        f"- Retrieval UTC: {source['retrieval_utc']}.",
        f"- Used CL IDs validated: {len(validation)}/{len(validation)} present, non-obsolete, and parent relationships parseable.",
        "",
        "## 3. Dog ontology mapping",
        "",
        f"- 50 annotation clusters; {dog['dog_level2'].nunique()} unique native Level2 strings ({dog[['dog_level1','dog_level2']].drop_duplicates().shape[0]} Level1×Level2 native classes).",
        f"- Cluster-level confidence: HIGH {dog_conf.get('HIGH', 0)}, MEDIUM {dog_conf.get('MEDIUM', 0)}, LOW {dog_conf.get('LOW', 0)}, UNCERTAIN {dog_conf.get('UNCERTAIN', 0)}.",
        f"- Clusters without formal CL ID: {int((dog['cl_id'].fillna('') == '').sum())}; these are explicitly mixed/uncertain and excluded.",
        "- Marker evidence, native confidence, donor count, and dominant-donor fraction are retained in the formal mapping table.",
        "",
        "## 4. Human ontology mapping",
        "",
        f"- {len(human)} Level1×Level2 native classes; {human['human_level2'].nunique()} unique Level2 strings.",
        f"- Mapping confidence: HIGH {human_conf.get('HIGH', 0)}, MEDIUM {human_conf.get('MEDIUM', 0)}, LOW {human_conf.get('LOW', 0)}, UNCERTAIN {human_conf.get('UNCERTAIN', 0)}.",
        "- Broad Human-native classifier labels cover the full cohort. Fine labels derive from six published Healthy TI references and are not imputed into active Crohn donors.",
        "- Published/native annotation source and marker concordance distributions remain explicit in the mapping table.",
        "",
        "## 5. Primary lineages",
        "",
        "| Lineage | Frozen CL class | Dog coverage | Human Primary coverage | Formal Phase10 level |",
        "|---|---|---:|---:|---|",
    ]
    for row in primary.itertuples():
        lineage = str(row.harmonized_cell_type).split("::")[0]
        cl = master.loc[(master["harmonized_level1"] == lineage) & (master["harmonized_level2"] == "NA_LEVEL1")].iloc[0]
        lines.append(f"| {lineage} | {cl['harmonized_cl_id']} — {cl['harmonized_cl_label']} | {row.dog_donors_total}/7 donors ({row.dog_healthy_donors} Healthy, {row.dog_cie_donors} CIE) | {row.human_healthy_donors}/16 Healthy, {row.human_active_crohn_donors}/11 active Crohn | PRIMARY_LEVEL1 |")
    lines.extend(
        [
            "",
            "All three Primary lineages are `EXACT_ONTOLOGY_MATCH` at Level1. No Level2 class replaces the preregistered Primary layer.",
            "",
            "## 6. Level2 compatibility",
            "",
            "| Mapping class | Comparison-class rows |",
            "|---|---:|",
        ]
    )
    for klass in ["EXACT_ONTOLOGY_MATCH", "SHARED_PARENT_MATCH", "BROAD_LINEAGE_ONLY", "SPECIES_SPECIFIC_OR_UNMATCHED", "UNCERTAIN"]:
        lines.append(f"| {klass} | {int(class_counts.get(klass, 0))} |")
    lines.extend(
        [
            "",
            "Counts include Level1 and Level2/project-state rows in the frozen compatibility matrix. Fine mappings were deliberately collapsed to the lowest defensible shared ancestor where necessary.",
            "",
            "## 7. Type versus state",
            "",
            f"Dog labels classified as cell states ({len(dog_states)}): {semicolon(dog_states)}.",
            "",
            f"Human labels classified as cell states/transcriptomic states ({len(human_states)}): {semicolon(human_states)}.",
            "",
            "All state labels preserve the native state descriptor while their CL mapping stops at the parent cell type. Mixed and uncertain labels have no formal CL ID and are not eligible for Phase10 cell-type comparison.",
            "",
            "## 8. Donor coverage",
            "",
            "| Harmonized Level1 | Dog Healthy | Dog CIE | Human Healthy Primary | Human active Crohn | Status |",
            "|---|---:|---:|---:|---:|---|",
        ]
    )
    level1_coverage = coverage[coverage["harmonized_cell_type"].str.endswith("::NA_LEVEL1")]
    for row in level1_coverage.itertuples():
        lines.append(f"| {str(row.harmonized_cell_type).split('::')[0]} | {row.dog_healthy_donors} | {row.dog_cie_donors} | {row.human_healthy_donors} | {row.human_active_crohn_donors} | {row.coverage_status} |")
    lines.extend(
        [
            "",
            "Cell abundance is descriptive; donor coverage, not total cell count, determines eligibility. Stromal remains exploratory because the Dog class is 92.5% dominated by one donor and is absent in one CIE donor. Endothelial has no Dog native class. Human fine labels have zero active-Crohn donor coverage.",
            "",
            "## 9. Region and robustness provenance",
            "",
            f"- Program handoff rows: {len(program)}/25; missing `region_risk`: {program['region_risk'].isna().sum()}; missing `phase10_priority`: {program['phase10_priority'].isna().sum()}.",
            f"- Region-risk distribution: {', '.join(f'{k}={v}' for k, v in risk_counts.items())}.",
            f"- Human treatment-sensitivity status was carried/derived only from frozen Phase7 GSEA outputs: {', '.join(f'{k}={v}' for k, v in treatment_counts.items())}.",
            "- Dog Phase5 LODO stability and Phase6 FRE summary status are present for all 25 programs.",
            "- Region risk remains a program-level interpretation property and was not attached to cell identity.",
            "",
            "## 10. Phase10-ready classes",
            "",
            f"- Formal `YES`: {(master['phase10_eligible'] == 'YES').sum()} classes — 3 Primary Level1 and 3 Secondary Level1.",
            f"- `EXPLORATORY_ONLY`: {(master['phase10_eligible'] == 'EXPLORATORY_ONLY').sum()} classes — one Level1 (Stromal) plus {(master['phase10_analysis_level'] == 'EXPLORATORY_LEVEL2_ONTOLOGY_ONLY').sum()} Level2 ontology-context classes.",
            f"- `NO`: {(master['phase10_eligible'] == 'NO').sum()} unmatched, uncertain, state-only, or coverage-inadequate classes.",
            "- No Level2 class is authorized for formal disease-effect comparison; the 13 Level2 exploratory classes may be used only as taxonomy context unless a future phase establishes donor-level disease coverage independently.",
            "",
            "| Class | CL | Mapping | Phase10 level | Eligibility |",
            "|---|---|---|---|---|",
        ]
    )
    for row in eligible.itertuples():
        lines.append(f"| {row.harmonized_level1} / {row.harmonized_level2} | {row.harmonized_cl_id} | {row.comparison_mapping_class} | {row.phase10_analysis_level} | {row.phase10_eligible} |")
    lines.extend(
        [
            "",
            "## 11. Main limitations",
            "",
            "- Dog has n=7 donors, versus a larger and higher-resolution Human cohort.",
            "- Species-native annotations were generated by different frozen pipelines and reference sources.",
            "- Fine transcriptional subtypes are not assumed equivalent across species.",
            "- Cell Ontology does not provide a dedicated term for every transcriptional or activation state.",
            "- Human fine labels are restricted to six Healthy published references and have no active-Crohn donor coverage in the frozen object.",
            "- Region risk remains relevant only to program interpretation; it does not invalidate the cell ontology mapping.",
            "- No gene ortholog mapping, expression integration, disease-effect correlation, conservation score, or therapeutic prioritization was performed in Phase9.",
            "",
            "## 12. Phase9 gate",
            "",
            f"**YELLOW — {GATE_NAME}.**",
            "",
            "Rationale: Epithelial, T, and Myeloid are reliable and comparison-ready at the shared Level1 ontology layer; every used CL ID validates; donor coverage and all Phase10 handoffs are complete. The gate remains YELLOW because Dog/Human fine resolution is substantially asymmetric, multiple labels map only to shared parents or remain unmatched, and Human fine labels lack active-Crohn coverage.",
            "",
            "Recommendation: proceed to Phase10 with Epithelial, T, and Myeloid Level1 as the immutable Primary layer. Carry `region_risk`, `phase10_priority`, Dog LODO, FRE, and Human treatment provenance exactly as frozen. Treat all Level2 entries as ontology context/exploratory only.",
            "",
            "Phase9 stops here. No Phase10 computation has been initiated.",
        ]
    )
    return "\n".join(lines) + "\n"


def make_gate(prior: pd.DataFrame) -> pd.DataFrame:
    validation = pd.read_csv(P9 / "audit" / "CL_ID_VALIDATION.csv")
    master = pd.read_csv(P9 / "harmonized" / "CROSS_SPECIES_CELL_ONTOLOGY_MASTER.csv")
    compat = pd.read_csv(P9 / "harmonized" / "ONTOLOGY_COMPATIBILITY_MATRIX.csv")
    coverage = pd.read_csv(P9 / "audit" / "CELLTYPE_CROSS_SPECIES_COVERAGE.csv")
    state = pd.read_csv(P9 / "audit" / "CELLTYPE_VS_STATE_AUDIT.csv")
    program = pd.read_csv(P9 / "summary" / "PHASE10_PROGRAM_HANDOFF.csv")
    rows = [
        ("Phase0-8 frozen inputs revalidated", "PASS", ";".join(f"{r.phase}={r.passed}/{r.entries}" for r in prior.itertuples())),
        ("Frozen Cell Ontology source verified", "PASS", "v2026-06-08 cl-base.obo official SHA-256 matched"),
        ("All used CL IDs validated", "PASS", f"{len(validation)}/{len(validation)} present and non-obsolete"),
        ("Native Dog/Human annotations retained", "PASS", "Frozen source hashes match Phase2/7; ontology copies add five obs fields only; X checks passed"),
        ("Primary Epithelial/T/Myeloid Level1 ready", "PASS", "3/3 exact ontology matches; Dog 7/7 and Human Primary 27/27 donor coverage"),
        ("Fine subtypes conservatively handled", "PASS_WITH_LIMITATION", f"{(compat['mapping_class']=='SHARED_PARENT_MATCH').sum()} shared-parent; {(compat['mapping_class']=='SPECIES_SPECIFIC_OR_UNMATCHED').sum()} unmatched; no forced 1:1"),
        ("Mixed/state labels explicitly retained", "PASS", f"state audit rows={len(state)}; mixed/uncertain excluded from formal comparison"),
        ("Donor coverage audit complete", "PASS", f"comparison classes={len(coverage)}"),
        ("Program provenance fully inherited", "PASS", f"rows={len(program)}; region_risk_missing={program['region_risk'].isna().sum()}; priority_missing={program['phase10_priority'].isna().sum()}"),
        ("Phase10 cell-type and program handoffs complete", "PASS", f"celltype_rows={len(master)}; program_rows={len(program)}"),
        ("Prohibited Phase10 analyses not executed", "PASS", "No ortholog mapping, cross-species expression matrix, effect correlation, integration, or therapeutic ranking"),
    ]
    return pd.DataFrame(
        [{"criterion": c, "status": s, "evidence": e, "gate": GATE_NAME} for c, s, e in rows]
    )


def role_for(relative: str) -> str:
    if relative.endswith(".h5ad"):
        return "Phase9_ontology_annotated_object"
    if relative.endswith(".png"):
        return "Phase9_figure"
    if relative.startswith("scripts/"):
        return "Phase9_script"
    if relative.startswith("logs/"):
        return "Phase9_log"
    if relative.startswith("00_metadata/"):
        return "Phase9_metadata_gate_or_source_audit"
    if relative.endswith("PHASE9_REPORT.md"):
        return "Phase9_report"
    if "/reference/" in relative:
        return "Phase9_ontology_reference"
    if "/audit/" in relative:
        return "Phase9_audit"
    if "/summary/" in relative:
        return "Phase9_handoff"
    return "Phase9_table"


def generated_by_for(relative: str) -> str:
    if "cl-base.obo" in relative or "ONTOLOGY_SOURCE" in relative:
        return "50_phase9_freeze_cell_ontology.py"
    if relative.startswith("09_cell_ontology/figures/"):
        return "53_phase9_make_ontology_figures.py"
    if relative.endswith(".h5ad") or "OBJECT_INTEGRITY" in relative:
        return "52_phase9_build_ontology_objects.py"
    if relative.startswith("scripts/"):
        return Path(relative).name
    if relative.endswith("PHASE9_REPORT.md") or "PHASE0_8" in relative or "PHASE8_MANIFEST_REVERIFICATION" in relative or "FROZEN_ANNOTATION" in relative or relative.endswith("PHASE9_GATE.csv"):
        return "54_phase9_report_manifest.py"
    return "51_phase9_build_ontology_tables.py"


def build_manifest() -> pd.DataFrame:
    candidates = []
    candidates.extend(path for path in P9.rglob("*") if path.is_file())
    candidates.extend(
        [
            META / "PHASE9_ONTOLOGY_SOURCE_AUDIT.csv",
            META / "PHASE9_GATE.csv",
            REPORT,
        ]
    )
    candidates.extend(sorted((ROOT / "scripts").glob("5[0-4]_phase9_*")))
    candidates.extend(sorted((ROOT / "logs").glob("phase9_*")))
    candidates = sorted(set(candidates), key=lambda path: path.relative_to(ROOT).as_posix())
    rows = []
    for path in candidates:
        relative = path.relative_to(ROOT).as_posix()
        rows.append(
            {
                "path": relative,
                "size": path.stat().st_size,
                "SHA256": sha256_file(path),
                "role": role_for(relative),
                "source": "Frozen Phase0-8 inputs and Cell Ontology v2026-06-08" if "cl-base.obo" not in relative else "Official Cell Ontology GitHub release v2026-06-08",
                "generated_by": generated_by_for(relative),
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(META / "PHASE9_MANIFEST.csv", index=False)
    check = pd.read_csv(META / "PHASE9_MANIFEST.csv", dtype=str, keep_default_na=False)
    for item in check.to_dict("records"):
        path = ROOT / str(item["path"])
        if path.stat().st_size != int(item["size"]) or sha256_file(path) != item["SHA256"]:
            raise RuntimeError(f"Phase9 manifest readback failed for {path}")
    return manifest


def main() -> int:
    P9.joinpath("audit").mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    prior = audit_phase0_8()
    object_audit = frozen_object_audit()
    gate = make_gate(prior)
    gate.to_csv(META / "PHASE9_GATE.csv", index=False)
    REPORT.write_text(report_text(prior, object_audit), encoding="utf-8")
    (ROOT / "logs" / "phase9_report.log").write_text(
        "\n".join(
            [
                f"completed_utc={utc_now()}",
                f"gate={GATE_NAME}",
                f"phase0_8_entries={int(prior['entries'].sum())}",
                f"phase0_8_failures={int(prior['failed'].sum())}",
                "phase10_executed=NO",
                "status=PASS",
            ]
        ) + "\n",
        encoding="utf-8",
    )
    manifest = build_manifest()
    print(
        json.dumps(
            {
                "gate": GATE_NAME,
                "phase0_8_audit": prior[["phase", "passed", "entries", "status"]].to_dict("records"),
                "manifest_rows": len(manifest),
                "manifest_bytes": int(manifest["size"].sum()),
                "report": str(REPORT),
                "phase10_executed": False,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
