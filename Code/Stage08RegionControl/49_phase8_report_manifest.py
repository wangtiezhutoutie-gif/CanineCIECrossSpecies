from __future__ import annotations

import csv
import hashlib
import os
import platform
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import anndata
import joblib
import matplotlib
import numpy as np
import pandas as pd
import scipy
import sklearn


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P8 = ROOT / "08_region_control"
META = P8 / "metadata"
SUMMARY = P8 / "summary"
RESULTS = ROOT / "08_results"
PROJECT_META = ROOT / "00_metadata"
REPORT = RESULTS / "PHASE8_REPORT.md"
GATE = PROJECT_META / "PHASE8_GATE.csv"
LOCAL_GATE = META / "PHASE8_GATE.csv"
MANIFEST = PROJECT_META / "PHASE8_MANIFEST.csv"
MANIFEST_VERIFICATION = PROJECT_META / "PHASE8_MANIFEST_VERIFICATION.csv"
SOFTWARE = PROJECT_META / "PHASE8_SOFTWARE_VERSIONS.csv"
REPORT_LOGS = (ROOT / "logs" / "phase8_report.log", P8 / "logs" / "phase8_report.log")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    temp = Path(name)
    try:
        frame.to_csv(temp, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_MINIMAL)
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def atomic_text(value: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    temp = Path(name)
    try:
        temp.write_text(value, encoding="utf-8")
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def all_pass(path: Path) -> bool:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    return bool(len(frame) and frame["status"].eq("PASS").all())


def software_versions() -> pd.DataFrame:
    rows = [
        ("Python", platform.python_version(), sys.executable),
        ("pandas", pd.__version__, "Phase8 Python environment"),
        ("numpy", np.__version__, "Phase8 Python environment"),
        ("scipy", scipy.__version__, "Phase8 Python environment"),
        ("anndata", anndata.__version__, "Phase8 Python environment"),
        ("scikit-learn", sklearn.__version__, "Frozen Phase7 classifier runtime"),
        ("joblib", joblib.__version__, "Frozen Phase7 classifier serialization"),
        ("matplotlib", matplotlib.__version__, "Phase8 figures"),
        ("R", "4.4.3", str(ROOT / "envs" / "phase3_r_portable" / "app" / "bin" / "x64" / "Rscript.exe")),
        ("edgeR", "4.4.2", "Portable Phase3/7 R library"),
        ("fgsea", "1.32.4", "Portable Phase3/7 R library"),
    ]
    frame = pd.DataFrame(rows, columns=["software", "version", "source"])
    atomic_csv(frame, SOFTWARE)
    return frame


def build_gate() -> pd.DataFrame:
    prior = pd.read_csv(META / "PHASE8_PRIOR_MANIFEST_AUDIT.csv", dtype=str, keep_default_na=False)
    audit = pd.read_csv(META / "HUMAN_AC_REGION_AUDIT.csv", dtype=str, keep_default_na=False)
    download = pd.read_csv(META / "HUMAN_AC_DOWNLOAD_MANIFEST.csv", dtype=str, keep_default_na=False)
    core = pd.read_csv(SUMMARY / "REGION_PROGRAM_STABILITY.csv", dtype=str, keep_default_na=False)
    matrix = pd.read_csv(SUMMARY / "REGION_INTERPRETATION_MATRIX.csv", dtype=str, keep_default_na=False)
    p8_files = [path for path in P8.rglob("*") if path.is_file()]
    prohibited = [
        str(path.relative_to(ROOT)).replace("\\", "/")
        for path in p8_files
        if path.name in {"GSE266616_RAW.tar", "scRNA_gene_counts.mtx", "Full_cohort.h5ad", "DGE.zip"}
        or path.suffix.lower() in {".fastq", ".sra", ".part", ".tmp"}
    ]
    execution_scopes = (P8, ROOT / "scripts", ROOT / "08_results", PROJECT_META)
    forbidden_execution = [
        str(path.relative_to(ROOT)).replace("\\", "/")
        for scope in execution_scopes
        if scope.exists()
        for path in scope.rglob("*")
        if path.is_file() and any(token in path.name.lower() for token in ("singist", "saturn", "concord", "phase9_ontology", "phase10_ortholog"))
    ]
    checks = [
        ("Phase0-7 frozen inputs revalidated", prior["status"].eq("PASS").all(), ";".join(f"{r.phase}={r.passed}/{r.entries}" for r in prior.itertuples())),
        ("AC metadata capability gate passed", all_pass(META / "PHASE8A_AC_CAPABILITY_GATE.csv"), "37 AC GSM / 31 donors; 16 unambiguous Healthy and 6 active Crohn analysis donors"),
        ("AC processed archive gate passed", len(download) == 22 and download["archive_valid"].eq("YES").all() and pd.to_numeric(download["actual_size"]).eq(pd.to_numeric(download["expected_size"])).all(), f"files={len(download)};bytes={pd.to_numeric(download['actual_size']).sum()}"),
        ("AC QC and frozen annotation gate passed", all_pass(META / "PHASE8A_AC_QC_ANNOTATION_GATE.csv"), "115556 raw; 90903 retained; 1053 Uncertain"),
        ("Donor-level region pseudobulk gate passed", all_pass(META / "PHASE8A_REGION_PSEUDOBULK_GATE.csv"), "9 donor-level lineage matrices"),
        ("edgeR QL region analysis gate passed", all_pass(META / "PHASE8A_REGION_EDGER_GATE.csv"), "9 preregistered fits; full tested genes retained"),
        ("Registered region-program gate passed", all_pass(META / "PHASE8A_REGION_PROGRAM_GATE.csv"), f"core={len(core)};Dog-forward={len(matrix)}"),
        ("Prohibited raw/partial routes absent", not prohibited, "none" if not prohibited else ";".join(prohibited)),
        ("Phase9/10 and cross-species algorithms not executed", not forbidden_execution and not (ROOT / "09_ontology").exists() and not (ROOT / "10_cross_species").exists(), "Stopped after Phase8A"),
    ]
    gate = pd.DataFrame([
        {"criterion": name, "status": "PASS" if okay else "FAIL", "evidence": evidence}
        for name, okay, evidence in checks
    ])
    gate["gate"] = "GREEN_PHASE8A_REGION_CONTROL_COMPLETE" if gate["status"].eq("PASS").all() else "YELLOW_OR_RED_PHASE8_REVIEW_REQUIRED"
    atomic_csv(gate, GATE)
    atomic_csv(gate, LOCAL_GATE)
    if not gate["status"].eq("PASS").all():
        raise RuntimeError("Phase8A final gate failed")
    return gate


def fmt(value: object, digits: int = 3) -> str:
    number = float(value)
    return f"{number:.{digits}f}"


def program_rows(core: pd.DataFrame, programs: list[str]) -> str:
    rows = ["| Lineage | Program | TI disease NES | AC disease NES | AC registered FDR | Status |", "|---|---|---:|---:|---:|---|"]
    for row in core.loc[core["program"].isin(programs)].itertuples(index=False):
        rows.append(f"| {row.lineage} | {row.program} | {float(row.TI_NES):.2f} | {float(row.AC_NES):.2f} | {float(row.AC_registered_FDR):.3g} | {row.region_status} |")
    return "\n".join(rows)


def build_report() -> str:
    audit = pd.read_csv(META / "HUMAN_AC_REGION_AUDIT.csv", dtype=str, keep_default_na=False)
    availability = pd.read_csv(META / "HUMAN_REGION_SAMPLE_AVAILABILITY.csv", dtype=str, keep_default_na=False)
    pairing = pd.read_csv(META / "HUMAN_REGION_PAIRING_AUDIT.csv", dtype=str, keep_default_na=False)
    download = pd.read_csv(META / "HUMAN_AC_DOWNLOAD_MANIFEST.csv")
    qc = pd.read_csv(META / "HUMAN_AC_QC_PER_SAMPLE.csv")
    coverage = pd.read_csv(META / "HUMAN_AC_LINEAGE_DONOR_COVERAGE.csv")
    de = pd.read_csv(SUMMARY / "HUMAN_REGION_DE_SUMMARY.csv")
    core = pd.read_csv(SUMMARY / "REGION_PROGRAM_STABILITY.csv")
    matrix = pd.read_csv(SUMMARY / "REGION_INTERPRETATION_MATRIX.csv")
    core_counts = core["region_status"].value_counts()
    risk_counts = matrix["region_risk"].value_counts()
    active_ac = audit.loc[audit["eligible"].eq("YES_ACTIVE_REGION_AND_AC_DISEASE")]
    treatment = active_ac["treatment"].value_counts()

    de_lines = ["| Analysis | Lineage | Design / n | Tested genes | FDR<0.05 |", "|---|---|---|---:|---:|"]
    for row in de.itertuples(index=False):
        de_lines.append(f"| {row.analysis} | {row.lineage} | `{row.design}`; {row.group_evidence} | {int(row.tested_genes):,} | {int(row.FDR05):,} |")

    report = f"""# Phase8A report — Human region control data audit

Final gate: **GREEN — PHASE8A REGION CONTROL COMPLETE**

Generated: {utc_now()}

## 1. Objective and scope

Phase8A assessed whether frozen Human Crohn TI disease programs are region-robust when secondary GSE266616 ascending-colon data are used. It did not perform ontology harmonization, ortholog mapping or formal Dog–Human conservation analysis. Phase7 TI disease effects remained frozen.

## 2. Frozen inputs and AC metadata

- Phase0–7 source/frozen manifests were revalidated by file size and SHA-256; all available formal manifests passed. Phase0/1, which do not have formal `PHASE*_MANIFEST.csv`, passed their preregistered audit/source-manifest checks.
- GSE266616 contains **37 AC GSM / 31 AC donors**: 17 Healthy donors and 14 Crohn donors. Crohn AC comprises 8 non-inflamed donors plus 6 donors with paired adjacent/active AC biopsies.
- Region-analysis scope after the preregistered identity rule: **16 unambiguous Healthy AC donors and 6 active Crohn AC donors**.
- OR95 was retained in metadata but excluded from every new region estimate because its AC and TI GSM have reciprocally swapped processed filenames. HA52 is Healthy AC-only. This leaves **15 unambiguous Healthy TI–AC pairs**.
- Active Crohn has **11 TI vs 6 AC donors and zero cross-region pairs**.
- Active AC treatment: NO_BIOLOGIC={int(treatment.get('NO_BIOLOGIC', 0))}, ANTI_TNF={int(treatment.get('ANTI_TNF', 0))}, ANTI_IL12_23={int(treatment.get('ANTI_IL12_23', 0))}. Treatment and cohort composition remain confounded with the unpaired active regional contrast.

## 3. Download, QC and annotation

- Downloaded **22/22** required per-GSM AC processed archives, **{int(download['actual_size'].sum()):,} bytes ({download['actual_size'].sum()/1024**3:.2f} GiB)**. Adjacent/non-inflamed AC, OR95, whole RAW, SRA and FASTQ were not downloaded.
- Every archive passed exact size, local SHA-256, tar path-safety and unique 10x matrix/features/barcodes checks.
- Phase7 Route B was reused unchanged: 200–6000 genes inclusive, >1000 UMI, mt%<50%, no new doublet filter.
- Raw cells: **{int(qc['raw_cells'].sum()):,}**; retained: **{int(qc['retained_cells'].sum()):,} ({qc['retained_cells'].sum()/qc['raw_cells'].sum():.1%})**; Uncertain: **{int(qc['uncertain_cells'].sum()):,} ({qc['uncertain_cells'].sum()/qc['retained_cells'].sum():.2%})**.
- The frozen Phase7 Human-native classifier was applied without retraining. All six lineage×condition groups remained analyzable. HA98 had only four Epithelial cells and was excluded from Epithelial pseudobulk only.

## 4. Preregistered donor-level analyses

{os.linesep.join(de_lines)}

Healthy AC–TI is a strong anatomical contrast, particularly in Epithelial cells. The active Crohn AC–TI comparison is exploratory because it is unpaired and treatment/cohort-confounded. AC active-Crohn-vs-Healthy is a secondary disease effect used only for region robustness; it does not replace Phase7 TI Primary.

## 5. Core program stability

Only 36 preregistered lineage×program tests plus frozen Dog-forward representatives were evaluated; no new all-pathway discovery search was performed. AC FDR is BH-adjusted within this registered family.

- `ROBUST_FDR`: {int(core_counts.get('ROBUST_FDR', 0))}/36.
- `ROBUST_DIRECTIONAL_LIMITED`: {int(core_counts.get('ROBUST_DIRECTIONAL_LIMITED', 0))}/36.
- `REGION_SENSITIVE_FDR`: {int(core_counts.get('REGION_SENSITIVE_FDR', 0))}/36.
- `REGION_UNKNOWN_OPPOSITE_WEAK`: {int(core_counts.get('REGION_UNKNOWN_OPPOSITE_WEAK', 0))}/36.

### Immune programs

{program_rows(core, ['TNF_NFkB', 'Cytokine_interleukin', 'IFN_alpha', 'IFN_gamma', 'JAK_STAT'])}

TNF/NF-kB, cytokine/interleukin, IFN and JAK-STAT disease directions were higher in active Crohn in both TI and AC for all three lineages. Their absolute regional baselines can still differ: for example, T TNF/NF-kB is Healthy-AC higher but active-Crohn-AC lower relative to TI, while Myeloid TNF/NF-kB is AC lower in both region contrasts. The defensible statement is disease-direction robustness, not equal anatomical baseline expression.

### Metabolic programs

{program_rows(core, ['Oxidative_phosphorylation', 'Fatty_acid_metabolism', 'Xenobiotic_metabolism'])}

Epithelial oxidative phosphorylation, fatty-acid metabolism and xenobiotic metabolism were Crohn-lower in both TI and AC with registered FDR<0.05. T metabolism was less stable: fatty-acid metabolism reversed with registered FDR support, oxidative phosphorylation reversed weakly, and xenobiotic metabolism stayed directionally positive but was not significant in AC. Myeloid metabolic signals were mostly direction-only or region-unknown and were already weak in TI.

### Stress and proliferation

Epithelial E2F and G2M disease effects reversed with registered FDR support in AC; Epithelial UPR reversed weakly. These programs require explicit region caution. T/Myeloid proliferation was more stable, but several results remain direction-only because active AC n is small.

## 6. Dog-forward region interpretation matrix

Of 25 frozen Dog-forward representatives:

- LOW_REGION_RISK={int(risk_counts.get('LOW_REGION_RISK', 0))}
- MODERATE_DIRECTION_ONLY={int(risk_counts.get('MODERATE_DIRECTION_ONLY', 0))}
- HIGH_REGION_RISK={int(risk_counts.get('HIGH_REGION_RISK', 0))}
- REGION_UNKNOWN={int(risk_counts.get('REGION_UNKNOWN', 0))}

High-risk representatives are Epithelial cell-cycle/proliferation, Epithelial ribosome/RNA biogenesis and Myeloid translation. Region-unknown representatives include T oxidative phosphorylation/mitochondrial respiration/mitochondrial translation and Myeloid respiratory/oxidative metabolism. These are Phase10 caution flags only; Phase8A does not declare or reject cross-species conservation.

## 7. Phase8A Gate and limitations

Final gate: **GREEN** because AC metadata are traceable, three donor-level comparisons were executable in all primary lineages, and the preregistered core programs received explicit region classifications.

Carried limitations: only six active AC donors; zero active cross-region pairs; biologic treatment and donor cohort differ between active TI and AC; the annotation transfer model was trained from published Healthy TI references; OR95 region identity remains technically anomalous; and AC registered FDR is not an all-pathway discovery FDR.

**Recommendation: YES — proceed to Phase9 ontology harmonization, carrying `region_risk` and `Phase10_priority` as mandatory provenance fields.** Do not discard high/unknown programs solely on Phase8; they require cautious interpretation or independent region control in Phase10.

Phase9, Phase10, singIST, SATURN and CONCORD were not executed.

## 8. Principal deliverables

- `08_region_control/PHASE8_PROTOCOL.md`
- `08_region_control/metadata/HUMAN_AC_REGION_AUDIT.csv`
- `08_region_control/metadata/HUMAN_AC_DOWNLOAD_MANIFEST.csv`
- `08_region_control/summary/HUMAN_REGION_DE_SUMMARY.csv`
- `08_region_control/summary/REGION_PROGRAM_STABILITY.csv`
- `08_region_control/summary/REGION_INTERPRETATION_MATRIX.csv`
- `08_region_control/figures/PHASE8_REGION_CORE_PROGRAM_HEATMAP.png`
- `00_metadata/PHASE8_GATE.csv`
- `00_metadata/PHASE8_MANIFEST.csv`
"""
    return report


def classify(path: Path) -> tuple[str, str, str]:
    relative = str(path.relative_to(ROOT)).replace("\\", "/")
    if relative.startswith("08_region_control/scripts/"):
        return "Phase8_script", "project-authored", path.name
    if "/raw_processed/" in relative:
        return "GSE266616_AC_processed_archive", "NCBI GEO per-GSM processed 10x", "44_phase8_download_ac"
    if relative.endswith(".h5ad"):
        return "Phase8_AC_per_sample_object", "GSE266616 AC processed counts", "45_phase8_ac_qc_annotation"
    if "/pseudobulk/matrices/" in relative:
        return "Phase8_region_pseudobulk_matrix", "Frozen TI and Phase8 AC counts", "46_phase8_region_pseudobulk"
    if "/de_results/" in relative or "/ranked_lists/" in relative:
        return "Phase8_region_edgeR_result", "Donor-level pseudobulk", "47_phase8_region_edger"
    if "/gsea/" in relative:
        return "Phase8_registered_GSEA_result", "Preregistered Phase7 Human-native pathways", "48_phase8_region_programs"
    if relative.endswith(".png"):
        return "Phase8_figure", "Phase8 region analysis", "Phase8 visualization"
    if "/metadata/" in relative or relative.startswith("00_metadata/"):
        return "Phase8_metadata_gate_or_manifest", "Phase7 frozen inputs / Phase8 audit", "Phase8 audit/report"
    if "/summary/" in relative or relative.startswith("08_results/"):
        return "Phase8_summary_or_report", "Phase8 frozen results", "Phase8 report"
    if "/logs/" in relative or relative.startswith("logs/"):
        return "Phase8_log", "Phase8 runtime", "Execution logging"
    return "Phase8_output", "Phase8 analysis", "Phase8"


def manifest_files() -> list[Path]:
    files = {path for path in P8.rglob("*") if path.is_file()}
    files.update(path for path in RESULTS.glob("PHASE8*") if path.is_file())
    files.update(path for path in (ROOT / "logs").glob("phase8*") if path.is_file())
    files.update({GATE, SOFTWARE})
    return sorted(
        path for path in files
        if path.is_file()
        and path not in {MANIFEST, MANIFEST_VERIFICATION}
        and not path.name.endswith((".tmp", ".part", ".pyc"))
        and "__pycache__" not in path.parts
    )


def build_manifest() -> pd.DataFrame:
    files = manifest_files()
    rows: list[dict[str, object]] = []
    for index, path in enumerate(files, start=1):
        role, source, generated_by = classify(path)
        rows.append({
            "path": str(path.relative_to(ROOT)).replace("\\", "/"),
            "size": path.stat().st_size,
            "SHA256": sha256(path),
            "role": role,
            "source": source,
            "generated_by": generated_by,
        })
        if index % 50 == 0:
            print(f"MANIFEST_HASH {index}/{len(files)}", flush=True)
    frame = pd.DataFrame(rows)
    atomic_csv(frame, MANIFEST)
    verification: list[dict[str, object]] = []
    for row in frame.itertuples(index=False):
        path = ROOT / row.path
        valid = path.is_file() and path.stat().st_size == int(row.size) and sha256(path).lower() == row.SHA256.lower()
        verification.append({"path": row.path, "status": "PASS" if valid else "FAIL"})
    verified = pd.DataFrame(verification)
    atomic_csv(verified, MANIFEST_VERIFICATION)
    if not verified["status"].eq("PASS").all():
        raise RuntimeError("Phase8 manifest verification failed")
    role, source, generated_by = classify(MANIFEST_VERIFICATION)
    verification_row = pd.DataFrame([{
        "path": str(MANIFEST_VERIFICATION.relative_to(ROOT)).replace("\\", "/"),
        "size": MANIFEST_VERIFICATION.stat().st_size,
        "SHA256": sha256(MANIFEST_VERIFICATION),
        "role": role,
        "source": source,
        "generated_by": generated_by,
    }])
    frame = pd.concat([frame, verification_row], ignore_index=True).sort_values("path").reset_index(drop=True)
    atomic_csv(frame, MANIFEST)
    return frame


def main() -> int:
    software_versions()
    gate = build_gate()
    atomic_text(build_report(), REPORT)
    for path in REPORT_LOGS:
        atomic_text(
            f"[{utc_now()}] Phase8A report/gate generated; gate={gate['gate'].iloc[0]}; criteria={len(gate)}\n",
            path,
        )
    manifest = build_manifest()
    print(
        f"PHASE8A_COMPLETE gate={gate['gate'].iloc[0]} manifest_entries={len(manifest)} "
        f"manifest_bytes={int(pd.to_numeric(manifest['size']).sum())}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
