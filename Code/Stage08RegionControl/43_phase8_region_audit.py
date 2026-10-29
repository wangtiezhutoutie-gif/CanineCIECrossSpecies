from __future__ import annotations

import hashlib
import os
import re
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
META = ROOT / "00_metadata"
P8 = ROOT / "08_region_control"
P8_META = P8 / "metadata"
P8_SUMMARY = P8 / "summary"
P8_LOGS = P8 / "logs"
PROJECT_LOG = ROOT / "logs" / "phase8_region_audit.log"
PROTOCOL = P8 / "PHASE8_PROTOCOL.md"

SAMPLE_MASTER = META / "GSE266616_HUMAN_SAMPLE_MASTER.csv"
PHASE7_GSEA = ROOT / "07_human_crohn" / "GSE266616" / "gsea"
DOG_FORWARD = PHASE7_GSEA / "DOG_FORWARD_IN_HUMAN_OBSERVATION.csv"


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
        frame.to_csv(temp, index=False, encoding="utf-8-sig")
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


def log(message: str) -> None:
    line = f"[{utc_now()}] {message}"
    print(line, flush=True)
    for path in (PROJECT_LOG, P8_LOGS / "phase8_region_audit.log"):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(line + "\n")


def validate_manifest(phase: int) -> dict[str, object]:
    manifest = META / f"PHASE{phase}_MANIFEST.csv"
    frame = pd.read_csv(manifest, dtype=str, keep_default_na=False)
    path_col = next(column for column in frame.columns if "path" in column.lower())
    sha_col = next(column for column in frame.columns if "sha256" in column.lower())
    size_col = next(column for column in frame.columns if column in {"bytes", "file_size", "size_bytes"})
    passed = 0
    failures: list[str] = []
    for row in frame.to_dict("records"):
        path = ROOT / str(row[path_col]).replace("\\", "/")
        valid = path.is_file() and path.stat().st_size == int(row[size_col])
        if valid:
            valid = sha256(path).lower() == str(row[sha_col]).lower()
        if valid:
            passed += 1
        else:
            failures.append(str(path))
    return {
        "phase": f"Phase{phase}",
        "audit_mode": "FULL_MANIFEST_SIZE_SHA256_RECOMPUTE",
        "manifest": str(manifest),
        "entries": len(frame),
        "passed": passed,
        "failed": len(failures),
        "status": "PASS" if not failures else "FAIL",
        "notes": ";".join(failures[:5]) if failures else "Frozen manifest unchanged",
    }


def audit_prior_freeze() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    phase0_required = [META / "PHASE0_AUDIT_REPORT.md", META / "MASTER_DATA_AUDIT.csv"]
    phase0_ok = all(path.is_file() and path.stat().st_size > 0 for path in phase0_required)
    rows.append({
        "phase": "Phase0",
        "audit_mode": "REQUIRED_OUTPUT_EXISTENCE_NO_FORMAL_PHASE_MANIFEST",
        "manifest": "NA_NO_FORMAL_PHASE0_MANIFEST",
        "entries": len(phase0_required),
        "passed": sum(path.is_file() and path.stat().st_size > 0 for path in phase0_required),
        "failed": 0 if phase0_ok else 1,
        "status": "PASS" if phase0_ok else "FAIL",
        "notes": "PHASE0_AUDIT_REPORT.md and MASTER_DATA_AUDIT.csv present",
    })
    phase1_manifest = META / "GSE254005_DOWNLOAD_MANIFEST.csv"
    phase1 = pd.read_csv(phase1_manifest, dtype=str, keep_default_na=False)
    if len(phase1) != 1:
        raise RuntimeError("Unexpected Phase1 GSE254005 download manifest cardinality")
    source = Path(phase1.loc[0, "local_path"])
    phase1_ok = (
        source.is_file()
        and source.stat().st_size == int(phase1.loc[0, "actual_bytes"])
        and sha256(source).lower() == phase1.loc[0, "sha256"].lower()
        and (META / "GSE254005_SAMPLE_MANIFEST.csv").is_file()
    )
    rows.append({
        "phase": "Phase1",
        "audit_mode": "SOURCE_DOWNLOAD_MANIFEST_SIZE_SHA256_RECOMPUTE",
        "manifest": str(phase1_manifest),
        "entries": 1,
        "passed": 1 if phase1_ok else 0,
        "failed": 0 if phase1_ok else 1,
        "status": "PASS" if phase1_ok else "FAIL",
        "notes": "No formal PHASE1_MANIFEST; GSE254005 source archive and sample manifest verified",
    })
    for phase in range(2, 8):
        result = validate_manifest(phase)
        rows.append(result)
        log(f"FREEZE_AUDIT {result['phase']} {result['passed']}/{result['entries']} {result['status']}")
    frame = pd.DataFrame(rows)
    atomic_csv(frame, P8_META / "PHASE8_PRIOR_MANIFEST_AUDIT.csv")
    if not frame["status"].eq("PASS").all():
        raise RuntimeError("Phase0-7 freeze audit failed")
    return frame


def build_ac_audit(sample: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ac = sample.loc[sample["region"].eq("AC")].copy()
    if len(ac) != 37 or ac["GSM"].duplicated().any():
        raise RuntimeError(f"Expected 37 unique explicit AC GSM, observed {len(ac)}")
    audit = pd.DataFrame({
        "GSM": ac["GSM"],
        "donor_id": ac["patient_id"],
        "disease": ac["healthy_or_crohn"],
        "region": ac["region"],
        "inflammation_state": ac["inflammation_state"],
        "treatment": ac["biologic_class"],
        "eligible": "NO_METADATA_ONLY",
        "exclusion_reason": "State not aligned to preregistered Healthy/active regional contrasts",
        "sample_title": ac["sample_title"],
        "processed_filename": ac["GEO_processed_filename"],
        "processed_file_size": ac["processed_file_size"],
        "processed_download_url": ac["processed_download_url"],
        "platform": ac["platform"],
        "metadata_source": ac["metadata_source"],
        "notes": ac["notes"],
    })
    healthy = audit["disease"].eq("Healthy") & audit["inflammation_state"].eq("Healthy")
    active = audit["disease"].eq("Crohn") & audit["inflammation_state"].eq("Active_inflamed")
    audit.loc[healthy, ["eligible", "exclusion_reason"]] = ["YES_HEALTHY_REGION_AND_AC_DISEASE", "NA"]
    audit.loc[active, ["eligible", "exclusion_reason"]] = ["YES_ACTIVE_REGION_AND_AC_DISEASE", "NA"]
    or95 = audit["donor_id"].eq("OR95")
    audit.loc[or95, ["eligible", "exclusion_reason"]] = [
        "NO_REGION_IDENTITY_ANOMALY",
        "GSM title/region and processed archive filename are reciprocally swapped between OR95 AC and TI; excluded from region estimates",
    ]
    audit = audit.sort_values(["disease", "inflammation_state", "donor_id", "GSM"]).reset_index(drop=True)
    atomic_csv(audit, P8_META / "HUMAN_AC_REGION_AUDIT.csv")

    availability_rows: list[dict[str, object]] = []
    for (region, disease, state), group in sample.groupby(["region", "healthy_or_crohn", "inflammation_state"], sort=True):
        availability_rows.append({
            "region": region,
            "disease": disease,
            "inflammation_state": state,
            "GSM_n": len(group),
            "donor_n": group["patient_id"].nunique(),
            "treatment_classes": ";".join(sorted(group["biologic_class"].unique())),
            "region_complete": "YES" if group["region"].ne("").all() else "NO",
            "donor_complete": "YES" if group["patient_id"].ne("").all() else "NO",
        })
    availability = pd.DataFrame(availability_rows)
    atomic_csv(availability, P8_META / "HUMAN_REGION_SAMPLE_AVAILABILITY.csv")

    def donor_set(region: str, disease: str, state: str) -> set[str]:
        return set(sample.loc[
            sample["region"].eq(region)
            & sample["healthy_or_crohn"].eq(disease)
            & sample["inflammation_state"].eq(state),
            "patient_id",
        ])

    healthy_nominal = donor_set("TI", "Healthy", "Healthy") & donor_set("AC", "Healthy", "Healthy")
    active_nominal = donor_set("TI", "Crohn", "Active_inflamed") & donor_set("AC", "Crohn", "Active_inflamed")
    pairing = pd.DataFrame([
        {
            "comparison": "Healthy_AC_vs_TI",
            "TI_donor_n": len(donor_set("TI", "Healthy", "Healthy")),
            "AC_donor_n": len(donor_set("AC", "Healthy", "Healthy")),
            "nominal_paired_donor_n": len(healthy_nominal),
            "analysis_paired_donor_n": len(healthy_nominal - {"OR95"}),
            "paired_donors": ";".join(sorted(healthy_nominal - {"OR95"})),
            "design": "~patient_id + region",
            "notes": "OR95 excluded for reciprocal processed-filename region anomaly",
        },
        {
            "comparison": "ActiveCrohn_AC_vs_TI",
            "TI_donor_n": len(donor_set("TI", "Crohn", "Active_inflamed")),
            "AC_donor_n": len(donor_set("AC", "Crohn", "Active_inflamed")),
            "nominal_paired_donor_n": len(active_nominal),
            "analysis_paired_donor_n": len(active_nominal),
            "paired_donors": ";".join(sorted(active_nominal)),
            "design": "~region",
            "notes": "Unpaired; treatment/cohort composition remains a limitation",
        },
        {
            "comparison": "AC_ActiveCrohn_vs_Healthy",
            "TI_donor_n": 0,
            "AC_donor_n": 22,
            "nominal_paired_donor_n": 0,
            "analysis_paired_donor_n": 0,
            "paired_donors": "NA",
            "design": "~condition",
            "notes": "16 Healthy AC (OR95 excluded) vs 6 active Crohn AC; secondary region-robustness effect only",
        },
    ])
    atomic_csv(pairing, P8_META / "HUMAN_REGION_PAIRING_AUDIT.csv")
    return audit, availability, pairing


CORE_PATHWAYS = [
    ("TNF_NFkB", "Hallmark", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
    ("Cytokine_interleukin", "Reactome", "R-HSA-449147"),
    ("IFN_alpha", "Hallmark", "HALLMARK_INTERFERON_ALPHA_RESPONSE"),
    ("IFN_gamma", "Hallmark", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
    ("JAK_STAT", "Hallmark", "HALLMARK_IL6_JAK_STAT3_SIGNALING"),
    ("Inflammatory_response", "Hallmark", "HALLMARK_INFLAMMATORY_RESPONSE"),
    ("Oxidative_phosphorylation", "Hallmark", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
    ("Fatty_acid_metabolism", "Hallmark", "HALLMARK_FATTY_ACID_METABOLISM"),
    ("Xenobiotic_metabolism", "Hallmark", "HALLMARK_XENOBIOTIC_METABOLISM"),
    ("UPR", "Hallmark", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
    ("Proliferation_E2F", "Hallmark", "HALLMARK_E2F_TARGETS"),
    ("Proliferation_G2M", "Hallmark", "HALLMARK_G2M_CHECKPOINT"),
]


def phase7_full_gsea(lineage: str, collection: str) -> pd.DataFrame:
    path = PHASE7_GSEA / "full" / f"Human_{lineage}_{collection}_GSEA.csv"
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    frame["NES"] = pd.to_numeric(frame["NES"])
    frame["padj"] = pd.to_numeric(frame["padj"])
    return frame


def freeze_programs() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    cache: dict[tuple[str, str], pd.DataFrame] = {}
    for lineage in ("Epithelial", "T", "Myeloid"):
        for concept, collection, pathway_id in CORE_PATHWAYS:
            frame = cache.setdefault((lineage, collection), phase7_full_gsea(lineage, collection))
            match = frame.loc[frame["pathway_id"].eq(pathway_id)]
            if len(match) != 1:
                raise RuntimeError(f"Core pathway lookup failed: {lineage} {collection} {pathway_id}")
            row = match.iloc[0]
            rows.append({
                "lineage": lineage,
                "program": concept,
                "collection": collection,
                "pathway_id": pathway_id,
                "pathway_name": row["pathway_name"],
                "phase7_TI_NES": row["NES"],
                "phase7_TI_FDR": row["padj"],
                "phase7_TI_direction": "ActiveCrohn_higher" if row["NES"] > 0 else "ActiveCrohn_lower",
                "registration_source": "Phase7 full Human-native GSEA; frozen before AC expression inspection",
            })
    core = pd.DataFrame(rows)
    atomic_csv(core, P8_META / "REGION_CORE_PROGRAM_REGISTRY.csv")

    dog = pd.read_csv(DOG_FORWARD, dtype=str, keep_default_na=False)
    dog_rows: list[dict[str, object]] = []
    evidence_pattern = re.compile(
        r"^(?P<collection>[^:]+):(?P<pathway>[^|]+)\|NES=(?P<NES>[^|]+)\|FDR=(?P<FDR>[^|]+)\|(?P<direction>.+)$"
    )
    for row in dog.itertuples(index=False):
        candidates = []
        for item in str(row.human_representative_evidence).split(";"):
            match = evidence_pattern.match(item)
            if match:
                parsed = match.groupdict()
                parsed["FDR_numeric"] = float(parsed["FDR"])
                parsed["NES_numeric"] = float(parsed["NES"])
                candidates.append(parsed)
        if not candidates:
            raise RuntimeError(f"No assessable frozen representative for {row.lineage} / {row.program}")
        best = sorted(candidates, key=lambda value: (value["FDR_numeric"], value["pathway"]))[0]
        dog_rows.append({
            "lineage": row.lineage,
            "program": row.program,
            "Dog_CIE_direction": row.dog_CIE_direction,
            "FRE_direction": row.FRE_direction or "NA_NOT_AVAILABLE",
            "representative_collection": best["collection"],
            "representative_pathway_id": best["pathway"],
            "phase7_TI_NES": best["NES_numeric"],
            "phase7_TI_FDR": best["FDR_numeric"],
            "phase7_TI_direction": best["direction"],
            "selection_rule": "Lowest Phase7 TI FDR among preregistered Phase6 Dog-forward Human representatives; tie by pathway ID",
            "freeze_time_utc": utc_now(),
        })
    dog_registry = pd.DataFrame(dog_rows)
    if len(dog_registry) != 25:
        raise RuntimeError("Expected 25 Dog-forward representative registrations")
    atomic_csv(dog_registry, P8_META / "DOG_FORWARD_REGION_REPRESENTATIVE_FREEZE.csv")
    return core, dog_registry


def write_protocol() -> None:
    value = f"""# Phase8A protocol — Region Control Data Audit

Frozen at: {utc_now()}

## Scope and questions

1. Are Human Crohn disease programs dependent on terminal ileum?
2. Could later Dog-duodenum/Human-TI observations be sensitive to anatomical region?
3. Which programs are region-robust, region-sensitive, or region-unknown?

Phase8A estimates Human region robustness only. It does not claim Dog–Human conservation and does not run ontology harmonization, ortholog mapping, singIST, SATURN or CONCORD.

## Frozen inputs

- Phase7 TI Primary remains 16 Healthy TI vs 11 active Crohn TI donors and is read-only.
- Dog, FRE, Human TI, IBDverse and Dog-forward Phase0–7 outputs remain read-only.
- AC metadata are selected only from the frozen 85-row GSE266616 sample master.

## Region identity rule

Direct GSM title/region metadata define AC/TI. OR95 is retained in the audit but excluded from every new region estimate because GSM8252541 and GSM8252542 have reciprocally swapped AC/TI processed filenames. This Phase8 exclusion does not modify Phase7.

## Preregistered analyses

- A — Healthy baseline region: 15 unambiguous paired donors, edgeR QL `~patient_id + region`; positive logFC/NES = AC higher.
- B — Active Crohn region: 11 TI vs 6 AC donors, zero paired donors, edgeR QL `~region`; positive = AC higher. This is non-paired and treatment/cohort-confounded.
- C — AC disease effect: 16 Healthy AC vs 6 active Crohn AC donors, edgeR QL `~condition`; positive = active Crohn higher. It is used only to assess the region robustness of frozen Phase7 TI programs and cannot replace the TI Primary result.

Every pseudobulk column is one donor/sample. Epithelial, T and Myeloid require at least 20 retained annotated cells per donor.

## QC and annotation

AC uses the Phase7 Human Route B thresholds unchanged: 200–6000 genes inclusive, >1000 UMI, mt% <50%, no added doublet filter. Broad-lineage labels use the frozen Phase7 Human-native classifier and frozen thresholds without retraining or Dog labels.

## Program testing and interpretation

No new all-pathway disease-program search is allowed. The 12 pathway concepts in `REGION_CORE_PROGRAM_REGISTRY.csv` and one deterministic frozen representative for each of the 25 Dog-forward programs are registered before AC expression inspection.

For TI-vs-AC disease effects:

- `ROBUST_FDR`: same NES sign and AC FDR <0.05.
- `ROBUST_DIRECTIONAL_LIMITED`: same sign but AC FDR >=0.05.
- `REGION_SENSITIVE_FDR`: opposite sign and AC FDR <0.05.
- `REGION_UNKNOWN_OPPOSITE_WEAK`: opposite sign but AC FDR >=0.05.
- `NOT_ASSESSABLE`: no valid statistic.

Region contrasts A/B are reported alongside disease-effect stability. With only six active AC donors and no active cross-region pairs, failure to reach FDR is not evidence of absence.

## Stop rule

After the Phase8 report, Gate and Manifest are generated, stop. Phase9/10 and all cross-species algorithms remain unexecuted.
"""
    atomic_text(value, PROTOCOL)


def build_capability_gate(audit: pd.DataFrame, prior: pd.DataFrame, core: pd.DataFrame, dog: pd.DataFrame) -> pd.DataFrame:
    healthy = audit.loc[audit["eligible"].eq("YES_HEALTHY_REGION_AND_AC_DISEASE")]
    active = audit.loc[audit["eligible"].eq("YES_ACTIVE_REGION_AND_AC_DISEASE")]
    checks = [
        ("Phase0-7 frozen inputs verified", prior["status"].eq("PASS").all(), ";".join(f"{r.phase}={r.status}" for r in prior.itertuples())),
        ("AC metadata cardinality frozen", len(audit) == 37 and audit["GSM"].nunique() == 37, f"GSM={len(audit)};donors={audit['donor_id'].nunique()}"),
        ("Healthy AC minimum met", healthy["donor_id"].nunique() >= 3, f"eligible_unambiguous_donors={healthy['donor_id'].nunique()}"),
        ("Active Crohn AC minimum met", active["donor_id"].nunique() >= 3, f"eligible_donors={active['donor_id'].nunique()}"),
        ("Region and donor identity complete", audit["region"].eq("AC").all() and audit["donor_id"].ne("").all(), "OR95 filename anomaly excluded, not imputed"),
        ("Per-GSM processed URLs available", audit.loc[audit["eligible"].str.startswith("YES"), "processed_download_url"].str.startswith("https://ftp.ncbi.nlm.nih.gov/geo/samples/").all(), "22 analysis-eligible per-GSM archives; OR95 excluded before download"),
        ("Program tests preregistered", len(core) == 36 and len(dog) == 25, f"core={len(core)};Dog-forward={len(dog)}"),
    ]
    frame = pd.DataFrame([
        {"criterion": name, "status": "PASS" if okay else "FAIL", "evidence": evidence}
        for name, okay, evidence in checks
    ])
    frame["gate"] = "GREEN_AC_REGION_ANALYSIS_CAPABLE" if frame["status"].eq("PASS").all() else "YELLOW_AC_REGION_CONTROL_LIMITED"
    atomic_csv(frame, P8_META / "PHASE8A_AC_CAPABILITY_GATE.csv")
    if not frame["status"].eq("PASS").all():
        raise RuntimeError("AC capability gate did not pass")
    return frame


def main() -> int:
    for directory in (
        P8_META,
        P8 / "human_AC_TI" / "raw_processed",
        P8 / "human_AC_TI" / "per_sample_qc",
        P8 / "human_AC_TI" / "per_sample_annotated",
        P8 / "human_AC_TI" / "pseudobulk",
        P8 / "human_AC_TI" / "de_results",
        P8 / "human_AC_TI" / "ranked_lists",
        P8 / "human_AC_TI" / "gsea",
        P8_SUMMARY,
        P8 / "figures",
        P8 / "scripts",
        P8_LOGS,
    ):
        directory.mkdir(parents=True, exist_ok=True)
    PROJECT_LOG.write_text("", encoding="utf-8")
    (P8_LOGS / "phase8_region_audit.log").write_text("", encoding="utf-8")
    write_protocol()
    sample = pd.read_csv(SAMPLE_MASTER, dtype=str, keep_default_na=False)
    if len(sample) != 85:
        raise RuntimeError("Frozen GSE266616 sample master changed")
    prior = audit_prior_freeze()
    audit, _, pairing = build_ac_audit(sample)
    core, dog = freeze_programs()
    gate = build_capability_gate(audit, prior, core, dog)
    log(
        f"PHASE8A_CAPABILITY {gate['gate'].iloc[0]} AC_GSM={len(audit)} AC_donors={audit['donor_id'].nunique()} "
        f"healthy_analysis={audit['eligible'].eq('YES_HEALTHY_REGION_AND_AC_DISEASE').sum()} "
        f"active_analysis={audit['eligible'].eq('YES_ACTIVE_REGION_AND_AC_DISEASE').sum()} "
        f"healthy_paired={int(pairing.loc[pairing['comparison'].eq('Healthy_AC_vs_TI'),'analysis_paired_donor_n'].iloc[0])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
