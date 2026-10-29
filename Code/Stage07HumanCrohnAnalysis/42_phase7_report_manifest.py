from __future__ import annotations

import csv
import hashlib
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
META = ROOT / "00_metadata"
RESULTS = ROOT / "08_results"
HUMAN = ROOT / "07_human_crohn" / "GSE266616"
IBD = ROOT / "07_human_crohn" / "IBDverse"
REPORT = RESULTS / "PHASE7_REPORT.md"
GATE = META / "PHASE7_GATE.csv"
MANIFEST = META / "PHASE7_MANIFEST.csv"
PRIOR_AUDIT = META / "PHASE7_PRIOR_MANIFEST_AUDIT.csv"
MANIFEST_VERIFICATION = META / "PHASE7_MANIFEST_VERIFICATION.csv"
REPORT_LOG = ROOT / "logs" / "phase7_report.log"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def atomic_text(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    tmp = Path(name)
    try:
        tmp.write_text(text, encoding="utf-8")
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def audit_prior_manifests() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for phase in range(2, 7):
        manifest = META / f"PHASE{phase}_MANIFEST.csv"
        frame = pd.read_csv(manifest)
        path_col = next(col for col in frame if "path" in col.lower())
        sha_col = next(col for col in frame if "sha256" in col.lower())
        size_col = next(col for col in frame if col in {"bytes", "file_size", "size_bytes"})
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
                failures.append(str(path.relative_to(ROOT)).replace("\\", "/"))
        rows.append({
            "phase": f"Phase{phase}", "manifest": str(manifest), "entries": len(frame),
            "passed": passed, "failed": len(failures),
            "status": "PASS" if not failures else "FAIL",
            "failure_examples": ";".join(failures[:5]),
        })
    result = pd.DataFrame(rows)
    atomic_csv(result, PRIOR_AUDIT)
    return result


def all_pass(path: Path, status_col: str = "status") -> bool:
    frame = pd.read_csv(path)
    return bool(len(frame) and frame[status_col].eq("PASS").all())


def exact_program(programs: pd.DataFrame, lineage: str, name: str) -> pd.Series | None:
    match = programs.loc[
        programs["lineage"].eq(lineage)
        & programs["program"].astype(str).str.casefold().eq(name.casefold())
    ].sort_values("FDR")
    return None if match.empty else match.iloc[0]


def program_table(programs: pd.DataFrame) -> str:
    focus = {
        "Epithelial": ["E2F Targets", "Interferon Alpha Response", "Tnfa Signaling Via Nfkb", "Oxidative Phosphorylation", "Fatty Acid Metabolism", "Xenobiotic Metabolism"],
        "T": ["Mtorc1 Signaling", "Interferon Gamma Response", "Tnfa Signaling Via Nfkb", "Inflammatory Response", "Il6 Jak Stat3 Signaling"],
        "Myeloid": ["Interferon Gamma Response", "Tnfa Signaling Via Nfkb", "Inflammatory Response", "Il6 Jak Stat3 Signaling", "Signaling by Interleukins", "Eukaryotic Translation Elongation"],
    }
    lines = ["| Lineage | Program | Collection | NES | FDR | Treatment sensitivity |", "|---|---|---:|---:|---:|---|"]
    for lineage, names in focus.items():
        for name in names:
            row = exact_program(programs, lineage, name)
            if row is None:
                continue
            lines.append(
                f"| {lineage} | {row['program']} | {row['collection']} | {row['NES']:.3f} | {row['FDR']:.3g} | {row['treatment_sensitivity']} |"
            )
    return "\n".join(lines)


def top_gene_text(lineage: str) -> str:
    path = HUMAN / "pseudobulk" / "de_results" / f"Human_{lineage}_ActiveCrohn_vs_Healthy.csv"
    frame = pd.read_csv(path)
    up = frame.loc[frame["logFC"] > 0].sort_values(["FDR", "PValue"]).head(3)
    down = frame.loc[frame["logFC"] < 0].sort_values(["FDR", "PValue"]).head(3)
    fmt = lambda x: ", ".join(f"{r.gene_symbol} (logFC {r.logFC:.2f}, FDR {r.FDR:.2g})" for r in x.itertuples())
    return f"higher: {fmt(up)}; lower: {fmt(down)}"


def build_gate(prior: pd.DataFrame) -> pd.DataFrame:
    download = pd.read_csv(META / "GSE266616_DOWNLOAD_MANIFEST.csv")
    programs = pd.read_csv(HUMAN / "gsea" / "HUMAN_DISEASE_PROGRAMS.csv")
    dog = pd.read_csv(HUMAN / "gsea" / "DOG_FORWARD_IN_HUMAN_OBSERVATION.csv")
    prohibited_names = {
        "GSE266616_RAW.tar", "scRNA_gene_counts.mtx", "Full_cohort.h5ad", "DGE.zip",
    }
    prohibited_found = [
        str(path.relative_to(ROOT)).replace("\\", "/")
        for path in (ROOT / "07_human_crohn").rglob("*")
        if path.is_file() and (path.name in prohibited_names or path.suffix.lower() in {".fastq", ".sra"})
    ]
    criteria = [
        ("Phase2-6 frozen manifests remain valid", prior["status"].eq("PASS").all(), ";".join(f"{r.phase}={r.passed}/{r.entries}" for r in prior.itertuples())),
        ("Phase7A metadata gate passed", all_pass(META / "PHASE7A_METADATA_GATE.csv"), "Primary sample/donor freeze and 79-vs-85 audit"),
        ("48 TI processed archives verified", len(download) == 48 and download["archive_valid"].eq("YES").all() and download["actual_size"].eq(download["expected_size"]).all(), f"n={len(download)};bytes={int(download['actual_size'].sum())}"),
        ("Prohibited large/raw routes avoided", not prohibited_found, "none" if not prohibited_found else ";".join(prohibited_found)),
        ("Human Route B QC gate passed", all_pass(HUMAN / "qc" / "PHASE7B_QC_GATE.csv"), "290814 raw; 237878 retained"),
        ("Human Level1 annotation gate passed", all_pass(HUMAN / "annotation" / "PHASE7B_ANNOTATION_GATE.csv"), "Epithelial/T/Myeloid donor-LOO F1 and donor coverage"),
        ("Donor-level pseudobulk gate passed", all_pass(HUMAN / "pseudobulk" / "metadata" / "PHASE7C_PSEUDOBULK_GATE.csv"), "27 donors per primary lineage"),
        ("Three complete Primary edgeR outputs exist", all((HUMAN / "pseudobulk" / "de_results" / f"Human_{x}_ActiveCrohn_vs_Healthy.csv").is_file() for x in ("Epithelial", "T", "Myeloid")), "edgeR QL ~condition; Healthy reference"),
        ("Human-native unbiased GSEA gate passed", all_pass(HUMAN / "gsea" / "metadata" / "PHASE7C_GSEA_GATE.csv"), f"frozen_programs={len(programs)}"),
        ("Dog-forward list read only after Human freeze", len(dog) == 25 and (HUMAN / "gsea" / "DOG_FORWARD_POSTHOC_PROVENANCE.csv").is_file(), "25/25; sequence assertion recorded"),
        ("IBDverse lightweight split-cohort replication passed", all_pass(IBD / "metadata" / "PHASE7D_IBDVERSE_GATE.csv"), "summary-only Discovery/Replication; Full_cohort.h5ad absent"),
        ("Phase8 and later analyses not executed", not (ROOT / "07_human_crohn" / "GSE185224").exists(), "Stopped before Phase8"),
    ]
    frame = pd.DataFrame([
        {"criterion": criterion, "status": "PASS" if passed else "FAIL", "evidence": evidence}
        for criterion, passed, evidence in criteria
    ])
    frame["gate"] = "GREEN_PHASE7_COMPLETE" if frame["status"].eq("PASS").all() else "YELLOW_OR_RED_REVIEW_REQUIRED"
    atomic_csv(frame, GATE)
    if not frame["status"].eq("PASS").all():
        raise RuntimeError("Phase7 final gate did not pass")
    return frame


def build_report() -> str:
    sample = pd.read_csv(META / "GSE266616_HUMAN_SAMPLE_MASTER.csv")
    donors = pd.read_csv(META / "GSE266616_HUMAN_DONOR_MASTER.csv")
    download = pd.read_csv(META / "GSE266616_DOWNLOAD_MANIFEST.csv")
    qc = pd.read_csv(HUMAN / "qc" / "HUMAN_QC_PER_SAMPLE.csv")
    cells = pd.read_csv(HUMAN / "annotation" / "HUMAN_LEVEL1_CELL_COUNTS.csv")
    programs = pd.read_csv(HUMAN / "gsea" / "HUMAN_DISEASE_PROGRAMS.csv")
    gene_sensitivity = pd.read_csv(HUMAN / "pseudobulk" / "sensitivity" / "HUMAN_TREATMENT_GENE_SENSITIVITY_SUMMARY.csv")
    paired = pd.read_csv(HUMAN / "pseudobulk" / "paired_de_results" / "HUMAN_PAIRED_DE_SUMMARY.csv")
    dog = pd.read_csv(HUMAN / "gsea" / "DOG_FORWARD_IN_HUMAN_OBSERVATION.csv")
    dog_summary = pd.read_csv(HUMAN / "gsea" / "DOG_FORWARD_IN_HUMAN_OBSERVATION_SUMMARY.csv")
    ibd_summary = pd.read_csv(IBD / "summary" / "IBDVERSE_HUMAN_PROGRAM_REPLICATION_SUMMARY.csv")
    ibd_programs = pd.read_csv(IBD / "summary" / "IBDVERSE_HUMAN_PROGRAM_REPLICATION.csv")

    active = sample.loc[sample["primary_contrast_eligible"].eq("YES_SELECTED") & sample["healthy_or_crohn"].eq("Crohn")]
    treatments = active["biologic_class"].value_counts()
    count_by_state = download["state"].value_counts()
    cell_totals = cells.groupby("human_level1")["n_cells"].sum().sort_values(ascending=False)
    de_rows = []
    for lineage in ("Epithelial", "T", "Myeloid"):
        frame = pd.read_csv(HUMAN / "pseudobulk" / "de_results" / f"Human_{lineage}_ActiveCrohn_vs_Healthy.csv")
        de_rows.append((lineage, len(frame), int((frame["FDR"] < 0.05).sum()), int(((frame["FDR"] < 0.05) & (frame["logFC"] > 0)).sum()), int(((frame["FDR"] < 0.05) & (frame["logFC"] < 0)).sum())))
    treatment_categories = programs["treatment_sensitivity"].value_counts()
    direction_change = int(treatment_categories[treatment_categories.index.str.startswith("DIRECTION_CHANGE")].sum())
    all_available_same = int(len(programs) - direction_change)
    ibd_counts = ibd_programs.loc[ibd_programs["ibdverse_mapped"].eq("YES"), "ibdverse_status"].value_counts()
    dog_counts = dict(zip(dog_summary["human_observation"], dog_summary["program_n"], strict=True))

    question_rows = []
    for row in dog.itertuples(index=False):
        if row.program in {"Mitochondrial respiration", "Metabolic programs", "TNF / NF-kappaB", "Cytokine / interleukin", "Cytokine / interleukin / JAK-STAT"}:
            question_rows.append(f"| {row.lineage} | {row.program} | {row.dog_CIE_direction} | {row.FRE_direction} | {row.human_observation} |")

    report = f"""# Phase7 report — Human Crohn terminal-ileum disease-program discovery

Final gate: **GREEN — PHASE7 HUMAN CROHN PROGRAM DISCOVERY COMPLETE**

Generated: {datetime.now(timezone.utc).isoformat()}

## 1. Metadata gate

- Current GEO explicit records: **{len(sample)} GSM / {donors['patient_id'].nunique()} donors** ({int((donors['healthy_or_crohn'] == 'Healthy').sum())} Healthy, {int((donors['healthy_or_crohn'] == 'Crohn').sum())} Crohn).
- TI donors: **{int((donors['n_TI_samples'] > 0).sum())}**.
- Primary was frozen before expression inspection at **16 Healthy TI vs 11 active inflamed Crohn TI donors**, one sample per donor.
- The original contiguous 79 records are GSM8252464–GSM8252542. Six later records, GSM8647246–GSM8647251, are three Crohn donors (HA98, HA102, HA104), each with an adjacent/inflamed ascending-colon pair. This explains **79 vs 85 GSM** as later-added AC records.
- A separate design-text discrepancy remains: the 85 explicit records resolve to 17 Healthy + 29 Crohn donors, not the stated 17 + 32. The three unrepresented Crohn donor identities were not guessed. This does not affect the explicit Primary TI freeze.
- One processed filename (`GSM8252542_OR95AC.tar.gz`) conflicts with its direct GEO sample title/region (`OR95TI`, TI); the formal region follows the explicit GEO region/title metadata and the anomaly is preserved in the audit.

## 2. Biologic treatment

Among 11 Primary active-TI Crohn donors: **NO_BIOLOGIC={int(treatments.get('NO_BIOLOGIC', 0))}**, **ANTI_TNF={int(treatments.get('ANTI_TNF', 0))}**, **ANTI_IL12_23={int(treatments.get('ANTI_IL12_23', 0))}**. All five preregistered sensitivity subsets met the minimum of three active Crohn donors. They are treatment-confounding sensitivities, not causal treatment-effect analyses.

## 3. Download

- **48/48** per-GSM processed TI archives, **{int(download['actual_size'].sum()):,} bytes ({download['actual_size'].sum() / 1024**3:.2f} GiB)**.
- State composition: Healthy={int(count_by_state.get('Healthy', 0))}, active={int(count_by_state.get('Active_inflamed', 0))}, adjacent={int(count_by_state.get('Adjacent', 0))}, non-inflamed={int(count_by_state.get('Non_inflamed', 0))}.
- Every archive has the 10x matrix/features/barcodes triplet and a local SHA-256. GEO did not provide repository checksums, recorded as `NA_NOT_PROVIDED`.
- Not downloaded: whole GSE266616 RAW tar, SRA/FASTQ, Zenodo 6.1-GB count matrix, IBDverse 25.48-GB Full_cohort.h5ad.

## 4. Human QC

- GEO per-GSM matrices are Cell Ranger filtered, not author-final QC: six published Healthy TI barcode sets were exact subsets of GEO, while GEO retained additional cells.
- Formal route: **Route B**, uniformly 200–6000 genes inclusive, >1000 UMI and mitochondrial fraction <50%; no additional doublet pipeline.
- Raw cells: **{int(qc['raw_cells'].sum()):,}**; retained cells: **{int(qc['retained_cells'].sum()):,} ({qc['retained_cells'].sum()/qc['raw_cells'].sum():.1%})**; all 48 samples retained cells.
- Retained by state: {', '.join(f'{state}={int(value):,}' for state, value in qc.groupby('state')['retained_cells'].sum().items())}.

## 5. Human Level1 annotation

- Published annotation provenance covers six Healthy TI libraries only. A Human-native sparse broad-lineage classifier was trained on 29,927 exact published reference cells and audited by leave-one-library-out CV.
- F1: Epithelial 0.9997, T 0.9990, Myeloid 0.9962; `Uncertain` **{int(cell_totals.get('Uncertain', 0)):,} / {int(cell_totals.sum()):,} ({cell_totals.get('Uncertain', 0)/cell_totals.sum():.2%})**.
- Final cells by Level1: {', '.join(f'{name}={int(value):,}' for name, value in cell_totals.items())}.
- All 16 Healthy and 11 active Crohn Primary donors had >=20 Epithelial, T and Myeloid cells. Counts remain in `layers['counts']`.

## 6. Primary donor-level pseudobulk and edgeR

Design: edgeR quasi-likelihood `~condition`, Healthy reference; positive logFC means active Crohn higher. No cell-level disease P values were used.

| Lineage | Tested genes | FDR<0.05 | Higher | Lower |
|---|---:|---:|---:|---:|
{os.linesep.join(f'| {a} | {b:,} | {c:,} | {d:,} | {e:,} |' for a,b,c,d,e in de_rows)}

- Epithelial top effects: {top_gene_text('Epithelial')}.
- T top effects: {top_gene_text('T')}.
- Myeloid top effects: {top_gene_text('Myeloid')}.
- Platform sensitivity was not fitted: NovaSeq dominates and NextSeq has one donor per condition, so no claim of platform adjustment is made.

## 7. Human-native disease programs

Full unbiased GSEA used Human-native Hallmark 2026.1.Hs, Reactome v97 and GO BP. It tested all eligible pathways before any Dog-forward file was read. **{len(programs):,}** Primary programs met FDR<0.05: Epithelial={int((programs['lineage']=='Epithelial').sum())}, T={int((programs['lineage']=='T').sum())}, Myeloid={int((programs['lineage']=='Myeloid').sum())}.

{program_table(programs)}

Direction is NES>0 active Crohn higher, NES<0 active Crohn lower. Major patterns are epithelial proliferation/IFN/inflammatory activation with oxidative, fatty-acid and xenobiotic programs lower; T activation/IFN/TNF/cytokine and mTOR programs higher; and Myeloid IFN/TNF/inflammatory/interleukin programs higher with selected translation programs lower.

## 8. Treatment sensitivity

- Gene-level Primary-vs-sensitivity Spearman logFC ranges: Epithelial {gene_sensitivity.loc[gene_sensitivity.lineage.eq('Epithelial'),'spearman_logFC'].min():.2f}–{gene_sensitivity.loc[gene_sensitivity.lineage.eq('Epithelial'),'spearman_logFC'].max():.2f}; T {gene_sensitivity.loc[gene_sensitivity.lineage.eq('T'),'spearman_logFC'].min():.2f}–{gene_sensitivity.loc[gene_sensitivity.lineage.eq('T'),'spearman_logFC'].max():.2f}; Myeloid {gene_sensitivity.loc[gene_sensitivity.lineage.eq('Myeloid'),'spearman_logFC'].min():.2f}–{gene_sensitivity.loc[gene_sensitivity.lineage.eq('Myeloid'),'spearman_logFC'].max():.2f}.
- Program level: **{all_available_same}/{len(programs)} ({all_available_same/len(programs):.1%})** kept the Primary direction across every assessable treatment sensitivity; **{direction_change}** changed direction in at least one sensitivity.
- Small n=3 treatment-only subsets make T/Myeloid significance less stable than direction; this remains a carried limitation.

## 9. Secondary paired within-Crohn analysis

- `active vs adjacent`: **11 paired donors**, donor-blocked `~patient_id + state`; FDR<0.05 genes were Epithelial={int(paired.loc[paired.lineage.eq('Epithelial'),'FDR05'].iloc[0])}, T={int(paired.loc[paired.lineage.eq('T'),'FDR05'].iloc[0])}, Myeloid={int(paired.loc[paired.lineage.eq('Myeloid'),'FDR05'].iloc[0])}.
- Paired GSEA recapitulated strong local active-inflammation-associated IFN/inflammatory programs in all three lineages; Epithelial also showed proliferation/UPR higher and compound-functionalization lower.
- `active vs non-inflamed`: zero complete donor pairs, so it was not fitted.
- These are cross-sectional within-Crohn active-inflammation associations, not a progression trajectory.

## 10. IBDverse lightweight replication

- Official BioStudies API exposed a 292,385,193-byte DGE zip despite the website's legacy directory links returning 404. HTTP Range extraction obtained only Discovery, Replication and Full inflamed-CD-vs-Healthy DE/GSEA summaries (**61,021,220 bytes**); the full 25,476,124,776-byte h5ad was not downloaded.
- Strict replication required the identical signed Reactome pathway and cell subtype to have FDR<0.05 in both Discovery and Replication with the same CD-vs-Healthy direction. Full-cohort results did not replace this split-cohort rule.
- Of 470 mapped significant GSE266616 Reactome programs: **same-direction={int(ibd_counts.get('SAME_DIRECTION_REPLICATED',0))}**, **mixed subtype={int(ibd_counts.get('MIXED_CELL_SUBTYPE_REPLICATION',0))}**, **opposite={int(ibd_counts.get('OPPOSITE_DIRECTION_REPLICATED',0))}**, **not replicated={int(ibd_counts.get('NOT_REPLICATED_SPLIT_COHORT',0))}**.
- Same-direction replicated by lineage: Epithelial=16, T=15, Myeloid=98. Examples include epithelial interferon/interleukin programs and Myeloid ITGA4+ macrophage IFN, interleukin, antigen-presentation and NF-kB programs.

## 11. Dog-forward post hoc observation

The Human program table was SHA-frozen before Dog inputs were opened. Among 25 Dog-forward programs: **same={int(dog_counts.get('SAME_DIRECTION',0))}**, **opposite={int(dog_counts.get('OPPOSITE_DIRECTION',0))}**, **weak={int(dog_counts.get('WEAK',0))}**, **not prominent={int(dog_counts.get('NOT_PROMINENT',0))}**, **not assessable={int(dog_counts.get('NOT_ASSESSABLE',0))}**.

| Lineage | Question program | Dog CIE | FRE summary | Human observation |
|---|---|---|---|---|
{os.linesep.join(question_rows)}

The strongest three-way relationship is that T/Myeloid TNF–NF-kB and cytokine/JAK programs are Dog-CIE higher, FRE whole-mucosa lower, but Human cell-type-specific higher. Metabolic/mitochondrial direction is lineage-specific: Epithelial tracks Dog-lower, whereas T respiratory/metabolic representatives are Human-opposite/higher and Myeloid representatives are mostly not prominent. These are program-level observations, not formal cross-species effect correlations.

## 12. Gate, limitations and recommendation

Final gate: **GREEN**.

Carried limitations: the GEO design-text Crohn donor count is unresolved; author annotation reference is six Healthy TI libraries with no published Crohn labels; treatment-only subsets are small; the minor platform has n=1 per condition; IBDverse replication is summary-level Reactome and not a donor-level reanalysis; 25-program Dog-Human calls are post hoc observations only.

**Recommendation: YES — proceed to Phase8 region filter with these limitations carried forward.** Phase8, ontology harmonization, ortholog mapping, singIST, SATURN and CONCORD were not run in Phase7.

## 13. Principal deliverables

- `00_metadata/GSE266616_HUMAN_SAMPLE_MASTER.csv`
- `00_metadata/GSE266616_HUMAN_DONOR_MASTER.csv`
- `00_metadata/GSE266616_DOWNLOAD_MANIFEST.csv`
- `07_human_crohn/GSE266616/metadata/GSE266616_PHASE7_METADATA_AUDIT.xlsx`
- `07_human_crohn/GSE266616/h5ad/human_GSE266616_TI_raw.h5ad`
- `07_human_crohn/GSE266616/h5ad/human_GSE266616_TI_annotated.h5ad`
- `07_human_crohn/GSE266616/annotation/HUMAN_CELLTYPE_MAPPING.csv`
- `07_human_crohn/GSE266616/pseudobulk/de_results/Human_*_ActiveCrohn_vs_Healthy.csv`
- `07_human_crohn/GSE266616/gsea/HUMAN_DISEASE_PROGRAMS.csv`
- `07_human_crohn/GSE266616/gsea/DOG_FORWARD_IN_HUMAN_OBSERVATION.csv`
- `07_human_crohn/IBDverse/summary/IBDVERSE_HUMAN_PROGRAM_REPLICATION.csv`
- `00_metadata/PHASE7_GATE.csv`
- `00_metadata/PHASE7_MANIFEST.csv`
"""
    return report


def classify_file(path: Path) -> tuple[str, str, str]:
    relative = str(path.relative_to(ROOT)).replace("\\", "/")
    if relative.startswith("scripts/"):
        return "phase7_script", "project-authored", "Phase7 implementation"
    if relative.startswith("logs/"):
        return "phase7_log", "project runtime", "execution logging"
    if "/raw_processed/" in relative:
        return "GSE266616_processed_archive", "NCBI GEO per-GSM processed 10x", "Phase7B download"
    if relative.endswith(".h5ad"):
        return "Human_AnnData_object", "GSE266616 processed counts", "Phase7B QC/annotation"
    if "IBDverse" in relative and ("raw_summary" in relative or "metadata" in relative):
        return "IBDverse_lightweight_source_or_audit", "EBI BioStudies S-BSST2944/E-MTAB-16999", "Phase7D lightweight replication"
    if "IBDverse" in relative:
        return "IBDverse_replication_result", "IBDverse summary statistics", "Phase7D replication"
    if "gsea/gene_sets" in relative:
        return "Human_native_gene_set", "MSigDB/Reactome/GOA", "Phase7C gene-set freeze"
    if relative.startswith("00_metadata/"):
        return "phase7_metadata_or_gate", "GEO/Zenodo/project audit", "Phase7 metadata/gate/report"
    if relative.startswith("08_results/"):
        return "phase7_report", "Phase7 frozen results", "Phase7 reporting"
    if relative.endswith(".png"):
        return "phase7_figure", "Phase7 analysis", "Phase7 visualization"
    return "phase7_result", "Phase7 analysis", "Phase7 analysis"


def collect_manifest_files() -> list[Path]:
    files: set[Path] = set()
    for base in (HUMAN, IBD):
        files.update(path for path in base.rglob("*") if path.is_file())
    for path in META.iterdir():
        if path.is_file() and (
            path.name.startswith(("GSE266616_", "PHASE7", "HUMAN_"))
            or path.name in {"GSE266616_PHASE7_METADATA_AUDIT.xlsx"}
        ):
            files.add(path)
    files.update(path for path in RESULTS.glob("PHASE7*") if path.is_file())
    files.update(path for path in (ROOT / "scripts").iterdir() if path.is_file() and path.name[:2].isdigit() and 28 <= int(path.name[:2]) <= 42)
    files.update(path for path in (ROOT / "logs").glob("phase7*") if path.is_file() and not path.name.endswith(".pid"))
    files.update([
        ROOT / "06_pseudobulk" / "gsea" / "gene_sets" / "raw" / "Homo_sapiens.gene_info.gz",
        ROOT / "06_pseudobulk" / "gsea" / "gene_sets" / "raw" / "goa_human.gaf.gz",
    ])
    return sorted(
        path for path in files
        if path.is_file()
        and path not in {MANIFEST, MANIFEST_VERIFICATION}
        and not path.name.endswith((".tmp", ".part"))
    )


def build_manifest() -> pd.DataFrame:
    rows = []
    files = collect_manifest_files()
    for index, path in enumerate(files, start=1):
        role, source, step = classify_file(path)
        rows.append({
            "path_relative_to_project": str(path.relative_to(ROOT)).replace("\\", "/"),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "file_role": role,
            "source": source,
            "generation_step": step,
        })
        if index % 100 == 0:
            print(f"MANIFEST_HASH {index}/{len(files)}")
    frame = pd.DataFrame(rows)
    atomic_csv(frame, MANIFEST)
    verification = []
    for row in frame.itertuples(index=False):
        path = ROOT / row.path_relative_to_project
        valid = path.is_file() and path.stat().st_size == row.bytes and sha256(path) == row.sha256
        verification.append({"path_relative_to_project": row.path_relative_to_project, "status": "PASS" if valid else "FAIL"})
    verification_frame = pd.DataFrame(verification)
    atomic_csv(verification_frame, MANIFEST_VERIFICATION)
    role, source, step = classify_file(MANIFEST_VERIFICATION)
    verification_row = pd.DataFrame([{
        "path_relative_to_project": str(MANIFEST_VERIFICATION.relative_to(ROOT)).replace("\\", "/"),
        "bytes": MANIFEST_VERIFICATION.stat().st_size,
        "sha256": sha256(MANIFEST_VERIFICATION),
        "file_role": role, "source": source, "generation_step": step,
    }])
    frame = pd.concat([frame, verification_row], ignore_index=True)
    frame = frame.sort_values("path_relative_to_project").reset_index(drop=True)
    atomic_csv(frame, MANIFEST)
    if not verification_frame["status"].eq("PASS").all():
        raise RuntimeError("Phase7 manifest verification failed")
    return frame


def main() -> int:
    prior = audit_prior_manifests()
    gate = build_gate(prior)
    atomic_text(build_report(), REPORT)
    atomic_text(
        f"Phase7 report/gate generated {datetime.now(timezone.utc).isoformat()}\n"
        f"Gate={gate['gate'].iloc[0]} criteria={len(gate)} pass={int(gate['status'].eq('PASS').sum())}\n",
        REPORT_LOG,
    )
    manifest = build_manifest()
    print(f"PHASE7_COMPLETE gate={gate['gate'].iloc[0]} manifest_entries={len(manifest)} manifest_bytes={int(manifest['bytes'].sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
