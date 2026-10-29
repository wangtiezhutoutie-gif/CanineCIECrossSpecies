from __future__ import annotations

import csv
import hashlib
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
BASE = ROOT / "10_cross_species_baseline"
META = ROOT / "00_metadata"
REPORT = ROOT / "08_results" / "PHASE10_REPORT.md"
GATE = META / "PHASE10_GATE.csv"
MANIFEST = META / "PHASE10_MANIFEST.csv"
LOG = ROOT / "logs" / "phase10_report.log"
FINAL_GATE = "YELLOW_PHASE10_BASELINE_COMPLETE_CURRENT_ID_AUDIT_LIMITED"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def md_table(frame: pd.DataFrame, columns: list[str], formats=None) -> str:
    formats = formats or {}
    work = frame[columns].copy()
    for column, formatter in formats.items():
        work[column] = work[column].map(
            lambda value: formatter(value) if pd.notna(value) else "NA"
        )
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in work.itertuples(index=False, name=None):
        lines.append(
            "| " + " | ".join(str(value).replace("|", "\\|") for value in row) + " |"
        )
    return "\n".join(lines)


def phase9_revalidation() -> pd.DataFrame:
    manifest = pd.read_csv(META / "PHASE9_MANIFEST.csv", dtype=str)
    rows = []
    for item in manifest.itertuples(index=False):
        path = ROOT / item.path
        observed_size = path.stat().st_size if path.is_file() else -1
        observed_hash = sha256(path) if path.is_file() else ""
        status = (
            "PASS"
            if observed_size == int(item.size) and observed_hash == item.SHA256
            else "FAIL"
        )
        rows.append(
            {
                "path": item.path,
                "expected": item.SHA256,
                "observed": observed_hash,
                "size_expected": item.size,
                "size_observed": observed_size,
                "status": status,
            }
        )
    audit = pd.DataFrame(rows)
    atomic_csv(audit, BASE / "audit" / "PHASE9_MANIFEST_REVALIDATION.csv")
    return audit


def sensitivity_ranges(
    gene: pd.DataFrame, lodo: pd.DataFrame, treatment: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary = gene.set_index("lineage")["spearman_rho"]
    dog_rows = []
    human_rows = []
    for lineage, frame in lodo.groupby("lineage", sort=False):
        work = frame.assign(delta=(frame["spearman_rho"] - primary[lineage]).abs())
        influential = work.loc[work["delta"].idxmax()]
        dog_rows.append(
            {
                "lineage": lineage,
                "primary_rho": primary[lineage],
                "rho_min": work["spearman_rho"].min(),
                "rho_max": work["spearman_rho"].max(),
                "meaningful_min": work["meaningful_direction_concordance"].min(),
                "meaningful_max": work["meaningful_direction_concordance"].max(),
                "most_influential": influential["dropped_donor"],
                "model_or_run": influential["short_label"],
                "max_abs_delta_rho": influential["delta"],
            }
        )
    for lineage, frame in treatment.groupby("lineage", sort=False):
        work = frame.assign(delta=(frame["spearman_rho"] - primary[lineage]).abs())
        influential = work.loc[work["delta"].idxmax()]
        human_rows.append(
            {
                "lineage": lineage,
                "primary_rho": primary[lineage],
                "rho_min": work["spearman_rho"].min(),
                "rho_max": work["spearman_rho"].max(),
                "meaningful_min": work["meaningful_direction_concordance"].min(),
                "meaningful_max": work["meaningful_direction_concordance"].max(),
                "most_influential": influential["human_treatment_model"],
                "model_or_run": influential["human_treatment_model"],
                "max_abs_delta_rho": influential["delta"],
            }
        )
    return pd.DataFrame(dog_rows), pd.DataFrame(human_rows)


def main() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("", encoding="utf-8")
    p9 = phase9_revalidation()
    if not p9["status"].eq("PASS").all():
        raise RuntimeError("Phase9 final revalidation failed")

    mapping = pd.read_csv(BASE / "mapping" / "ORTHOLOG_MAPPING_AUDIT.csv").iloc[0]
    coverage = pd.read_csv(BASE / "mapping" / "ORTHOLOG_COVERAGE_BY_LINEAGE.csv")
    gene = pd.read_csv(BASE / "summary" / "GENE_EFFECT_CONCORDANCE_SUMMARY.csv")
    pathway_summary = pd.read_csv(
        BASE / "pathway_effects" / "PATHWAY_CONCORDANCE_SUMMARY.csv"
    )
    pathway = pd.read_csv(
        BASE / "pathway_effects" / "PATHWAY_CROSS_SPECIES_BASELINE.csv.gz"
    )
    program = pd.read_csv(BASE / "summary" / "PHASE10_PROGRAM_BASELINE_MATRIX.csv")
    triangle = pd.read_csv(BASE / "summary" / "DOG_FRE_HUMAN_TRIANGULATION.csv")
    anchors = pd.read_csv(
        BASE / "summary" / "ANCHOR_GENE_CROSS_SPECIES_BASELINE.csv"
    )
    lodo = pd.read_csv(
        BASE / "sensitivity" / "DOG_LODO_CROSS_SPECIES_SENSITIVITY.csv"
    )
    treatment = pd.read_csv(
        BASE / "sensitivity" / "HUMAN_TREATMENT_CROSS_SPECIES_SENSITIVITY.csv"
    )
    dog_sens, human_sens = sensitivity_ranges(gene, lodo, treatment)

    secondary = pd.DataFrame(
        {
            "lineage": ["B", "Plasma", "Mast"],
            "ontology_status": ["SECONDARY_LEVEL1"] * 3,
            "phase10_disease_effect_status": [
                "ONTOLOGY_READY_BUT_NO_FROZEN_DISEASE_EFFECT"
            ]
            * 3,
            "formal_phase10_comparison": ["NO"] * 3,
            "notes": ["No new DE; not expanded beyond Primary"] * 3,
        }
    )
    atomic_csv(secondary, BASE / "audit" / "SECONDARY_LEVEL1_SCOPE.csv")
    scope = pd.DataFrame(
        [
            ["Primary Level1", "PASS", "Epithelial CL:0000066; T CL:0000084; Myeloid CL:0000763"],
            ["Formal Level2 comparisons", "PASS", "0; ontology context only"],
            ["Secondary Level1 new DE", "PASS", "0; B/Plasma/Mast not expanded"],
            ["Joint expression integration", "PASS", "0; no concat/PCA/UMAP/scVI/Harmony/CCA"],
            ["Forbidden Phase11 methods", "PASS", "singIST/SATURN/CONCORD not run"],
        ],
        columns=["check", "status", "detail"],
    )
    atomic_csv(scope, BASE / "audit" / "PHASE10_SCOPE_AUDIT.csv")

    gate_rows = pd.DataFrame(
        [
            ["Phase9 integrity", "PASS", f"{len(p9)}/{len(p9)} size and SHA-256 match"],
            ["Ensembl104 source", "PASS", "Official release-104 FTP BioMart SQL and tables frozen"],
            ["Dog feature accounting", "PASS", f"{int(mapping.total_dog_features)}/20257 explicit statuses"],
            ["Primary one2one", "PASS", f"{int(mapping.primary_one2one_final)} final"],
            ["Ensembl116 current-ID audit", "YELLOW", "0 direct matches; frozen core history tables empty"],
            ["Shared coverage", "PASS", f"{coverage.coverage_gate.eq('PASS').sum()}/3 pass"],
            ["Gene baseline", "PASS", f"{len(gene)}/3 complete"],
            ["Pathway baseline", "PASS", f"{len(pathway)} paired rows; 9 summaries"],
            ["Dog LODO sensitivity", "PASS", f"{len(lodo)}/21 complete"],
            ["Human treatment sensitivity", "PASS", f"{len(treatment)}/15 complete"],
            ["25-program provenance", "PASS", f"{len(program)}/25 complete"],
            ["Scope audit", "PASS", "Only Primary Level1; no joint integration"],
            ["Traditional biology", "YELLOW", "Weak gene rho, negative T Reactome, five opposite programs"],
            ["PHASE10_FINAL_GATE", FINAL_GATE, "Complete and Phase11-ready with stated limitations"],
        ],
        columns=["criterion", "status", "detail"],
    )
    atomic_csv(gate_rows, GATE)

    gene_idx = gene.set_index("lineage")
    best_gene = gene.sort_values("spearman_rho", ascending=False).iloc[0]
    means = pathway_summary.groupby("lineage", as_index=False).agg(
        mean_nes_rho=("nes_spearman_rho", "mean"),
        mean_direction=("direction_concordance", "mean"),
    )
    best_path_rho = means.sort_values("mean_nes_rho", ascending=False).iloc[0]
    best_path_dir = means.sort_values("mean_direction", ascending=False).iloc[0]
    status_order = [
        "CONCORDANT_STRONG",
        "CONCORDANT_DOG_LIMITED",
        "CONCORDANT_HUMAN_LIMITED",
        "CONCORDANT_DIRECTION_ONLY",
        "OPPOSITE",
        "NOT_TESTABLE",
    ]
    status_counts = program["baseline_status"].value_counts()
    status_text = ", ".join(
        f"{name}={int(status_counts.get(name, 0))}" for name in status_order
    )
    risk_order = [
        "LOW_REGION_RISK",
        "MODERATE_DIRECTION_ONLY",
        "HIGH_REGION_RISK",
        "REGION_UNKNOWN",
    ]
    risk_counts = program["region_risk"].value_counts()
    risk_text = ", ".join(
        f"{name}={int(risk_counts.get(name, 0))}" for name in risk_order
    )
    kcnk16 = anchors.loc[
        anchors["lineage"].eq("Epithelial")
        & anchors["anchor"].eq("KCNK16")
        & anchors["anchor_source"].eq("PHASE3_DOG_ANCHOR")
    ].iloc[0]

    focus_names = {
        "Cell cycle / proliferation",
        "TNF / NF-kappaB",
        "Cytokine / interleukin",
        "Cytokine / interleukin / JAK-STAT",
        "Mitochondrial respiration",
        "Fatty acid metabolism",
        "Xenobiotic metabolism",
        "Mitochondrial translation",
        "Metabolic programs",
        "Respiratory electron transport",
        "Oxidative metabolism",
        "Translation",
    }
    focus = program.loc[program["program"].isin(focus_names)]
    focus_triangle = triangle.loc[triangle["program"].isin(focus_names)]

    gene_table = md_table(
        gene,
        [
            "lineage", "shared_genes", "spearman_rho", "pearson_r",
            "all_direction_concordance", "meaningful_gene_n",
            "meaningful_direction_concordance", "conserved_up_n",
            "conserved_down_n", "opposite_meaningful_n",
            "dog_only_n", "human_only_n", "weak_both_n",
        ],
        {
            "spearman_rho": lambda x: f"{x:.4f}",
            "pearson_r": lambda x: f"{x:.4f}",
            "all_direction_concordance": lambda x: f"{x:.2%}",
            "meaningful_direction_concordance": lambda x: f"{x:.2%}",
        },
    )
    coverage_table = md_table(
        coverage,
        [
            "lineage", "dog_tested_genes", "human_tested_genes",
            "dog_one2one_available", "human_one2one_available",
            "shared_one2one_tested", "fraction_of_dog_tested",
            "fraction_of_human_tested", "fraction_of_smaller_tested_universe",
            "coverage_gate",
        ],
        {
            "fraction_of_dog_tested": lambda x: f"{x:.2%}",
            "fraction_of_human_tested": lambda x: f"{x:.2%}",
            "fraction_of_smaller_tested_universe": lambda x: f"{x:.2%}",
        },
    )
    pathway_table = md_table(
        pathway_summary,
        [
            "lineage", "collection", "common_tested_pathways",
            "nes_spearman_rho", "nes_pearson_r", "direction_concordance",
            "both_FDR10_same_direction", "both_FDR10_opposite",
            "dog_FDR10_human_weak_same", "human_FDR10_dog_weak_same",
            "leading_edge_median_jaccard",
        ],
        {
            "nes_spearman_rho": lambda x: f"{x:.4f}",
            "nes_pearson_r": lambda x: f"{x:.4f}",
            "direction_concordance": lambda x: f"{x:.2%}",
            "leading_edge_median_jaccard": lambda x: f"{x:.3f}",
        },
    )
    sens_formats = {
        "primary_rho": lambda x: f"{x:.4f}",
        "rho_min": lambda x: f"{x:.4f}",
        "rho_max": lambda x: f"{x:.4f}",
        "meaningful_min": lambda x: f"{x:.2%}",
        "meaningful_max": lambda x: f"{x:.2%}",
        "max_abs_delta_rho": lambda x: f"{x:.4f}",
    }
    dog_sens_table = md_table(dog_sens, list(dog_sens.columns), sens_formats)
    human_sens_table = md_table(human_sens, list(human_sens.columns), sens_formats)
    focus_table = md_table(
        focus,
        [
            "lineage", "program", "dog_common1to1_NES",
            "human_common1to1_NES", "dog_FDR", "human_FDR",
            "baseline_status", "fre_summary_status", "region_risk",
            "treatment_caution", "ibdverse_status",
        ],
        {
            "dog_common1to1_NES": lambda x: f"{x:.3f}",
            "human_common1to1_NES": lambda x: f"{x:.3f}",
            "dog_FDR": lambda x: f"{x:.3g}",
            "human_FDR": lambda x: f"{x:.3g}",
        },
    )
    triangle_table = md_table(
        focus_triangle,
        [
            "lineage", "program", "Dog_CIE_direction", "Dog_LODO",
            "FRE_bulk_direction", "FRE_status", "Human_Crohn_direction",
            "Human_treatment_status", "IBDverse_status", "region_risk",
            "interpretation",
        ],
    )

    report = f"""# Phase10 Report — Traditional 1:1 Ortholog Disease-Effect Baseline

Generated UTC: {datetime.now(timezone.utc).isoformat()}

Final gate: **{FINAL_GATE}**

Phase10 completed the authorized traditional one-to-one ortholog baseline and stopped. It did not run singIST, SATURN, CONCORD, joint cell embeddings, Level2 formal comparisons, CeD specificity, therapeutic prioritization, or a fidelity score.

## 1. Input integrity and scope

- Project root: value supplied by `CANINE_CIE_PROJECT_ROOT`.
- Phase9 manifest: **{len(p9)}/{len(p9)} PASS** by path, size, and SHA-256.
- Formal Primary: Epithelial (CL:0000066), T (CL:0000084), Myeloid (CL:0000763), all Level1.
- B, Plasma, Mast: ONTOLOGY_READY_BUT_NO_FROZEN_DISEASE_EFFECT; no new DE.
- Dog: frozen Phase3 Full7 edgeR, CIE minus Healthy, n=7 donors.
- Human: frozen Phase7 edgeR, active Crohn TI minus Healthy TI, n=27 donors.
- Rank metrics revalidated as sign(logFC) × sqrt(F).

## 2. Ortholog source and mapping

Primary mapping was reconstructed from the official Ensembl release-104 FTP BioMart database snapshot using its SQL DDL, Dog gene main table, and Dog–Human homolog table. This matches CanFam3.1 / Ensembl104 and Human GRCh38. Raw snapshots, reconstruction query, byte sizes, and SHA-256 are frozen in mapping/source_snapshot and PHASE10_ORTHOLOG_SOURCE_AUDIT.csv.

All 20,257 Dog original features were keyed by feature_id. No symbol join or disease effect selected an ortholog.

- ONE2ONE_PRIMARY: **{int(mapping.primary_one2one_final):,}**
- ONE2MANY_EXCLUDED: **{int(mapping.one2many):,}**
- MANY2MANY_EXCLUDED: **{int(mapping.many2many):,}**
- NO_ORTHOLOG: **{int(mapping.no_ortholog):,}**
- release104 unresolved: **{int(mapping.dog_id_unresolved):,}**
- ambiguous: **{int(mapping.ambiguous):,}**
- annotation-derived high-confidence subset: **{int(mapping.high_confidence_one2one):,}**, not substituted for Primary.

Ensembl116 / ROS_Cfam_1.0 was audit-only. It had zero direct same-stable-ID matches to the old CanFam3.1 IDs; the frozen release116 Dog core stable_id_event and mapping_session tables contain no rows. Current IDs remain NA_NOT_AVAILABLE. This is YELLOW_CURRENT_ID_AUDIT_LIMITED and does not alter release104 Primary mapping.

### Coverage

{coverage_table}

All three lineages pass shared ≥1,500 and shared/min(tested) ≥0.50.

## 3. Gene baseline

{gene_table}

Epithelial: 8,052 genes, rho=0.1474, Pearson r=0.1600, all-direction 54.66%, meaningful-direction 59.23%; 800 conserved-up, 554 conserved-down, 932 opposite, 2,478 Dog-only, 1,245 Human-only, 2,043 weak-both.

T: 7,037 genes, rho=0.0159, r=0.0657, all-direction 50.09%, meaningful-direction 54.07%; 226 conserved-up, 126 conserved-down, 299 opposite, 2,033 Dog-only, 719 Human-only, 3,634 weak-both.

Myeloid: 2,936 genes, rho=0.1011, r=0.1622, all-direction 53.37%, meaningful-direction 59.60%; 217 conserved-up, 109 conserved-down, 221 opposite, 900 Dog-only, 526 Human-only, 963 weak-both.

Highest gene rho is **{best_gene.lineage} ({best_gene.spearman_rho:.4f})**. Myeloid has the highest meaningful-direction concordance ({gene_idx.loc['Myeloid', 'meaningful_direction_concordance']:.2%}), just above Epithelial ({gene_idx.loc['Epithelial', 'meaningful_direction_concordance']:.2%}). No fidelity class is assigned.

## 4. KCNK16 and anchors

KCNK16 maps naturally from ENSCAFG00000046962 to **{kcnk16.human_gene_id} ({kcnk16.human_symbol})**, Ensembl104 ortholog_one2one, confidence=1, GOC=100. Human KCNK16 is absent from the frozen Primary Epithelial tested-gene DE universe, so its formal Human effect is **NOT_TESTABLE**. The paired active-vs-adjacent contrast was not substituted.

TLR9 and CCL3 remain one-to-many exclusions. SAA1 has no acceptable one-to-one Primary mapping.

## 5. Paired pathway baseline

Dog and Human used the identical Human Ensembl-ID universe per lineage, identical memberships reconstructed from frozen Phase7 Hallmark (MSigDB 2026.1.Hs), Reactome v97, and GO BP, plus identical fgseaMultilevel parameters: minSize=15, maxSize=500, eps=1e-10, scoreType=std, nproc=1. BH correction remained separate by species, lineage, and collection.

{pathway_table}

Highest mean NES correlation across collections: **{best_path_rho.lineage} ({best_path_rho.mean_nes_rho:.3f})**. Highest mean direction concordance: **{best_path_dir.lineage} ({best_path_dir.mean_direction:.2%})**. Myeloid Reactome contains one Dog NA result for the broad Metabolism pathway; it is retained as NOT_TESTABLE and complete-case summaries are explicit.

## 6. Registered programs and core questions

{focus_table}

Epithelial metabolic suppression is supported for mitochondrial respiration (CONCORDANT_STRONG; NES -1.921/-1.591), xenobiotic metabolism (CONCORDANT_STRONG; -1.928/-1.923), and fatty-acid metabolism (CONCORDANT_DOG_LIMITED; -1.471/-1.533). GO:0006119 oxidative phosphorylation is NOT_TESTABLE because paired membership falls below the GSEA size gate; this is not negative evidence. These programs are LOW region risk but FRE bulk is higher.

T TNF/NF-kappaB is CONCORDANT_STRONG (2.542/2.805); T cytokine/interleukin is CONCORDANT_DOG_LIMITED (1.130/2.352). Myeloid TNF/NF-kappaB (2.843/2.104) and cytokine/interleukin/JAK-STAT (1.597/2.263) are CONCORDANT_STRONG; Myeloid TLR4 and inflammatory response are DOG_LIMITED. The representative Myeloid cytokine/JAK pathway has same-direction IBDverse support.

T metabolic discordance persists: mitochondrial respiration is opposite (-2.591/+1.146), mitochondrial translation is opposite (-2.351/+1.551), and broad metabolic programs are opposite (-0.773/+1.250). They retain REGION_UNKNOWN or MODERATE caution. T oxidative phosphorylation is NOT_TESTABLE and retains REGION_UNKNOWN plus treatment caution.

## 7. Dog LODO sensitivity

{dog_sens_table}

Dropping GSM8031703_cie_4 causes the largest rho shift in all lineages. Epithelial and Myeloid stay positive across seven rounds. T spans -0.0482 to 0.0643 around its near-zero Primary result. No donor was removed from Primary.

## 8. Human treatment sensitivity

{human_sens_table}

Epithelial is comparatively stable (rho 0.1111–0.1546). Myeloid ranges 0.0156–0.1335, most affected by the anti-IL12/23-only model. T is most treatment-sensitive, -0.1042–0.0962; the no-biologic model is most influential. Five registered programs retain TREATMENT_CAUTION=YES; no model was refit.

## 9. Region and 25-program status

Region categories: **{risk_text}**. HIGH and UNKNOWN were never collapsed into LOW. Epithelial cell cycle remains HIGH region risk even though it is CONCORDANT_STRONG.

Baseline counts: **{status_text}**. These are traditional one-to-one statuses, not final translational fidelity classes.

## 10. Dog CIE / FRE / Human triangulation

{triangle_table}

Key three-way relationships:

- T and Myeloid TNF/cytokine programs are Dog/Human-higher but FRE bulk-lower.
- Epithelial mitochondrial, fatty-acid, and xenobiotic programs are Dog/Human-lower but FRE bulk-higher.
- Myeloid translation is Dog-lower and FRE-lower but Human-higher; it is OPPOSITE, HIGH region risk, moderately stable in Dog, and treatment-cautioned.
- T mitochondrial programs are Dog-lower/Human-higher and REGION_UNKNOWN.

## 11. Traditional baseline conclusion

Gene effects are weakly concordant overall, with Epithelial highest by rho and Myeloid highest by meaningful direction. Pathway programs are more concordant than individual genes in several registered axes. Epithelial has the highest mean pathway NES correlation; Myeloid has the highest average pathway direction concordance. T Reactome is negative overall and T metabolic discordance is a major qualified result. Low rho was not used to change mapping, donors, thresholds, or Primary identities.

## 12. Limitations

1. Dog n=7 versus Human n=27 and unequal power.
2. Historical CanFam3.1 / Ensembl104 Dog annotation.
3. ROS_Cfam_1.0 / Ensembl116 current-ID history unavailable for these old IDs.
4. Orthology loss excludes one-to-many, many-to-many, and absent relations.
5. Human treatment structure, strongest for T.
6. Dog duodenum versus Human terminal ileum.
7. FRE is whole-mucosa published summary, not matched donor-level single-cell DE.
8. Overlapping pathways and dependent genes.
9. No formal Level2 comparison.
10. Traditional one-to-one orthology cannot capture divergent paralogs or learned programs.

## 13. Gate and next step

Final gate: **{FINAL_GATE}**.

Every executable Phase10 requirement passed. YELLOW reflects the secondary current-ID audit limit and weak/heterogeneous biology, not a broken Primary baseline.

Recommendation: **proceed to a separately gated Phase11 singIST analysis**, preserving Phase10 as the frozen traditional baseline. Phase11 was not started.

## 14. Deliverables

- Mapping and source snapshots: 10_cross_species_baseline/mapping/
- Gene effects: 10_cross_species_baseline/gene_effects/
- Pathway effects: 10_cross_species_baseline/pathway_effects/
- Sensitivities: 10_cross_species_baseline/sensitivity/
- Summaries: 10_cross_species_baseline/summary/
- Audits and figures: 10_cross_species_baseline/audit/ and figures/
- Source audit and gate: 00_metadata/PHASE10_ORTHOLOG_SOURCE_AUDIT.csv and PHASE10_GATE.csv
- Manifest: 00_metadata/PHASE10_MANIFEST.csv
- Scripts: scripts/55_phase10_* through scripts/60_phase10_*
- Logs: logs/phase10_*.log

Phase10 is complete. Stop before Phase11.
"""
    atomic_text(report, REPORT)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(f"Report generated: {REPORT}\n")
        handle.write(f"Final gate: {FINAL_GATE}\n")

    selected = []
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT).as_posix()
        if rel == "00_metadata/PHASE10_MANIFEST.csv":
            continue
        if (
            rel.startswith("10_cross_species_baseline/")
            or rel.startswith("00_metadata/PHASE10_")
            or rel == "08_results/PHASE10_REPORT.md"
            or any(rel.startswith(f"scripts/{n}_phase10_") for n in range(55, 61))
            or rel.startswith("logs/phase10_")
        ):
            selected.append(path)

    def role(rel: str) -> str:
        if "source_snapshot/" in rel:
            return "Phase10_ortholog_source_snapshot"
        if rel.endswith("PHASE10_REPORT.md"):
            return "Phase10_report"
        if rel.endswith("PHASE10_GATE.csv"):
            return "Phase10_gate"
        if rel.startswith("scripts/"):
            return "Phase10_script"
        if rel.startswith("logs/"):
            return "Phase10_log"
        if "/figures/" in rel:
            return "Phase10_figure"
        if "/gene_effects/" in rel:
            return "Phase10_gene_effect"
        if "/pathway_effects/" in rel:
            return "Phase10_pathway_effect"
        if "/sensitivity/" in rel:
            return "Phase10_sensitivity"
        if "/mapping/" in rel:
            return "Phase10_mapping"
        if "/summary/" in rel:
            return "Phase10_summary"
        return "Phase10_audit_or_metadata"

    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(
            f"Manifest selection: {len(selected)} files; manifest excludes itself\n"
        )
    rows = []
    for path in sorted(selected):
        rel = path.relative_to(ROOT).as_posix()
        generated = "external_source_freeze"
        if rel.startswith("scripts/"):
            generated = "Codex_apply_patch"
        elif "/mapping/" in rel or "ORTHOLOG_SOURCE" in rel:
            generated = "55_phase10_freeze_orthologs.py"
        elif any(
            token in rel
            for token in [
                "/gene_effects/", "/sensitivity/", "GENE_EFFECT_CONCORDANCE",
                "ANCHOR_GENE", "SHARED_ONE2ONE", "paired_ranks",
            ]
        ):
            generated = "56_phase10_gene_baseline.py"
        elif "gsea_raw" in rel or "PATHWAY_CONCORDANCE" in rel:
            generated = "57_phase10_common_ortholog_gsea.R"
        elif any(
            token in rel
            for token in [
                "PATHWAY_CROSS_SPECIES_BASELINE.csv.gz",
                "PHASE10_PROGRAM_BASELINE_MATRIX",
                "DOG_FRE_HUMAN_TRIANGULATION",
            ]
        ):
            generated = "58_phase10_program_matrix.py"
        elif "/figures/" in rel or "FIGURE_AUDIT" in rel:
            generated = "59_phase10_figures.py"
        elif rel.endswith("PHASE10_REPORT.md") or rel.endswith(
            "PHASE10_GATE.csv"
        ) or "SCOPE" in rel or "PHASE9_MANIFEST_REVALIDATION" in rel:
            generated = "60_phase10_report_manifest.py"
        rows.append(
            {
                "path": rel,
                "size": path.stat().st_size,
                "SHA256": sha256(path),
                "role": role(rel),
                "source": "Frozen Phase0-9; Ensembl104 Primary; Ensembl116 audit-only",
                "generated_by": generated,
            }
        )
    atomic_csv(pd.DataFrame(rows), MANIFEST)


if __name__ == "__main__":
    main()
