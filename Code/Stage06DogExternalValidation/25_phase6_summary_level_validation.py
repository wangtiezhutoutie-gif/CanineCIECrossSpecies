#!/usr/bin/env python3
"""Phase6A published summary-level external support assessment.

This is not an independent whole-transcriptome reanalysis.  It maps the frozen
25 Phase5 programs only to explicitly preregistered close/related pathways in
the published 25-row FRE KEGG summary.  It never infers a missing gene-level
rank, Hallmark/Reactome/GO result, or leading edge.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
FORWARD = ROOT / "06_pseudobulk" / "lodo" / "summary" / "PHASE5_PROGRAM_FORWARD_LIST.csv"
PARSED = ROOT / "01_raw_processed" / "dog_FRE_2026" / "reference" / "parsed_tables"
SUMMARY = ROOT / "07_external_validation" / "FRE_2026" / "summary"
LOG = ROOT / "logs" / "phase6_summary_level_validation.log"

MATRIX = SUMMARY / "PHASE6_SUMMARY_LEVEL_VALIDATION_MATRIX.csv"
HUMAN_FORWARD = SUMMARY / "PHASE6_HUMAN_FORWARD_LIST.csv"
CONTEXTUAL = SUMMARY / "FRE_CONTEXTUAL_FINDINGS.csv"
INTERPRETATION = SUMMARY / "DISCOVERY_VS_FRE_INTERPRETATION_MATRIX.csv"
COUNTS = SUMMARY / "PHASE6_SUMMARY_SUPPORT_COUNTS.csv"

EVIDENCE_LEVEL = "PUBLISHED_SUMMARY_ONLY"
BULK_WARNING = (
    "Whole-duodenal-mucosa bulk cannot identify the source lineage or validate "
    "cell-state abundance."
)
SUBTYPE_WARNING = "FRE is a food-responsive treatment subtype and is not equivalent to all canine CIE."


@dataclass(frozen=True)
class MappingSpec:
    pathways: tuple[str, ...]
    relationship: str
    support_class: str
    priority: str
    notes: str


def spec(
    pathways: tuple[str, ...],
    relationship: str,
    support_class: str,
    priority: str,
    notes: str,
) -> MappingSpec:
    return MappingSpec(pathways, relationship, support_class, priority, notes)


# Frozen before matrix construction.  The keys must match the Phase5 list
# exactly.  Empty pathway tuples are intentional NOT_ASSESSABLE decisions.
SPECS: dict[tuple[str, str], MappingSpec] = {
    ("Epithelial", "Cell cycle / proliferation"): spec(
        (), "NO_PUBLIC_MATCH", "NOT_ASSESSABLE", "MEDIUM",
        "No cell-cycle or proliferation term is present in the published 25-row KEGG summary.",
    ),
    ("Epithelial", "Ribosome / RNA biogenesis"): spec(
        ("Ribosome, eukaryotes",), "EXACT_OR_CLOSE_PATHWAY", "SUMMARY_DISCORDANT", "MEDIUM",
        "The ribosome component is directly represented and is FRE-lower; RNA biogenesis is not separately testable.",
    ),
    ("Epithelial", "Stress / UPR / MAPK"): spec(
        ("Protein processing in endoplasmic reticulum",), "BIOLOGICALLY_RELATED", "SUMMARY_PARTIAL_SUPPORT", "MEDIUM",
        "The UPR/ER-processing component is directionally concordant; MAPK is not directly represented.",
    ),
    ("Epithelial", "TNF / NF-kappaB"): spec(
        ("NF-kappa B signaling pathway",), "EXACT_OR_CLOSE_PATHWAY", "SUMMARY_DISCORDANT", "HIGH",
        "Direct bulk-level NF-kappaB discordance; the external signal cannot be assigned to epithelium.",
    ),
    ("Epithelial", "Mitochondrial respiration"): spec(
        ("Citrate cycle (TCA cycle)", "AMPK signaling pathway", "Mitophagy - animal"),
        "BIOLOGICALLY_RELATED", "SUMMARY_DISCORDANT", "HIGH",
        "CIE lineage-specific mitochondrial suppression contrasts with FRE bulk metabolic/mitochondrial enrichment.",
    ),
    ("Epithelial", "Oxidative phosphorylation"): spec(
        ("Citrate cycle (TCA cycle)", "AMPK signaling pathway", "Mitophagy - animal"),
        "BIOLOGICALLY_RELATED", "SUMMARY_DISCORDANT", "HIGH",
        "No exact oxidative-phosphorylation term is public; closely related mitochondrial programs are FRE-higher.",
    ),
    ("Epithelial", "Fatty acid metabolism"): spec(
        ("PPAR signaling pathway", "Peroxisome"), "BIOLOGICALLY_RELATED", "SUMMARY_DISCORDANT", "HIGH",
        "Related lipid oxidative programs are FRE-higher, opposite to the Phase5 CIE-lower program.",
    ),
    ("Epithelial", "Xenobiotic metabolism"): spec(
        (), "NO_PUBLIC_MATCH", "NOT_ASSESSABLE", "MEDIUM",
        "No sufficiently specific xenobiotic-metabolism term is present; PPAR/peroxisome was not used as a surrogate.",
    ),
    ("T", "TNF / NF-kappaB"): spec(
        ("NF-kappa B signaling pathway",), "EXACT_OR_CLOSE_PATHWAY", "SUMMARY_DISCORDANT", "HIGH",
        "Direct bulk-level NF-kappaB discordance; the external signal cannot be assigned to T cells.",
    ),
    ("T", "Inflammatory response"): spec(
        ("Inflammatory bowel disease (IBD)", "NF-kappa B signaling pathway"),
        "BIOLOGICALLY_RELATED", "SUMMARY_DISCORDANT", "HIGH",
        "Multiple related published immune pathways are FRE-lower; this remains whole-mucosa evidence.",
    ),
    ("T", "Cytokine / interleukin"): spec(
        ("Cytokine-cytokine receptor interaction",), "EXACT_OR_CLOSE_PATHWAY", "SUMMARY_DISCORDANT", "HIGH",
        "Published bulk cytokine-receptor enrichment is FRE-lower and cannot establish T-cell source.",
    ),
    ("T", "Stress / UPR / MAPK"): spec(
        ("Protein processing in endoplasmic reticulum",), "BIOLOGICALLY_RELATED", "SUMMARY_PARTIAL_SUPPORT", "MEDIUM",
        "The UPR/ER-processing component is directionally concordant; MAPK is not directly represented.",
    ),
    ("T", "Oxidative phosphorylation"): spec(
        ("Citrate cycle (TCA cycle)", "AMPK signaling pathway", "Mitophagy - animal"),
        "BIOLOGICALLY_RELATED", "SUMMARY_DISCORDANT", "HIGH",
        "Related mitochondrial programs are FRE-higher, opposite to the Phase5 CIE-lower T-cell program.",
    ),
    ("T", "Mitochondrial respiration"): spec(
        ("Citrate cycle (TCA cycle)", "AMPK signaling pathway", "Mitophagy - animal"),
        "BIOLOGICALLY_RELATED", "SUMMARY_DISCORDANT", "HIGH",
        "Related mitochondrial programs are FRE-higher, opposite to the Phase5 CIE-lower T-cell program.",
    ),
    ("T", "Mitochondrial translation"): spec(
        (), "NO_PUBLIC_MATCH", "NOT_ASSESSABLE", "HIGH",
        "The general eukaryotic ribosome term was not treated as mitochondrial translation.",
    ),
    ("T", "Metabolic programs"): spec(
        ("Citrate cycle (TCA cycle)", "Carbon metabolism", "AMPK signaling pathway", "Peroxisome"),
        "BIOLOGICALLY_RELATED", "SUMMARY_DISCORDANT", "HIGH",
        "A coherent set of published FRE bulk metabolic pathways is enriched in the opposite direction.",
    ),
    ("Myeloid", "TNF / NF-kappaB"): spec(
        ("NF-kappa B signaling pathway",), "EXACT_OR_CLOSE_PATHWAY", "SUMMARY_DISCORDANT", "HIGH",
        "Direct bulk-level NF-kappaB discordance; the external signal cannot be assigned to myeloid cells.",
    ),
    ("Myeloid", "Inflammatory response"): spec(
        ("Inflammatory bowel disease (IBD)", "NF-kappa B signaling pathway"),
        "BIOLOGICALLY_RELATED", "SUMMARY_DISCORDANT", "HIGH",
        "Multiple related published immune pathways are FRE-lower; this remains whole-mucosa evidence.",
    ),
    ("Myeloid", "TLR4"): spec(
        (), "NO_PUBLIC_MATCH", "NOT_ASSESSABLE", "HIGH",
        "No TLR4-specific published term is available; general immune pathways were not used as substitutes.",
    ),
    ("Myeloid", "Cytokine / interleukin / JAK-STAT"): spec(
        ("Cytokine-cytokine receptor interaction",), "EXACT_OR_CLOSE_PATHWAY", "SUMMARY_DISCORDANT", "HIGH",
        "The cytokine/interleukin component is bulk-discordant; JAK-STAT itself remains NOT_ASSESSABLE.",
    ),
    ("Myeloid", "Stress response"): spec(
        ("Protein processing in endoplasmic reticulum",), "BIOLOGICALLY_RELATED", "SUMMARY_PARTIAL_SUPPORT", "MEDIUM",
        "ER protein-processing provides related same-direction evidence but does not cover the full stress program.",
    ),
    ("Myeloid", "Apoptosis"): spec(
        (), "NO_PUBLIC_MATCH", "NOT_ASSESSABLE", "MEDIUM",
        "Autophagy was not treated as apoptosis; no sufficiently specific public match is available.",
    ),
    ("Myeloid", "Respiratory electron transport"): spec(
        ("Citrate cycle (TCA cycle)", "AMPK signaling pathway", "Mitophagy - animal"),
        "BIOLOGICALLY_RELATED", "SUMMARY_DISCORDANT", "HIGH",
        "No exact electron-transport term is public; related mitochondrial programs are FRE-higher.",
    ),
    ("Myeloid", "Oxidative metabolism"): spec(
        ("Citrate cycle (TCA cycle)", "Carbon metabolism", "AMPK signaling pathway", "Peroxisome"),
        "BIOLOGICALLY_RELATED", "SUMMARY_DISCORDANT", "HIGH",
        "Related oxidative/metabolic programs are FRE-higher, opposite to the Phase5 CIE-lower program.",
    ),
    ("Myeloid", "Translation"): spec(
        ("Ribosome, eukaryotes",), "EXACT_OR_CLOSE_PATHWAY", "SUMMARY_DIRECTIONAL_SUPPORT", "MEDIUM",
        "The published ribosome term is FRE-lower, directionally concordant with CIE-lower translation.",
    ),
}


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def log(message: str) -> None:
    line = f"[{now_utc()}] {message}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(line + "\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def human_status(stability: str, support_class: str) -> str:
    if stability == "ROBUST":
        return "YES"
    if support_class == "SUMMARY_DISCORDANT":
        return "EXPLORATORY_OR_QUALIFIED"
    return "YES_WITH_QUALIFICATION"


def external_status(support_class: str) -> str:
    return {
        "SUMMARY_DIRECTIONAL_SUPPORT": "DIRECTIONALLY_SUPPORTIVE_FRE_BULK",
        "SUMMARY_PARTIAL_SUPPORT": "PARTIAL_FRE_BULK_SUPPORT",
        "SUMMARY_DISCORDANT": "DISCORDANT_FRE_BULK",
        "NOT_ASSESSABLE": "NOT_ASSESSABLE_FRE_SUMMARY",
    }[support_class]


def build_validation_matrix() -> list[dict[str, object]]:
    forward = read_csv(FORWARD)
    pathway_rows = read_csv(PARSED / "Supplementary_Table_S4_KEGG_GSEA.csv")
    pathway_by_name = {row["Pathway"]: row for row in pathway_rows}
    observed_keys = {(row["lineage"], row["program"]) for row in forward}
    if len(forward) != 25 or observed_keys != set(SPECS):
        missing = sorted(set(SPECS) - observed_keys)
        unexpected = sorted(observed_keys - set(SPECS))
        raise RuntimeError(
            f"Frozen target mismatch: rows={len(forward)}, missing={missing}, unexpected={unexpected}"
        )

    matrix: list[dict[str, object]] = []
    for row in forward:
        key = (row["lineage"], row["program"])
        mapping = SPECS[key]
        published = [pathway_by_name[name] for name in mapping.pathways]
        directions = {p["direction"] for p in published}
        if len(directions) > 1:
            raise RuntimeError(f"Mixed FRE directions for {key}: {directions}")
        published_direction = next(iter(directions)) if directions else "NA"
        if mapping.support_class == "SUMMARY_DISCORDANT":
            expected_opposite = (
                (row["full7_direction"] == "CIE_higher" and published_direction == "FRE_lower")
                or (row["full7_direction"] == "CIE_lower" and published_direction == "FRE_higher")
            )
            if not expected_opposite:
                raise RuntimeError(f"Discordance direction assertion failed for {key}")
        if mapping.support_class in {"SUMMARY_DIRECTIONAL_SUPPORT", "SUMMARY_PARTIAL_SUPPORT"}:
            expected_same = (
                (row["full7_direction"] == "CIE_higher" and published_direction == "FRE_higher")
                or (row["full7_direction"] == "CIE_lower" and published_direction == "FRE_lower")
            )
            if not expected_same:
                raise RuntimeError(f"Support direction assertion failed for {key}")

        matrix.append(
            {
                "lineage_origin": row["lineage"],
                "program": row["program"],
                "phase5_direction": row["full7_direction"],
                "phase5_stability": row["stability_class"],
                "phase5_forward_status": row["forward_to_phase6"],
                "published_fre_pathway": ";".join(p["Pathway"] for p in published) or "NA",
                "published_pathway_id": ";".join(p["Pathway_ID"] for p in published) or "NA",
                "published_database": "KEGG" if published else "NA",
                "published_NES": ";".join(p["NES"] for p in published) or "NA",
                "published_padj": ";".join(p["padj"] for p in published) or "NA",
                "published_direction": published_direction,
                "relationship_type": mapping.relationship,
                "summary_support_class": mapping.support_class,
                "evidence_level": EVIDENCE_LEVEL,
                "external_canine_status": external_status(mapping.support_class),
                "bulk_resolution_warning": BULK_WARNING,
                "subtype_warning": SUBTYPE_WARNING,
                "forward_to_human": human_status(row["stability_class"], mapping.support_class),
                "human_priority": mapping.priority,
                "notes": mapping.notes,
            }
        )

    fields = [
        "lineage_origin",
        "program",
        "phase5_direction",
        "phase5_stability",
        "phase5_forward_status",
        "published_fre_pathway",
        "published_pathway_id",
        "published_database",
        "published_NES",
        "published_padj",
        "published_direction",
        "relationship_type",
        "summary_support_class",
        "evidence_level",
        "external_canine_status",
        "bulk_resolution_warning",
        "subtype_warning",
        "forward_to_human",
        "human_priority",
        "notes",
    ]
    write_csv(MATRIX, fields, matrix)
    return matrix


def build_counts(matrix: list[dict[str, object]]) -> dict[str, int]:
    order = [
        "SUMMARY_DIRECTIONAL_SUPPORT",
        "SUMMARY_PARTIAL_SUPPORT",
        "SUMMARY_DISCORDANT",
        "NOT_ASSESSABLE",
    ]
    counts = {name: sum(row["summary_support_class"] == name for row in matrix) for name in order}
    expected = {
        "SUMMARY_DIRECTIONAL_SUPPORT": 1,
        "SUMMARY_PARTIAL_SUPPORT": 3,
        "SUMMARY_DISCORDANT": 16,
        "NOT_ASSESSABLE": 5,
    }
    if counts != expected or sum(counts.values()) != 25:
        raise RuntimeError(f"Unexpected support-class counts: {counts}")
    rows = [
        {
            "summary_support_class": name,
            "program_n": counts[name],
            "denominator": 25,
            "percent": f"{100 * counts[name] / 25:.1f}",
            "evidence_level": EVIDENCE_LEVEL,
        }
        for name in order
    ]
    write_csv(
        COUNTS,
        ["summary_support_class", "program_n", "denominator", "percent", "evidence_level"],
        rows,
    )
    return counts


def build_human_forward(matrix: list[dict[str, object]]) -> None:
    rows = []
    for row in matrix:
        support = str(row["summary_support_class"])
        if support == "SUMMARY_DISCORDANT":
            reason = (
                "Phase5 stability is retained; FRE whole-mucosa direction differs, making independent Human testing high-value."
            )
        elif support == "NOT_ASSESSABLE":
            reason = "Phase5 evidence is retained because the public FRE summary cannot test this program."
        elif support == "SUMMARY_PARTIAL_SUPPORT":
            reason = "Related same-direction FRE bulk evidence is limited to part of the frozen discovery program."
        else:
            reason = "Published FRE bulk summary is directionally concordant, with summary-only evidence limitations."
        rows.append(
            {
                "lineage": row["lineage_origin"],
                "program": row["program"],
                "phase5_stability": row["phase5_stability"],
                "CIE_direction": row["phase5_direction"],
                "FRE_summary_status": support,
                "FRE_direction": row["published_direction"],
                "evidence_level": EVIDENCE_LEVEL,
                "external_canine_status": row["external_canine_status"],
                "forward_to_human": row["forward_to_human"],
                "priority": row["human_priority"],
                "reason": reason,
                "notes": row["notes"],
            }
        )
    write_csv(
        HUMAN_FORWARD,
        [
            "lineage",
            "program",
            "phase5_stability",
            "CIE_direction",
            "FRE_summary_status",
            "FRE_direction",
            "evidence_level",
            "external_canine_status",
            "forward_to_human",
            "priority",
            "reason",
            "notes",
        ],
        rows,
    )


def build_contextual_findings() -> None:
    pathways = read_csv(PARSED / "Supplementary_Table_S4_KEGG_GSEA.csv")
    pathway = next(row for row in pathways if row["Pathway"] == "Bile secretion")
    de_rows = read_csv(PARSED / "Main_Table1_targeted_bile_acid_DE.csv")
    means = read_csv(PARSED / "Supplementary_Table_S6_targeted_bile_acid_group_means.csv")
    means_by_symbol = {row["Gene symbol"]: row for row in means}
    output: list[dict[str, object]] = [
        {
            "finding_type": "PATHWAY",
            "feature": pathway["Pathway"],
            "feature_id": pathway["Pathway_ID"],
            "published_source": "Supplementary Table S4",
            "published_database": "KEGG",
            "published_NES_or_log2FC": pathway["NES"],
            "published_padj": pathway["padj"],
            "published_direction": pathway["direction"],
            "mean_raw_counts_CTRL": "NA",
            "mean_raw_counts_FRE": "NA",
            "mean_TPM_CTRL": "NA",
            "mean_TPM_FRE": "NA",
            "context_class": "CONTEXTUAL_ONLY",
            "evidence_level": EVIDENCE_LEVEL,
            "relationship_to_phase5": "Not used to select or reclassify Phase5 targets",
            "notes": "FRE subtype-specific bile-acid context; not a Phase5 preregistered program.",
        }
    ]
    for row in de_rows:
        gene = row["Gene"]
        mean = means_by_symbol[gene]
        direction = "FRE_higher" if float(row["log2FC"]) > 0 else "FRE_lower"
        output.append(
            {
                "finding_type": "SELECTED_GENE",
                "feature": gene,
                "feature_id": mean["Ensembl ID"],
                "published_source": "Main Table 1 and Supplementary Table S6",
                "published_database": "NA",
                "published_NES_or_log2FC": row["log2FC"],
                "published_padj": row["padj"],
                "published_direction": direction,
                "mean_raw_counts_CTRL": mean["Mean raw counts (CTRL)"],
                "mean_raw_counts_FRE": mean["Mean raw counts (FRE)"],
                "mean_TPM_CTRL": mean["Mean TPM (CTRL)"],
                "mean_TPM_FRE": mean["Mean TPM (FRE)"],
                "context_class": "CONTEXTUAL_ONLY",
                "evidence_level": EVIDENCE_LEVEL,
                "relationship_to_phase5": "Not used to select or reclassify Phase5 targets",
                "notes": "Hypothesis-selected bile-acid gene; not whole-transcriptome external validation.",
            }
        )
    write_csv(
        CONTEXTUAL,
        [
            "finding_type",
            "feature",
            "feature_id",
            "published_source",
            "published_database",
            "published_NES_or_log2FC",
            "published_padj",
            "published_direction",
            "mean_raw_counts_CTRL",
            "mean_raw_counts_FRE",
            "mean_TPM_CTRL",
            "mean_TPM_FRE",
            "context_class",
            "evidence_level",
            "relationship_to_phase5",
            "notes",
        ],
        output,
    )


def build_interpretation_matrix() -> None:
    rows = [
        {
            "finding": "Cross-lineage mitochondrial/metabolic programs",
            "discovery_result": "CIE-lower across Epithelial, T and Myeloid lineages",
            "fre_published_result": "FRE-higher TCA, carbon metabolism, AMPK, mitophagy and peroxisome",
            "relationship": "SUMMARY_DISCORDANT",
            "possible_explanation": "FRE subtype biology; whole-mucosa composition; cell resolution; cohort/source/breed; pathway/rank differences",
            "can_disprove_discovery": "NO",
            "human_followup_priority": "HIGH",
            "evidence_level": EVIDENCE_LEVEL,
        },
        {
            "finding": "TNF/NF-kappaB",
            "discovery_result": "CIE-higher in T and Myeloid; Epithelial moderately stable",
            "fre_published_result": "Bulk FRE-lower NF-kappa B signaling (NES=-1.886; padj=0.005)",
            "relationship": "BULK_SUMMARY_DISCORDANT",
            "possible_explanation": "Lineage-specific activation can be diluted or reversed in mixed whole-mucosal expression",
            "can_disprove_discovery": "NO",
            "human_followup_priority": "HIGH",
            "evidence_level": EVIDENCE_LEVEL,
        },
        {
            "finding": "Cytokine/interleukin",
            "discovery_result": "CIE-higher T and Myeloid cytokine/interleukin programs",
            "fre_published_result": "Bulk FRE-lower cytokine-cytokine receptor interaction (NES=-2.372; padj=0.005)",
            "relationship": "BULK_SUMMARY_DISCORDANT",
            "possible_explanation": "Subtype and whole-mucosa composition differences; external result has no lineage assignment",
            "can_disprove_discovery": "NO",
            "human_followup_priority": "HIGH",
            "evidence_level": EVIDENCE_LEVEL,
        },
        {
            "finding": "Myeloid TLR4",
            "discovery_result": "ROBUST CIE-higher",
            "fre_published_result": "No TLR4-specific published summary term",
            "relationship": "NOT_ASSESSABLE",
            "possible_explanation": "Insufficient public pathway granularity",
            "can_disprove_discovery": "NO",
            "human_followup_priority": "HIGH",
            "evidence_level": EVIDENCE_LEVEL,
        },
        {
            "finding": "Myeloid JAK-STAT / STAT3 arm",
            "discovery_result": "ROBUST CIE-higher cytokine/interleukin/JAK-STAT; STAT3 arm only, not IL-23-specific",
            "fre_published_result": "General cytokine-receptor pathway is FRE-lower; no direct JAK-STAT or IL-23-specific term",
            "relationship": "CYTOKINE_COMPONENT_DISCORDANT_JAK_STAT_NOT_ASSESSABLE",
            "possible_explanation": "Published summary covers the upstream cytokine category but cannot resolve JAK-STAT components or lineage",
            "can_disprove_discovery": "NO",
            "human_followup_priority": "HIGH",
            "evidence_level": EVIDENCE_LEVEL,
        },
        {
            "finding": "Epithelial KCNK16",
            "discovery_result": "Full7 logFC=-3.772; 7/7 LODO same direction; ROBUST",
            "fre_published_result": "Absent from selected-gene table; volcano has no reusable exact statistic",
            "relationship": "NOT_ASSESSABLE",
            "possible_explanation": "Complete external gene-level statistics are not public",
            "can_disprove_discovery": "NO",
            "human_followup_priority": "HIGH",
            "evidence_level": EVIDENCE_LEVEL,
        },
        {
            "finding": "Stress/UPR",
            "discovery_result": "CIE-higher stress/UPR programs across lineages",
            "fre_published_result": "FRE-higher protein processing in endoplasmic reticulum",
            "relationship": "SUMMARY_PARTIAL_SUPPORT",
            "possible_explanation": "Published KEGG term captures the UPR/ER-processing component but not the full MAPK/stress programs",
            "can_disprove_discovery": "NO",
            "human_followup_priority": "MEDIUM",
            "evidence_level": EVIDENCE_LEVEL,
        },
        {
            "finding": "Ribosome/translation",
            "discovery_result": "Epithelial ribosome/RNA biogenesis CIE-higher; Myeloid translation CIE-lower",
            "fre_published_result": "FRE-lower eukaryotic ribosome",
            "relationship": "LINEAGE_DEPENDENT_MIXED",
            "possible_explanation": "The same bulk term is discordant with Epithelial but directionally supportive of Myeloid translation",
            "can_disprove_discovery": "NO",
            "human_followup_priority": "MEDIUM",
            "evidence_level": EVIDENCE_LEVEL,
        },
        {
            "finding": "Bile-acid axis",
            "discovery_result": "Not a Phase5 preregistered core program",
            "fre_published_result": "Bile secretion FRE-higher plus coordinated selected-gene changes",
            "relationship": "CONTEXTUAL_ONLY",
            "possible_explanation": "Hypothesis-driven FRE subtype-specific transcriptional context",
            "can_disprove_discovery": "NO",
            "human_followup_priority": "CONTEXTUAL",
            "evidence_level": EVIDENCE_LEVEL,
        },
    ]
    write_csv(
        INTERPRETATION,
        [
            "finding",
            "discovery_result",
            "fre_published_result",
            "relationship",
            "possible_explanation",
            "can_disprove_discovery",
            "human_followup_priority",
            "evidence_level",
        ],
        rows,
    )


def main() -> int:
    SUMMARY.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("", encoding="utf-8")
    log("Phase6A Option A started: published summary-level evidence only")
    matrix = build_validation_matrix()
    counts = build_counts(matrix)
    build_human_forward(matrix)
    build_contextual_findings()
    build_interpretation_matrix()
    if sum(1 for row in matrix if row["forward_to_human"] in {"YES", "YES_WITH_QUALIFICATION", "EXPLORATORY_OR_QUALIFIED"}) != 25:
        raise RuntimeError("All 25 frozen targets must remain represented in the Human forward list")
    log(f"Completed 25-program classification: {counts}")
    log("All 25 programs retained for Human testing with explicit status/priority; no Human data accessed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        log(f"FATAL: {type(exc).__name__}: {exc}")
        raise
