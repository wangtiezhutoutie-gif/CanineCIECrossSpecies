from __future__ import annotations

import csv
import os
import re
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
BASE = ROOT / "07_human_crohn" / "IBDverse"
RAW = BASE / "raw_summary"
META = BASE / "metadata"
SUMMARY = BASE / "summary"
FIGURES = BASE / "figures"
HUMAN_GSEA = ROOT / "07_human_crohn" / "GSE266616" / "gsea"
HUMAN_PROGRAMS = HUMAN_GSEA / "HUMAN_DISEASE_PROGRAMS.csv"

COHORT_FILES = {
    cohort: RAW / f"{cohort}_inflamed_vs_healthy_gsea_results.tsv.gz"
    for cohort in ("Discovery", "Replication", "Full")
}


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


def normalize_pathway(value: str) -> str:
    text = str(value).upper().replace("REACTOME_", "")
    replacements = {
        "NF-KAPPA-B": "NFKB", "NF-KAPPA B": "NFKB", "NFKAPPAB": "NFKB",
        "BETA": "BETA", "ALPHA": "ALPHA",
    }
    for source, target in replacements.items():
        text = text.replace(source, target)
    return re.sub(r"[^A-Z0-9]", "", text)


def level1(label: str) -> str:
    if label.startswith(("Entero", "Goblet_", "Stem_", "Tuft_")):
        return "Epithelial"
    if label.startswith("T_cells_") or label == "ILCs":
        return "T"
    if label.startswith(("Macrophages_", "Monocytes_")) or label in {"cDC1", "pDCscDC2"}:
        return "Myeloid"
    return "Other"


def load_cohort(cohort: str) -> pd.DataFrame:
    path = COHORT_FILES[cohort]
    if not path.is_file():
        raise FileNotFoundError(path)
    columns = [
        "annot_id", "qvalue_bh", "annot_coef", "leadingEdge", "reference_value", "coef_value",
        "cell_label", "signed_ranking", "gset_database", "n_gene_contained",
    ]
    frame = pd.read_csv(path, sep="\t", usecols=columns, low_memory=False)
    frame = frame.loc[frame["signed_ranking"].eq(True) & frame["gset_database"].eq("c2.cp.reactome")].copy()
    if frame.duplicated(["annot_id", "cell_label"]).any():
        raise RuntimeError(f"Duplicated signed pathway/cell identity in {cohort}")
    if not frame["reference_value"].eq("disease_statusHealthy").all() or not frame["coef_value"].eq("disease_statusCD").all():
        raise RuntimeError(f"IBDverse contrast direction changed in {cohort}")
    frame["human_level1"] = frame["cell_label"].astype(str).map(level1)
    frame["pathway_norm"] = frame["annot_id"].astype(str).map(normalize_pathway)
    return frame


def split_replication(discovery: pd.DataFrame, replication: pd.DataFrame) -> pd.DataFrame:
    keys = ["annot_id", "pathway_norm", "cell_label", "human_level1"]
    keep = keys + ["qvalue_bh", "annot_coef", "leadingEdge", "n_gene_contained"]
    merged = discovery[keep].merge(
        replication[keep], on=keys, how="inner", suffixes=("_discovery", "_replication"), validate="one_to_one"
    )
    merged["direction_discovery"] = np.where(merged["annot_coef_discovery"] > 0, "CD_higher", "CD_lower")
    merged["direction_replication"] = np.where(merged["annot_coef_replication"] > 0, "CD_higher", "CD_lower")
    merged["same_direction_between_cohorts"] = merged["direction_discovery"].eq(merged["direction_replication"])
    merged["split_cohort_replicated"] = (
        merged["qvalue_bh_discovery"].lt(0.05)
        & merged["qvalue_bh_replication"].lt(0.05)
        & merged["same_direction_between_cohorts"]
    )
    merged["replicated_direction"] = np.where(
        merged["split_cohort_replicated"], merged["direction_discovery"], "NOT_REPLICATED"
    )
    return merged


def annotate_human_programs(programs: pd.DataFrame, replicated: pd.DataFrame, ibdverse_norms: set[str]) -> pd.DataFrame:
    relevant = replicated.loc[replicated["human_level1"].isin(["Epithelial", "T", "Myeloid"])].copy()
    rows: list[dict[str, object]] = []
    for row in programs.itertuples(index=False):
        record = row._asdict()
        if row.collection != "Reactome":
            record.update({
                "ibdverse_status": "NOT_ASSESSABLE_COLLECTION_NOT_REACTOME",
                "ibdverse_mapped": "NO", "replicated_same_subtypes": 0, "replicated_opposite_subtypes": 0,
                "ibdverse_best_discovery_FDR": np.nan, "ibdverse_best_replication_FDR": np.nan,
                "ibdverse_supporting_cell_labels": "", "ibdverse_opposite_cell_labels": "",
            })
            rows.append(record)
            continue
        pathway_norm = normalize_pathway(row.program)
        if pathway_norm not in ibdverse_norms:
            record.update({
                "ibdverse_status": "NOT_ASSESSABLE_PATHWAY_NAME_MAPPING",
                "ibdverse_mapped": "NO", "replicated_same_subtypes": 0, "replicated_opposite_subtypes": 0,
                "ibdverse_best_discovery_FDR": np.nan, "ibdverse_best_replication_FDR": np.nan,
                "ibdverse_supporting_cell_labels": "", "ibdverse_opposite_cell_labels": "",
            })
            rows.append(record)
            continue
        match = relevant.loc[
            relevant["human_level1"].eq(row.lineage)
            & relevant["pathway_norm"].eq(pathway_norm)
            & relevant["split_cohort_replicated"]
        ].copy()
        expected = "CD_higher" if row.NES > 0 else "CD_lower"
        same = match.loc[match["replicated_direction"].eq(expected)]
        opposite = match.loc[~match["replicated_direction"].eq(expected)]
        if len(same) and len(opposite):
            status = "MIXED_CELL_SUBTYPE_REPLICATION"
        elif len(same):
            status = "SAME_DIRECTION_REPLICATED"
        elif len(opposite):
            status = "OPPOSITE_DIRECTION_REPLICATED"
        else:
            status = "NOT_REPLICATED_SPLIT_COHORT"
        best_pool = match if len(match) else relevant.loc[
            relevant["human_level1"].eq(row.lineage) & relevant["pathway_norm"].eq(pathway_norm)
        ]
        record.update({
            "ibdverse_status": status,
            "ibdverse_mapped": "YES",
            "replicated_same_subtypes": int(len(same)),
            "replicated_opposite_subtypes": int(len(opposite)),
            "ibdverse_best_discovery_FDR": best_pool["qvalue_bh_discovery"].min() if len(best_pool) else np.nan,
            "ibdverse_best_replication_FDR": best_pool["qvalue_bh_replication"].min() if len(best_pool) else np.nan,
            "ibdverse_supporting_cell_labels": ";".join(sorted(same["cell_label"].astype(str).unique())),
            "ibdverse_opposite_cell_labels": ";".join(sorted(opposite["cell_label"].astype(str).unique())),
        })
        rows.append(record)
    return pd.DataFrame(rows)


def plot_summary(summary: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    status_order = [
        "SAME_DIRECTION_REPLICATED", "MIXED_CELL_SUBTYPE_REPLICATION",
        "OPPOSITE_DIRECTION_REPLICATED", "NOT_REPLICATED_SPLIT_COHORT",
    ]
    colours = ["#2F855A", "#D69E2E", "#C53030", "#A0AEC0"]
    pivot = summary.pivot(index="lineage", columns="ibdverse_status", values="program_n").fillna(0)
    pivot = pivot.reindex(index=["Epithelial", "T", "Myeloid"], columns=status_order, fill_value=0)
    fig, ax = plt.subplots(figsize=(9, 5.2))
    left = np.zeros(len(pivot))
    for status, colour in zip(status_order, colours, strict=True):
        values = pivot[status].to_numpy()
        ax.barh(pivot.index, values, left=left, label=status.replace("_", " ").title(), color=colour)
        left += values
    ax.set_xlabel("Mapped significant GSE266616 Reactome programs")
    ax.set_ylabel("")
    ax.set_title("IBDverse split-cohort program replication")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.38), ncol=2, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGURES / "IBDVERSE_SPLIT_COHORT_PROGRAM_REPLICATION.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    SUMMARY.mkdir(parents=True, exist_ok=True)
    discovery = load_cohort("Discovery")
    replication = load_cohort("Replication")
    full = load_cohort("Full")
    merged = split_replication(discovery, replication)
    replicated = merged.loc[merged["split_cohort_replicated"]].copy()
    atomic_csv(replicated, SUMMARY / "IBDVERSE_SPLIT_COHORT_REPLICATED_PATHWAYS.csv")

    programs = pd.read_csv(HUMAN_PROGRAMS)
    if not programs["primary_or_secondary"].eq("PRIMARY").all():
        raise RuntimeError("Human program freeze contains non-primary rows")
    annotated = annotate_human_programs(programs, merged, set(discovery["pathway_norm"]))
    atomic_csv(annotated, SUMMARY / "IBDVERSE_HUMAN_PROGRAM_REPLICATION.csv")

    mapped = annotated.loc[annotated["ibdverse_mapped"].eq("YES")].copy()
    summary = mapped.groupby(["lineage", "ibdverse_status"], as_index=False).size().rename(columns={"size": "program_n"})
    atomic_csv(summary, SUMMARY / "IBDVERSE_HUMAN_PROGRAM_REPLICATION_SUMMARY.csv")
    plot_summary(summary)

    focus = mapped.loc[
        ~mapped["treatment_sensitivity"].str.startswith("DIRECTION_CHANGE", na=False)
    ].sort_values(
        ["ibdverse_status", "lineage", "FDR", "replicated_same_subtypes"],
        ascending=[True, True, True, False],
    )
    atomic_csv(focus, SUMMARY / "IBDVERSE_TREATMENT_STABLE_HUMAN_PROGRAMS.csv")

    same_counts = mapped.loc[mapped["ibdverse_status"].eq("SAME_DIRECTION_REPLICATED")].groupby("lineage").size()
    gate_pass = all(int(same_counts.get(lineage, 0)) > 0 for lineage in ("Epithelial", "T", "Myeloid"))
    gate = pd.DataFrame([
        {"criterion": "Discovery and Replication summary files readable", "status": "PASS", "evidence": f"rows={len(discovery)}/{len(replication)}"},
        {"criterion": "Contrast direction is CD vs Healthy", "status": "PASS", "evidence": "reference=disease_statusHealthy; coefficient=disease_statusCD"},
        {"criterion": "Split-cohort criterion frozen", "status": "PASS", "evidence": "signed ranking; FDR<0.05 in both; same sign for identical pathway and cell subtype"},
        {"criterion": "At least one same-direction replicated Human program in each primary lineage", "status": "PASS" if gate_pass else "FAIL", "evidence": ";".join(f"{x}={int(same_counts.get(x,0))}" for x in ("Epithelial", "T", "Myeloid"))},
        {"criterion": "Full 25.48-GB h5ad not downloaded", "status": "PASS", "evidence": "summary-only HTTP Range route"},
    ])
    gate["gate"] = "GREEN_LIGHTWEIGHT_REPLICATION" if gate["status"].eq("PASS").all() else "YELLOW_LIGHTWEIGHT_REPLICATION_LIMITED"
    atomic_csv(gate, META / "PHASE7D_IBDVERSE_GATE.csv")

    status_counts = mapped["ibdverse_status"].value_counts()
    report = f"""# Phase7D IBDverse lightweight replication report

Generated: {datetime.now(timezone.utc).isoformat()}

## Availability

- BioStudies S-BSST2944 API lists `DGE.zip` ({292_385_193:,} bytes); the website's legacy directory links returned HTTP 404, while the `fire` endpoint returned HTTP 200 with byte ranges.
- Only six nested summary files were extracted by HTTP Range ({sum(path.stat().st_size for path in RAW.glob('*_results.tsv.gz')):,} bytes): Discovery, Replication and Full inflamed-CD-vs-Healthy DE/GSEA.
- `Full_cohort.h5ad` ({25_476_124_776:,} bytes in the E-MTAB-16999 API) was not downloaded.

## Replication definition

IBDverse signed Reactome GSEA was considered replicated only when the identical pathway and cell subtype had FDR < 0.05 in both Discovery and Replication and the CD-vs-Healthy coefficient had the same sign. GSE266616 Human programs were then compared only within the matching broad lineage. Full-cohort results were not used to replace split-cohort replication.

## Result

- Mapped significant GSE266616 Reactome programs: {len(mapped)}.
- Same-direction replicated: {int(status_counts.get('SAME_DIRECTION_REPLICATED', 0))}.
- Mixed cell-subtype directions: {int(status_counts.get('MIXED_CELL_SUBTYPE_REPLICATION', 0))}.
- Opposite-direction replicated: {int(status_counts.get('OPPOSITE_DIRECTION_REPLICATED', 0))}.
- Not replicated under the strict split-cohort rule: {int(status_counts.get('NOT_REPLICATED_SPLIT_COHORT', 0))}.
- Same-direction replicated by lineage: {', '.join(f'{x}={int(same_counts.get(x, 0))}' for x in ('Epithelial', 'T', 'Myeloid'))}.

This is program-level external Human support, not a donor-level reanalysis of IBDverse and not formal Dog-Human fidelity statistics.
"""
    atomic_text(report, SUMMARY / "PHASE7D_IBDVERSE_REPORT.md")
    if not gate["status"].eq("PASS").all():
        raise RuntimeError("Phase7D lightweight replication gate did not fully pass")


if __name__ == "__main__":
    main()
