#!/usr/bin/env python3
"""Generate Phase9 annotation-structure and donor-coverage figures."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P9 = ROOT / "09_cell_ontology"
FIG = P9 / "figures"
COLORS = {
    "Epithelial": "#4C78A8",
    "T": "#E45756",
    "Myeloid": "#F2CF5B",
    "B": "#72B7B2",
    "Plasma": "#54A24B",
    "Mast": "#B279A2",
    "Stromal": "#FF9DA6",
    "Endothelial": "#9D755D",
    "Uncertain": "#BAB0AC",
}


def save(fig, name: str) -> None:
    fig.savefig(FIG / name, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def mapping_plot(table: pd.DataFrame, species: str, level1: str, level2: str, count: str, name: str) -> None:
    grouped = (
        table.groupby([level1, level2], as_index=False)[count].sum()
        .sort_values(count, ascending=True)
    )
    fig, ax = plt.subplots(figsize=(10, max(6, len(grouped) * 0.30)))
    y = np.arange(len(grouped))
    ax.barh(y, grouped[count], color=[COLORS.get(x, "#888888") for x in grouped[level1]])
    ax.set_yticks(y, grouped[level2])
    ax.set_xlabel("Cells")
    ax.set_title(f"{species} native labels mapped to frozen Cell Ontology")
    ax.grid(axis="x", alpha=0.2)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.text(0.01, 0.005, "Bar length is abundance, not cross-species similarity.", fontsize=8, color="#555555")
    fig.tight_layout()
    save(fig, name)


def hierarchy_plot(master: pd.DataFrame) -> None:
    level1_order = ["Epithelial", "T", "Myeloid", "B", "Plasma", "Mast", "Stromal", "Endothelial", "Uncertain"]
    level1_order = [x for x in level1_order if x in set(master["harmonized_level1"])]
    fig, ax = plt.subplots(figsize=(15, 9))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    root_x, root_y = 0.5, 0.93
    ax.add_patch(FancyBboxPatch((root_x - 0.08, root_y - 0.035), 0.16, 0.07, boxstyle="round,pad=0.01", facecolor="#1F2937", edgecolor="none"))
    ax.text(root_x, root_y, "Frozen CL\ncomparison layer", ha="center", va="center", color="white", weight="bold")
    xs = np.linspace(0.06, 0.94, len(level1_order))
    for x, lineage in zip(xs, level1_order):
        ax.plot([root_x, x], [root_y - 0.04, 0.76], color="#CBD5E1", lw=1.2, zorder=0)
        ax.add_patch(FancyBboxPatch((x - 0.048, 0.70), 0.096, 0.08, boxstyle="round,pad=0.008", facecolor=COLORS[lineage], edgecolor="white"))
        ax.text(x, 0.74, lineage, ha="center", va="center", fontsize=9, color="white", weight="bold")
        fine = master[(master["harmonized_level1"] == lineage) & (master["harmonized_level2"] != "NA_LEVEL1")]
        fine = fine[~fine["harmonized_level2"].str.contains("state|uncertain", case=False, na=False)].head(5)
        for idx, (_, row) in enumerate(fine.iterrows()):
            yy = 0.61 - idx * 0.105
            ax.plot([x, x], [0.70, yy + 0.03], color=COLORS[lineage], lw=1, alpha=0.6)
            ax.add_patch(FancyBboxPatch((x - 0.052, yy - 0.027), 0.104, 0.054, boxstyle="round,pad=0.005", facecolor="white", edgecolor=COLORS[lineage], linewidth=1))
            label = str(row["harmonized_level2"])
            if len(label) > 22:
                label = label[:20] + "…"
            ax.text(x, yy, label, ha="center", va="center", fontsize=7)
    ax.set_title("Cross-species shared ontology hierarchy (annotation structure only)", fontsize=16, pad=18)
    fig.text(0.5, 0.02, "Edges represent ontology/project hierarchy; they do not encode expression similarity.", ha="center", fontsize=9, color="#555555")
    save(fig, "CROSS_SPECIES_SHARED_ONTOLOGY_HIERARCHY.png")


def primary_overview(master: pd.DataFrame) -> None:
    primary = master[(master["harmonized_level2"] == "NA_LEVEL1") & (master["harmonized_level1"].isin(["Epithelial", "T", "Myeloid"]))]
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, (_, row) in zip(axes, primary.sort_values("harmonized_level1").iterrows()):
        lineage = row["harmonized_level1"]
        ax.axis("off")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        for x, title, subtitle in [
            (0.17, "Dog", f"{int(row['dog_cells']):,} cells\n{int(row['dog_donors'])}/7 donors"),
            (0.50, lineage, f"{row['harmonized_cl_id']}\n{row['harmonized_cl_label']}"),
            (0.83, "Human TI", f"{int(row['human_cells']):,} cells\n{int(row['human_primary_donors'])}/27 Primary donors"),
        ]:
            ax.add_patch(FancyBboxPatch((x - 0.14, 0.37), 0.28, 0.27, boxstyle="round,pad=0.02", facecolor=COLORS[lineage] if x == 0.5 else "#F8FAFC", edgecolor=COLORS[lineage], linewidth=1.5))
            ax.text(x, 0.54, title, ha="center", va="center", weight="bold", color="white" if x == 0.5 else "#111827")
            ax.text(x, 0.43, subtitle, ha="center", va="center", fontsize=7, color="white" if x == 0.5 else "#475569")
        ax.annotate("", xy=(0.36, 0.505), xytext=(0.31, 0.505), arrowprops={"arrowstyle": "->", "color": "#64748B"})
        ax.annotate("", xy=(0.64, 0.505), xytext=(0.69, 0.505), arrowprops={"arrowstyle": "<-", "color": "#64748B"})
        ax.text(0.5, 0.20, "PRIMARY_LEVEL1 • EXACT_ONTOLOGY_MATCH", ha="center", fontsize=8, weight="bold")
    fig.suptitle("Primary Phase10 ontology layer: independently mapped to the same CL class", fontsize=15)
    save(fig, "PRIMARY_ETM_MAPPING_OVERVIEW.png")


def compatibility_counts(compat: pd.DataFrame) -> None:
    order = ["EXACT_ONTOLOGY_MATCH", "SHARED_PARENT_MATCH", "BROAD_LINEAGE_ONLY", "SPECIES_SPECIFIC_OR_UNMATCHED", "UNCERTAIN"]
    counts = compat["mapping_class"].value_counts().reindex(order, fill_value=0)
    colors = ["#2E8B57", "#4C78A8", "#F2CF5B", "#E45756", "#8C8C8C"]
    fig, ax = plt.subplots(figsize=(10, 5.5))
    bars = ax.bar(np.arange(len(counts)), counts.values, color=colors)
    ax.set_xticks(np.arange(len(counts)), [x.replace("_", "\n") for x in counts.index], fontsize=8)
    ax.set_ylabel("Compatibility classes")
    ax.set_title("Comparison-ready versus unmatched ontology classes")
    ax.bar_label(bars, padding=3)
    ax.grid(axis="y", alpha=0.2)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    save(fig, "COMPARISON_READY_VS_UNMATCHED_COUNTS.png")


def donor_coverage(coverage: pd.DataFrame) -> None:
    x = coverage[coverage["harmonized_cell_type"].str.endswith("::NA_LEVEL1")].copy()
    x["lineage"] = x["harmonized_cell_type"].str.split("::").str[0]
    order = ["Epithelial", "T", "Myeloid", "B", "Plasma", "Mast", "Stromal", "Endothelial", "Uncertain"]
    x["lineage"] = pd.Categorical(x["lineage"], order, ordered=True)
    x = x.sort_values("lineage")
    columns = ["dog_healthy_donors", "dog_cie_donors", "human_healthy_donors", "human_active_crohn_donors"]
    labels = ["Dog Healthy (max 3)", "Dog CIE (max 4)", "Human Healthy Primary (max 16)", "Human Active Crohn Primary (max 11)"]
    fig, ax = plt.subplots(figsize=(13, 6))
    positions = np.arange(len(x))
    width = 0.19
    palette = ["#72B7B2", "#E45756", "#4C78A8", "#F58518"]
    for i, (column, label, color) in enumerate(zip(columns, labels, palette)):
        ax.bar(positions + (i - 1.5) * width, x[column], width=width, label=label, color=color)
    ax.set_xticks(positions, x["lineage"].astype(str), rotation=30, ha="right")
    ax.set_ylabel("Donors with cells")
    ax.set_title("Condition-level donor coverage of Harmonized Level1 classes")
    ax.legend(frameon=False, ncol=2)
    ax.grid(axis="y", alpha=0.2)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    save(fig, "CELLTYPE_DONOR_COVERAGE.png")


def main() -> int:
    FIG.mkdir(parents=True, exist_ok=True)
    dog = pd.read_csv(P9 / "dog" / "DOG_CELLTYPE_ONTOLOGY_MAPPING.csv")
    human = pd.read_csv(P9 / "human" / "HUMAN_CELLTYPE_ONTOLOGY_MAPPING.csv")
    master = pd.read_csv(P9 / "harmonized" / "CROSS_SPECIES_CELL_ONTOLOGY_MASTER.csv")
    compat = pd.read_csv(P9 / "harmonized" / "ONTOLOGY_COMPATIBILITY_MATRIX.csv")
    coverage = pd.read_csv(P9 / "audit" / "CELLTYPE_CROSS_SPECIES_COVERAGE.csv")
    mapping_plot(dog, "Dog", "dog_level1", "dog_level2", "n_cells", "DOG_NATIVE_TO_CL_MAPPING.png")
    mapping_plot(human, "Human TI", "human_level1", "human_level2", "n_cells", "HUMAN_NATIVE_TO_CL_MAPPING.png")
    hierarchy_plot(master)
    primary_overview(master)
    compatibility_counts(compat)
    donor_coverage(coverage)
    files = sorted(path.name for path in FIG.glob("*.png"))
    (ROOT / "logs" / "phase9_figures.log").write_text(
        "\n".join(f"generated={name}" for name in files) + "\nstatus=PASS\n",
        encoding="utf-8",
    )
    print("\n".join(files))
    return 0


if __name__ == "__main__":
    sys.exit(main())
