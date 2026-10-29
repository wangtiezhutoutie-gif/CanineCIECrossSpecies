from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
BASE = ROOT / "10_cross_species_baseline"
GENE = BASE / "gene_effects"
PATHWAY = BASE / "pathway_effects"
SUMMARY = BASE / "summary"
SENS = BASE / "sensitivity"
FIG = BASE / "figures"
LOG = ROOT / "logs" / "phase10_figures.log"
LINEAGES = ("Epithelial", "T", "Myeloid")
COLLECTIONS = ("Hallmark", "Reactome", "GO_BP")
COLORS = {
    "CONSERVED_UP": "#B33A3A",
    "CONSERVED_DOWN": "#2E63A3",
    "OPPOSITE_MEANINGFUL": "#7B4FA3",
    "DOG_ONLY_EFFECT": "#D9902F",
    "HUMAN_ONLY_EFFECT": "#3A9D8F",
    "WEAK_BOTH": "#B8B8B8",
}
ANCHOR_SYMBOLS = {
    "Epithelial": {"KCNK16", "MMP7", "GAL3ST2", "CYP1B1"},
    "T": {"GAL3ST2", "TLR9", "TMEM255A", "PTPRT"},
    "Myeloid": {"SERPINE1", "CCL3", "ADAMDEC1", "MMP9", "SAA1"},
}


def save(fig: plt.Figure, path: Path, width: float, height: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.stem + ".", suffix=".png", dir=path.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        fig.set_size_inches(width, height)
        fig.savefig(tmp, dpi=220, bbox_inches="tight", facecolor="white")
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)
    plt.close(fig)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")
    outputs: list[dict[str, object]] = []
    effects: dict[str, pd.DataFrame] = {}
    for lineage in LINEAGES:
        path = GENE / f"{lineage.upper()}_DOG_HUMAN_GENE_EFFECTS.csv.gz"
        frame = pd.read_csv(path)
        effects[lineage] = frame
        fig, ax = plt.subplots()
        hb = ax.hexbin(
            frame["dog_logFC"],
            frame["human_logFC"],
            gridsize=65,
            mincnt=1,
            bins="log",
            cmap="viridis",
            linewidths=0,
        )
        fig.colorbar(hb, ax=ax, label="log10 bin count")
        for value, style, alpha in [
            (0, "-", 0.75),
            (0.25, "--", 0.45),
            (-0.25, "--", 0.45),
        ]:
            ax.axvline(value, color="#444444", linestyle=style, linewidth=0.8, alpha=alpha)
            ax.axhline(value, color="#444444", linestyle=style, linewidth=0.8, alpha=alpha)
        labels = frame.loc[frame["dog_symbol"].isin(ANCHOR_SYMBOLS[lineage])]
        for row in labels.itertuples(index=False):
            ax.scatter(row.dog_logFC, row.human_logFC, s=24, color="#E34A33", zorder=4)
            ax.annotate(
                row.dog_symbol,
                (row.dog_logFC, row.human_logFC),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=7,
            )
        rho = frame[["dog_logFC", "human_logFC"]].corr(method="spearman").iloc[0, 1]
        ax.set(
            title=f"{lineage}: Dog CIE vs Human active Crohn gene effects",
            xlabel="Dog logFC (CIE − Healthy)",
            ylabel="Human logFC (active Crohn TI − Healthy TI)",
        )
        ax.text(
            0.02,
            0.98,
            f"shared 1:1 genes = {len(frame):,}\nSpearman ρ = {rho:.3f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "#CCCCCC"},
        )
        out = FIG / f"{lineage}_GENE_EFFECT_SCATTER.png"
        save(fig, out, 7.2, 6.3)
        outputs.append({"figure": out.name, "role": "gene_effect_scatter", "lineage": lineage})

    classes = [
        "CONSERVED_UP",
        "CONSERVED_DOWN",
        "OPPOSITE_MEANINGFUL",
        "DOG_ONLY_EFFECT",
        "HUMAN_ONLY_EFFECT",
        "WEAK_BOTH",
    ]
    counts = []
    for lineage, frame in effects.items():
        for label in classes:
            counts.append(
                {
                    "lineage": lineage,
                    "effect_class": label,
                    "genes": int(frame["effect_class"].eq(label).sum()),
                }
            )
    count_frame = pd.DataFrame(counts)
    pivot = count_frame.pivot(index="lineage", columns="effect_class", values="genes").loc[
        list(LINEAGES), classes
    ]
    fig, ax = plt.subplots()
    bottom = np.zeros(len(pivot))
    for label in classes:
        values = pivot[label].to_numpy()
        ax.bar(pivot.index, values, bottom=bottom, label=label, color=COLORS[label])
        bottom += values
    ax.set(
        title="Meaningful-effect quadrants and species-skewed effects",
        xlabel=None,
        ylabel="Shared one-to-one genes",
    )
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)
    out = FIG / "MEANINGFUL_EFFECT_QUADRANT_COUNTS.png"
    save(fig, out, 9.0, 6.2)
    outputs.append({"figure": out.name, "role": "effect_class_counts", "lineage": "ALL"})

    pathway = pd.read_csv(PATHWAY / "PATHWAY_CROSS_SPECIES_BASELINE.csv.gz")
    for collection in COLLECTIONS:
        sub = pathway.loc[pathway["collection"].eq(collection)]
        fig, axes = plt.subplots(1, 3, sharex=False, sharey=False)
        for ax, lineage in zip(axes, LINEAGES, strict=True):
            frame = sub.loc[sub["lineage"].eq(lineage)].dropna(
                subset=["dog_NES", "human_NES"]
            )
            both = (frame["dog_FDR"] < 0.10) & (frame["human_FDR"] < 0.10)
            ax.scatter(
                frame.loc[~both, "dog_NES"],
                frame.loc[~both, "human_NES"],
                s=11,
                alpha=0.35,
                color="#7B8DA3",
                edgecolors="none",
                label="Other tested",
            )
            ax.scatter(
                frame.loc[both, "dog_NES"],
                frame.loc[both, "human_NES"],
                s=24,
                alpha=0.85,
                color="#C9473A",
                edgecolors="white",
                linewidths=0.25,
                label="Both FDR < 0.10",
            )
            limit = max(
                1.5,
                np.nanmax(np.abs(frame[["dog_NES", "human_NES"]].to_numpy())) * 1.05,
            )
            ax.plot([-limit, limit], [-limit, limit], color="#444444", linewidth=0.7)
            ax.axhline(0, color="#999999", linewidth=0.5)
            ax.axvline(0, color="#999999", linewidth=0.5)
            rho = frame[["dog_NES", "human_NES"]].corr(method="spearman").iloc[0, 1]
            ax.set(
                title=f"{lineage}\nρ={rho:.3f}; n={len(frame)}",
                xlabel="Dog NES",
                ylabel="Human NES",
                xlim=(-limit, limit),
                ylim=(-limit, limit),
            )
        axes[-1].legend(loc="lower right", fontsize=7, frameon=True)
        fig.suptitle(f"{collection}: paired common-one2one pathway NES", y=1.02)
        out = FIG / f"{collection}_PATHWAY_NES_SCATTER.png"
        save(fig, out, 13.5, 4.6)
        outputs.append({"figure": out.name, "role": "pathway_nes_scatter", "lineage": "ALL"})

    matrix = pd.read_csv(SUMMARY / "PHASE10_PROGRAM_BASELINE_MATRIX.csv")
    matrix["row_label"] = matrix["lineage"] + " · " + matrix["program"]
    numeric = matrix.set_index("row_label")[
        ["dog_common1to1_NES", "human_common1to1_NES"]
    ]
    fig, (ax_heat, ax_text) = plt.subplots(
        1,
        2,
        gridspec_kw={"width_ratios": [1.0, 3.0]},
    )
    sns.heatmap(
        numeric,
        cmap="RdBu_r",
        center=0,
        vmin=-3,
        vmax=3,
        annot=True,
        fmt=".2f",
        annot_kws={"fontsize": 6},
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "NES", "shrink": 0.55},
        ax=ax_heat,
    )
    ax_heat.set(xlabel=None, ylabel=None, title="Common 1:1 NES")
    ax_heat.set_xticklabels(["Dog", "Human"], rotation=0)
    ax_heat.set_yticks(np.arange(len(numeric)) + 0.5)
    ax_heat.set_yticklabels(numeric.index, rotation=0, fontsize=6.2)
    ax_text.set_xlim(0, 5)
    ax_text.set_ylim(len(matrix), 0)
    ax_text.axis("off")
    headers = ["FRE", "Region", "Dog LODO", "Treatment", "Baseline"]
    for j, header in enumerate(headers):
        ax_text.text(j + 0.45, -0.5, header, ha="center", va="bottom", weight="bold", fontsize=8)
    for i, row in matrix.iterrows():
        values = [
            row["fre_summary_status"],
            row["region_caution"],
            row["dog_lodo_stability"],
            row["treatment_caution"],
            row["baseline_status"],
        ]
        for j, value in enumerate(values):
            ax_text.text(
                j + 0.45,
                i + 0.5,
                str(value).replace("CONCORDANT_", "C_").replace("_DIRECTIONALLY_STABLE", ""),
                ha="center",
                va="center",
                fontsize=5.5,
                color="#222222",
            )
            ax_text.add_patch(
                plt.Rectangle(
                    (j, i),
                    0.9,
                    1,
                    facecolor="#F4F4F4" if i % 2 == 0 else "#EAEAEA",
                    edgecolor="white",
                    linewidth=0.3,
                    zorder=-1,
                )
            )
    fig.suptitle("Phase10 25-program baseline and frozen provenance", y=1.01)
    out = FIG / "PHASE10_25_PROGRAM_HEATMAP.png"
    save(fig, out, 17.5, 12.0)
    outputs.append({"figure": out.name, "role": "program_heatmap", "lineage": "ALL"})

    gene_summary = pd.read_csv(SUMMARY / "GENE_EFFECT_CONCORDANCE_SUMMARY.csv")
    lodo = pd.read_csv(SENS / "DOG_LODO_CROSS_SPECIES_SENSITIVITY.csv")
    treatment = pd.read_csv(SENS / "HUMAN_TREATMENT_CROSS_SPECIES_SENSITIVITY.csv")
    fig, axes = plt.subplots(1, 3, sharex=True)
    for ax, lineage in zip(axes, LINEAGES, strict=True):
        primary = gene_summary.loc[gene_summary["lineage"].eq(lineage)].iloc[0]
        rows = [
            {
                "label": "Primary Full7",
                "rho": primary["spearman_rho"],
                "axis": "Primary",
            }
        ]
        for row in lodo.loc[lodo["lineage"].eq(lineage)].itertuples(index=False):
            rows.append({"label": f"LODO {row.short_label}", "rho": row.spearman_rho, "axis": "Dog donor"})
        for row in treatment.loc[treatment["lineage"].eq(lineage)].itertuples(index=False):
            label = str(row.human_treatment_model).replace("_active_Crohn_TI_vs_Healthy_TI", "")
            rows.append({"label": label, "rho": row.spearman_rho, "axis": "Human treatment"})
        plot = pd.DataFrame(rows)
        colors = plot["axis"].map(
            {"Primary": "#1B1B1B", "Dog donor": "#2E63A3", "Human treatment": "#B35A3A"}
        )
        y = np.arange(len(plot))
        ax.scatter(plot["rho"], y, c=colors, s=35, zorder=3)
        ax.axvline(0, color="#777777", linewidth=0.7)
        ax.axvline(float(primary["spearman_rho"]), color="#1B1B1B", linestyle="--", linewidth=0.8)
        ax.set_yticks(y, plot["label"], fontsize=6.5)
        ax.set_title(lineage)
        ax.set_xlabel("Spearman ρ")
        ax.invert_yaxis()
    fig.suptitle("Cross-species gene-effect sensitivity (separate axes; no pooled P value)", y=1.02)
    out = FIG / "CROSS_SPECIES_SENSITIVITY_FOREST.png"
    save(fig, out, 15.5, 7.4)
    outputs.append({"figure": out.name, "role": "sensitivity_forest", "lineage": "ALL"})

    pd.DataFrame(outputs).to_csv(
        BASE / "audit" / "PHASE10_FIGURE_AUDIT.csv", index=False, quoting=csv.QUOTE_MINIMAL
    )
    with LOG.open("w", encoding="utf-8") as handle:
        handle.write(f"Generated {len(outputs)} Phase10 figures\n")
        for item in outputs:
            handle.write(f"{item['figure']}\t{item['role']}\t{item['lineage']}\n")


if __name__ == "__main__":
    main()
