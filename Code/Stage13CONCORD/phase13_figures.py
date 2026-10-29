#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import umap

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P13 = ROOT / "13_concord"
FIG = P13 / "figures"
FIG.mkdir(parents=True, exist_ok=True)
SEED = 20260827
PALETTE = {"Epithelial": "#2E86AB", "T": "#E07A5F", "Myeloid": "#6A994E"}


def finish(name):
    plt.tight_layout()
    plt.savefig(FIG / name, dpi=220, bbox_inches="tight")
    plt.close()


def primary_umap():
    path = P13 / "primary" / "CONCORD_UMAP_COORDINATES.csv.gz"
    emb = ad.read_h5ad(P13 / "training" / f"seed{SEED}" / "concord_embedding.h5ad")
    if path.exists():
        coords = pd.read_csv(path)
        if len(coords) != emb.n_obs:
            raise RuntimeError("Stored UMAP coordinate count mismatch")
    else:
        z = np.asarray(emb.obsm["X_concord"], dtype=np.float32)
        xy = umap.UMAP(
            n_neighbors=30, min_dist=0.3, n_components=2, metric="cosine",
            random_state=SEED, n_jobs=1, low_memory=True,
        ).fit_transform(z)
        coords = pd.DataFrame({
            "cell_id": emb.obs_names.astype(str),
            "UMAP1": xy[:, 0], "UMAP2": xy[:, 1],
            "species": emb.obs["species"].astype(str).to_numpy(),
            "lineage": emb.obs["harmonized_level1"].astype(str).to_numpy(),
        })
        coords.to_csv(path, index=False, compression="gzip")
    side = pd.read_csv(P13 / "input" / "CONCORD_CELL_CONDITION_SIDECAR.csv", low_memory=False)
    coords = coords.merge(side[["cell_id", "condition"]], on="cell_id", how="left", validate="one_to_one")
    coords["species_condition"] = coords["species"] + " " + coords["condition"]

    for column, name, palette in [
        ("lineage", "01_CONCORD_UMAP_LEVEL1.png", PALETTE),
        ("species", "02_CONCORD_UMAP_SPECIES.png", {"dog": "#7B2CBF", "human": "#F4A261"}),
        ("species_condition", "03_CONCORD_UMAP_CONDITION_POSTHOC.png", None),
    ]:
        plt.figure(figsize=(7.2, 6))
        for value in sorted(coords[column].unique()):
            block = coords.loc[coords[column] == value]
            color = palette.get(value) if palette else None
            plt.scatter(block["UMAP1"], block["UMAP2"], s=2, alpha=0.55, label=value, c=color, rasterized=True)
        plt.legend(markerscale=4, frameon=False, fontsize=8)
        plt.xlabel("UMAP1 (visualization only)")
        plt.ylabel("UMAP2 (visualization only)")
        plt.title(f"CONCORD Primary seed {SEED}: {column}")
        finish(name)


def integration_figures():
    transfer = pd.read_csv(P13 / "primary" / "CONCORD_CELLTYPE_TRANSFER.csv")
    mixing = pd.read_csv(P13 / "primary" / "CONCORD_SPECIES_MIXING.csv")
    rows = []
    for seed in sorted(transfer["seed"].unique()):
        for lineage in PALETTE:
            recall = transfer.loc[
                (transfer["seed"] == seed) & (transfer["metric"] == f"recall_{lineage}"), "observed"
            ].median()
            cross = mixing.loc[
                (mixing["seed"] == seed) & (mixing["lineage"] == lineage)
                & (mixing["species_summary"] == "ALL"), "mean_cross_species_neighbor_fraction"
            ].iloc[0]
            rows.append({"seed": seed, "lineage": lineage, "recall": recall, "cross30": cross})
    q = pd.DataFrame(rows)
    plt.figure(figsize=(6.5, 5.2))
    sns.scatterplot(data=q, x="cross30", y="recall", hue="lineage", style="seed", palette=PALETTE, s=90)
    plt.axvline(0.5, color="0.5", linestyle="--", linewidth=1, label="balanced-mixing reference")
    plt.xlabel("Mean cross-species fraction among 30NN")
    plt.ylabel("Median cross-species lineage recall")
    plt.title("Mixing versus cell-type preservation")
    finish("04_MIXING_VS_CELLTYPE_PRESERVATION.png")

    comp = pd.read_csv(P13 / "primary" / "CONCORD_VS_SATURN_INTEGRATION_QC.csv")
    med = comp.groupby("lineage").agg(
        CONCORD=("mean_cross_species_neighbor_fraction", "median"),
        SATURN=("SATURN_cross_species30NN_fraction", "median"),
    ).reindex(PALETTE).reset_index()
    plot = med.melt("lineage", var_name="method", value_name="cross30")
    plt.figure(figsize=(7, 5))
    sns.barplot(data=plot, x="lineage", y="cross30", hue="method", palette={"CONCORD": "#3A86FF", "SATURN": "#8338EC"})
    plt.ylabel("Median cross-species 30NN fraction")
    plt.title("CONCORD improves local mixing but species remains predictable")
    finish("12_SATURN_VS_CONCORD_INTEGRATION.png")


def disease_figures():
    metrics = pd.read_csv(P13 / "primary" / "CONCORD_DISEASE_VECTOR_METRICS.csv")
    metrics = metrics.loc[metrics["seed"] == SEED].set_index("lineage").reindex(PALETTE)
    perm = pd.read_csv(P13 / "primary" / "CONCORD_DISEASE_VECTOR_PERMUTATION.csv")
    perm = perm.loc[(perm["seed"] == SEED) & (perm["test"] == "joint_disease_vector_cosine")].set_index("lineage").reindex(PALETTE)
    x = np.arange(3)
    obs = metrics["cosine_similarity"].to_numpy()
    lo, hi = perm["null_q025"].to_numpy(), perm["null_q975"].to_numpy()
    plt.figure(figsize=(7, 5))
    plt.vlines(x, lo, hi, color="#8D99AE", linewidth=4, label="joint-null 95% interval")
    plt.scatter(x, obs, color="#1D3557", s=55, zorder=3, label="observed")
    plt.legend(frameon=False)
    plt.axhline(0, color="0.5", linewidth=1)
    plt.xticks(x, list(PALETTE))
    plt.ylabel("Dog–Human disease-vector cosine")
    plt.title("Primary disease alignment with joint-null 95% interval")
    for i, lineage in enumerate(PALETTE):
        plt.text(i, obs[i] + 0.05, f'p={perm.loc[lineage, "empirical_p_two_sided"]:.3f}', ha="center", fontsize=8)
    finish("05_DISEASE_VECTOR_COSINE_NULL.png")

    proj = pd.read_csv(P13 / "primary" / "CONCORD_DISEASE_AXIS_PROJECTION.csv")
    proj = proj.loc[proj["seed"] == SEED].copy()
    plt.figure(figsize=(10, 5.5))
    sns.stripplot(data=proj, x="lineage", y="projection", hue="projected_species", dodge=True, jitter=0.25, alpha=0.7)
    plt.axhline(0, color="0.5", linewidth=1)
    plt.ylabel("Projection on other-species disease axis")
    plt.title("Donor-level disease-axis projections")
    finish("06_DONOR_DISEASE_AXIS_PROJECTION.png")

    stable = pd.read_csv(P13 / "sensitivity" / "CONCORD_SEED_STABILITY.csv").set_index("lineage").reindex(PALETTE)
    plt.figure(figsize=(7, 5))
    med = stable["median_3seed"].to_numpy()
    plt.errorbar(x, med, yerr=np.vstack([med - stable["min_3seed"], stable["max_3seed"] - med]), fmt="o", capsize=5)
    plt.axhline(0, color="0.5", linewidth=1)
    plt.xticks(x, list(PALETTE))
    plt.ylabel("Disease-vector cosine")
    plt.title("Three-seed stability (min / median / max)")
    finish("07_SEED_STABILITY.png")

    lodo = pd.read_csv(P13 / "sensitivity" / "CONCORD_DOG_LODO_DISEASE_VECTOR.csv")
    lodo["donor_short"] = lodo["dropped_donor"].str.replace("GSM803", "", regex=False)
    plt.figure(figsize=(10, 5.5))
    sns.lineplot(data=lodo, x="donor_short", y="cosine_similarity", hue="lineage", marker="o", palette=PALETTE)
    plt.ylabel("Disease-vector cosine after dropping donor")
    plt.xlabel("Dropped Dog donor")
    plt.title("Post-hoc donor-centroid LODO")
    finish("08_DOG_LODO.png")

    tx = pd.read_csv(P13 / "sensitivity" / "CONCORD_HUMAN_TREATMENT_DISEASE_VECTOR.csv")
    plt.figure(figsize=(10, 5.5))
    sns.lineplot(data=tx, x="subset", y="cosine_similarity", hue="lineage", marker="o", palette=PALETTE)
    plt.xticks(rotation=30, ha="right")
    plt.ylabel("Disease-vector cosine")
    plt.title("Human treatment sensitivity (post-hoc; no drug conclusion)")
    finish("09_HUMAN_TREATMENT_SENSITIVITY.png")


def attribution_figures():
    theme = pd.read_csv(P13 / "attribution" / "CONCORD_ATTRIBUTION_THEME_SUMMARY.csv")
    sig = theme.loc[(theme["tail_relation_to_disease"] == "DISEASE_ALIGNED") & (theme["theme"] != "OTHER")].copy()
    sig["context"] = sig["lineage"] + " " + sig["species"]
    pivot = sig.pivot_table(index="context", columns="theme", values="significant_pathway_dimension_hits", aggfunc="sum", fill_value=0)
    plt.figure(figsize=(11, 5.5))
    sns.heatmap(np.log1p(pivot), cmap="mako", annot=pivot.astype(int), fmt="d", cbar_kws={"label": "log(1 + hits)"})
    plt.xlabel("Attribution theme")
    plt.ylabel("Lineage × species")
    plt.title("Disease-aligned attribution enrichment (FDR ≤ 0.10)")
    finish("10_DISEASE_ALIGNED_ATTRIBUTION_HEATMAP.png")

    t = theme.loc[
        (theme["lineage"] == "T")
        & theme["theme"].isin(["IMMUNE_TNF_IFN_JAK", "OXPHOS_METABOLISM", "TRANSLATION"])
    ].copy()
    if t.empty:
        t = pd.DataFrame({"species": [], "tail_relation_to_disease": [], "theme": [], "significant_pathway_dimension_hits": []})
    plt.figure(figsize=(9, 5.5))
    sns.barplot(
        data=t, x="theme", y="significant_pathway_dimension_hits",
        hue="species",
        palette={"dog": "#7B2CBF", "human": "#F4A261"},
        errorbar=None,
    )
    for relation, linestyle in [("DISEASE_ALIGNED", "-"), ("DISEASE_OPPOSED", "--")]:
        subset = t.loc[t["tail_relation_to_disease"] == relation]
        if not subset.empty:
            totals = subset.groupby("theme")["significant_pathway_dimension_hits"].sum()
            for i, label in enumerate(["IMMUNE_TNF_IFN_JAK", "OXPHOS_METABOLISM", "TRANSLATION"]):
                if label in totals:
                    plt.text(i, totals[label] + 1, f"{relation.replace('DISEASE_', '')}:{int(totals[label])}", ha="center", fontsize=7)
    plt.ylabel("Significant pathway-dimension hits")
    plt.xlabel("")
    plt.title("T attribution: immune structure dominates; metabolic/translation is limited")
    finish("11_T_IMMUNE_METABOLIC_ATTRIBUTION.png")


def main():
    sns.set_theme(style="whitegrid", context="notebook")
    primary_umap()
    integration_figures()
    disease_figures()
    attribution_figures()
    expected = [FIG / f"{i:02d}_" for i in range(1, 13)]
    found = list(FIG.glob("*.png"))
    if len(found) < 12:
        raise RuntimeError(f"Expected >=12 figures, found {len(found)}")
    print(f"figures={len(found)}")


if __name__ == "__main__":
    main()
