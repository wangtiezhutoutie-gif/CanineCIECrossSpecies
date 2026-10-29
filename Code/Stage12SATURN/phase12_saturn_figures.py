#!/usr/bin/env python
"""Create the preregistered Phase12 SATURN figures from frozen outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA


LINEAGES = ("Epithelial", "T", "Myeloid")
COLORS = {
    "Epithelial": "#D55E00",
    "T": "#0072B2",
    "Myeloid": "#009E73",
    "dog": "#E69F00",
    "human": "#56B4E9",
    "dog_Healthy": "#F0C36E",
    "dog_CIE": "#C44E52",
    "human_Healthy": "#8EC5E8",
    "human_Crohn": "#4C72B0",
    "INFLAMMATORY_LIKE": "#C44E52",
    "METABOLIC_LIKE": "#4C72B0",
    "MIXED_INFLAMMATORY_METABOLIC": "#8172B2",
    "OTHER_ANNOTATION": "#999999",
}


def discover_final(root: Path) -> Path:
    directory = root / "training" / "seed20260827" / "run" / "saturn_results"
    paths = [
        path for path in directory.glob("*_seed_20260827.h5ad")
        if "_pretrain" not in path.name and "_ep_" not in path.name
    ]
    if len(paths) != 1:
        raise RuntimeError(f"Expected one Primary H5AD, found {paths}")
    return paths[0]


def save(fig, path: Path) -> None:
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def categorical_scatter(coords, labels, title, path, order=None):
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    order = order or sorted(pd.unique(labels))
    for label in order:
        mask = np.asarray(labels == label)
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=2,
            alpha=0.45,
            linewidths=0,
            rasterized=True,
            color=COLORS.get(label, None),
            label=str(label),
        )
    ax.set(title=title, xlabel="UMAP1", ylabel="UMAP2")
    ax.legend(markerscale=4, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.text(
        0.01, 0.01, "Visualization only; no quantitative metric uses UMAP",
        transform=ax.transAxes, fontsize=8, color="#555555"
    )
    save(fig, path)


def centered_donor_pca(centroids):
    dims = [column for column in centroids if column.startswith("latent_")]
    records = []
    pca_by_lineage = {}
    for lineage in LINEAGES:
        block = centroids.loc[centroids["lineage"].eq(lineage)].copy()
        matrix = block[dims].to_numpy(float)
        centered = matrix.copy()
        for species in ("dog", "human"):
            species_mask = block["species"].eq(species).to_numpy()
            healthy_mask = species_mask & block["condition"].eq("Healthy").to_numpy()
            centered[species_mask] -= matrix[healthy_mask].mean(axis=0)
        pca = PCA(n_components=2, random_state=20260827)
        coordinates = pca.fit_transform(centered)
        pca_by_lineage[lineage] = (pca, block, centered)
        for row, xy in zip(block.itertuples(), coordinates):
            records.append(
                {
                    "lineage": lineage,
                    "species": row.species,
                    "donor_id": row.donor_id,
                    "condition": row.condition,
                    "PC1": xy[0],
                    "PC2": xy[1],
                }
            )
    return pd.DataFrame(records), pca_by_lineage


def theme_from_text(text):
    upper = str(text).upper()
    inflammatory = any(
        token in upper for token in (
            "TNF", "NFK", "NF-K", "INTERFERON", "CYTOKINE", "INTERLEUKIN",
            "JAK", "TOLL", "INFLAM", "IMMUNE", "CHEMOKINE"
        )
    )
    metabolic = any(
        token in upper for token in (
            "OXIDATIVE", "RESPIR", "MITOCHON", "METABOL", "GLYCOL",
            "TRICARBOXY", "TCA", "FATTY ACID", "ELECTRON TRANSPORT"
        )
    )
    if inflammatory and metabolic:
        return "MIXED_INFLAMMATORY_METABOLIC"
    if inflammatory:
        return "INFLAMMATORY_LIKE"
    if metabolic:
        return "METABOLIC_LIKE"
    return "OTHER_ANNOTATION"


def thematic_table(root, lineage):
    priority = pd.read_csv(root / "macrogenes" / "PRIORITY_MACROGENES.csv")
    annotation = pd.read_csv(
        root / "macrogenes" / "MACROGENE_FUNCTIONAL_ANNOTATION.csv"
    )
    selected = priority.loc[priority["lineage"].eq(lineage)].copy()
    text = (
        annotation.loc[annotation["lineage"].eq(lineage)]
        .groupby("macrogene")["pathway_name"]
        .agg(lambda values: " | ".join(map(str, values)))
    )
    selected["annotation_text"] = selected["macrogene"].map(text).fillna("")
    selected["functional_theme"] = selected["annotation_text"].map(theme_from_text)
    selected["interpretation_scope"] = "post-hoc annotation; no SATURN retraining"
    return selected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root", type=Path, default=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
    )
    args = parser.parse_args()
    project = args.project_root.resolve()
    root = project / "12_saturn"
    figures = root / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    final = discover_final(root)
    adata = sc.read_h5ad(final)
    sidecar = pd.read_csv(root / "input" / "CELL_CONDITION_SIDECAR.csv")
    obs = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "species": adata.obs["species"].astype(str).to_numpy(),
            "lineage": adata.obs["labels2"].astype(str).to_numpy(),
        }
    ).merge(sidecar, on=["cell_id", "species", "lineage"], validate="one_to_one")
    if len(obs) != adata.n_obs:
        raise RuntimeError("UMAP sidecar join lost cells")
    obs["species_condition"] = obs["species"] + "_" + obs["condition"]

    umap_path = root / "primary" / "SATURN_UMAP_COORDINATES.csv.gz"
    if umap_path.is_file():
        umap = pd.read_csv(umap_path)
        if not umap["cell_id"].astype(str).eq(obs["cell_id"].astype(str)).all():
            raise RuntimeError("Cached UMAP cell order mismatch")
        coords = umap[["UMAP1", "UMAP2"]].to_numpy(float)
    else:
        viz = sc.AnnData(X=np.asarray(adata.X, dtype=np.float32), obs=obs.copy())
        sc.pp.neighbors(
            viz, n_neighbors=30, metric="cosine", random_state=20260827
        )
        sc.tl.umap(viz, min_dist=0.3, random_state=20260827)
        coords = np.asarray(viz.obsm["X_umap"])
        pd.DataFrame(
            {"cell_id": obs["cell_id"], "UMAP1": coords[:, 0], "UMAP2": coords[:, 1]}
        ).to_csv(umap_path, index=False, compression="gzip")

    categorical_scatter(
        coords, obs["lineage"].to_numpy(),
        "SATURN joint embedding by Primary Level1 lineage",
        figures / "FIG01_SATURN_JOINT_BY_CELLTYPE.png",
        list(LINEAGES),
    )
    categorical_scatter(
        coords, obs["species"].to_numpy(),
        "SATURN joint embedding by species",
        figures / "FIG02_SATURN_JOINT_BY_SPECIES.png",
        ["dog", "human"],
    )
    categorical_scatter(
        coords, obs["species_condition"].to_numpy(),
        "SATURN joint embedding by post-hoc condition",
        figures / "FIG03_SATURN_JOINT_BY_CONDITION_POSTHOC.png",
        ["dog_Healthy", "dog_CIE", "human_Healthy", "human_Crohn"],
    )

    centroids = pd.read_csv(root / "disease_state" / "DONOR_CENTROIDS_SATURN.csv.gz")
    donor_pca, pca_by_lineage = centered_donor_pca(centroids)
    donor_pca.to_csv(
        root / "primary" / "SATURN_DONOR_CENTROID_PCA_VISUALIZATION.csv", index=False
    )
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, lineage in zip(axes, LINEAGES):
        block = donor_pca.loc[donor_pca["lineage"].eq(lineage)]
        for species, marker in (("dog", "o"), ("human", "^")):
            for condition in sorted(block.loc[block["species"].eq(species), "condition"].unique()):
                selected = block["species"].eq(species) & block["condition"].eq(condition)
                label = f"{species}_{condition}"
                ax.scatter(
                    block.loc[selected, "PC1"], block.loc[selected, "PC2"],
                    marker=marker, s=38, color=COLORS.get(label, COLORS[species]),
                    label=label, alpha=0.85
                )
        ax.set(title=lineage, xlabel="PC1", ylabel="PC2")
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="upper center", ncol=4)
    fig.suptitle("Species-centered SATURN donor centroids")
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    save(fig, figures / "FIG04_SATURN_DONOR_CENTROIDS.png")

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    for ax, lineage in zip(axes, LINEAGES):
        pca, block, centered = pca_by_lineage[lineage]
        for species in ("dog", "human"):
            disease = "CIE" if species == "dog" else "Crohn"
            species_mask = block["species"].eq(species).to_numpy()
            disease_mask = species_mask & block["condition"].eq(disease).to_numpy()
            vector = centered[disease_mask].mean(axis=0)
            vector2 = vector @ pca.components_.T
            ax.arrow(
                0, 0, vector2[0], vector2[1], width=0.01,
                length_includes_head=True, color=COLORS[species], label=species
            )
            ax.text(vector2[0], vector2[1], species, color=COLORS[species], fontsize=9)
        ax.axhline(0, color="#CCCCCC", lw=0.6)
        ax.axvline(0, color="#CCCCCC", lw=0.6)
        ax.set(title=lineage, xlabel="PC1", ylabel="PC2")
    fig.suptitle("Dog CIE and Human Crohn disease-vector arrows")
    fig.tight_layout()
    save(fig, figures / "FIG05_SATURN_DISEASE_VECTOR_ARROWS.png")

    disease = pd.read_csv(root / "primary" / "SATURN_DISEASE_VECTOR_METRICS.csv")
    permutation = pd.read_csv(
        root / "primary" / "SATURN_DISEASE_VECTOR_PERMUTATION.csv"
    )
    joint = permutation.loc[permutation["test"].eq("joint_disease_vector_cosine")]
    plot = disease.merge(
        joint[["lineage", "null_q025", "null_q975", "empirical_percentile"]],
        on="lineage",
    ).set_index("lineage").loc[list(LINEAGES)].reset_index()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    x = np.arange(len(plot))
    y = plot["cosine_similarity"].to_numpy()
    null_low = plot["null_q025"].to_numpy()
    null_high = plot["null_q975"].to_numpy()
    ax.bar(x, y, color=[COLORS[value] for value in plot["lineage"]])
    ax.vlines(x, null_low, null_high, color="black", linewidth=2, label="95% null interval")
    ax.scatter(x, null_low, marker="_", color="black", s=70)
    ax.scatter(x, null_high, marker="_", color="black", s=70)
    ax.axhline(0, color="black", lw=0.8)
    ax.set(xticks=x, xticklabels=plot["lineage"], ylabel="Disease-vector cosine")
    ax.set_title("Observed cosine with donor-permutation 95% null interval")
    save(fig, figures / "FIG06_SATURN_DISEASE_COSINE_PERMUTATION.png")

    effects = pd.read_csv(root / "macrogenes" / "MACROGENE_DISEASE_EFFECTS.csv.gz")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, lineage in zip(axes, LINEAGES):
        block = effects.loc[effects["lineage"].eq(lineage)]
        colors = np.where(block["same_direction"], "#009E73", "#C44E52")
        ax.scatter(block["dog_effect"], block["human_effect"], s=7, c=colors, alpha=0.5)
        ax.axhline(0, color="#999999", lw=0.6)
        ax.axvline(0, color="#999999", lw=0.6)
        ax.set(title=lineage, xlabel="Dog macrogene effect", ylabel="Human effect")
    fig.tight_layout()
    save(fig, figures / "FIG07_MACROGENE_EFFECT_CONCORDANCE.png")

    priority = pd.read_csv(root / "macrogenes" / "PRIORITY_MACROGENES.csv")
    fig, axes = plt.subplots(1, 3, figsize=(15, 7))
    for ax, lineage in zip(axes, LINEAGES):
        block = priority.loc[priority["lineage"].eq(lineage)].copy()
        block["label"] = block["macrogene"] + " " + block["priority_class"].str.replace("TOP_", "")
        matrix = block[["dog_effect", "human_effect"]].to_numpy()
        bound = np.nanmax(np.abs(matrix))
        image = ax.imshow(matrix, aspect="auto", cmap="coolwarm", vmin=-bound, vmax=bound)
        ax.set(
            title=lineage, xticks=[0, 1], xticklabels=["Dog", "Human"],
            yticks=np.arange(len(block)), yticklabels=block["label"], xlabel="Effect"
        )
    fig.colorbar(image, ax=axes, shrink=0.6, label="Donor-centroid disease effect")
    save(fig, figures / "FIG08_PRIORITY_MACROGENES.png")

    stability = pd.read_csv(root / "sensitivity" / "SATURN_SEED_STABILITY.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for ax, metric, title in (
        (axes[0], "disease_vector_cosine", "Disease-vector cosine"),
        (axes[1], "macrogene_spearman_rho", "Macrogene effect Spearman"),
    ):
        block = stability.loc[stability["metric"].eq(metric)]
        for lineage in LINEAGES:
            selected = block.loc[block["lineage"].eq(lineage)].sort_values("seed")
            ax.plot(selected["seed"].astype(str), selected["value"], marker="o", label=lineage)
        ax.axhline(0, color="#999999", lw=0.6)
        ax.set(title=title, xlabel="Seed", ylabel="Value")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    save(fig, figures / "FIG09_SATURN_SEED_STABILITY.png")

    lodo = pd.read_csv(root / "sensitivity" / "SATURN_DOG_LODO_DISEASE_VECTOR.csv")
    fig, ax = plt.subplots(figsize=(10, 5))
    for lineage in LINEAGES:
        block = lodo.loc[lodo["lineage"].eq(lineage)]
        ax.plot(
            block["dropped_donor"].str.replace("GSM803", "", regex=False),
            block["cosine_similarity"], marker="o", label=lineage
        )
    ax.axhline(0, color="#999999", lw=0.6)
    ax.set(xlabel="Dropped Dog donor", ylabel="Cosine", title="Post-hoc donor-centroid LODO")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, figures / "FIG10_SATURN_DOG_LODO.png")

    t_table = thematic_table(root, "T")
    t_table.to_csv(root / "macrogenes" / "T_FUNCTIONAL_ASYMMETRY.csv", index=False)
    myeloid_table = thematic_table(root, "Myeloid")
    myeloid_table.to_csv(
        root / "macrogenes" / "MYELOID_FUNCTIONAL_STRUCTURE.csv", index=False
    )
    fig, ax = plt.subplots(figsize=(7, 6))
    for theme, block in t_table.groupby("functional_theme"):
        ax.scatter(
            block["dog_effect"], block["human_effect"], s=45,
            color=COLORS.get(theme, "#999999"), label=theme, alpha=0.85
        )
        for row in block.itertuples():
            ax.text(row.dog_effect, row.human_effect, row.macrogene, fontsize=6)
    ax.axhline(0, color="#999999", lw=0.6)
    ax.axvline(0, color="#999999", lw=0.6)
    ax.set(
        xlabel="Dog T macrogene effect", ylabel="Human T macrogene effect",
        title="T functional asymmetry among priority macrogenes"
    )
    ax.legend(frameon=False, fontsize=8)
    save(fig, figures / "FIG11_T_FUNCTIONAL_ASYMMETRY.png")
    print("PASS: 11 Phase12 figures and thematic audit tables")


if __name__ == "__main__":
    main()

