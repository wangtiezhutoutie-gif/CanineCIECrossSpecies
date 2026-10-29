#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P13 = ROOT / "13_concord"
sys.path.insert(0, str(ROOT / "scripts"))
import phase13_integration_qc as iqc
import phase13_disease_vectors as dv

SEED = 20260827
MODES = {
    "PRIMARY_HCL_DOMAIN_SPECIES": P13 / "training" / f"seed{SEED}" / "concord_embedding.h5ad",
    "HCL_DOMAIN_SPECIES_DONOR": P13 / "sensitivity" / "domain_species_donor" / "concord_embedding.h5ad",
    "KNN_DOMAIN_SPECIES": P13 / "sensitivity" / "knn_mode" / "concord_embedding.h5ad",
}


def load(path: Path):
    a = ad.read_h5ad(path)
    z = np.asarray(a.obsm["X_concord"], dtype=np.float32)
    if z.shape != (28046, 300) or not np.isfinite(z).all():
        raise RuntimeError(f"Bad embedding: {path}")
    return a, z, a.obs.copy()


def transfer_rows(mode: str, z: np.ndarray, obs: pd.DataFrame):
    rows = []
    species = obs["species"].astype(str).to_numpy()
    labels = obs["harmonized_level1"].astype(str).to_numpy()
    for source, target in [("dog", "human"), ("human", "dog")]:
        src, tgt = species == source, species == target
        xs, xt, ys, yt = z[src], z[tgt], labels[src], labels[tgt]
        nn = NearestNeighbors(n_neighbors=iqc.K, metric="cosine", algorithm="brute", n_jobs=-1)
        nn.fit(xs)
        pred_knn = iqc.vote(ys[nn.kneighbors(xt, return_distance=False)])
        pred_linear = iqc.ridge_probe_predictions(xs, ys, xt, [])[0]
        for probe, pred in [("knn", pred_knn), ("linear", pred_linear)]:
            for metric, value in iqc.metric_dict(yt, pred).items():
                rows.append({
                    "mode": mode, "direction": f"{source}_to_{target}",
                    "probe": probe, "metric": metric, "value": value,
                })
    return rows


def disease_rows(mode: str, a: ad.AnnData, z: np.ndarray):
    training = pd.DataFrame({
        "cell_id": a.obs_names.astype(str),
        "train_species": a.obs["species"].astype(str).to_numpy(),
        "train_lineage": a.obs["harmonized_level1"].astype(str).to_numpy(),
    })
    sidecar = pd.read_csv(P13 / "input" / "CONCORD_CELL_CONDITION_SIDECAR.csv", low_memory=False)
    meta = training.merge(sidecar, on="cell_id", how="left", validate="one_to_one")
    if meta["species"].isna().any() or not meta["train_species"].eq(meta["species"]).all():
        raise RuntimeError("Posthoc sidecar join failed")
    centroids = dv.centroid_table(SEED, meta.reset_index(drop=True), z.astype(np.float64))
    rows = []
    for lineage in iqc.LINEAGES:
        frame = centroids.loc[centroids["lineage"] == lineage].reset_index(drop=True)
        dog, _, _ = dv.group_vector(frame, "dog")
        human, _, _ = dv.group_vector(frame, "human")
        metrics = dv.vector_metrics(dog, human)
        rows.append({"mode": mode, "lineage": lineage, **metrics})
    return rows


def mode_summary(mode: str, transfer: pd.DataFrame, mixing: pd.DataFrame, disease: pd.DataFrame):
    rows = []
    balance = transfer.loc[transfer["metric"] == "balanced_accuracy", "value"].median()
    for lineage in iqc.LINEAGES:
        recall = transfer.loc[transfer["metric"] == f"recall_{lineage}", "value"].median()
        mix = mixing.loc[
            (mixing["lineage"] == lineage) & (mixing["species_summary"] == "ALL")
        ].iloc[0]
        dis = disease.loc[disease["lineage"] == lineage].iloc[0]
        rows.append({
            "mode": mode, "lineage": lineage,
            "celltype_transfer_balanced_accuracy_median": float(balance),
            "celltype_transfer_lineage_recall_median": float(recall),
            "mean_cross_species_30nn_fraction": float(mix["mean_cross_species_neighbor_fraction"]),
            "species_silhouette": float(mix["species_silhouette"]),
            "species_linear_probe_accuracy": float(mix["species_linear_probe_accuracy"]),
            "species_knn_probe_accuracy": float(mix["species_knn_probe_accuracy"]),
            "disease_vector_cosine": float(dis["cosine_similarity"]),
            "disease_vector_spearman": float(dis["spearman_rho"]),
            "metric_space": "FULL_CONCORD_LATENT_300_NOT_UMAP",
        })
    return pd.DataFrame(rows)


def comparison(summary: pd.DataFrame, secondary: str):
    primary = summary.loc[summary["mode"] == "PRIMARY_HCL_DOMAIN_SPECIES"].set_index("lineage")
    other = summary.loc[summary["mode"] == secondary].set_index("lineage")
    rows = []
    for lineage in iqc.LINEAGES:
        p, s = primary.loc[lineage], other.loc[lineage]
        severe = (
            np.sign(p["disease_vector_cosine"]) != np.sign(s["disease_vector_cosine"])
            or abs(s["disease_vector_cosine"]) < 0.5 * abs(p["disease_vector_cosine"])
        )
        rows.append({
            "lineage": lineage,
            "primary_mode": "PRIMARY_HCL_DOMAIN_SPECIES",
            "sensitivity_mode": secondary,
            **{f"primary_{k}": p[k] for k in summary.columns if k not in {"mode", "lineage", "metric_space"}},
            **{f"sensitivity_{k}": s[k] for k in summary.columns if k not in {"mode", "lineage", "metric_space"}},
            "delta_cross_species_30nn_fraction": s["mean_cross_species_30nn_fraction"] - p["mean_cross_species_30nn_fraction"],
            "delta_disease_vector_cosine": s["disease_vector_cosine"] - p["disease_vector_cosine"],
            "possible_overintegration_of_donor_disease_signal": bool(severe),
            "interpretation": (
                "POSSIBLE_OVERINTEGRATION_OF_DONOR_DISEASE_SIGNAL"
                if severe and secondary == "HCL_DOMAIN_SPECIES_DONOR"
                else "SENSITIVITY_DIFFERENCE_NO_PRIMARY_REPLACEMENT"
            ),
        })
    return pd.DataFrame(rows)


def main():
    all_transfer, all_mixing, all_disease, summaries = [], [], [], []
    for mode, path in MODES.items():
        a, z, obs = load(path)
        tr = pd.DataFrame(transfer_rows(mode, z, obs))
        mix = pd.DataFrame(iqc.mixing(SEED, z, obs))
        mix.insert(0, "mode", mode)
        dis = pd.DataFrame(disease_rows(mode, a, z))
        all_transfer.append(tr)
        all_mixing.append(mix)
        all_disease.append(dis)
        summaries.append(mode_summary(mode, tr, mix, dis))
    summary = pd.concat(summaries, ignore_index=True)
    comparison(summary, "HCL_DOMAIN_SPECIES_DONOR").to_csv(
        P13 / "sensitivity" / "CONCORD_DOMAIN_KEY_SENSITIVITY.csv", index=False
    )
    comparison(summary, "KNN_DOMAIN_SPECIES").to_csv(
        P13 / "sensitivity" / "CONCORD_KNN_MODE_SENSITIVITY.csv", index=False
    )
    pd.concat(all_transfer, ignore_index=True).to_csv(
        P13 / "sensitivity" / "CONCORD_SENSITIVITY_CELLTYPE_TRANSFER.csv", index=False
    )
    pd.concat(all_mixing, ignore_index=True).to_csv(
        P13 / "sensitivity" / "CONCORD_SENSITIVITY_SPECIES_MIXING.csv", index=False
    )
    pd.concat(all_disease, ignore_index=True).to_csv(
        P13 / "sensitivity" / "CONCORD_SENSITIVITY_DISEASE_VECTOR.csv", index=False
    )
    summary.to_csv(P13 / "sensitivity" / "CONCORD_SENSITIVITY_SUMMARY.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
