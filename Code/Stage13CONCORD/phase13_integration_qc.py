#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score, recall_score, silhouette_score
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P13 = ROOT / "13_concord"
SEEDS = [20260827, 20260828, 20260829]
LINEAGES = ["Epithelial", "T", "Myeloid"]
CLASSES = np.array(LINEAGES)
N_SHUFFLE = 100
K = 30


def load_seed(seed: int) -> tuple[np.ndarray, pd.DataFrame]:
    a = ad.read_h5ad(P13 / "training" / f"seed{seed}" / "concord_embedding.h5ad")
    z = np.asarray(a.obsm["X_concord"], dtype=np.float32)
    obs = a.obs.copy()
    if z.shape != (28_046, 300) or not np.isfinite(z).all():
        raise RuntimeError(f"Invalid seed {seed} embedding")
    return z, obs


def vote(labels: np.ndarray) -> np.ndarray:
    counts = np.stack([(labels == c).sum(axis=1) for c in CLASSES], axis=1)
    return CLASSES[counts.argmax(axis=1)]


def metric_dict(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    recalls = recall_score(y_true, y_pred, labels=CLASSES, average=None, zero_division=0)
    out = {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
    }
    out.update({f"recall_{c}": v for c, v in zip(CLASSES, recalls)})
    return out


def ridge_probe_predictions(
    xs: np.ndarray, ys: np.ndarray, xt: np.ndarray, shuffled_sets: list[np.ndarray], alpha: float = 1.0
) -> list[np.ndarray]:
    """Batched ridge least-squares linear probes for observed and shuffled labels."""
    xmean = xs.mean(axis=0, keepdims=True)
    xsd = xs.std(axis=0, keepdims=True)
    xsd[xsd < 1e-6] = 1.0
    xa = np.column_stack([((xs - xmean) / xsd).astype(np.float32), np.ones(len(xs), dtype=np.float32)])
    xta = np.column_stack([((xt - xmean) / xsd).astype(np.float32), np.ones(len(xt), dtype=np.float32)])
    label_sets = [ys] + shuffled_sets
    targets = np.zeros((len(xs), len(label_sets) * len(CLASSES)), dtype=np.float32)
    for rep, labels in enumerate(label_sets):
        for ci, cls in enumerate(CLASSES):
            targets[:, rep * len(CLASSES) + ci] = labels == cls
    gram = xa.T @ xa
    penalty = np.eye(gram.shape[0], dtype=np.float32) * alpha
    penalty[-1, -1] = 0.0
    weights = np.linalg.solve(gram + penalty, xa.T @ targets)
    scores = (xta @ weights).reshape(len(xt), len(label_sets), len(CLASSES))
    return [CLASSES[scores[:, rep, :].argmax(axis=1)] for rep in range(len(label_sets))]

def transfer(seed: int, z: np.ndarray, obs: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    rows: list[dict] = []
    conf_rows: list[dict] = []
    species = obs["species"].astype(str).to_numpy()
    y = obs["harmonized_level1"].astype(str).to_numpy()
    for source, target in [("dog", "human"), ("human", "dog")]:
        src = species == source
        tgt = species == target
        xs, xt, ys, yt = z[src], z[tgt], y[src], y[tgt]
        nn = NearestNeighbors(n_neighbors=K, metric="cosine", algorithm="brute", n_jobs=-1)
        nn.fit(xs)
        indices = nn.kneighbors(xt, return_distance=False)
        observed_knn = vote(ys[indices])

        rng = np.random.default_rng(seed + (11 if source == "dog" else 29))
        shuffled_sets = [rng.permutation(ys) for _ in range(N_SHUFFLE)]
        linear_predictions = ridge_probe_predictions(xs, ys, xt, shuffled_sets)
        observed_linear = linear_predictions[0]
        null = {"knn": {m: [] for m in metric_dict(yt, observed_knn)},
                "linear": {m: [] for m in metric_dict(yt, observed_linear)}}
        for rep, shuffled in enumerate(shuffled_sets):
            pred_knn = vote(shuffled[indices])
            for m, v in metric_dict(yt, pred_knn).items():
                null["knn"][m].append(v)
            pred_linear = linear_predictions[rep + 1]
            for m, v in metric_dict(yt, pred_linear).items():
                null["linear"][m].append(v)

        for probe, pred in [("knn", observed_knn), ("linear", observed_linear)]:
            observed = metric_dict(yt, pred)
            for metric, value in observed.items():
                arr = np.asarray(null[probe][metric], dtype=float)
                rows.append({
                    "seed": seed,
                    "direction": f"{source}_to_{target}",
                    "probe": probe,
                    "metric": metric,
                    "observed": value,
                    "null_n": N_SHUFFLE,
                    "null_median": float(np.median(arr)),
                    "null_q95": float(np.quantile(arr, 0.95)),
                    "empirical_percentile": float((arr <= value).mean()),
                    "empirical_p_upper": float((1 + (arr >= value).sum()) / (N_SHUFFLE + 1)),
                    "null_definition": "source Level1 labels permuted; fixed embedding and target; no retraining of CONCORD",
                })
            for truth in CLASSES:
                for prediction in CLASSES:
                    conf_rows.append({
                        "seed": seed,
                        "direction": f"{source}_to_{target}",
                        "probe": probe,
                        "true_lineage": truth,
                        "predicted_lineage": prediction,
                        "n_cells": int(((yt == truth) & (pred == prediction)).sum()),
                    })
    return rows, conf_rows


def mixing(seed: int, z: np.ndarray, obs: pd.DataFrame) -> list[dict]:
    rows: list[dict] = []
    species_all = obs["species"].astype(str).to_numpy()
    lineage_all = obs["harmonized_level1"].astype(str).to_numpy()
    for lineage in LINEAGES:
        take = lineage_all == lineage
        zl = z[take]
        species = species_all[take]
        nbr = NearestNeighbors(n_neighbors=K + 1, metric="cosine", algorithm="brute", n_jobs=-1)
        nbr.fit(zl)
        indices = nbr.kneighbors(zl, return_distance=False)[:, 1:]
        cross = (species[indices] != species[:, None]).mean(axis=1)
        sil = silhouette_score(
            zl, species, metric="cosine", sample_size=min(4000, len(zl)), random_state=seed
        )
        train, test = train_test_split(
            np.arange(len(zl)), test_size=0.2, stratify=species, random_state=seed
        )
        lin = SGDClassifier(loss="log_loss", alpha=1e-4, max_iter=1000, tol=1e-4, random_state=seed)
        lin.fit(zl[train], species[train])
        lin_acc = accuracy_score(species[test], lin.predict(zl[test]))
        knn = KNeighborsClassifier(n_neighbors=K, metric="cosine", algorithm="brute", n_jobs=-1)
        knn.fit(zl[train], species[train])
        knn_acc = accuracy_score(species[test], knn.predict(zl[test]))
        for species_summary in ["dog", "human", "ALL"]:
            subset = np.ones(len(zl), dtype=bool) if species_summary == "ALL" else species == species_summary
            rows.append({
                "seed": seed,
                "lineage": lineage,
                "species_summary": species_summary,
                "n_cells": int(subset.sum()),
                "species_silhouette": float(sil),
                "silhouette_n_sampled": min(4000, len(zl)),
                "knn_k": K,
                "mean_same_species_neighbor_fraction": float(1 - cross[subset].mean()),
                "mean_cross_species_neighbor_fraction": float(cross[subset].mean()),
                "expected_cross_species_fraction_balanced": 0.5,
                "species_linear_probe_accuracy": float(lin_acc) if species_summary == "ALL" else np.nan,
                "species_linear_probe_error": float(1 - lin_acc) if species_summary == "ALL" else np.nan,
                "species_knn_probe_accuracy": float(knn_acc) if species_summary == "ALL" else np.nan,
                "species_knn_probe_error": float(1 - knn_acc) if species_summary == "ALL" else np.nan,
                "probe_split": "stratified_80_20_posthoc" if species_summary == "ALL" else "NA",
                "metric_space": "CONCORD_FULL_LATENT_300_COSINE_NOT_UMAP",
            })
    return rows


def main() -> None:
    transfer_rows, conf_rows, mixing_rows = [], [], []
    for seed in SEEDS:
        z, obs = load_seed(seed)
        t, c = transfer(seed, z, obs)
        transfer_rows.extend(t)
        conf_rows.extend(c)
        mixing_rows.extend(mixing(seed, z, obs))

    primary = P13 / "primary"
    method = P13 / "method_gate"
    primary.mkdir(parents=True, exist_ok=True)
    method.mkdir(parents=True, exist_ok=True)
    transfer_df = pd.DataFrame(transfer_rows)
    conf_df = pd.DataFrame(conf_rows)
    mixing_df = pd.DataFrame(mixing_rows)
    transfer_df.to_csv(primary / "CONCORD_CELLTYPE_TRANSFER.csv", index=False)
    conf_df.to_csv(primary / "CONCORD_CELLTYPE_TRANSFER_CONFUSION.csv", index=False)
    mixing_df.to_csv(primary / "CONCORD_SPECIES_MIXING.csv", index=False)

    key = transfer_df["metric"].eq("balanced_accuracy")
    transfer_pass = bool((transfer_df.loc[key, "observed"] > transfer_df.loc[key, "null_q95"]).all())
    recall_key = transfer_df["metric"].str.startswith("recall_")
    recall_pass = bool((transfer_df.loc[recall_key, "observed"] > 0.20).all())
    prediction_coverage = conf_df.groupby(["seed", "direction", "probe", "predicted_lineage"])["n_cells"].sum()
    collapse_pass = bool((prediction_coverage > 0).all() and len(prediction_coverage) == len(SEEDS) * 2 * 2 * 3)
    gate = pd.DataFrame([
        {"criterion": "celltype_transfer_balanced_accuracy_above_100shuffle_q95", "status": "PASS" if transfer_pass else "FAIL"},
        {"criterion": "all_lineage_recall_above_0.20", "status": "PASS" if recall_pass else "FAIL"},
        {"criterion": "all_three_predicted_lineages_present", "status": "PASS" if collapse_pass else "FAIL"},
        {"criterion": "embedding_finite_all_seeds", "status": "PASS"},
        {"criterion": "condition_not_joined_before_gate", "status": "PASS"},
        {"criterion": "CELLTYPE_INTEGRATION_GATE", "status": "PASS" if all([transfer_pass, recall_pass, collapse_pass]) else "FAIL"},
    ])
    gate.to_csv(method / "CELLTYPE_INTEGRATION_GATE.csv", index=False)
    if gate.iloc[-1]["status"] != "PASS":
        raise RuntimeError("Cell-type Integration Gate failed; stop before disease-state analysis")

    sat_mixing = pd.read_csv(ROOT / "12_saturn" / "primary" / "SATURN_SPECIES_MIXING.csv")
    sat_all = sat_mixing.loc[sat_mixing["species_summary"].eq("ALL"), [
        "seed", "lineage", "species_silhouette", "mean_cross_species_neighbor_fraction"
    ]].rename(columns={
        "species_silhouette": "SATURN_species_silhouette",
        "mean_cross_species_neighbor_fraction": "SATURN_cross_species30NN_fraction",
    })
    con_all = mixing_df.loc[mixing_df["species_summary"].eq("ALL")].copy()
    compare = con_all.merge(sat_all, on=["seed", "lineage"], how="left")
    compare["CONCORD_minus_SATURN_cross_species30NN"] = (
        compare["mean_cross_species_neighbor_fraction"] - compare["SATURN_cross_species30NN_fraction"]
    )
    sat_transfer = pd.read_csv(ROOT / "12_saturn" / "primary" / "SATURN_CELLTYPE_ALIGNMENT.csv")
    sat_transfer = sat_transfer.loc[
        sat_transfer["metric"].eq("balanced_accuracy"), ["seed", "direction", "observed"]
    ].rename(columns={"observed": "SATURN_celltype_transfer_balanced_accuracy"})
    con_transfer = transfer_df.loc[
        transfer_df["metric"].eq("balanced_accuracy") & transfer_df["probe"].eq("knn"),
        ["seed", "direction", "observed"],
    ].rename(columns={"observed": "CONCORD_celltype_transfer_balanced_accuracy"})
    transfer_compare = con_transfer.merge(sat_transfer, on=["seed", "direction"], how="left")
    compare.to_csv(primary / "CONCORD_VS_SATURN_INTEGRATION_QC.csv", index=False)
    transfer_compare.to_csv(primary / "CONCORD_VS_SATURN_CELLTYPE_TRANSFER.csv", index=False)
    print(gate.to_string(index=False))
    print(compare.groupby("lineage")[["mean_cross_species_neighbor_fraction", "SATURN_cross_species30NN_fraction"]].median())


if __name__ == "__main__":
    main()

