#!/usr/bin/env python
"""Donor-aware post-hoc analysis for frozen Phase12 SATURN runs."""

from __future__ import annotations

import argparse
import itertools
import pickle
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, pearsonr, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, silhouette_score
from sklearn.model_selection import train_test_split
from sklearn.neighbors import NearestNeighbors


SEEDS = (20260827, 20260828, 20260829)
LINEAGES = ("Epithelial", "T", "Myeloid")
RNG_SEED = 20260827


def discover_final(root: Path, seed: int) -> Path:
    result_dir = root / "training" / f"seed{seed}" / "run" / "saturn_results"
    paths = [
        path for path in result_dir.glob(f"*_seed_{seed}.h5ad")
        if "_pretrain" not in path.name and "_ep_" not in path.name
    ]
    if len(paths) != 1:
        raise RuntimeError(f"Expected one formal seed {seed} H5AD, found {paths}")
    return paths[0]


def safe_cosine(a: np.ndarray, b: np.ndarray) -> float:
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / denominator) if denominator > 0 else np.nan


def safe_pearson(a: np.ndarray, b: np.ndarray) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(pearsonr(a, b).statistic)


def safe_spearman(a: np.ndarray, b: np.ndarray) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(spearmanr(a, b).statistic)


def vector_metrics(dog: np.ndarray, human: np.ndarray) -> dict[str, float]:
    cosine = safe_cosine(dog, human)
    dog_norm = float(np.linalg.norm(dog))
    human_norm = float(np.linalg.norm(human))
    return {
        "cosine_similarity": cosine,
        "pearson_r": safe_pearson(dog, human),
        "spearman_rho": safe_spearman(dog, human),
        "angle_degrees": (
            float(np.degrees(np.arccos(np.clip(cosine, -1, 1))))
            if np.isfinite(cosine) else np.nan
        ),
        "dog_vector_norm": dog_norm,
        "human_vector_norm": human_norm,
        "norm_ratio_dog_over_human": (
            dog_norm / human_norm if human_norm > 0 else np.nan
        ),
    }


def load_joined(root: Path, seed: int):
    path = discover_final(root, seed)
    adata = ad.read_h5ad(path)
    latent = np.asarray(adata.X, dtype=np.float64)
    macro = np.asarray(adata.obsm["macrogenes"], dtype=np.float64)
    if latent.shape != (28046, 256) or macro.shape != (28046, 2000):
        raise RuntimeError(f"Unexpected formal shape for {path}: {latent.shape}, {macro.shape}")
    if not np.isfinite(latent).all() or not np.isfinite(macro).all():
        raise RuntimeError(f"NaN/Inf in formal SATURN output: {path}")
    obs = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "saturn_species": adata.obs["species"].astype(str).to_numpy(),
            "saturn_lineage": adata.obs["labels2"].astype(str).to_numpy(),
        }
    )
    sidecar = pd.read_csv(root / "input" / "CELL_CONDITION_SIDECAR.csv")
    if sidecar["cell_id"].duplicated().any():
        raise RuntimeError("Condition sidecar has duplicate cell IDs")
    meta = obs.merge(sidecar, on="cell_id", how="left", validate="one_to_one")
    if meta["species"].isna().any():
        raise RuntimeError("At least one SATURN cell is absent from the condition sidecar")
    if not meta["saturn_species"].eq(meta["species"]).all():
        raise RuntimeError("Species mismatch after condition-blind sidecar join")
    if not meta["saturn_lineage"].eq(meta["lineage"]).all():
        raise RuntimeError("Lineage mismatch after condition-blind sidecar join")
    if set(meta["lineage"]) != set(LINEAGES):
        raise RuntimeError("One or more Primary lineages missing after join")
    return path, adata, meta, latent, macro


def centroid_table(meta: pd.DataFrame, matrix: np.ndarray, prefix: str) -> pd.DataFrame:
    records = []
    dimensions = [f"{prefix}_{i:04d}" for i in range(matrix.shape[1])]
    group_cols = ["species", "donor_id", "lineage"]
    for keys, indices in meta.groupby(group_cols, sort=True).indices.items():
        block = meta.iloc[indices]
        record = dict(zip(group_cols, keys))
        for column in (
            "condition", "state", "treatment", "treatment_class", "region",
            "primary_contrast_eligible"
        ):
            values = block[column].dropna().astype(str).unique()
            if len(values) != 1:
                raise RuntimeError(f"Non-unique donor metadata for {keys}, {column}: {values}")
            record[column] = values[0]
        record["n_cells"] = len(indices)
        record.update(zip(dimensions, matrix[indices].mean(axis=0)))
        records.append(record)
    frame = pd.DataFrame(records)
    expected = 34 * 3
    if len(frame) != expected:
        raise RuntimeError(f"Expected {expected} species-donor-lineage centroids, got {len(frame)}")
    return frame


def disease_masks(frame: pd.DataFrame, species: str) -> tuple[np.ndarray, np.ndarray]:
    selected = frame["species"].eq(species).to_numpy()
    if species == "dog":
        healthy = selected & frame["condition"].eq("Healthy").to_numpy()
        disease = selected & frame["condition"].eq("CIE").to_numpy()
    else:
        healthy = selected & frame["condition"].eq("Healthy").to_numpy()
        disease = selected & frame["condition"].eq("Crohn").to_numpy()
    return healthy, disease


def group_vector(
    frame: pd.DataFrame, dimensions: list[str], species: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    healthy, disease = disease_masks(frame, species)
    matrix = frame[dimensions].to_numpy(dtype=np.float64)
    healthy_mean = matrix[healthy].mean(axis=0)
    disease_mean = matrix[disease].mean(axis=0)
    return disease_mean - healthy_mean, healthy_mean, disease_mean


def official_transfer_scores(
    root: Path, final_path: Path
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    sys.path.insert(0, str(root / "software" / "SATURN"))
    import score_adata as official_score

    original_lr = official_score.LogisticRegression
    official_score.LogisticRegression = lambda **kwargs: original_lr(
        **{key: value for key, value in kwargs.items() if key != "multi_class"}
    )
    mapping = pd.read_csv(root / "reference" / "PHASE9_CELLTYPE_EVALUATION_MAP.csv")
    adata = ad.read_h5ad(final_path)
    result_rows = []
    audit_rows = []
    for source, target in (("dog", "human"), ("human", "dog")):
        result = official_score.classify_cell_types(
            adata,
            mapping,
            cell_type_column="labels2",
            species_column="species",
            species_1=source,
            species_2=target,
            C=0.005,
            verbose=False,
        )
        result_rows.append(
            {
                "direction": f"{source}_to_{target}",
                "accuracy": result["species_2_logreg_accuracy"],
                "balanced_accuracy": result["species_2_logreg_accuracy_balanced"],
            }
        )
        audit_rows.append(
            {
                "direction": f"{source}_to_{target}",
                "status": "PASS",
                "frozen_function": "score_adata.classify_cell_types",
                "compatibility_shim": (
                    "removed obsolete sklearn multi_class keyword only; "
                    "lbfgs multiclass remains multinomial"
                ),
                "training_supervision": "NONE_EVALUATION_ONLY",
                "ct_map_use": "PHASE9_EVALUATION_ONLY",
            }
        )
    official_score.LogisticRegression = original_lr
    return result_rows, audit_rows


def transfer_with_null(
    latent: np.ndarray,
    meta: pd.DataFrame,
    source: str,
    target: str,
    seed: int,
    permutation_seed: int,
    permutations: int,
) -> list[dict[str, object]]:
    source_mask = meta["species"].eq(source).to_numpy()
    target_mask = meta["species"].eq(target).to_numpy()
    source_x = latent[source_mask]
    source_y = meta.loc[source_mask, "lineage"].astype(str).to_numpy()
    target_x = latent[target_mask]
    target_y = meta.loc[target_mask, "lineage"].astype(str).to_numpy()
    train_idx, _ = train_test_split(
        np.arange(len(source_x)), test_size=0.2, random_state=0
    )
    model = LogisticRegression(C=0.005, random_state=0, n_jobs=-1)
    model.fit(source_x[train_idx], source_y[train_idx])
    predictions = model.predict(target_x)
    metrics = {
        "accuracy": accuracy_score(target_y, predictions),
        "balanced_accuracy": balanced_accuracy_score(target_y, predictions),
    }
    for lineage in LINEAGES:
        lineage_mask = target_y == lineage
        metrics[f"recall_{lineage}"] = float(
            np.mean(predictions[lineage_mask] == target_y[lineage_mask])
        )
    rng = np.random.default_rng(permutation_seed)
    null = {metric: [] for metric in metrics}
    for _ in range(permutations):
        permuted = rng.permutation(source_y)
        null_model = LogisticRegression(C=0.005, random_state=0, n_jobs=-1)
        null_model.fit(source_x[train_idx], permuted[train_idx])
        null_prediction = null_model.predict(target_x)
        null["accuracy"].append(accuracy_score(target_y, null_prediction))
        null["balanced_accuracy"].append(
            balanced_accuracy_score(target_y, null_prediction)
        )
        for lineage in LINEAGES:
            lineage_mask = target_y == lineage
            null[f"recall_{lineage}"].append(
                float(np.mean(null_prediction[lineage_mask] == target_y[lineage_mask]))
            )
    rows = []
    for metric, observed in metrics.items():
        values = np.asarray(null[metric], dtype=float)
        rows.append(
            {
                "seed": seed,
                "direction": f"{source}_to_{target}",
                "metric": metric,
                "observed": observed,
                "null_n": permutations,
                "null_median": float(np.median(values)) if len(values) else np.nan,
                "null_q95": float(np.quantile(values, 0.95)) if len(values) else np.nan,
                "empirical_percentile": (
                    float((np.sum(values < observed) + 0.5 * np.sum(values == observed)) / len(values))
                    if len(values) else np.nan
                ),
                "empirical_p_upper": (
                    float((1 + np.sum(values >= observed)) / (len(values) + 1))
                    if len(values) else np.nan
                ),
                "null_definition": "source cell-type labels permuted; fixed source split; target unchanged",
            }
        )
    return rows


def species_mixing(
    latent: np.ndarray, meta: pd.DataFrame, seed: int
) -> pd.DataFrame:
    rows = []
    for lineage_index, lineage in enumerate(LINEAGES):
        mask = meta["lineage"].eq(lineage).to_numpy()
        values = latent[mask]
        species = meta.loc[mask, "species"].astype(str).to_numpy()
        sample_n = min(4000, len(values))
        silhouette = silhouette_score(
            values,
            species,
            metric="euclidean",
            sample_size=sample_n,
            random_state=seed + lineage_index,
        )
        neighbors = min(30, len(values) - 1)
        nn = NearestNeighbors(n_neighbors=neighbors + 1, metric="euclidean", n_jobs=-1)
        indices = nn.fit(values).kneighbors(values, return_distance=False)[:, 1:]
        same_fraction = np.mean(species[indices] == species[:, None], axis=1)
        for sp in ("dog", "human", "ALL"):
            selected = np.ones(len(species), dtype=bool) if sp == "ALL" else species == sp
            rows.append(
                {
                    "seed": seed,
                    "lineage": lineage,
                    "species_summary": sp,
                    "n_cells": int(selected.sum()),
                    "species_silhouette": float(silhouette),
                    "silhouette_n_sampled": sample_n,
                    "knn_k": neighbors,
                    "mean_same_species_neighbor_fraction": float(same_fraction[selected].mean()),
                    "mean_cross_species_neighbor_fraction": float(1 - same_fraction[selected].mean()),
                    "expected_cross_species_fraction_balanced": 0.5,
                    "metric_space": "SATURN_LATENT_256_NOT_UMAP",
                }
            )
    return pd.DataFrame(rows)


def unique_combinations(
    n: int, k: int, count: int, rng: np.random.Generator
) -> list[tuple[int, ...]]:
    selected: set[tuple[int, ...]] = set()
    while len(selected) < count:
        selected.add(tuple(sorted(rng.choice(n, size=k, replace=False).tolist())))
    return sorted(selected)


def effect_for_assignment(matrix: np.ndarray, disease_indices: tuple[int, ...]) -> np.ndarray:
    disease_mask = np.zeros(len(matrix), dtype=bool)
    disease_mask[list(disease_indices)] = True
    return matrix[disease_mask].mean(axis=0) - matrix[~disease_mask].mean(axis=0)


def primary_disease_analysis(
    centroids: pd.DataFrame,
    macro_profiles: pd.DataFrame,
):
    latent_dims = [column for column in centroids if column.startswith("latent_")]
    macro_dims = [column for column in macro_profiles if column.startswith("MG_")]
    metric_rows = []
    projection_rows = []
    permutation_rows = []
    joint_null_rows = []
    projection_null_rows = []
    effect_rows = []
    macro_summary_rows = []
    macro_null_rows = []
    rng = np.random.default_rng(RNG_SEED)

    for lineage in LINEAGES:
        latent_frame = centroids.loc[
            centroids["lineage"].eq(lineage)
        ].reset_index(drop=True)
        macro_frame = macro_profiles.loc[
            macro_profiles["lineage"].eq(lineage)
        ].reset_index(drop=True)
        dog_vector, dog_healthy, _ = group_vector(latent_frame, latent_dims, "dog")
        human_vector, human_healthy, _ = group_vector(latent_frame, latent_dims, "human")
        metrics = vector_metrics(dog_vector, human_vector)
        metric_rows.append({"seed": 20260827, "lineage": lineage, **metrics})

        dog_unit = dog_vector / np.linalg.norm(dog_vector)
        human_unit = human_vector / np.linalg.norm(human_vector)
        for axis_species, target_species, unit, center in (
            ("human", "dog", human_unit, dog_healthy),
            ("dog", "human", dog_unit, human_healthy),
        ):
            target = latent_frame.loc[
                latent_frame["species"].eq(target_species)
            ].copy()
            target_matrix = target[latent_dims].to_numpy(dtype=float)
            projections = (target_matrix - center) @ unit
            for row, projection in zip(target.itertuples(), projections):
                projection_rows.append(
                    {
                        "seed": 20260827,
                        "lineage": lineage,
                        "axis_species": axis_species,
                        "projected_species": target_species,
                        "donor_id": row.donor_id,
                        "condition": row.condition,
                        "state": row.state,
                        "treatment_class": row.treatment_class,
                        "projection": float(projection),
                        "centering": f"{target_species}_healthy_donor_mean",
                    }
                )

        dog_frame = latent_frame.loc[
            latent_frame["species"].eq("dog")
        ].reset_index(drop=True)
        human_frame = latent_frame.loc[
            latent_frame["species"].eq("human")
        ].reset_index(drop=True)
        dog_matrix = dog_frame[latent_dims].to_numpy(float)
        human_matrix = human_frame[latent_dims].to_numpy(float)
        dog_observed = tuple(np.flatnonzero(dog_frame["condition"].eq("CIE")))
        human_observed = tuple(np.flatnonzero(human_frame["condition"].eq("Crohn")))
        dog_assignments = list(
            itertools.combinations(range(len(dog_frame)), len(dog_observed))
        )
        human_projection_assignments = unique_combinations(
            len(human_frame), len(human_observed), 10000, rng
        )

        dog_projection_values = np.array(
            [
                float(effect_for_assignment(dog_matrix, assignment) @ human_unit)
                for assignment in dog_assignments
            ]
        )
        human_projection_values = np.array(
            [
                float(effect_for_assignment(human_matrix, assignment) @ dog_unit)
                for assignment in human_projection_assignments
            ]
        )
        observed_dog_projection = float(dog_vector @ human_unit)
        observed_human_projection = float(human_vector @ dog_unit)
        for species, observed, values, n_perms in (
            ("dog", observed_dog_projection, dog_projection_values, 35),
            ("human", observed_human_projection, human_projection_values, 10000),
        ):
            projection_null_rows.extend(
                {
                    "lineage": lineage,
                    "species": species,
                    "permutation_index": index,
                    "projection_group_difference": float(value),
                }
                for index, value in enumerate(values)
            )
            permutation_rows.append(
                {
                    "seed": 20260827,
                    "lineage": lineage,
                    "test": f"{species}_disease_axis_projection",
                    "observed": observed,
                    "null_n": n_perms,
                    "null_median": float(np.median(values)),
                    "null_q025": float(np.quantile(values, 0.025)),
                    "null_q975": float(np.quantile(values, 0.975)),
                    "empirical_percentile": float(
                        (
                            np.sum(values < observed)
                            + 0.5 * np.sum(values == observed)
                        )
                        / len(values)
                    ),
                    "empirical_p_two_sided": float(
                        min(
                            1.0,
                            2
                            * min(
                                (1 + np.sum(values >= observed))
                                / (len(values) + 1),
                                (1 + np.sum(values <= observed))
                                / (len(values) + 1),
                            ),
                        )
                    ),
                    "permutation_unit": "donor",
                }
            )

        human_joint_assignments = unique_combinations(
            len(human_frame),
            len(human_observed),
            len(dog_assignments) * 200,
            rng,
        )
        latent_cosines = []
        for index, dog_assignment in enumerate(dog_assignments):
            dog_effect = effect_for_assignment(dog_matrix, dog_assignment)
            human_slice = human_joint_assignments[index * 200:(index + 1) * 200]
            for human_assignment in human_slice:
                human_effect = effect_for_assignment(human_matrix, human_assignment)
                latent_cosines.append(safe_cosine(dog_effect, human_effect))
        latent_cosines = np.asarray(latent_cosines)
        observed_cosine = safe_cosine(dog_vector, human_vector)
        joint_null_rows.extend(
            {
                "lineage": lineage,
                "permutation_index": index,
                "cosine_similarity": float(value),
            }
            for index, value in enumerate(latent_cosines)
        )
        permutation_rows.append(
            {
                "seed": 20260827,
                "lineage": lineage,
                "test": "joint_disease_vector_cosine",
                "observed": observed_cosine,
                "null_n": len(latent_cosines),
                "null_median": float(np.median(latent_cosines)),
                "null_q025": float(np.quantile(latent_cosines, 0.025)),
                "null_q975": float(np.quantile(latent_cosines, 0.975)),
                "empirical_percentile": float(
                    (
                        np.sum(latent_cosines < observed_cosine)
                        + 0.5 * np.sum(latent_cosines == observed_cosine)
                    )
                    / len(latent_cosines)
                ),
                "empirical_p_two_sided": float(
                    min(
                        1.0,
                        2
                        * min(
                            (1 + np.sum(latent_cosines >= observed_cosine))
                            / (len(latent_cosines) + 1),
                            (1 + np.sum(latent_cosines <= observed_cosine))
                            / (len(latent_cosines) + 1),
                        ),
                    )
                ),
                "permutation_unit": (
                    "donor; 35 Dog assignments x 200 Human assignments"
                ),
            }
        )

        dog_macro_vector, _, _ = group_vector(macro_frame, macro_dims, "dog")
        human_macro_vector, _, _ = group_vector(macro_frame, macro_dims, "human")
        for index, (dog_effect, human_effect) in enumerate(
            zip(dog_macro_vector, human_macro_vector)
        ):
            same = bool(np.sign(dog_effect) == np.sign(human_effect))
            effect_rows.append(
                {
                    "seed": 20260827,
                    "lineage": lineage,
                    "macrogene": f"MG_{index:04d}",
                    "macrogene_index": index,
                    "dog_effect": float(dog_effect),
                    "human_effect": float(human_effect),
                    "same_direction": same,
                    "min_abs_effect": float(
                        min(abs(dog_effect), abs(human_effect))
                    ),
                    "effect_class": "CONSERVED" if same else "OPPOSITE",
                    "statistical_unit": "donor_centroid",
                }
            )
        macro_metrics = vector_metrics(dog_macro_vector, human_macro_vector)
        direction = float(
            np.mean(np.sign(dog_macro_vector) == np.sign(human_macro_vector))
        )
        macro_summary_rows.append(
            {
                "seed": 20260827,
                "lineage": lineage,
                "macrogene_n": len(macro_dims),
                "pearson_r": macro_metrics["pearson_r"],
                "spearman_rho": macro_metrics["spearman_rho"],
                "cosine_similarity": macro_metrics["cosine_similarity"],
                "direction_concordance": direction,
            }
        )

        dog_macro_frame = macro_frame.loc[
            macro_frame["species"].eq("dog")
        ].reset_index(drop=True)
        human_macro_frame = macro_frame.loc[
            macro_frame["species"].eq("human")
        ].reset_index(drop=True)
        dog_macro_matrix = dog_macro_frame[macro_dims].to_numpy(float)
        human_macro_matrix = human_macro_frame[macro_dims].to_numpy(float)
        macro_cosines = []
        macro_spearmans = []
        macro_directions = []
        for index, dog_assignment in enumerate(dog_assignments):
            dog_effect = effect_for_assignment(dog_macro_matrix, dog_assignment)
            human_slice = human_joint_assignments[index * 200:(index + 1) * 200]
            for human_assignment in human_slice:
                human_effect = effect_for_assignment(
                    human_macro_matrix, human_assignment
                )
                macro_cosines.append(safe_cosine(dog_effect, human_effect))
                macro_spearmans.append(safe_spearman(dog_effect, human_effect))
                macro_directions.append(
                    float(np.mean(np.sign(dog_effect) == np.sign(human_effect)))
                )
        for index, values in enumerate(
            zip(macro_cosines, macro_spearmans, macro_directions)
        ):
            macro_null_rows.append(
                {
                    "lineage": lineage,
                    "permutation_index": index,
                    "cosine_similarity": values[0],
                    "spearman_rho": values[1],
                    "direction_concordance": values[2],
                }
            )
        observed_map = macro_summary_rows[-1]
        lineage_null = [
            row for row in macro_null_rows if row["lineage"] == lineage
        ]
        for metric in (
            "cosine_similarity",
            "spearman_rho",
            "direction_concordance",
        ):
            values = np.asarray([row[metric] for row in lineage_null])
            observed = float(observed_map[metric])
            observed_map[f"{metric}_null_median"] = float(np.median(values))
            observed_map[f"{metric}_null_q025"] = float(
                np.quantile(values, 0.025)
            )
            observed_map[f"{metric}_null_q975"] = float(
                np.quantile(values, 0.975)
            )
            observed_map[f"{metric}_empirical_percentile"] = float(
                (
                    np.sum(values < observed)
                    + 0.5 * np.sum(values == observed)
                )
                / len(values)
            )

    return (
        pd.DataFrame(metric_rows),
        pd.DataFrame(projection_rows),
        pd.DataFrame(permutation_rows),
        pd.DataFrame(joint_null_rows),
        pd.DataFrame(projection_null_rows),
        pd.DataFrame(effect_rows),
        pd.DataFrame(macro_summary_rows),
        pd.DataFrame(macro_null_rows),
    )


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def orthology_maps(project_root: Path):
    path = (
        project_root
        / "10_cross_species_baseline"
        / "mapping"
        / "DOG_HUMAN_ORTHOLOG_MASTER.csv"
    )
    frame = pd.read_csv(path, low_memory=False)

    def normalize(value) -> str:
        text = str(value).lower()
        if "one2one" in text:
            return "ONE2ONE"
        if "one2many" in text:
            return "ONE2MANY"
        if "many2many" in text:
            return "MANY2MANY"
        return "NO_ONE2ONE"

    dog_map = {}
    human_map = {}
    for row in frame.itertuples():
        cls = normalize(row.homology_type_release104)
        dog_id = str(row.dog_original_feature_id)
        human_id = str(row.human_ortholog_release104_id)
        old = dog_map.get(dog_id, "NO_ONE2ONE")
        if old == "NO_ONE2ONE" or cls == "ONE2ONE":
            dog_map[dog_id] = cls
        if human_id and human_id != "nan":
            old = human_map.get(human_id, "NO_ONE2ONE")
            if old == "NO_ONE2ONE" or cls == "ONE2ONE":
                human_map[human_id] = cls
    return dog_map, human_map


def gene_symbols(root: Path, species: str) -> dict[str, str]:
    frame = pd.read_csv(
        root
        / "protein_sequences"
        / f"{species}_stable_id_protein_map.csv.gz",
        usecols=["ensembl_gene_id", "gene_symbol"],
    ).dropna(subset=["gene_symbol"])
    frame = frame.drop_duplicates("ensembl_gene_id")
    return dict(
        zip(
            frame["ensembl_gene_id"].astype(str),
            frame["gene_symbol"].astype(str),
        )
    )


def macrogene_composition(
    root: Path,
    project_root: Path,
    effects: pd.DataFrame,
    final_path: Path,
):
    weights_path = final_path.with_name(
        final_path.stem + "_genes_to_macrogenes.pkl"
    )
    with weights_path.open("rb") as handle:
        weights = pickle.load(handle)
    by_species = {"dog": {}, "human": {}}
    for key, value in weights.items():
        species, gene_id = key.split("_", 1)
        by_species[species][gene_id] = np.asarray(value, dtype=float)
    dog_orth, human_orth = orthology_maps(project_root)
    orth_maps = {"dog": dog_orth, "human": human_orth}
    symbols = {
        species: gene_symbols(root, species)
        for species in ("dog", "human")
    }

    priority_rows = []
    for lineage in LINEAGES:
        block = effects.loc[effects["lineage"].eq(lineage)]
        conserved = block.loc[
            block["effect_class"].eq("CONSERVED")
        ].nlargest(10, "min_abs_effect")
        opposite = block.loc[
            block["effect_class"].eq("OPPOSITE")
        ].nlargest(10, "min_abs_effect")
        for priority_class, selected in (
            ("TOP_CONSERVED", conserved),
            ("TOP_OPPOSITE", opposite),
        ):
            for rank, row in enumerate(selected.itertuples(), start=1):
                priority_rows.append(
                    {
                        "lineage": lineage,
                        "macrogene": row.macrogene,
                        "macrogene_index": row.macrogene_index,
                        "priority_class": priority_class,
                        "priority_rank": rank,
                        "dog_effect": row.dog_effect,
                        "human_effect": row.human_effect,
                        "same_direction": row.same_direction,
                        "min_abs_effect": row.min_abs_effect,
                    }
                )
    priority = pd.DataFrame(priority_rows)

    weight_rows = []
    composition_rows = []
    for priority_row in priority.itertuples():
        index = int(priority_row.macrogene_index)
        for species in ("dog", "human"):
            gene_ids = np.array(list(by_species[species]))
            gene_weights = np.array(
                [by_species[species][gene][index] for gene in gene_ids]
            )
            order = np.argsort(gene_weights)[::-1]
            for rank, position in enumerate(order[:100], start=1):
                gene_id = str(gene_ids[position])
                orthology = orth_maps[species].get(
                    gene_id, "NO_ONE2ONE"
                )
                weight_rows.append(
                    {
                        "lineage": priority_row.lineage,
                        "macrogene": priority_row.macrogene,
                        "priority_class": priority_row.priority_class,
                        "priority_rank": priority_row.priority_rank,
                        "species": species,
                        "gene_rank": rank,
                        "ensembl_gene_id": gene_id,
                        "gene_symbol": symbols[species].get(
                            gene_id, "NA_NOT_AVAILABLE"
                        ),
                        "weight": float(gene_weights[position]),
                        "orthology_class": orthology,
                    }
                )
            top50 = order[:50]
            classes = [
                orth_maps[species].get(
                    str(gene_ids[position]), "NO_ONE2ONE"
                )
                for position in top50
            ]
            top_weights = gene_weights[top50]
            total_weight = float(top_weights.sum())
            for orthology_class in (
                "ONE2ONE",
                "ONE2MANY",
                "MANY2MANY",
                "NO_ONE2ONE",
            ):
                class_mask = np.array(
                    [value == orthology_class for value in classes]
                )
                composition_rows.append(
                    {
                        "lineage": priority_row.lineage,
                        "macrogene": priority_row.macrogene,
                        "priority_class": priority_row.priority_class,
                        "species": species,
                        "top_gene_scope": 50,
                        "orthology_class": orthology_class,
                        "gene_n": int(class_mask.sum()),
                        "gene_fraction": float(class_mask.mean()),
                        "weight_fraction": (
                            float(
                                top_weights[class_mask].sum()
                                / total_weight
                            )
                            if total_weight > 0
                            else np.nan
                        ),
                        "contains_non_one2one_top20": bool(
                            any(
                                value != "ONE2ONE"
                                for value in classes[:20]
                            )
                        ),
                    }
                )
    return (
        priority,
        pd.DataFrame(weight_rows),
        pd.DataFrame(composition_rows),
    )


def functional_annotation(
    project_root: Path,
    root: Path,
    priority: pd.DataFrame,
    weights: pd.DataFrame,
) -> pd.DataFrame:
    gene_set_dir = (
        project_root
        / "07_human_crohn"
        / "GSE266616"
        / "gsea"
        / "gene_sets"
    )
    sets = pd.concat(
        [
            pd.read_csv(path, sep="\t")
            for path in sorted(
                gene_set_dir.glob("HUMAN_*_GENE_SETS.tsv.gz")
            )
        ],
        ignore_index=True,
    )
    human_symbols = set(gene_symbols(root, "human").values())
    sets = sets.loc[sets["gene_symbol"].isin(human_symbols)].copy()
    grouped_sets = {
        keys: set(block["gene_symbol"].astype(str))
        for keys, block in sets.groupby(
            ["collection", "pathway_id", "pathway_name"]
        )
    }
    rows = []
    N = len(human_symbols)
    for row in priority.itertuples():
        selected = weights.loc[
            weights["lineage"].eq(row.lineage)
            & weights["macrogene"].eq(row.macrogene)
            & weights["species"].eq("human")
            & weights["gene_rank"].le(100)
        ]
        query = set(selected["gene_symbol"]) - {"NA_NOT_AVAILABLE"}
        n = len(query)
        for (
            collection,
            pathway_id,
            pathway_name,
        ), members in grouped_sets.items():
            K = len(members)
            if K < 5 or K > 500:
                continue
            overlap = query & members
            k = len(overlap)
            if k == 0:
                continue
            p_value = float(hypergeom.sf(k - 1, N, K, n))
            rows.append(
                {
                    "lineage": row.lineage,
                    "macrogene": row.macrogene,
                    "priority_class": row.priority_class,
                    "collection": collection,
                    "pathway_id": pathway_id,
                    "pathway_name": pathway_name,
                    "query_top_human_genes_n": n,
                    "universe_genes_n": N,
                    "pathway_universe_genes_n": K,
                    "overlap_n": k,
                    "overlap_genes": "|".join(sorted(overlap)),
                    "hypergeometric_p": p_value,
                    "annotation_only": "YES_POSTHOC_NO_RETRAINING",
                }
            )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame["FDR"] = np.nan
    groups = frame.groupby(
        ["lineage", "macrogene", "collection"]
    ).groups
    for _, indices in groups.items():
        frame.loc[indices, "FDR"] = bh_adjust(
            frame.loc[indices, "hypergeometric_p"].to_numpy()
        )
    frame = frame.sort_values(
        [
            "lineage",
            "macrogene",
            "collection",
            "FDR",
            "hypergeometric_p",
        ]
    )
    return frame.groupby(
        ["lineage", "macrogene", "collection"],
        group_keys=False,
    ).head(10).reset_index(drop=True)


def sensitivity(
    root: Path,
    per_seed,
    primary_centroids: pd.DataFrame,
):
    seed_rows = []
    for seed, bundle in per_seed.items():
        alignment, mixing, disease, macro_summary = bundle
        for row in alignment.itertuples():
            if row.metric in {"accuracy", "balanced_accuracy"}:
                seed_rows.append(
                    {
                        "metric": f"celltype_transfer_{row.metric}",
                        "lineage": "ALL",
                        "direction": row.direction,
                        "seed": seed,
                        "value": row.observed,
                    }
                )
        all_mixing = mixing.loc[
            mixing["species_summary"].eq("ALL")
        ]
        for row in all_mixing.itertuples():
            seed_rows.extend(
                [
                    {
                        "metric": "species_silhouette",
                        "lineage": row.lineage,
                        "direction": "NA",
                        "seed": seed,
                        "value": row.species_silhouette,
                    },
                    {
                        "metric": "cross_species_neighbor_fraction",
                        "lineage": row.lineage,
                        "direction": "NA",
                        "seed": seed,
                        "value": row.mean_cross_species_neighbor_fraction,
                    },
                ]
            )
        for row in disease.itertuples():
            seed_rows.append(
                {
                    "metric": "disease_vector_cosine",
                    "lineage": row.lineage,
                    "direction": "dog_vs_human",
                    "seed": seed,
                    "value": row.cosine_similarity,
                }
            )
        for row in macro_summary.itertuples():
            for metric in (
                "spearman_rho",
                "pearson_r",
                "cosine_similarity",
                "direction_concordance",
            ):
                seed_rows.append(
                    {
                        "metric": f"macrogene_{metric}",
                        "lineage": row.lineage,
                        "direction": "dog_vs_human",
                        "seed": seed,
                        "value": getattr(row, metric),
                    }
                )
    seed_frame = pd.DataFrame(seed_rows)
    summaries = []
    for keys, block in seed_frame.groupby(
        ["metric", "lineage", "direction"]
    ):
        values = block["value"].to_numpy(float)
        summaries.append(
            {
                "metric": keys[0],
                "lineage": keys[1],
                "direction": keys[2],
                "median_3seed": float(np.median(values)),
                "min_3seed": float(np.min(values)),
                "max_3seed": float(np.max(values)),
                "sign_stability": (
                    "ALL_POSITIVE"
                    if np.all(values > 0)
                    else "ALL_NEGATIVE"
                    if np.all(values < 0)
                    else "MIXED_SIGN"
                ),
            }
        )
    seed_frame = seed_frame.merge(
        pd.DataFrame(summaries),
        on=["metric", "lineage", "direction"],
        how="left",
    )

    dimensions = [
        column
        for column in primary_centroids
        if column.startswith("latent_")
    ]
    lodo_rows = []
    definitions = pd.read_csv(
        root / "sensitivity" / "DOG_LODO_DEFINITIONS.csv"
    )
    for lineage in LINEAGES:
        frame = primary_centroids.loc[
            primary_centroids["lineage"].eq(lineage)
        ].copy()
        full_dog, _, _ = group_vector(frame, dimensions, "dog")
        human, _, _ = group_vector(frame, dimensions, "human")
        human_unit = human / np.linalg.norm(human)
        full_metrics = vector_metrics(full_dog, human)
        for definition in definitions.itertuples():
            reduced = frame.loc[
                ~frame["donor_id"].eq(definition.dropped_donor)
            ].copy()
            dog, _, _ = group_vector(reduced, dimensions, "dog")
            metrics = vector_metrics(dog, human)
            lodo_rows.append(
                {
                    "lineage": lineage,
                    "dropped_donor": definition.dropped_donor,
                    "dropped_condition": definition.dropped_condition,
                    "analysis_name": definition.analysis_name,
                    "embedding_retrained": "NO",
                    "human_vector": "FIXED_FULL_PRIMARY",
                    "cosine_similarity": metrics["cosine_similarity"],
                    "pearson_r": metrics["pearson_r"],
                    "projection_on_human_axis": float(
                        dog @ human_unit
                    ),
                    "full7_cosine_similarity": full_metrics[
                        "cosine_similarity"
                    ],
                    "absolute_cosine_change": abs(
                        metrics["cosine_similarity"]
                        - full_metrics["cosine_similarity"]
                    ),
                }
            )
    lodo = pd.DataFrame(lodo_rows)
    lodo["largest_influence_within_lineage"] = False
    for _, indices in lodo.groupby("lineage").groups.items():
        best = lodo.loc[
            indices, "absolute_cosine_change"
        ].idxmax()
        lodo.loc[best, "largest_influence_within_lineage"] = True

    treatment_rows = []
    treatment_defs = pd.read_csv(
        root
        / "sensitivity"
        / "HUMAN_TREATMENT_SUBSET_DEFINITIONS.csv"
    )
    for lineage in LINEAGES:
        frame = primary_centroids.loc[
            primary_centroids["lineage"].eq(lineage)
        ].copy()
        dog, _, _ = group_vector(frame, dimensions, "dog")
        dog_unit = dog / np.linalg.norm(dog)
        human_frame = frame.loc[
            frame["species"].eq("human")
        ].copy()
        healthy = human_frame.loc[
            human_frame["condition"].eq("Healthy")
        ]
        healthy_mean = healthy[dimensions].to_numpy(float).mean(axis=0)
        full_human, _, _ = group_vector(frame, dimensions, "human")
        full_cosine = safe_cosine(dog, full_human)
        for definition in treatment_defs.itertuples():
            donors = str(definition.active_donors).split("|")
            active = human_frame.loc[
                human_frame["donor_id"].isin(donors)
            ]
            if len(active) != int(definition.n_active):
                raise RuntimeError(
                    f"Treatment subset donor mismatch: "
                    f"{definition.subset}, {lineage}"
                )
            human = (
                active[dimensions].to_numpy(float).mean(axis=0)
                - healthy_mean
            )
            metrics = vector_metrics(dog, human)
            treatment_rows.append(
                {
                    "lineage": lineage,
                    "subset": definition.subset,
                    "active_donors": definition.active_donors,
                    "n_active": definition.n_active,
                    "n_healthy": len(healthy),
                    "rule": definition.rule,
                    "embedding_retrained": "NO",
                    "dog_vector": "FIXED_FULL7_PRIMARY",
                    "cosine_similarity": metrics["cosine_similarity"],
                    "pearson_r": metrics["pearson_r"],
                    "spearman_rho": metrics["spearman_rho"],
                    "human_projection_on_dog_axis": float(
                        human @ dog_unit
                    ),
                    "full_active_cosine": full_cosine,
                    "absolute_cosine_change_vs_full": abs(
                        metrics["cosine_similarity"] - full_cosine
                    ),
                }
            )
    return seed_frame, lodo, pd.DataFrame(treatment_rows)


def phase10_11_comparison(
    project_root: Path,
    disease: pd.DataFrame,
    macro_summary: pd.DataFrame,
    treatment: pd.DataFrame,
):
    gene = pd.read_csv(
        project_root
        / "10_cross_species_baseline"
        / "summary"
        / "GENE_EFFECT_CONCORDANCE_SUMMARY.csv"
    )
    pathway = pd.read_csv(
        project_root
        / "10_cross_species_baseline"
        / "pathway_effects"
        / "PATHWAY_CONCORDANCE_SUMMARY.csv"
    )
    program = pd.read_csv(
        project_root
        / "10_cross_species_baseline"
        / "summary"
        / "PHASE10_PROGRAM_BASELINE_MATRIX.csv"
    )
    pathway_mean = pathway.groupby(
        "lineage", as_index=False
    ).agg(
        phase10_mean_pathway_rho=("nes_spearman_rho", "mean")
    )
    region = program.groupby("lineage", as_index=False).agg(
        region_risk_summary=(
            "region_risk",
            lambda x: "|".join(sorted(set(map(str, x)))),
        ),
        phase10_treatment_sensitivity=(
            "human_treatment_stability",
            lambda x: "|".join(sorted(set(map(str, x)))),
        ),
    )
    phase12 = disease.merge(
        macro_summary[
            ["lineage", "spearman_rho", "direction_concordance"]
        ].rename(
            columns={
                "spearman_rho": "SATURN_macrogene_rho",
                "direction_concordance": "SATURN_macrogene_direction",
            }
        ),
        on="lineage",
    )
    comparison10 = (
        gene[["lineage", "spearman_rho"]]
        .rename(columns={"spearman_rho": "phase10_gene_rho"})
        .merge(pathway_mean, on="lineage")
        .merge(phase12, on="lineage")
        .merge(region, on="lineage")
    )
    comparison10 = comparison10.rename(
        columns={
            "cosine_similarity": "SATURN_embedding_cosine",
            "spearman_rho": "SATURN_embedding_spearman",
        }
    )
    treatment_summary = treatment.groupby(
        "lineage"
    )["cosine_similarity"].agg(["min", "max"])
    comparison10["treatment_sensitivity"] = comparison10[
        "lineage"
    ].map(
        lambda lineage: (
            "SATURN treatment-subset cosine range "
            f"{treatment_summary.loc[lineage, 'min']:.6g} to "
            f"{treatment_summary.loc[lineage, 'max']:.6g}"
        )
    )
    comparison10["interpretation"] = (
        "parallel continuous metrics; no weighting or composite fidelity score"
    )

    phase11 = pd.read_csv(
        project_root
        / "11_singIST"
        / "summary"
        / "PHASE11_VS_PHASE10_BASELINE.csv"
    )
    comparison11 = phase11.merge(
        comparison10[
            [
                "lineage",
                "phase10_gene_rho",
                "phase10_mean_pathway_rho",
                "SATURN_embedding_cosine",
                "SATURN_embedding_spearman",
                "SATURN_macrogene_rho",
                "SATURN_macrogene_direction",
            ]
        ],
        on="lineage",
        how="left",
    )
    comparison11["comparison_rule"] = (
        "Phase10, singIST, and SATURN retained side-by-side; "
        "no averaging and no adjudication"
    )
    return comparison10, comparison11


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"]),
    )
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    root = project_root / "12_saturn"
    for directory in (
        "primary",
        "disease_state",
        "macrogenes",
        "sensitivity",
        "summary",
        "audit",
    ):
        (root / directory).mkdir(parents=True, exist_ok=True)

    per_seed = {}
    primary_objects = None
    all_alignment = []
    all_mixing = []
    official_rows, official_audit = official_transfer_scores(
        root, discover_final(root, 20260827)
    )
    pd.DataFrame(official_rows).to_csv(
        root / "primary" / "SATURN_OFFICIAL_SCORE_ADATA.csv",
        index=False,
    )
    pd.DataFrame(official_audit).to_csv(
        root
        / "audit"
        / "SATURN_OFFICIAL_SCORE_ADATA_COMPATIBILITY.csv",
        index=False,
    )

    for seed in SEEDS:
        path, adata, meta, latent, macro = load_joined(root, seed)
        alignment = []
        for direction_index, (source, target) in enumerate(
            (("dog", "human"), ("human", "dog"))
        ):
            alignment.extend(
                transfer_with_null(
                    latent,
                    meta,
                    source,
                    target,
                    seed,
                    seed + direction_index,
                    100 if seed == 20260827 else 0,
                )
            )
        alignment = pd.DataFrame(alignment)
        if seed == 20260827:
            for official in official_rows:
                block = alignment.loc[
                    alignment["direction"].eq(
                        official["direction"]
                    )
                    & alignment["metric"].isin(
                        ["accuracy", "balanced_accuracy"]
                    )
                ]
                for row in block.itertuples():
                    difference = abs(
                        row.observed - official[row.metric]
                    )
                    if difference > 1e-12:
                        raise RuntimeError(
                            "Official score_adata and audited "
                            "implementation differ"
                        )
        mixing = species_mixing(latent, meta, seed)
        centroids = centroid_table(meta, latent, "latent")
        macro_profiles = centroid_table(meta, macro, "MG")
        latent_dims = [
            column
            for column in centroids
            if column.startswith("latent_")
        ]
        macro_dims = [
            column
            for column in macro_profiles
            if column.startswith("MG_")
        ]
        disease_rows = []
        macro_rows = []
        for lineage in LINEAGES:
            latent_frame = centroids.loc[
                centroids["lineage"].eq(lineage)
            ]
            dog, _, _ = group_vector(
                latent_frame, latent_dims, "dog"
            )
            human, _, _ = group_vector(
                latent_frame, latent_dims, "human"
            )
            disease_rows.append(
                {
                    "seed": seed,
                    "lineage": lineage,
                    **vector_metrics(dog, human),
                }
            )
            macro_frame = macro_profiles.loc[
                macro_profiles["lineage"].eq(lineage)
            ]
            dog_macro, _, _ = group_vector(
                macro_frame, macro_dims, "dog"
            )
            human_macro, _, _ = group_vector(
                macro_frame, macro_dims, "human"
            )
            metrics = vector_metrics(dog_macro, human_macro)
            macro_rows.append(
                {
                    "seed": seed,
                    "lineage": lineage,
                    "pearson_r": metrics["pearson_r"],
                    "spearman_rho": metrics["spearman_rho"],
                    "cosine_similarity": metrics[
                        "cosine_similarity"
                    ],
                    "direction_concordance": float(
                        np.mean(
                            np.sign(dog_macro)
                            == np.sign(human_macro)
                        )
                    ),
                }
            )
        disease = pd.DataFrame(disease_rows)
        macro_summary = pd.DataFrame(macro_rows)
        per_seed[seed] = (
            alignment,
            mixing,
            disease,
            macro_summary,
        )
        all_alignment.append(alignment)
        all_mixing.append(mixing)
        if seed == 20260827:
            primary_objects = (
                path,
                meta,
                centroids,
                macro_profiles,
            )

    alignment_frame = pd.concat(
        all_alignment, ignore_index=True
    )
    mixing_frame = pd.concat(all_mixing, ignore_index=True)
    alignment_frame.to_csv(
        root / "primary" / "SATURN_CELLTYPE_ALIGNMENT.csv",
        index=False,
    )
    mixing_frame.to_csv(
        root / "primary" / "SATURN_SPECIES_MIXING.csv",
        index=False,
    )

    final_path, primary_meta, centroids, macro_profiles = (
        primary_objects
    )
    centroids.to_csv(
        root
        / "disease_state"
        / "DONOR_CENTROIDS_SATURN.csv.gz",
        index=False,
        compression="gzip",
    )
    centroids.to_csv(
        root / "primary" / "SATURN_DONOR_CENTROIDS.csv.gz",
        index=False,
        compression="gzip",
    )
    macro_profiles.to_csv(
        root
        / "disease_state"
        / "DONOR_MACROGENE_PROFILES.csv.gz",
        index=False,
        compression="gzip",
    )

    (
        disease_metrics,
        projections,
        permutations,
        joint_null,
        projection_null,
        effects,
        macro_summary,
        macro_null,
    ) = primary_disease_analysis(centroids, macro_profiles)
    disease_metrics.to_csv(
        root
        / "primary"
        / "SATURN_DISEASE_VECTOR_METRICS.csv",
        index=False,
    )
    projections.to_csv(
        root
        / "primary"
        / "SATURN_DISEASE_AXIS_PROJECTION.csv",
        index=False,
    )
    permutations.to_csv(
        root
        / "primary"
        / "SATURN_DISEASE_VECTOR_PERMUTATION.csv",
        index=False,
    )
    joint_null.to_csv(
        root
        / "disease_state"
        / "SATURN_DISEASE_VECTOR_JOINT_NULL.csv.gz",
        index=False,
        compression="gzip",
    )
    projection_null.to_csv(
        root
        / "disease_state"
        / "SATURN_DISEASE_AXIS_PROJECTION_NULL.csv.gz",
        index=False,
        compression="gzip",
    )
    effects.to_csv(
        root
        / "macrogenes"
        / "MACROGENE_DISEASE_EFFECTS.csv.gz",
        index=False,
        compression="gzip",
    )
    macro_summary.to_csv(
        root
        / "macrogenes"
        / "MACROGENE_CONCORDANCE_SUMMARY.csv",
        index=False,
    )
    macro_null.to_csv(
        root
        / "macrogenes"
        / "MACROGENE_CONCORDANCE_PERMUTATION.csv.gz",
        index=False,
        compression="gzip",
    )

    priority, weights, composition = macrogene_composition(
        root, project_root, effects, final_path
    )
    priority.to_csv(
        root / "macrogenes" / "PRIORITY_MACROGENES.csv",
        index=False,
    )
    weights.to_csv(
        root
        / "macrogenes"
        / "MACROGENE_GENE_WEIGHTS.csv.gz",
        index=False,
        compression="gzip",
    )
    composition.to_csv(
        root
        / "macrogenes"
        / "MACROGENE_ORTHOLOGY_COMPOSITION.csv",
        index=False,
    )
    annotation = functional_annotation(
        project_root, root, priority, weights
    )
    annotation.to_csv(
        root
        / "macrogenes"
        / "MACROGENE_FUNCTIONAL_ANNOTATION.csv",
        index=False,
    )

    seed_stability, lodo, treatment = sensitivity(
        root, per_seed, centroids
    )
    seed_stability.to_csv(
        root
        / "sensitivity"
        / "SATURN_SEED_STABILITY.csv",
        index=False,
    )
    lodo.to_csv(
        root
        / "sensitivity"
        / "SATURN_DOG_LODO_DISEASE_VECTOR.csv",
        index=False,
    )
    treatment.to_csv(
        root
        / "sensitivity"
        / "SATURN_HUMAN_TREATMENT_DISEASE_VECTOR.csv",
        index=False,
    )

    comparison10, comparison11 = phase10_11_comparison(
        project_root,
        disease_metrics,
        macro_summary,
        treatment,
    )
    comparison10.to_csv(
        root
        / "summary"
        / "PHASE12_VS_PHASE10_BASELINE.csv",
        index=False,
    )
    comparison11.to_csv(
        root
        / "summary"
        / "PHASE12_VS_PHASE11_SINGIST.csv",
        index=False,
    )
    provenance = comparison11[
        [
            "lineage",
            "program",
            "region_risk",
            "fre_status",
            "dog_lodo",
            "human_treatment",
            "phase10_priority",
            "phase10_baseline_status",
            "reference_model_status",
            "SATURN_embedding_cosine",
            "SATURN_macrogene_rho",
        ]
    ].copy()
    provenance[
        "saturn_training_condition_blind"
    ] = "YES"
    provenance["saturn_training_level"] = "LEVEL1_ONLY"
    provenance["phase9_ct_map_use"] = "EVALUATION_ONLY"
    provenance[
        "FRE_in_SATURN"
    ] = "NO_SUMMARY_CONTEXT_ONLY"
    provenance[
        "orthology_in_SATURN_training"
    ] = "NO_PROTEIN_EMBEDDING_MACROGENES"
    provenance[
        "region_caution_carried_forward"
    ] = "YES"
    provenance.to_csv(
        root
        / "summary"
        / "PHASE12_PROVENANCE_MATRIX.csv",
        index=False,
    )
    print(
        "PASS: Phase12 donor-aware post-hoc statistics complete"
    )


if __name__ == "__main__":
    main()

