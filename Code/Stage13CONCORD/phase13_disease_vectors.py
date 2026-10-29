#!/usr/bin/env python
from __future__ import annotations

import itertools
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P13 = ROOT / "13_concord"
SEEDS = [20260827, 20260828, 20260829]
LINEAGES = ["Epithelial", "T", "Myeloid"]
DIMS = [f"latent_{i:03d}" for i in range(300)]
RNG_SEED = 20260827


def safe_cosine(a: np.ndarray, b: np.ndarray) -> float:
    den = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / den) if den > 0 else np.nan


def safe_corr(a: np.ndarray, b: np.ndarray, kind: str) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float((pearsonr if kind == "pearson" else spearmanr)(a, b).statistic)


def vector_metrics(dog: np.ndarray, human: np.ndarray) -> dict[str, float]:
    cosine = safe_cosine(dog, human)
    dog_norm = float(np.linalg.norm(dog))
    human_norm = float(np.linalg.norm(human))
    return {
        "cosine_similarity": cosine,
        "pearson_r": safe_corr(dog, human, "pearson"),
        "spearman_rho": safe_corr(dog, human, "spearman"),
        "angle_degrees": float(np.degrees(np.arccos(np.clip(cosine, -1, 1)))),
        "dog_vector_norm": dog_norm,
        "human_vector_norm": human_norm,
        "norm_ratio_dog_over_human": dog_norm / human_norm if human_norm > 0 else np.nan,
    }


def load_joined(seed: int) -> tuple[pd.DataFrame, np.ndarray]:
    a = ad.read_h5ad(P13 / "training" / f"seed{seed}" / "concord_embedding.h5ad")
    z = np.asarray(a.obsm["X_concord"], dtype=np.float64)
    if z.shape != (28_046, 300) or not np.isfinite(z).all():
        raise RuntimeError(f"Invalid formal seed {seed}")
    training = pd.DataFrame({
        "cell_id": a.obs_names.astype(str),
        "train_species": a.obs["species"].astype(str).to_numpy(),
        "train_lineage": a.obs["harmonized_level1"].astype(str).to_numpy(),
    })
    sidecar = pd.read_csv(P13 / "input" / "CONCORD_CELL_CONDITION_SIDECAR.csv", low_memory=False)
    meta = training.merge(sidecar, on="cell_id", how="left", validate="one_to_one")
    if meta["species"].isna().any() or len(meta) != 28_046:
        raise RuntimeError("Post-hoc condition sidecar join failed")
    if not meta["train_species"].eq(meta["species"]).all():
        raise RuntimeError("Species mismatch after post-hoc join")
    if not meta["train_lineage"].eq(meta["lineage"]).all():
        raise RuntimeError("Lineage mismatch after post-hoc join")
    return meta.reset_index(drop=True), z


def centroid_table(seed: int, meta: pd.DataFrame, z: np.ndarray) -> pd.DataFrame:
    records = []
    group_cols = ["species", "donor_id", "lineage"]
    for keys, positions in meta.groupby(group_cols, sort=True).indices.items():
        block = meta.iloc[positions]
        record = {"seed": seed, **dict(zip(group_cols, keys))}
        for col in ["condition", "state", "treatment", "treatment_class", "region", "primary_contrast_eligible"]:
            values = block[col].dropna().astype(str).unique()
            if len(values) != 1:
                raise RuntimeError(f"Non-unique donor metadata {keys} {col}: {values}")
            record[col] = values[0]
        record["n_cells"] = len(positions)
        record.update(zip(DIMS, z[positions].mean(axis=0)))
        records.append(record)
    frame = pd.DataFrame(records)
    if len(frame) != 102:
        raise RuntimeError(f"Expected 102 centroids for seed {seed}, found {len(frame)}")
    return frame


def masks(frame: pd.DataFrame, species: str) -> tuple[np.ndarray, np.ndarray]:
    selected = frame["species"].eq(species).to_numpy()
    healthy = selected & frame["condition"].eq("Healthy").to_numpy()
    disease_label = "CIE" if species == "dog" else "Crohn"
    disease = selected & frame["condition"].eq(disease_label).to_numpy()
    return healthy, disease


def group_vector(frame: pd.DataFrame, species: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    healthy, disease = masks(frame, species)
    matrix = frame[DIMS].to_numpy(float)
    h = matrix[healthy].mean(axis=0)
    d = matrix[disease].mean(axis=0)
    return d - h, h, d


def effect_for_assignment(matrix: np.ndarray, disease_indices: tuple[int, ...]) -> np.ndarray:
    disease = np.zeros(len(matrix), dtype=bool)
    disease[list(disease_indices)] = True
    return matrix[disease].mean(axis=0) - matrix[~disease].mean(axis=0)


def unique_combinations(n: int, k: int, count: int, rng: np.random.Generator) -> list[tuple[int, ...]]:
    chosen: set[tuple[int, ...]] = set()
    while len(chosen) < count:
        chosen.add(tuple(sorted(rng.choice(n, size=k, replace=False).tolist())))
    return sorted(chosen)


def assignments_by_lineage() -> dict[str, dict[str, list[tuple[int, ...]]]]:
    rng = np.random.default_rng(RNG_SEED)
    result = {}
    for lineage in LINEAGES:
        result[lineage] = {
            "human_projection": unique_combinations(27, 11, 10_000, rng),
            "human_joint": unique_combinations(27, 11, 7_000, rng),
        }
    return result


def empirical_summary(values: np.ndarray, observed: float) -> dict[str, float]:
    return {
        "null_n": len(values),
        "null_median": float(np.median(values)),
        "null_q025": float(np.quantile(values, 0.025)),
        "null_q975": float(np.quantile(values, 0.975)),
        "empirical_percentile": float((np.sum(values < observed) + 0.5 * np.sum(values == observed)) / len(values)),
        "empirical_p_two_sided": float(min(
            1.0,
            2 * min(
                (1 + np.sum(values >= observed)) / (len(values) + 1),
                (1 + np.sum(values <= observed)) / (len(values) + 1),
            ),
        )),
    }


def analyze_seed(seed: int, centroids: pd.DataFrame, assignment_map):
    metrics_rows, projection_rows, permutation_rows = [], [], []
    joint_raw, projection_raw, effect_rows = [], [], []
    for lineage in LINEAGES:
        frame = centroids.loc[centroids["lineage"].eq(lineage)].reset_index(drop=True)
        dog, dog_h, _ = group_vector(frame, "dog")
        human, human_h, _ = group_vector(frame, "human")
        metrics_rows.append({"seed": seed, "lineage": lineage, **vector_metrics(dog, human)})
        dog_unit = dog / np.linalg.norm(dog)
        human_unit = human / np.linalg.norm(human)

        for axis_species, target_species, unit, center in [
            ("human", "dog", human_unit, dog_h),
            ("dog", "human", dog_unit, human_h),
        ]:
            target = frame.loc[frame["species"].eq(target_species)].copy()
            values = (target[DIMS].to_numpy(float) - center) @ unit
            for row, value in zip(target.itertuples(), values):
                projection_rows.append({
                    "seed": seed, "lineage": lineage, "axis_species": axis_species,
                    "projected_species": target_species, "donor_id": row.donor_id,
                    "condition": row.condition, "state": row.state,
                    "treatment_class": row.treatment_class, "projection": float(value),
                    "centering": f"{target_species}_healthy_donor_mean",
                    "statistical_unit": "donor_centroid",
                })

        dog_frame = frame.loc[frame["species"].eq("dog")].reset_index(drop=True)
        human_frame = frame.loc[frame["species"].eq("human")].reset_index(drop=True)
        dog_matrix = dog_frame[DIMS].to_numpy(float)
        human_matrix = human_frame[DIMS].to_numpy(float)
        dog_obs = tuple(np.flatnonzero(dog_frame["condition"].eq("CIE")))
        human_obs = tuple(np.flatnonzero(human_frame["condition"].eq("Crohn")))
        if len(dog_obs) != 4 or len(human_obs) != 11:
            raise RuntimeError("Disease donor counts changed")
        dog_assignments = list(itertools.combinations(range(7), 4))
        human_projection = assignment_map[lineage]["human_projection"]
        human_joint = assignment_map[lineage]["human_joint"]

        dog_proj = np.asarray([
            effect_for_assignment(dog_matrix, assignment) @ human_unit
            for assignment in dog_assignments
        ])
        human_proj = np.asarray([
            effect_for_assignment(human_matrix, assignment) @ dog_unit
            for assignment in human_projection
        ])
        for species_name, observed, values, definition in [
            ("dog", float(dog @ human_unit), dog_proj, "all C(7,4)=35 donor assignments"),
            ("human", float(human @ dog_unit), human_proj, "10000 unique 11-of-27 donor assignments"),
        ]:
            permutation_rows.append({
                "seed": seed, "lineage": lineage, "test": f"{species_name}_disease_axis_projection",
                "observed": observed, **empirical_summary(values, observed),
                "permutation_unit": "donor", "null_definition": definition,
            })
            projection_raw.extend({
                "seed": seed, "lineage": lineage, "species": species_name,
                "permutation_index": i, "projection_group_difference": float(v),
            } for i, v in enumerate(values))

        joint_values = []
        for i, dog_assignment in enumerate(dog_assignments):
            dog_effect = effect_for_assignment(dog_matrix, dog_assignment)
            for human_assignment in human_joint[i * 200:(i + 1) * 200]:
                human_effect = effect_for_assignment(human_matrix, human_assignment)
                joint_values.append(safe_cosine(dog_effect, human_effect))
        joint_values = np.asarray(joint_values)
        observed_cosine = safe_cosine(dog, human)
        permutation_rows.append({
            "seed": seed, "lineage": lineage, "test": "joint_disease_vector_cosine",
            "observed": observed_cosine, **empirical_summary(joint_values, observed_cosine),
            "permutation_unit": "donor; 35 Dog assignments x 200 Human assignments",
            "null_definition": "approximately 7000 joint donor-label assignments; exact Phase12 logic",
        })
        joint_raw.extend({
            "seed": seed, "lineage": lineage, "permutation_index": i,
            "cosine_similarity": float(v),
        } for i, v in enumerate(joint_values))

        for i, (dv, hv) in enumerate(zip(dog, human)):
            effect_rows.append({
                "seed": seed, "lineage": lineage, "latent_dimension": f"latent_{i:03d}",
                "latent_index": i, "dog_disease_loading": float(dv),
                "human_disease_loading": float(hv),
                "same_direction": bool(np.sign(dv) == np.sign(hv)),
                "min_abs_loading": float(min(abs(dv), abs(hv))),
                "abs_human_loading": float(abs(hv)),
            })
    return tuple(pd.DataFrame(x) for x in [
        metrics_rows, projection_rows, permutation_rows, joint_raw, projection_raw, effect_rows
    ])


def lodo_and_treatment(primary: pd.DataFrame):
    lodo_rows, treatment_rows = [], []
    definitions = pd.read_csv(ROOT / "12_saturn" / "sensitivity" / "DOG_LODO_DEFINITIONS.csv")
    treatment_defs = pd.read_csv(ROOT / "12_saturn" / "sensitivity" / "HUMAN_TREATMENT_SUBSET_DEFINITIONS.csv")
    for lineage in LINEAGES:
        frame = primary.loc[primary["lineage"].eq(lineage)].copy()
        full_dog, _, _ = group_vector(frame, "dog")
        human, _, _ = group_vector(frame, "human")
        human_unit = human / np.linalg.norm(human)
        full_cosine = safe_cosine(full_dog, human)
        for d in definitions.itertuples():
            reduced = frame.loc[~frame["donor_id"].eq(d.dropped_donor)].copy()
            dog, _, _ = group_vector(reduced, "dog")
            m = vector_metrics(dog, human)
            lodo_rows.append({
                "lineage": lineage, "dropped_donor": d.dropped_donor,
                "dropped_condition": d.dropped_condition,
                "analysis_name": "POSTHOC_DONOR_CENTROID_LODO",
                "embedding_retrained": "NO", "human_vector": "FIXED_FULL_PRIMARY",
                "cosine_similarity": m["cosine_similarity"], "pearson_r": m["pearson_r"],
                "projection_on_human_axis": float(dog @ human_unit),
                "full7_cosine_similarity": full_cosine,
                "absolute_cosine_change": abs(m["cosine_similarity"] - full_cosine),
            })

        dog_unit = full_dog / np.linalg.norm(full_dog)
        human_frame = frame.loc[frame["species"].eq("human")].copy()
        healthy = human_frame.loc[human_frame["condition"].eq("Healthy")]
        healthy_mean = healthy[DIMS].to_numpy(float).mean(axis=0)
        for d in treatment_defs.itertuples():
            donors = str(d.active_donors).split("|")
            active = human_frame.loc[human_frame["donor_id"].isin(donors)]
            if len(active) != int(d.n_active):
                raise RuntimeError(f"Treatment subset mismatch {d.subset} {lineage}")
            human_subset = active[DIMS].to_numpy(float).mean(axis=0) - healthy_mean
            m = vector_metrics(full_dog, human_subset)
            treatment_rows.append({
                "lineage": lineage, "subset": d.subset, "active_donors": d.active_donors,
                "n_active": d.n_active, "n_healthy": len(healthy), "rule": d.rule,
                "embedding_retrained": "NO", "dog_vector": "FIXED_FULL7_PRIMARY",
                "cosine_similarity": m["cosine_similarity"], "pearson_r": m["pearson_r"],
                "spearman_rho": m["spearman_rho"],
                "human_projection_on_dog_axis": float(human_subset @ dog_unit),
                "full_active_cosine": full_cosine,
                "absolute_cosine_change_vs_full": abs(m["cosine_similarity"] - full_cosine),
            })
    lodo = pd.DataFrame(lodo_rows)
    lodo["largest_influence_within_lineage"] = False
    for _, indices in lodo.groupby("lineage").groups.items():
        lodo.loc[lodo.loc[indices, "absolute_cosine_change"].idxmax(), "largest_influence_within_lineage"] = True
    return lodo, pd.DataFrame(treatment_rows)


def main() -> None:
    assignment_map = assignments_by_lineage()
    all_centroids, bundles = [], {}
    for seed in SEEDS:
        meta, z = load_joined(seed)
        centroids = centroid_table(seed, meta, z)
        all_centroids.append(centroids)
        bundles[seed] = analyze_seed(seed, centroids, assignment_map)

    centroids = pd.concat(all_centroids, ignore_index=True)
    metrics = pd.concat([bundles[s][0] for s in SEEDS], ignore_index=True)
    projections = pd.concat([bundles[s][1] for s in SEEDS], ignore_index=True)
    permutations = pd.concat([bundles[s][2] for s in SEEDS], ignore_index=True)
    joint_raw = pd.concat([bundles[s][3] for s in SEEDS], ignore_index=True)
    projection_raw = pd.concat([bundles[s][4] for s in SEEDS], ignore_index=True)
    effects = pd.concat([bundles[s][5] for s in SEEDS], ignore_index=True)

    primary_centroids = centroids.loc[centroids["seed"].eq(20260827)].copy()
    lodo, treatment = lodo_and_treatment(primary_centroids)
    stability_rows = []
    for lineage in LINEAGES:
        vals = metrics.loc[metrics["lineage"].eq(lineage), "cosine_similarity"].to_numpy()
        pvals = permutations.loc[
            permutations["lineage"].eq(lineage) & permutations["test"].eq("joint_disease_vector_cosine"),
            "empirical_p_two_sided",
        ].to_numpy()
        stability_rows.append({
            "metric": "disease_vector_cosine", "lineage": lineage, "direction": "dog_vs_human",
            "median_3seed": float(np.median(vals)), "min_3seed": float(vals.min()),
            "max_3seed": float(vals.max()),
            "sign_stability": "ALL_POSITIVE" if np.all(vals > 0) else "ALL_NEGATIVE" if np.all(vals < 0) else "MIXED_SIGN",
            "joint_null_significance_stability": "ALL_FDR05" if np.all(pvals < 0.05) else "NONE_FDR05" if np.all(pvals >= 0.05) else "MIXED",
            "joint_null_p_min": float(pvals.min()), "joint_null_p_max": float(pvals.max()),
        })
    stability = pd.DataFrame(stability_rows)

    primary_effects = effects.loc[effects["seed"].eq(20260827)].copy()
    primary_effects["rank_abs_human"] = primary_effects.groupby("lineage")["abs_human_loading"].rank(method="first", ascending=False)
    same_score = primary_effects["min_abs_loading"].where(primary_effects["same_direction"], -1)
    primary_effects["same_direction_score"] = same_score
    primary_effects["rank_same_direction_min_abs"] = primary_effects.groupby("lineage")["same_direction_score"].rank(method="first", ascending=False)
    primary_effects["selection_class"] = "NOT_SELECTED"
    primary_effects.loc[primary_effects["rank_abs_human"].le(10), "selection_class"] = "TOP10_ABS_HUMAN"
    both = primary_effects["rank_same_direction_min_abs"].le(10)
    primary_effects.loc[both & primary_effects["selection_class"].eq("NOT_SELECTED"), "selection_class"] = "TOP10_SAME_DIRECTION_MIN_ABS"
    primary_effects.loc[both & primary_effects["selection_class"].eq("TOP10_ABS_HUMAN"), "selection_class"] = "TOP10_BOTH_RULES"
    primary_effects["selection_rule_frozen_before_attribution"] = "YES"

    disease = P13 / "disease_state"
    primary = P13 / "primary"
    sensitivity = P13 / "sensitivity"
    attribution = P13 / "attribution"
    for p in [disease, primary, sensitivity, attribution]:
        p.mkdir(parents=True, exist_ok=True)
    centroids.to_csv(disease / "CONCORD_DONOR_CENTROIDS.csv.gz", index=False, compression="gzip")
    primary_centroids.to_csv(primary / "CONCORD_DONOR_CENTROIDS.csv.gz", index=False, compression="gzip")
    metrics.to_csv(primary / "CONCORD_DISEASE_VECTOR_METRICS.csv", index=False)
    projections.to_csv(primary / "CONCORD_DISEASE_AXIS_PROJECTION.csv", index=False)
    permutations.to_csv(primary / "CONCORD_DISEASE_VECTOR_PERMUTATION.csv", index=False)
    joint_raw.to_csv(disease / "CONCORD_JOINT_NULL_VALUES.csv.gz", index=False, compression="gzip")
    projection_raw.to_csv(disease / "CONCORD_PROJECTION_NULL_VALUES.csv.gz", index=False, compression="gzip")
    effects.to_csv(disease / "CONCORD_DISEASE_VECTOR_LOADINGS_ALL_SEEDS.csv.gz", index=False, compression="gzip")
    stability.to_csv(sensitivity / "CONCORD_SEED_STABILITY.csv", index=False)
    lodo.to_csv(sensitivity / "CONCORD_DOG_LODO_DISEASE_VECTOR.csv", index=False)
    treatment.to_csv(sensitivity / "CONCORD_HUMAN_TREATMENT_DISEASE_VECTOR.csv", index=False)
    primary_effects.to_csv(attribution / "CONCORD_DISEASE_ALIGNED_DIMENSIONS.csv", index=False)
    print(metrics.to_string(index=False))
    print(permutations.loc[permutations["test"].eq("joint_disease_vector_cosine")].to_string(index=False))
    print(lodo.loc[lodo["largest_influence_within_lineage"]].to_string(index=False))


if __name__ == "__main__":
    main()

