from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
OUT = ROOT / "WORKING" / "12_SELECTION_AWARE_NULL"
SELECTED_PATH = ROOT / "10_cross_species_baseline" / "summary" / "PHASE10_PROGRAM_BASELINE_MATRIX.csv"
BASELINE_PATH = ROOT / "10_cross_species_baseline" / "pathway_effects" / "PATHWAY_CROSS_SPECIES_BASELINE.csv.gz"
LODO_PATH = ROOT / "06_pseudobulk" / "lodo" / "summary" / "CORE_PROGRAM_LODO_SUMMARY.csv"

PRIMARY_K = 20
PRIMARY_SEED = 20261014
SENSITIVITY = ((10, 20261015), (30, 20261016))
N_DRAWS = 5000
MAX_ATTEMPTS_PER_DRAW = 100


def finite(series: pd.Series) -> pd.Series:
    return np.isfinite(pd.to_numeric(series, errors="coerce"))


def robust_center_scale(values: pd.Series) -> tuple[float, float]:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(float)
    center = float(np.nanmedian(arr))
    mad = float(np.nanmedian(np.abs(arr - center))) * 1.4826
    if not np.isfinite(mad) or mad <= 0:
        mad = float(np.nanstd(arr, ddof=1))
    if not np.isfinite(mad) or mad <= 0:
        mad = 1.0
    return center, mad


def safe_spearman(x: pd.Series, y: pd.Series) -> float:
    xv = pd.to_numeric(x, errors="coerce").to_numpy(float)
    yv = pd.to_numeric(y, errors="coerce").to_numpy(float)
    keep = np.isfinite(xv) & np.isfinite(yv)
    if keep.sum() < 3:
        return math.nan
    xr = pd.Series(xv[keep]).rank(method="average").to_numpy(float)
    yr = pd.Series(yv[keep]).rank(method="average").to_numpy(float)
    return float(np.corrcoef(xr, yr)[0, 1])


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected = pd.read_csv(SELECTED_PATH)
    baseline = pd.read_csv(BASELINE_PATH)
    lodo = pd.read_csv(LODO_PATH)
    if len(selected) != 25:
        raise RuntimeError(f"Expected 25 selected units, found {len(selected)}")
    if selected[["lineage", "program"]].duplicated().any():
        raise RuntimeError("Selected lineage-program identity is not unique")
    if baseline[["lineage", "collection", "pathway_id"]].duplicated().any():
        raise RuntimeError("Baseline lineage-collection-pathway identity is not unique")
    return selected, baseline, lodo


def prepare_selected(selected: pd.DataFrame, baseline: pd.DataFrame, lodo: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    keys = selected[["lineage", "representative_collection", "representative_pathway"]].rename(
        columns={"representative_collection": "collection", "representative_pathway": "pathway_id"}
    )
    selected_keys = set(map(tuple, keys.to_numpy()))
    paired_cols = [
        "lineage", "collection", "pathway_id", "pathway_name", "dog_NES", "dog_FDR",
        "human_NES", "human_FDR", "same_direction", "dog_size", "human_size",
        "dog_pvalue", "human_pvalue", "source_version", "comparison_scope"
    ]
    paired = baseline[paired_cols].copy()
    merged = selected.merge(
        paired,
        left_on=["lineage", "representative_collection", "representative_pathway"],
        right_on=["lineage", "collection", "pathway_id"],
        how="left",
        validate="one_to_one",
    )
    merged["testable"] = finite(merged["dog_NES"]) & finite(merged["human_NES"])
    merged["observed_concordant"] = np.where(
        merged["testable"],
        np.sign(pd.to_numeric(merged["dog_NES"], errors="coerce"))
        == np.sign(pd.to_numeric(merged["human_NES"], errors="coerce")),
        pd.NA,
    )
    lodo_small = lodo[["lineage", "program", "stability_class", "minimum_tested_lodo_n", "minimum_same_direction_n"]]
    merged = merged.merge(lodo_small, on=["lineage", "program"], how="left", validate="one_to_one")
    if merged["stability_class"].isna().any():
        raise RuntimeError("LODO join failed for selected program units")

    baseline = baseline.copy()
    baseline["selected_key"] = [
        (lineage, collection, pathway_id) in selected_keys
        for lineage, collection, pathway_id in zip(
            baseline["lineage"], baseline["collection"], baseline["pathway_id"]
        )
    ]
    background = baseline.loc[
        (~baseline["selected_key"]) & finite(baseline["dog_NES"]) & finite(baseline["human_NES"])
    ].copy()
    background["candidate_key"] = (
        background["lineage"].astype(str)
        + "::" + background["collection"].astype(str)
        + "::" + background["pathway_id"].astype(str)
    )
    background["same_direction_bool"] = (
        np.sign(background["dog_NES"].to_numpy(float))
        == np.sign(background["human_NES"].to_numpy(float))
    )
    return merged, background


def build_neighbors(selected_testable: pd.DataFrame, background: pd.DataFrame) -> tuple[dict[int, pd.DataFrame], pd.DataFrame]:
    neighbors: dict[int, pd.DataFrame] = {}
    diagnostics: list[dict[str, object]] = []
    for idx, row in selected_testable.iterrows():
        candidates = background.loc[
            background["lineage"].eq(row["lineage"])
            & background["collection"].eq(row["representative_collection"])
        ].copy()
        if candidates.empty:
            raise RuntimeError(f"No background candidates for selected row {idx}")
        log_size = np.log1p(candidates["dog_size"].astype(float))
        abs_nes = candidates["dog_NES"].abs().astype(float)
        size_center, size_scale = robust_center_scale(log_size)
        nes_center, nes_scale = robust_center_scale(abs_nes)
        selected_log_size = math.log1p(float(row["dog_size"]))
        selected_abs_nes = abs(float(row["dog_NES"]))
        candidates["match_distance"] = np.sqrt(
            ((log_size - selected_log_size) / size_scale) ** 2
            + ((abs_nes - selected_abs_nes) / nes_scale) ** 2
        )
        candidates["size_log2_ratio_abs"] = np.abs(
            np.log2((candidates["dog_size"].astype(float) + 1.0) / (float(row["dog_size"]) + 1.0))
        )
        candidates["abs_dog_NES_difference"] = np.abs(abs_nes - selected_abs_nes)
        candidates = candidates.sort_values(
            ["match_distance", "pathway_id"], kind="stable"
        ).reset_index(drop=True)
        neighbors[idx] = candidates
        k20_index = min(PRIMARY_K, len(candidates)) - 1
        k20_distance = float(candidates.loc[k20_index, "match_distance"])
        diagnostics.append(
            {
                "selected_index": idx,
                "lineage": row["lineage"],
                "program": row["program"],
                "collection": row["representative_collection"],
                "pathway_id": row["representative_pathway"],
                "candidate_n": len(candidates),
                "nearest_distance": float(candidates.loc[0, "match_distance"]),
                "k20_distance": k20_distance,
                "nearest_pathway_id": candidates.loc[0, "pathway_id"],
                "nearest_pathway_name": candidates.loc[0, "pathway_name"],
                "nearest_size_log2_ratio_abs": float(candidates.loc[0, "size_log2_ratio_abs"]),
                "nearest_abs_dog_NES_difference": float(candidates.loc[0, "abs_dog_NES_difference"]),
                "match_status": (
                    "MATCH_LIMITED" if len(candidates) < 10 or k20_distance > 2.5 else "ADEQUATE"
                ),
            }
        )
    return neighbors, pd.DataFrame(diagnostics)


def draw_matched_sets(
    selected_testable: pd.DataFrame,
    neighbors: dict[int, pd.DataFrame],
    k: int,
    seed: int,
    n_draws: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    rng = np.random.default_rng(seed)
    selected_indices = selected_testable.index.to_numpy()
    rows: list[dict[str, object]] = []
    failed_draws = 0
    for draw_id in range(1, n_draws + 1):
        success = False
        chosen_rows = None
        for _ in range(MAX_ATTEMPTS_PER_DRAW):
            used: set[str] = set()
            chosen: list[pd.Series] = []
            for idx in rng.permutation(selected_indices):
                pool = neighbors[int(idx)].iloc[: min(k, len(neighbors[int(idx)]))]
                available = pool.loc[~pool["candidate_key"].isin(used)]
                if available.empty:
                    break
                pick_position = int(rng.integers(0, len(available)))
                picked = available.iloc[pick_position]
                used.add(str(picked["candidate_key"]))
                chosen.append(picked)
            if len(chosen) == len(selected_indices):
                chosen_rows = pd.DataFrame(chosen)
                success = True
                break
        if not success or chosen_rows is None:
            failed_draws += 1
            continue
        concordance = float(chosen_rows["same_direction_bool"].mean())
        set_hash = hashlib.sha256(
            "|".join(sorted(chosen_rows["candidate_key"].astype(str))).encode("utf-8")
        ).hexdigest()
        rows.append(
            {
                "matching_k": k,
                "seed": seed,
                "draw_id": draw_id,
                "concordant_n": int(chosen_rows["same_direction_bool"].sum()),
                "denominator": len(chosen_rows),
                "concordance": concordance,
                "nes_spearman_rho": safe_spearman(chosen_rows["dog_NES"], chosen_rows["human_NES"]),
                "set_hash": set_hash,
            }
        )
    result = pd.DataFrame(rows)
    metadata = {
        "matching_k": k,
        "seed": seed,
        "target_draws": n_draws,
        "valid_draws": len(result),
        "failed_draws": failed_draws,
        "unique_set_hashes": int(result["set_hash"].nunique()) if len(result) else 0,
    }
    return result, metadata


def summarize_draws(draws: pd.DataFrame, observed_s: float) -> dict[str, object]:
    if draws.empty:
        return {"valid_draws": 0}
    values = draws["concordance"].to_numpy(float)
    return {
        "valid_draws": int(len(draws)),
        "background_mean": float(np.mean(values)),
        "background_sd": float(np.std(values, ddof=1)),
        "q025": float(np.quantile(values, 0.025)),
        "q05": float(np.quantile(values, 0.05)),
        "q50": float(np.quantile(values, 0.50)),
        "q95": float(np.quantile(values, 0.95)),
        "q975": float(np.quantile(values, 0.975)),
        "observed_minus_background_mean": float(observed_s - np.mean(values)),
        "upper_tail_count": int(np.sum(values >= observed_s)),
        "empirical_upper_tail": float((1 + np.sum(values >= observed_s)) / (1 + len(values))),
        "unique_set_hashes": int(draws["set_hash"].nunique()),
    }


def subgroup_summary(frame: pd.DataFrame, group_col: str, value_col: str) -> list[dict[str, object]]:
    output = []
    for group, part in frame.groupby(group_col, sort=False):
        valid = part.loc[part["testable"]].copy()
        output.append(
            {
                "group": str(group),
                "testable_n": int(len(valid)),
                "concordant_n": int(valid[value_col].astype(bool).sum()) if len(valid) else 0,
                "concordance": float(valid[value_col].astype(bool).mean()) if len(valid) else math.nan,
                "nes_spearman_rho": safe_spearman(valid["dog_NES"], valid["human_NES"]),
            }
        )
    return output


def main() -> None:
    started = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    selected, baseline, lodo = load_inputs()
    selected_merged, background = prepare_selected(selected, baseline, lodo)
    selected_testable = selected_merged.loc[selected_merged["testable"]].copy()
    if len(selected_testable) != 23:
        raise RuntimeError(f"Expected 23 testable selected units, found {len(selected_testable)}")
    observed_s = float(selected_testable["observed_concordant"].astype(bool).mean())
    observed_n = int(selected_testable["observed_concordant"].astype(bool).sum())

    neighbors, diagnostics = build_neighbors(selected_testable, background)
    all_draws = []
    run_metadata = []
    primary_draws, primary_meta = draw_matched_sets(
        selected_testable, neighbors, PRIMARY_K, PRIMARY_SEED, N_DRAWS
    )
    all_draws.append(primary_draws)
    run_metadata.append(primary_meta)
    sensitivity_summaries = []
    for k, seed in SENSITIVITY:
        draws, meta = draw_matched_sets(selected_testable, neighbors, k, seed, N_DRAWS)
        all_draws.append(draws)
        run_metadata.append(meta)
        sensitivity_summaries.append(
            {"matching_k": k, **summarize_draws(draws, observed_s)}
        )
    draws_all = pd.concat(all_draws, ignore_index=True)

    diagnostics_for_merge = diagnostics.drop(columns=["program"]).rename(
        columns={"collection": "representative_collection", "pathway_id": "representative_pathway"}
    )
    program_table = selected_merged.merge(
        diagnostics_for_merge,
        on=["lineage", "representative_collection", "representative_pathway"],
        how="left",
        validate="many_to_one",
    )
    program_table["selected_unit_id"] = (
        program_table["lineage"].astype(str)
        + "::" + program_table["program"].astype(str)
    )
    program_table["directional_concordance"] = program_table["observed_concordant"]
    program_table["background_null_role"] = np.where(
        program_table["testable"], "PRIMARY_23_TESTABLE_UNITS", "NOT_TESTABLE_IN_COMMON_UNIVERSE"
    )
    ordered = [
        "selected_unit_id", "lineage", "program", "representative_collection",
        "representative_pathway", "dog_direction", "dog_common_direction",
        "human_direction", "dog_lodo_stability", "stability_class",
        "dog_common1to1_NES", "human_common1to1_NES", "dog_NES", "human_NES",
        "dog_FDR", "human_FDR", "dog_size", "human_size", "testable",
        "directional_concordance", "baseline_status", "candidate_n",
        "nearest_distance", "k20_distance", "match_status", "nearest_pathway_id",
        "nearest_pathway_name", "nearest_size_log2_ratio_abs",
        "nearest_abs_dog_NES_difference", "background_null_role"
    ]
    program_table = program_table[[c for c in ordered if c in program_table.columns]]
    program_table.to_csv(OUT / "02_PROGRAM_LEVEL_COMPARISON.tsv", sep="\t", index=False, na_rep="NA")
    diagnostics.to_csv(OUT / "NULL_A_MATCHING_DIAGNOSTICS.tsv", sep="\t", index=False, na_rep="NA")
    draws_all.to_csv(OUT / "NULL_A_MATCHED_SET_DRAWS.tsv.gz", sep="\t", index=False, compression="gzip")

    primary_summary = summarize_draws(primary_draws, observed_s)
    summary = {
        "selected_total_n": int(len(selected_merged)),
        "selected_testable_n": int(len(selected_testable)),
        "observed_concordant_n": observed_n,
        "observed_concordance": observed_s,
        "observed_nes_spearman_rho": safe_spearman(
            selected_testable["dog_NES"], selected_testable["human_NES"]
        ),
        "by_lineage": subgroup_summary(selected_testable, "lineage", "observed_concordant"),
        "by_collection": subgroup_summary(
            selected_testable.rename(columns={"representative_collection": "collection_group"}),
            "collection_group",
            "observed_concordant",
        ),
        "primary_matching_k": PRIMARY_K,
        "primary": primary_summary,
        "sensitivity": sensitivity_summaries,
        "matching_run_metadata": run_metadata,
        "match_limited_n": int((diagnostics["match_status"] == "MATCH_LIMITED").sum()),
        "runtime_seconds": float(time.time() - started),
        "python": sys.version,
        "platform": platform.platform(),
    }
    with (OUT / "NULL_A_SUMMARY.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, allow_nan=True)
    print(json.dumps(summary, indent=2, allow_nan=True))


if __name__ == "__main__":
    main()
