from __future__ import annotations

import csv
import hashlib
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
BASE = ROOT / "10_cross_species_baseline"
MAPPING = BASE / "mapping"
GENE_OUT = BASE / "gene_effects"
SENS_OUT = BASE / "sensitivity"
SUMMARY = BASE / "summary"
AUDIT = BASE / "audit"
RANKS = BASE / "pathway_effects" / "paired_ranks"
LOG_SHARED = ROOT / "logs" / "phase10_shared_universe.log"
LOG_GENE = ROOT / "logs" / "phase10_gene_baseline.log"
LOG_SENS = ROOT / "logs" / "phase10_sensitivity.log"

LINEAGES = ("Epithelial", "T", "Myeloid")
DOG_DE = {
    lineage: ROOT / "06_pseudobulk" / "de_results" / f"{lineage}_CIE_vs_Healthy_edger.csv"
    for lineage in LINEAGES
}
HUMAN_DE = {
    lineage: ROOT / "07_human_crohn" / "GSE266616" / "pseudobulk" / "de_results"
    / f"Human_{lineage}_ActiveCrohn_vs_Healthy.csv"
    for lineage in LINEAGES
}
HUMAN_SENS_DIR = ROOT / "07_human_crohn" / "GSE266616" / "pseudobulk" / "sensitivity"
LODO_ROOT = ROOT / "06_pseudobulk" / "lodo" / "gene_level"
LODO_STABILITY = ROOT / "06_pseudobulk" / "lodo" / "summary" / "GENE_LODO_STABILITY.csv.gz"
MASTER_PATH = MAPPING / "DOG_HUMAN_ORTHOLOG_MASTER.csv"
PHASE9_MANIFEST = ROOT / "00_metadata" / "PHASE9_MANIFEST.csv"
LOCAL_REVALIDATION = Path(
    os.environ.get(
        "CANINE_CIE_PHASE10_REVALIDATION",
        ROOT / "work" / "phase10_phase9_revalidation.csv",
    )
)

ANCHORS = {
    "Epithelial": ("KCNK16", "MMP7", "GAL3ST2", "CYP1B1"),
    "T": ("GAL3ST2", "TLR9", "TMEM255A", "PTPRT"),
    "Myeloid": ("SERPINE1", "CCL3", "ADAMDEC1", "MMP9", "SAA1"),
}


def write_log(path: Path, message: str) -> None:
    stamp = datetime.now(timezone.utc).isoformat()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="") as handle:
        handle.write(f"[{stamp}] {message}\n")
    print(message, flush=True)


def atomic_csv(frame: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        frame.to_csv(tmp, index=False, quoting=csv.QUOTE_MINIMAL, compression=compression)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_rank_formula(frame: pd.DataFrame, path: Path) -> None:
    expected = np.sign(frame["logFC"].to_numpy(float)) * np.sqrt(frame["F"].to_numpy(float))
    observed = frame["rank_metric"].to_numpy(float)
    if not np.allclose(expected, observed, rtol=1e-9, atol=1e-10, equal_nan=False):
        raise RuntimeError(f"Frozen rank formula mismatch: {path}")


def correlations(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    valid = np.isfinite(x.to_numpy(float)) & np.isfinite(y.to_numpy(float))
    if valid.sum() < 3:
        return np.nan, np.nan
    rho = float(spearmanr(x.to_numpy(float)[valid], y.to_numpy(float)[valid]).statistic)
    pearson = float(pearsonr(x.to_numpy(float)[valid], y.to_numpy(float)[valid]).statistic)
    return rho, pearson


def effect_class(dog: np.ndarray, human: np.ndarray) -> np.ndarray:
    result = np.full(len(dog), "WEAK_BOTH", dtype=object)
    dog_meaningful = np.abs(dog) >= 0.25
    human_meaningful = np.abs(human) >= 0.25
    result[dog_meaningful & ~human_meaningful] = "DOG_ONLY_EFFECT"
    result[~dog_meaningful & human_meaningful] = "HUMAN_ONLY_EFFECT"
    both = dog_meaningful & human_meaningful
    result[both & (dog > 0) & (human > 0)] = "CONSERVED_UP"
    result[both & (dog < 0) & (human < 0)] = "CONSERVED_DOWN"
    result[both & (np.sign(dog) != np.sign(human))] = "OPPOSITE_MEANINGFUL"
    return result


def treatment_files(lineage: str) -> list[Path]:
    files = sorted(HUMAN_SENS_DIR.glob(f"Human_{lineage}_*.csv"))
    return [p for p in files if p.name != "HUMAN_TREATMENT_GENE_SENSITIVITY_SUMMARY.csv"]


def human_treatment_gene_labels(lineage: str, primary: pd.DataFrame) -> dict[str, str]:
    files = treatment_files(lineage)
    if len(files) != 5:
        raise RuntimeError(f"Expected five Human treatment gene-effect files for {lineage}, got {len(files)}")
    primary_sign = primary.set_index("gene_id")["logFC"].astype(float).map(np.sign).to_dict()
    observed: dict[str, list[bool]] = {}
    for path in files:
        frame = pd.read_csv(path)
        validate_rank_formula(frame, path)
        for row in frame[["gene_id", "logFC"]].itertuples(index=False):
            if row.gene_id in primary_sign:
                observed.setdefault(row.gene_id, []).append(
                    np.sign(float(row.logFC)) == primary_sign[row.gene_id]
                )
    labels: dict[str, str] = {}
    for gene_id, values in observed.items():
        changed = len(values) - sum(values)
        if len(values) == 5 and changed == 0:
            labels[gene_id] = "ALL_5_SAME_DIRECTION"
        elif changed == 0:
            labels[gene_id] = f"TESTED_{len(values)}_OF_5_ALL_SAME_DIRECTION"
        else:
            labels[gene_id] = f"DIRECTION_CHANGE_IN_{changed}_OF_{len(values)}"
    return labels


def load_primary_pair(
    lineage: str, master: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    dog = pd.read_csv(DOG_DE[lineage])
    human = pd.read_csv(HUMAN_DE[lineage])
    validate_rank_formula(dog, DOG_DE[lineage])
    validate_rank_formula(human, HUMAN_DE[lineage])
    if dog["feature_id"].duplicated().any() or human["gene_id"].duplicated().any():
        raise RuntimeError(f"DE identity duplication: {lineage}")
    primary = master.loc[master["primary_one2one"].astype(str).str.lower().eq("true")].copy()
    pair = primary.merge(
        dog,
        left_on="dog_original_feature_id",
        right_on="feature_id",
        how="inner",
        validate="one_to_one",
    )
    pair = pair.merge(
        human,
        left_on="human_ortholog_release104_id",
        right_on="gene_id",
        how="inner",
        validate="one_to_one",
        suffixes=("_dog", "_human"),
    )
    if pair["human_ortholog_release104_id"].duplicated().any():
        raise RuntimeError(f"Shared Human ID collision: {lineage}")
    return dog, human, pair


def phase9_revalidation_audit() -> pd.DataFrame:
    manifest = pd.read_csv(PHASE9_MANIFEST, dtype=str)
    if LOCAL_REVALIDATION.is_file():
        prior = pd.read_csv(LOCAL_REVALIDATION, dtype=str)
        if len(prior) == len(manifest) and set(prior["status"]) == {"PASS"}:
            return prior
    rows = []
    for item in manifest.itertuples(index=False):
        path = ROOT / item.path
        observed_size = path.stat().st_size if path.is_file() else -1
        observed_hash = hash_file(path) if path.is_file() else ""
        status = (
            "PASS"
            if observed_size == int(item.size) and observed_hash == item.SHA256
            else "FAIL"
        )
        rows.append(
            {
                "path": item.path,
                "expected": item.SHA256,
                "observed": observed_hash,
                "size_expected": item.size,
                "size_observed": observed_size,
                "status": status,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    for path in (LOG_SHARED, LOG_GENE, LOG_SENS):
        path.unlink(missing_ok=True)
    for directory in (GENE_OUT, SENS_OUT, SUMMARY, AUDIT, RANKS):
        directory.mkdir(parents=True, exist_ok=True)

    revalidation = phase9_revalidation_audit()
    atomic_csv(revalidation, AUDIT / "PHASE9_MANIFEST_REVALIDATION.csv")
    if not revalidation["status"].eq("PASS").all():
        raise RuntimeError("Phase9 manifest revalidation failed")
    write_log(
        LOG_SHARED,
        f"Phase9 manifest revalidation PASS: {len(revalidation)}/{len(revalidation)}",
    )

    master = pd.read_csv(MASTER_PATH, dtype=str, keep_default_na=False)
    primary = master.loc[master["primary_one2one"].str.lower().eq("true")].copy()
    if len(primary) != 14_945:
        raise RuntimeError(f"Unexpected frozen Primary one2one count: {len(primary)}")
    lodo_stability = pd.read_csv(LODO_STABILITY)

    coverage_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    effects_by_lineage: dict[str, pd.DataFrame] = {}
    dog_by_lineage: dict[str, pd.DataFrame] = {}
    human_by_lineage: dict[str, pd.DataFrame] = {}

    for lineage in LINEAGES:
        dog, human, pair = load_primary_pair(lineage, master)
        dog_by_lineage[lineage] = dog
        human_by_lineage[lineage] = human
        dog_available = int(dog["feature_id"].isin(primary["dog_original_feature_id"]).sum())
        human_available = int(
            human["gene_id"].isin(primary["human_ortholog_release104_id"]).sum()
        )
        shared = len(pair)
        fraction_dog = shared / len(dog)
        fraction_human = shared / len(human)
        fraction_smaller = shared / min(len(dog), len(human))
        coverage_pass = shared >= 1500 and fraction_smaller >= 0.50
        coverage_rows.append(
            {
                "lineage": lineage,
                "dog_tested_genes": len(dog),
                "human_tested_genes": len(human),
                "dog_one2one_available": dog_available,
                "human_one2one_available": human_available,
                "shared_one2one_tested": shared,
                "fraction_of_dog_tested": fraction_dog,
                "fraction_of_human_tested": fraction_human,
                "fraction_of_smaller_tested_universe": fraction_smaller,
                "minimum_shared_n_pass": shared >= 1500,
                "minimum_fraction_pass": fraction_smaller >= 0.50,
                "coverage_gate": (
                    "PASS" if coverage_pass else "YELLOW_MAPPING_COVERAGE_LIMITED"
                ),
            }
        )

        universe = pair[
            [
                "dog_original_feature_id",
                "dog_original_symbol",
                "dog_unique_var_name",
                "human_ortholog_release104_id",
                "gene_symbol_human",
                "mapping_confidence",
                "orthology_confidence",
                "dog_percent_identity",
                "human_percent_identity",
                "goc_score_if_available",
            ]
        ].rename(
            columns={
                "dog_original_feature_id": "dog_feature_id",
                "dog_original_symbol": "dog_symbol",
                "human_ortholog_release104_id": "human_gene_id",
                "gene_symbol_human": "human_symbol",
            }
        ).sort_values("human_gene_id", kind="stable")
        atomic_csv(universe, MAPPING / f"{lineage}_SHARED_ONE2ONE_UNIVERSE.csv")

        treatment_labels = human_treatment_gene_labels(lineage, human)
        stability = lodo_stability.loc[
            lodo_stability["lineage"].eq(lineage), ["feature_id", "stability_class"]
        ]
        pair = pair.merge(stability, on="feature_id", how="left", validate="one_to_one")
        dog_logfc = pair["logFC_dog"].to_numpy(float)
        human_logfc = pair["logFC_human"].to_numpy(float)
        classes = effect_class(dog_logfc, human_logfc)
        same = np.sign(dog_logfc) == np.sign(human_logfc)
        pair["same_direction"] = same
        pair["effect_class"] = classes
        pair["human_treatment_stability_if_available"] = (
            pair["gene_id"].map(treatment_labels).fillna("NOT_AVAILABLE")
        )
        pair["common_effect_strength"] = np.minimum(
            np.abs(dog_logfc), np.abs(human_logfc)
        )

        effects = pd.DataFrame(
            {
                "human_gene_id": pair["gene_id"],
                "human_symbol": pair["gene_symbol_human"],
                "dog_feature_id": pair["feature_id"],
                "dog_symbol": pair["gene_symbol_dog"],
                "dog_logFC": pair["logFC_dog"],
                "dog_PValue": pair["PValue_dog"],
                "dog_FDR": pair["FDR_dog"],
                "dog_rank_metric": pair["rank_metric_dog"],
                "human_logFC": pair["logFC_human"],
                "human_PValue": pair["PValue_human"],
                "human_FDR": pair["FDR_human"],
                "human_rank_metric": pair["rank_metric_human"],
                "same_direction": pair["same_direction"],
                "effect_class": pair["effect_class"],
                "common_effect_strength": pair["common_effect_strength"],
                "dog_lodo_stability_if_available": pair["stability_class"].fillna(
                    "NOT_AVAILABLE"
                ),
                "human_treatment_stability_if_available": pair[
                    "human_treatment_stability_if_available"
                ],
                "mapping_confidence": pair["mapping_confidence"],
                "orthology_confidence": pair["orthology_confidence"],
                "notes": (
                    "Effect class uses |logFC|=0.25 independently of FDR; "
                    "mapping frozen before effects"
                ),
            }
        ).sort_values("human_gene_id", kind="stable")
        effects_by_lineage[lineage] = effects
        atomic_csv(
            effects,
            GENE_OUT / f"{lineage.upper()}_DOG_HUMAN_GENE_EFFECTS.csv.gz",
            compression="gzip",
        )

        dog_rank = effects[
            ["human_gene_id", "human_symbol", "dog_rank_metric", "dog_logFC"]
        ].rename(columns={"dog_rank_metric": "rank_metric", "dog_logFC": "logFC"})
        human_rank = effects[
            ["human_gene_id", "human_symbol", "human_rank_metric", "human_logFC"]
        ].rename(columns={"human_rank_metric": "rank_metric", "human_logFC": "logFC"})
        atomic_csv(
            dog_rank.sort_values("human_gene_id"),
            RANKS / f"{lineage}_Dog_common_one2one_rank.csv",
        )
        atomic_csv(
            human_rank.sort_values("human_gene_id"),
            RANKS / f"{lineage}_Human_common_one2one_rank.csv",
        )

        rho, pearson = correlations(effects["dog_logFC"], effects["human_logFC"])
        meaningful = effects["effect_class"].isin(
            ["CONSERVED_UP", "CONSERVED_DOWN", "OPPOSITE_MEANINGFUL"]
        )
        meaningful_n = int(meaningful.sum())
        meaningful_concordance = (
            float(effects.loc[meaningful, "same_direction"].mean())
            if meaningful_n
            else np.nan
        )
        summary_rows.append(
            {
                "lineage": lineage,
                "shared_genes": len(effects),
                "spearman_rho": rho,
                "pearson_r": pearson,
                "all_direction_concordance": float(effects["same_direction"].mean()),
                "meaningful_gene_n": meaningful_n,
                "meaningful_direction_concordance": meaningful_concordance,
                "conserved_up_n": int(
                    (effects["effect_class"] == "CONSERVED_UP").sum()
                ),
                "conserved_down_n": int(
                    (effects["effect_class"] == "CONSERVED_DOWN").sum()
                ),
                "opposite_meaningful_n": int(
                    (effects["effect_class"] == "OPPOSITE_MEANINGFUL").sum()
                ),
                "dog_only_n": int(
                    (effects["effect_class"] == "DOG_ONLY_EFFECT").sum()
                ),
                "human_only_n": int(
                    (effects["effect_class"] == "HUMAN_ONLY_EFFECT").sum()
                ),
                "weak_both_n": int((effects["effect_class"] == "WEAK_BOTH").sum()),
                "notes": (
                    "All shared tested one2one genes; P values are not the "
                    "primary concordance evidence"
                ),
            }
        )

        top_specs = {
            "TOP_CONSERVED_UP_GENES.csv": "CONSERVED_UP",
            "TOP_CONSERVED_DOWN_GENES.csv": "CONSERVED_DOWN",
            "TOP_OPPOSITE_GENES.csv": "OPPOSITE_MEANINGFUL",
        }
        for suffix, class_name in top_specs.items():
            top = effects.loc[effects["effect_class"].eq(class_name)].sort_values(
                ["common_effect_strength", "human_gene_id"],
                ascending=[False, True],
                kind="stable",
            ).head(25)
            atomic_csv(top, SUMMARY / f"{lineage}_{suffix}")
        write_log(
            LOG_GENE,
            f"{lineage}: shared={shared}; rho={rho:.6f}; "
            f"meaningful concordance={meaningful_concordance:.6f}",
        )

    coverage = pd.DataFrame(coverage_rows)
    gene_summary = pd.DataFrame(summary_rows)
    atomic_csv(coverage, MAPPING / "ORTHOLOG_COVERAGE_BY_LINEAGE.csv")
    atomic_csv(gene_summary, SUMMARY / "GENE_EFFECT_CONCORDANCE_SUMMARY.csv")
    if not coverage["coverage_gate"].eq("PASS").all():
        write_log(
            LOG_SHARED,
            "At least one lineage is below preregistered coverage; retained with YELLOW flag",
        )
    else:
        write_log(
            LOG_SHARED,
            "All three Primary lineages pass preregistered shared-one2one coverage",
        )

    lodo_rows: list[dict[str, object]] = []
    for lineage in LINEAGES:
        human = human_by_lineage[lineage]
        for path in sorted((LODO_ROOT / lineage).glob("*.csv")):
            dog_lodo = pd.read_csv(path)
            validate_rank_formula(dog_lodo, path)
            pair = primary.merge(
                dog_lodo,
                left_on="dog_original_feature_id",
                right_on="feature_id",
                how="inner",
            )
            pair = pair.merge(
                human,
                left_on="human_ortholog_release104_id",
                right_on="gene_id",
                how="inner",
                suffixes=("_dog", "_human"),
            )
            rho, pearson = correlations(pair["logFC_dog"], pair["logFC_human"])
            same = np.sign(pair["logFC_dog"].to_numpy(float)) == np.sign(
                pair["logFC_human"].to_numpy(float)
            )
            meaningful = (pair["logFC_dog"].abs() >= 0.25) & (
                pair["logFC_human"].abs() >= 0.25
            )
            lodo_rows.append(
                {
                    "lineage": lineage,
                    "lodo_run": pair["lodo_run"].iloc[0],
                    "short_label": pair["short_label"].iloc[0],
                    "dropped_donor": pair["dropped_donor"].iloc[0],
                    "dropped_condition": pair["dropped_condition"].iloc[0],
                    "shared_one2one_tested": len(pair),
                    "spearman_rho": rho,
                    "pearson_r": pearson,
                    "direction_concordance": float(same.mean()),
                    "meaningful_gene_n": int(meaningful.sum()),
                    "meaningful_direction_concordance": (
                        float(same[meaningful.to_numpy()].mean())
                        if meaningful.any()
                        else np.nan
                    ),
                    "comparison": "Frozen Dog LODO vs Human Primary",
                }
            )
    lodo_summary = pd.DataFrame(lodo_rows)
    if len(lodo_summary) != 21:
        raise RuntimeError(
            f"Expected 21 Dog LODO cross-species rows, got {len(lodo_summary)}"
        )
    atomic_csv(
        lodo_summary, SENS_OUT / "DOG_LODO_CROSS_SPECIES_SENSITIVITY.csv"
    )
    write_log(
        LOG_SENS, f"Dog LODO cross-species sensitivities complete: {len(lodo_summary)}"
    )

    treatment_rows: list[dict[str, object]] = []
    for lineage in LINEAGES:
        dog = dog_by_lineage[lineage]
        files = treatment_files(lineage)
        if len(files) != 5:
            raise RuntimeError(f"Expected 5 treatment models for {lineage}, got {len(files)}")
        for path in files:
            human_sens = pd.read_csv(path)
            validate_rank_formula(human_sens, path)
            pair = primary.merge(
                dog,
                left_on="dog_original_feature_id",
                right_on="feature_id",
                how="inner",
            )
            pair = pair.merge(
                human_sens,
                left_on="human_ortholog_release104_id",
                right_on="gene_id",
                how="inner",
                suffixes=("_dog", "_human"),
            )
            rho, pearson = correlations(pair["logFC_dog"], pair["logFC_human"])
            same = np.sign(pair["logFC_dog"].to_numpy(float)) == np.sign(
                pair["logFC_human"].to_numpy(float)
            )
            meaningful = (pair["logFC_dog"].abs() >= 0.25) & (
                pair["logFC_human"].abs() >= 0.25
            )
            model = str(human_sens["analysis"].iloc[0])
            treatment_rows.append(
                {
                    "lineage": lineage,
                    "human_treatment_model": model,
                    "source_file": str(path),
                    "human_n_reference": int(human_sens["n_reference"].iloc[0]),
                    "human_n_test": int(human_sens["n_test"].iloc[0]),
                    "shared_one2one_tested": len(pair),
                    "spearman_rho": rho,
                    "pearson_r": pearson,
                    "direction_concordance": float(same.mean()),
                    "meaningful_gene_n": int(meaningful.sum()),
                    "meaningful_direction_concordance": (
                        float(same[meaningful.to_numpy()].mean())
                        if meaningful.any()
                        else np.nan
                    ),
                    "comparison": (
                        "Dog Full7 vs frozen Human treatment-sensitivity model"
                    ),
                }
            )
    treatment_summary = pd.DataFrame(treatment_rows)
    if len(treatment_summary) != 15:
        raise RuntimeError(
            f"Expected 15 Human treatment cross-species rows, got {len(treatment_summary)}"
        )
    atomic_csv(
        treatment_summary,
        SENS_OUT / "HUMAN_TREATMENT_CROSS_SPECIES_SENSITIVITY.csv",
    )
    write_log(
        LOG_SENS,
        f"Human treatment cross-species sensitivities complete: {len(treatment_summary)}",
    )

    anchor_rows: list[dict[str, object]] = []
    for lineage, symbols in ANCHORS.items():
        effects = effects_by_lineage[lineage]
        for symbol in symbols:
            candidates = master.loc[master["dog_original_symbol"].eq(symbol)]
            if candidates.empty:
                anchor_rows.append(
                    {
                        "lineage": lineage,
                        "anchor": symbol,
                        "anchor_source": "PHASE3_DOG_ANCHOR",
                        "mapping_status": "DOG_ID_UNRESOLVED",
                    }
                )
                continue
            for candidate in candidates.itertuples(index=False):
                effect = effects.loc[
                    effects["dog_feature_id"].eq(candidate.dog_original_feature_id)
                ]
                values = effect.iloc[0].to_dict() if len(effect) == 1 else {}
                anchor_rows.append(
                    {
                        "lineage": lineage,
                        "anchor": symbol,
                        "anchor_source": "PHASE3_DOG_ANCHOR",
                        "dog_feature_id": candidate.dog_original_feature_id,
                        "human_gene_id": candidate.human_ortholog_release104_id,
                        "human_symbol": candidate.human_ortholog_symbol,
                        "mapping_status": candidate.mapping_status,
                        "mapping_confidence": candidate.mapping_confidence,
                        "dog_logFC": values.get("dog_logFC", np.nan),
                        "dog_FDR": values.get("dog_FDR", np.nan),
                        "human_logFC": values.get("human_logFC", np.nan),
                        "human_FDR": values.get("human_FDR", np.nan),
                        "effect_class": values.get("effect_class", "NOT_TESTABLE"),
                        "dog_lodo_stability": values.get(
                            "dog_lodo_stability_if_available", "NOT_AVAILABLE"
                        ),
                        "human_treatment_stability": values.get(
                            "human_treatment_stability_if_available", "NOT_AVAILABLE"
                        ),
                        "notes": (
                            "Anchor observed after mapping freeze; never used to select ortholog"
                        ),
                    }
                )
        human_top = effects.assign(
            abs_human_logFC=effects["human_logFC"].abs()
        ).sort_values(
            ["abs_human_logFC", "human_gene_id"],
            ascending=[False, True],
            kind="stable",
        ).head(5)
        for row in human_top.itertuples(index=False):
            anchor_rows.append(
                {
                    "lineage": lineage,
                    "anchor": row.human_symbol,
                    "anchor_source": "PHASE7_HUMAN_TOP_EFFECT_REFERENCE",
                    "dog_feature_id": row.dog_feature_id,
                    "human_gene_id": row.human_gene_id,
                    "human_symbol": row.human_symbol,
                    "mapping_status": "ONE2ONE_PRIMARY",
                    "mapping_confidence": row.mapping_confidence,
                    "dog_logFC": row.dog_logFC,
                    "dog_FDR": row.dog_FDR,
                    "human_logFC": row.human_logFC,
                    "human_FDR": row.human_FDR,
                    "effect_class": row.effect_class,
                    "dog_lodo_stability": row.dog_lodo_stability_if_available,
                    "human_treatment_stability": (
                        row.human_treatment_stability_if_available
                    ),
                    "notes": (
                        "Reference selected after unbiased full gene baseline; "
                        "not a mapping criterion"
                    ),
                }
            )
    atomic_csv(
        pd.DataFrame(anchor_rows),
        SUMMARY / "ANCHOR_GENE_CROSS_SPECIES_BASELINE.csv",
    )
    write_log(LOG_GENE, "Anchor and Human top-effect reference table complete")


if __name__ == "__main__":
    main()
