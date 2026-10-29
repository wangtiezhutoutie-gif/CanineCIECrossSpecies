#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P13 = ROOT / "13_concord"
SUMMARY = P13 / "summary"
SUMMARY.mkdir(parents=True, exist_ok=True)
LINEAGES = ["Epithelial", "T", "Myeloid"]
SEED = 20260827


def joined_inputs():
    p10 = pd.read_csv(ROOT / "10_cross_species_baseline" / "summary" / "PHASE10_PROGRAM_BASELINE_MATRIX.csv")
    p11 = pd.read_csv(ROOT / "11_singIST" / "primary" / "method_faithful" / "CELLTYPE_RECAPITULATION.csv")
    p12 = pd.read_csv(ROOT / "12_saturn" / "summary" / "PHASE12_VS_PHASE10_BASELINE.csv")
    cmet = pd.read_csv(P13 / "primary" / "CONCORD_DISEASE_VECTOR_METRICS.csv")
    cmet = cmet.loc[cmet["seed"] == SEED].copy()
    perm = pd.read_csv(P13 / "primary" / "CONCORD_DISEASE_VECTOR_PERMUTATION.csv")
    perm = perm.loc[(perm["seed"] == SEED) & (perm["test"] == "joint_disease_vector_cosine")].copy()
    stable = pd.read_csv(P13 / "sensitivity" / "CONCORD_SEED_STABILITY.csv")
    lodo = pd.read_csv(P13 / "sensitivity" / "CONCORD_DOG_LODO_DISEASE_VECTOR.csv")
    treatment = pd.read_csv(P13 / "sensitivity" / "CONCORD_HUMAN_TREATMENT_DISEASE_VECTOR.csv")
    return p10, p11, p12, cmet, perm, stable, lodo, treatment


def aggregate_provenance(p10, lineage):
    x = p10.loc[p10["lineage"] == lineage]
    def counts(col):
        return "|".join(f"{k}:{v}" for k, v in x[col].value_counts().sort_index().items())
    return {
        "region_risk_provenance": counts("region_risk"),
        "fre_provenance": counts("fre_summary_status"),
        "phase10_treatment_provenance": counts("human_treatment_stability"),
        "phase10_priority_provenance": counts("phase10_priority"),
    }


def main():
    p10, p11, p12, cmet, perm, stable, lodo, treatment = joined_inputs()
    rows10, rows11, rows12, cross = [], [], [], []
    for lineage in LINEAGES:
        b = p12.loc[p12["lineage"] == lineage].iloc[0]
        c = cmet.loc[cmet["lineage"] == lineage].iloc[0]
        n = perm.loc[perm["lineage"] == lineage].iloc[0]
        s = stable.loc[stable["lineage"] == lineage].iloc[0]
        ct = p11.loc[p11["celltype"] == lineage]
        lod = lodo.loc[(lodo["lineage"] == lineage) & (lodo["largest_influence_within_lineage"].astype(str) == "True")].iloc[0]
        tx = treatment.loc[treatment["lineage"] == lineage]
        tx_nonfull = tx.loc[tx["subset"] != "FULL_ACTIVE_CROHN"]
        tx_max = tx_nonfull.sort_values("absolute_cosine_change_vs_full", ascending=False).iloc[0]
        prov = aggregate_provenance(p10, lineage)

        if lineage == "Epithelial":
            interp10 = "Traditional gene/pathway baseline strongest; CONCORD is moderate positive with mixed null significance across seeds."
            interp12 = "CONCORD exceeds SATURN disease cosine and markedly improves 30NN mixing, but neither removes species separability."
        elif lineage == "T":
            interp10 = "Traditional baseline is near zero; CONCORD Primary is weak positive and nonsignificant, despite stronger immune attribution dimensions."
            interp12 = "SATURN showed high T cosine, whereas independent CONCORD Primary is near zero; model representation materially changes the T conclusion."
        else:
            interp10 = "Traditional baseline is modest; CONCORD gives the highest and consistently null-significant latent alignment."
            interp12 = "CONCORD is substantially stronger than SATURN for Myeloid and recovers shared inflammatory attribution structure."

        rows10.append({
            "lineage": lineage,
            "phase10_gene_rho": b["phase10_gene_rho"],
            "phase10_pathway_mean_rho": b["phase10_mean_pathway_rho"],
            "CONCORD_seed": SEED,
            "CONCORD_disease_cosine": c["cosine_similarity"],
            "CONCORD_disease_spearman": c["spearman_rho"],
            "CONCORD_joint_null_p": n["empirical_p_two_sided"],
            "CONCORD_seed_min": s["min_3seed"],
            "CONCORD_seed_median": s["median_3seed"],
            "CONCORD_seed_max": s["max_3seed"],
            "interpretation": interp10,
        })
        rows11.append({
            "lineage": lineage,
            "singIST_celltype_recapitulation_median": ct["recapitulation"].median(),
            "singIST_celltype_recapitulation_min": ct["recapitulation"].min(),
            "singIST_celltype_recapitulation_max": ct["recapitulation"].max(),
            "singIST_reference_models_validated": int((ct["reference_model_status"] != "WEAK_REFERENCE_MODEL").sum()),
            "singIST_reference_model_status": "46_OF_46_WEAK_REFERENCE_MODELS",
            "CONCORD_whole_state_cosine": c["cosine_similarity"],
            "CONCORD_joint_null_p": n["empirical_p_two_sided"],
            "comparison_scope": "singIST pathway/reference transfer versus CONCORD whole-state condition-blind latent; not a winner ranking",
        })
        rows12.append({
            "lineage": lineage,
            "SATURN_seed": b["seed"],
            "SATURN_disease_cosine": b["SATURN_embedding_cosine"],
            "SATURN_disease_spearman": b["SATURN_embedding_spearman"],
            "CONCORD_seed": SEED,
            "CONCORD_disease_cosine": c["cosine_similarity"],
            "CONCORD_disease_spearman": c["spearman_rho"],
            "CONCORD_joint_null_p": n["empirical_p_two_sided"],
            "CONCORD_seed_range": f'{s["min_3seed"]:.6f}..{s["max_3seed"]:.6f}',
            "CONCORD_largest_LODO_donor": lod["dropped_donor"],
            "CONCORD_largest_LODO_abs_change": lod["absolute_cosine_change"],
            "CONCORD_most_sensitive_treatment_subset": tx_max["subset"],
            "CONCORD_treatment_max_abs_change": tx_max["absolute_cosine_change_vs_full"],
            "representation_difference": "SATURN protein/macrogene latent; CONCORD shared-1:1 transcript coexpression latent",
            "interpretation": interp12,
        })
        cross.append({
            "lineage": lineage,
            "Phase10_gene_rho": b["phase10_gene_rho"],
            "Phase10_pathway_mean_rho": b["phase10_mean_pathway_rho"],
            "Phase11_singIST_summary": f'median={ct["recapitulation"].median():.3f}; min={ct["recapitulation"].min():.3f}; max={ct["recapitulation"].max():.3f}; all reference models weak',
            "Phase12_SATURN_cosine": b["SATURN_embedding_cosine"],
            "Phase13_CONCORD_cosine": c["cosine_similarity"],
            "Phase13_CONCORD_null_p": n["empirical_p_two_sided"],
            "Phase13_CONCORD_seed_range": f'{s["min_3seed"]:.6f}..{s["max_3seed"]:.6f}',
            "region": prov["region_risk_provenance"],
            "treatment": f'CONCORD max change {tx_max["absolute_cosine_change_vs_full"]:.6f} ({tx_max["subset"]}); {prov["phase10_treatment_provenance"]}',
            "LODO": f'CONCORD max {lod["dropped_donor"]} abs_change={lod["absolute_cosine_change"]:.6f}',
            "FRE": prov["fre_provenance"],
            "no_composite_score": "YES",
        })

    pd.DataFrame(rows10).to_csv(SUMMARY / "PHASE13_VS_PHASE10_BASELINE.csv", index=False)
    pd.DataFrame(rows11).to_csv(SUMMARY / "PHASE13_VS_PHASE11_SINGIST.csv", index=False)
    pd.DataFrame(rows12).to_csv(SUMMARY / "PHASE13_VS_PHASE12_SATURN.csv", index=False)
    pd.DataFrame(cross).to_csv(SUMMARY / "CROSS_METHOD_LINEAGE_SUMMARY.csv", index=False)

    prov = p10.copy()
    primary = cmet[["lineage", "cosine_similarity", "spearman_rho"]].rename(columns={
        "cosine_similarity": "phase13_CONCORD_cosine",
        "spearman_rho": "phase13_CONCORD_spearman",
    })
    nulls = perm[["lineage", "empirical_p_two_sided"]].rename(columns={"empirical_p_two_sided": "phase13_joint_null_p"})
    prov = prov.merge(primary, on="lineage", how="left").merge(nulls, on="lineage", how="left")
    prov["phase13_data_role"] = "POSTHOC_PROVENANCE_ONLY_NOT_MODEL_INPUT"
    prov["phase13_level"] = "LEVEL1_ONLY"
    prov["phase13_external_validation"] = "NO_INTERNAL_CROSS_DATASET_ALIGNMENT"
    prov.to_csv(SUMMARY / "PHASE13_PROVENANCE_MATRIX.csv", index=False)
    print(pd.DataFrame(cross).to_string(index=False))


if __name__ == "__main__":
    main()
