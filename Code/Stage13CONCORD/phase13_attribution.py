#!/usr/bin/env python
from __future__ import annotations

import gzip
import sys
import time
import types
from pathlib import Path

sys.modules.setdefault("resource", types.ModuleType("resource"))

import anndata as ad
import numpy as np
import pandas as pd
import torch
import concord as ccd

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P13 = ROOT / "13_concord"
INPUT = P13 / "input" / "concord_primary_input.h5ad"
MAPPING = P13 / "input" / "CONCORD_SHARED_FEATURE_MAPPING.csv"
DIMS = P13 / "attribution" / "CONCORD_DISEASE_ALIGNED_DIMENSIONS.csv"
OUT_FULL = P13 / "attribution" / "CONCORD_GENE_ATTRIBUTION.csv.gz"
OUT_TOP = P13 / "attribution" / "CONCORD_DISEASE_ALIGNED_FEATURES.csv"
OUT_AUDIT = P13 / "audit" / "CONCORD_ATTRIBUTION_API_AUDIT.csv"
SEED = 20260827
LAYER_INDEX = 6
TOP_N_PER_TAIL = 50


def checkpoint() -> Path:
    paths = sorted((P13 / "training" / f"seed{SEED}" / "model").glob("final_model_*.pt"))
    if len(paths) != 1:
        raise RuntimeError(f"Expected one checkpoint; found {paths}")
    return paths[0]


def donor_profiles(a: ad.AnnData, species: str, lineage: str):
    mask = (
        (a.obs["species"].astype(str).to_numpy() == species)
        & (a.obs["harmonized_level1"].astype(str).to_numpy() == lineage)
    )
    donor_values = a.obs["donor_id"].astype(str).to_numpy()
    donors = sorted(np.unique(donor_values[mask]))
    profiles = []
    for donor in donors:
        idx = np.flatnonzero(mask & (donor_values == donor))
        profiles.append(np.asarray(a.X[idx].mean(axis=0), dtype=np.float32).ravel())
    result = np.vstack(profiles).astype(np.float32, copy=False)
    if result.shape != (len(donors), a.n_vars) or not np.isfinite(result).all():
        raise RuntimeError(f"Invalid profiles: {species}/{lineage}")
    return result, donors


def make_model(a: ad.AnnData):
    wrapper = ccd.Concord(
        a, save_dir=None, copy_adata=False, verbose=False,
        input_feature=a.var_names.astype(str).tolist(),
        domain_key="dataset_domain", latent_dim=300, encoder_dims=[1000],
        use_decoder=False, use_classifier=False, dropout_prob=0.0,
        norm_type="layer_norm", seed=SEED, device=torch.device("cuda:0"),
    )
    wrapper.init_model()
    state = torch.load(checkpoint(), map_location="cuda:0", weights_only=True)
    if not isinstance(state, dict) or not all(torch.isfinite(v).all() for v in state.values()):
        raise RuntimeError("Checkpoint is not a finite state_dict")
    wrapper.model.load_state_dict(state)
    wrapper.model.eval()
    return wrapper


def main() -> None:
    started = time.perf_counter()
    a = ad.read_h5ad(INPUT)
    if a.shape != (28046, 10000) or "condition" in a.obs.columns:
        raise RuntimeError("Input shape or leakage audit failed")
    dims = pd.read_csv(DIMS)
    selected = dims.loc[dims["selection_class"] != "NOT_SELECTED"].copy()
    expected_counts = {"Epithelial": 17, "Myeloid": 15, "T": 18}
    if selected.groupby("lineage")["latent_index"].nunique().to_dict() != expected_counts:
        raise RuntimeError("Frozen dimension counts changed")
    mapping = pd.read_csv(MAPPING, dtype=str)
    mapping = mapping.set_index("human_feature_id", drop=False).loc[a.var_names.astype(str)].reset_index(drop=True)
    if len(mapping) != a.n_vars or mapping["human_feature_id"].duplicated().any():
        raise RuntimeError("Mapping does not match the 10k input")

    wrapper = make_model(a)
    feature_arrays = {
        "human_feature_id": mapping["human_feature_id"].astype(str).to_numpy(),
        "dog_feature_id": mapping["dog_feature_id"].astype(str).to_numpy(),
        "human_symbol": mapping["human_symbol"].fillna("").astype(str).to_numpy(),
        "dog_symbol": mapping["dog_symbol"].fillna("").astype(str).to_numpy(),
    }
    top_frames, audit_rows = [], []
    wrote_header = False
    with gzip.open(OUT_FULL, "wt", encoding="utf-8", newline="") as gz:
        for lineage in ["Epithelial", "T", "Myeloid"]:
            dim_meta = selected.loc[selected["lineage"] == lineage].sort_values("latent_index")
            for species in ["dog", "human"]:
                profiles, donors = donor_profiles(a, species, lineage)
                x = torch.as_tensor(profiles, dtype=torch.float32, device="cuda:0")
                torch.cuda.empty_cache()
                t0 = time.perf_counter()
                importance = ccd.ul.compute_feature_importance(wrapper.model, x, LAYER_INDEX)
                elapsed = time.perf_counter() - t0
                finite = bool(np.isfinite(importance).all())
                if importance.shape != (10000, 300) or not finite:
                    raise RuntimeError(f"Official attribution failed: {species}/{lineage}")
                audit_rows.append({
                    "seed": SEED, "lineage": lineage, "species": species,
                    "n_donor_profiles": len(donors), "donors": "|".join(donors),
                    "n_features": 10000, "n_latent_dimensions": 300,
                    "layer_index": LAYER_INDEX,
                    "official_api": "concord.ul.compute_feature_importance",
                    "input_summary": "DONOR_EQUAL_WEIGHT_MEAN_NORMALIZED_EXPRESSION_PROFILE",
                    "condition_used": "NO", "finite": finite,
                    "elapsed_seconds": elapsed,
                    "checkpoint_load": "torch.load(weights_only=True)",
                })
                for row in dim_meta.itertuples(index=False):
                    values = importance[:, int(row.latent_index)].astype(np.float64)
                    order = np.argsort(-values, kind="mergesort")
                    rank = np.empty(10000, dtype=np.int32)
                    rank[order] = np.arange(1, 10001, dtype=np.int32)
                    tail = np.full(10000, "MIDDLE", dtype=object)
                    tail[order[:TOP_N_PER_TAIL]] = "TOP_POSITIVE"
                    tail[order[-TOP_N_PER_TAIL:]] = "TOP_NEGATIVE"
                    gene = feature_arrays["dog_symbol"] if species == "dog" else feature_arrays["human_symbol"]
                    species_id = feature_arrays["dog_feature_id"] if species == "dog" else feature_arrays["human_feature_id"]
                    frame = pd.DataFrame({
                        "seed": SEED, "lineage": lineage, "species": species,
                        "latent_dimension": row.latent_dimension, "latent_index": int(row.latent_index),
                        "dog_disease_loading": float(row.dog_disease_loading),
                        "human_disease_loading": float(row.human_disease_loading),
                        "same_direction": bool(row.same_direction), "selection_class": row.selection_class,
                        **feature_arrays, "species_gene_id": species_id, "gene": gene,
                        "attribution": values, "attribution_rank": rank,
                        "attribution_tail": tail,
                        "phase10_orthology_status": "ortholog_one2one",
                        "input_summary": "DONOR_EQUAL_WEIGHT_MEAN_NORMALIZED_EXPRESSION_PROFILE",
                    })
                    frame.to_csv(gz, index=False, header=not wrote_header)
                    wrote_header = True
                    top_frames.append(frame.loc[frame["attribution_tail"] != "MIDDLE"].copy())
                del importance, x

    top = pd.concat(top_frames, ignore_index=True)
    top["pathway_annotation"] = "PENDING_FROZEN_GENESET_ENRICHMENT"
    top["notes"] = "Official CONCORD gradient attribution; representation contributor, not differential expression"
    top.to_csv(OUT_TOP, index=False)
    pd.DataFrame(audit_rows).to_csv(OUT_AUDIT, index=False)
    print(f"contexts={len(audit_rows)} full_rows={2 * 10000 * sum(expected_counts.values())}")
    print(f"top_rows={len(top)} elapsed_seconds={time.perf_counter()-started:.2f}")


if __name__ == "__main__":
    main()
