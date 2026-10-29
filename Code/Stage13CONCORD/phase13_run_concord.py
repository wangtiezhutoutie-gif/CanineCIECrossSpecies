#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import types
from datetime import datetime, timezone
from pathlib import Path

sys.modules.setdefault("resource", types.ModuleType("resource"))

import anndata as ad
import numpy as np
import pandas as pd
import torch
from scipy import sparse
import concord as ccd

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P13 = ROOT / "13_concord"
INPUT = P13 / "input" / "concord_primary_input.h5ad"
ALLOWED_OBS = ["cell_id", "species", "donor_id", "harmonized_level1", "dataset_domain"]

BASE = {
    "n_features": 10_000,
    "batch_size": 256,
    "latent_dim": 300,
    "encoder_dims": [1000],
    "n_epochs": 15,
    "element_mask_prob": 0.4,
    "feature_mask_prob": 0.2,
    "clr_temperature": 0.3,
    "clr_beta": 1.0,
    "p_intra_domain": 1.0,
    "p_intra_knn": 0.0,
    "dropout_prob": 0.0,
    "norm_type": "layer_norm",
    "lr": 1e-2,
    "schedule_ratio": 0.97,
    "use_decoder": False,
    "use_classifier": False,
    "use_faiss": True,
    "use_ivf": True,
    "load_data_into_memory": True,
    "domain_key": "dataset_domain",
    "num_workers_requested": 0,
    "num_workers_effective": 0,
    "device": "cuda:0",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def get_out(mode: str, seed: int) -> Path:
    if mode == "primary":
        return P13 / "training" / f"seed{seed}"
    if mode == "pilot":
        return P13 / "method_gate" / "project_pilot"
    if mode == "species_donor":
        return P13 / "sensitivity" / "domain_species_donor"
    if mode == "knn":
        return P13 / "sensitivity" / "knn_mode"
    raise ValueError(mode)


def prepare_adata(mode: str, seed: int) -> ad.AnnData:
    a = ad.read_h5ad(INPUT)
    if list(a.obs.columns) != ALLOWED_OBS:
        raise RuntimeError(f"Unauthorized training obs schema: {list(a.obs.columns)}")
    if a.shape != (28_046, 10_000):
        raise RuntimeError(f"Unexpected primary input shape {a.shape}")
    if not sparse.issparse(a.X) or not np.isfinite(a.X.data).all():
        raise RuntimeError("Input matrix is not finite sparse normalized expression")
    if mode == "pilot":
        rng = np.random.default_rng(seed)
        positions = []
        for _, idx in a.obs.groupby(["species", "harmonized_level1"], observed=True).indices.items():
            idx = np.asarray(idx)
            positions.extend(rng.choice(idx, size=min(500, len(idx)), replace=False).tolist())
        positions = np.asarray(sorted(positions))
        a = a[positions].copy()
    if mode == "species_donor":
        a.obs["dataset_domain"] = (
            a.obs["species"].astype(str) + "__" + a.obs["donor_id"].astype(str)
        ).astype("category")
    else:
        a.obs["dataset_domain"] = a.obs["species"].astype("category")
    a.obs["species"] = a.obs["species"].astype("category")
    a.obs["harmonized_level1"] = a.obs["harmonized_level1"].astype("category")
    return a


def run(mode: str, seed: int) -> None:
    out = get_out(mode, seed)
    out.mkdir(parents=True, exist_ok=True)
    a = prepare_adata(mode, seed)
    cfg = dict(BASE)
    cfg["seed"] = seed
    cfg["mode"] = mode
    cfg["n_cells"] = a.n_obs
    if mode == "pilot":
        cfg["n_epochs"] = 3
    if mode == "knn":
        cfg.update({
            "clr_beta": 0.0,
            "p_intra_knn": 0.3,
            "sampler_knn": 3000,
            "knn_warmup_epochs": 2,
        })
    else:
        cfg.update({
            "sampler_knn": None,
            "knn_warmup_epochs": 2,
        })

    if mode == "pilot":
        a.obs.reset_index(drop=True).to_csv(out / "PROJECT_PILOT_CELL_SAMPLING.csv", index=False)

    kwargs = {
        "input_feature": a.var_names.astype(str).tolist(),
        "domain_key": cfg["domain_key"],
        "batch_size": cfg["batch_size"],
        "n_epochs": cfg["n_epochs"],
        "lr": cfg["lr"],
        "schedule_ratio": cfg["schedule_ratio"],
        "latent_dim": cfg["latent_dim"],
        "encoder_dims": cfg["encoder_dims"],
        "element_mask_prob": cfg["element_mask_prob"],
        "feature_mask_prob": cfg["feature_mask_prob"],
        "clr_temperature": cfg["clr_temperature"],
        "clr_beta": cfg["clr_beta"],
        "p_intra_domain": cfg["p_intra_domain"],
        "p_intra_knn": cfg["p_intra_knn"],
        "sampler_knn": cfg["sampler_knn"],
        "knn_warmup_epochs": cfg["knn_warmup_epochs"],
        "dropout_prob": cfg["dropout_prob"],
        "norm_type": cfg["norm_type"],
        "use_decoder": False,
        "use_classifier": False,
        "use_faiss": cfg["use_faiss"],
        "use_ivf": cfg["use_ivf"],
        "load_data_into_memory": cfg["load_data_into_memory"],
        "num_workers": 0,
        "seed": seed,
        "device": torch.device(cfg["device"]),
    }

    with (out / "PHASE13_RUN_CONFIG.json").open("w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=2)
    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats()
    started = datetime.now(timezone.utc)
    t0 = time.perf_counter()
    model = ccd.Concord(
        a,
        save_dir=out / "model",
        copy_adata=False,
        verbose=True,
        **kwargs,
    )
    model.fit_transform(
        output_key="X_concord",
        return_decoded=False,
        return_class=False,
        return_class_prob=False,
        save_model=True,
    )
    elapsed = time.perf_counter() - t0
    z = np.asarray(model.adata.obsm["X_concord"], dtype=np.float32)
    if z.shape != (a.n_obs, cfg["latent_dim"]) or not np.isfinite(z).all():
        raise RuntimeError(f"Invalid embedding shape/values: {z.shape}")

    emb = ad.AnnData(
        X=sparse.csr_matrix((a.n_obs, 0), dtype=np.float32),
        obs=a.obs[ALLOWED_OBS].copy(),
    )
    emb.obs_names = a.obs_names.astype(str)
    emb.obsm["X_concord"] = z
    emb.uns["phase13_run_config"] = cfg
    emb_path = out / "concord_embedding.h5ad"
    emb.write_h5ad(emb_path, compression="gzip")

    summary = pd.DataFrame([{
        "mode": mode,
        "seed": seed,
        "n_cells": a.n_obs,
        "n_features": a.n_vars,
        "latent_dim": z.shape[1],
        "finite": bool(np.isfinite(z).all()),
        "latent_norm_mean": float(np.linalg.norm(z, axis=1).mean()),
        "latent_norm_sd": float(np.linalg.norm(z, axis=1).std()),
        "elapsed_seconds": elapsed,
        "cuda_peak_allocated_bytes": int(torch.cuda.max_memory_allocated()),
        "cuda_peak_reserved_bytes": int(torch.cuda.max_memory_reserved()),
        "resource_fallback": "NO",
        "embedding_sha256": sha256(emb_path),
        "started_utc": started.isoformat(),
        "completed_utc": datetime.now(timezone.utc).isoformat(),
    }])
    summary.to_csv(out / "TRAINING_SUMMARY.csv", index=False)
    print(summary.to_string(index=False), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["primary", "pilot", "species_donor", "knn"], required=True)
    parser.add_argument("--seed", type=int, default=20260827)
    args = parser.parse_args()
    run(args.mode, args.seed)


if __name__ == "__main__":
    main()

