#!/usr/bin/env python
"""Run the condition-blind Phase12 SATURN project pilot and three frozen seeds."""

from __future__ import annotations

import argparse
import csv
import os
import pickle
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


SEEDS = (20260827, 20260828, 20260829)
LINEAGES = {"Epithelial", "T", "Myeloid"}
FORMAL_CELLS = 28046
PILOT_CELLS = 300


def command_for(
    root: Path,
    work: Path,
    run_csv: Path,
    seed: int,
    pretrain_epochs: int,
    metric_epochs: int,
    batch_size: int,
    org: str,
) -> list[str]:
    return [
        sys.executable,
        str(root / "software" / "SATURN" / "train-saturn.py"),
        "--in_data", run_csv.name,
        "--device", "cuda",
        "--device_num", "0",
        "--org", org,
        "--log_dir", "tboard/",
        "--work_dir", str(work.name) + "/",
        "--seed", str(seed),
        "--in_label_col", "harmonized_level1",
        "--ref_label_col", "harmonized_level1",
        "--hv_genes", "8000",
        "--hv_span", "0.3",
        "--num_macrogenes", "2000",
        "--centroids_init_path", f"centroids_seed{seed}.pkl",
        "--embedding_model", "ESM1b",
        "--centroid_score_func", "default",
        "--hidden_dim", "256",
        "--model_dim", "256",
        "--pretrain_lr", "0.0005",
        "--pretrain_batch_size", str(batch_size),
        "--l1_penalty", "0.0",
        "--pe_sim_penalty", "0.2",
        "--pretrain_epochs", str(pretrain_epochs),
        "--batch_size", str(batch_size),
        "--metric_lr", "0.001",
        "--epochs", str(metric_epochs),
        "--polling_freq", "25",
    ]


def run_once(
    root: Path,
    run_dir: Path,
    run_csv: Path,
    seed: int,
    pretrain_epochs: int,
    metric_epochs: int,
    batch_size: int,
    org: str,
    log_path: Path,
) -> int:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "software" / "SATURN")
    env["TORCH_HOME"] = str(root / "cache" / "torch")
    env["HF_HOME"] = str(root / "cache" / "huggingface")
    env["XDG_CACHE_HOME"] = str(root / "cache" / "xdg")
    env["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
    command = command_for(
        root, run_dir, run_csv, seed, pretrain_epochs, metric_epochs,
        batch_size, org
    )
    with log_path.open("a", encoding="utf-8") as log:
        log.write(
            f"\nUTC_START={datetime.now(timezone.utc).isoformat()}\n"
            f"BATCH_SIZE={batch_size}\nCOMMAND={' '.join(command)}\n"
        )
        log.flush()
        result = subprocess.run(
            command,
            cwd=run_csv.parent,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        log.write(
            f"\nUTC_END={datetime.now(timezone.utc).isoformat()}\n"
            f"EXIT_CODE={result.returncode}\n"
        )
    return result.returncode


def discover_final(run_dir: Path, seed: int) -> Path:
    results = run_dir / "saturn_results"
    candidates = [
        path for path in results.glob(f"*_seed_{seed}.h5ad")
        if not path.name.endswith("_pretrain.h5ad") and "_ep_" not in path.name
    ]
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one final H5AD for seed {seed}; found {candidates}")
    return candidates[0]


def validate_output(
    final_path: Path, expected_cells: int, expected_macrogenes: int
) -> list[dict[str, object]]:
    adata = ad.read_h5ad(final_path)
    latent = np.asarray(adata.X)
    macro = np.asarray(adata.obsm["macrogenes"])
    labels = set(adata.obs["labels2"].astype(str))
    rows = [
        {
            "check": "final_cell_count",
            "observed": adata.n_obs,
            "expected": expected_cells,
            "status": "PASS" if adata.n_obs == expected_cells else "FAIL",
        },
        {
            "check": "latent_dimension",
            "observed": adata.n_vars,
            "expected": 256,
            "status": "PASS" if adata.n_vars == 256 else "FAIL",
        },
        {
            "check": "latent_finite",
            "observed": bool(np.isfinite(latent).all()),
            "expected": True,
            "status": "PASS" if np.isfinite(latent).all() else "FAIL",
        },
        {
            "check": "macrogene_shape",
            "observed": f"{macro.shape[0]}x{macro.shape[1]}",
            "expected": f"{expected_cells}x{expected_macrogenes}",
            "status": "PASS"
            if macro.shape == (expected_cells, expected_macrogenes)
            else "FAIL",
        },
        {
            "check": "macrogene_finite",
            "observed": bool(np.isfinite(macro).all()),
            "expected": True,
            "status": "PASS" if np.isfinite(macro).all() else "FAIL",
        },
        {
            "check": "three_primary_lineages",
            "observed": "|".join(sorted(labels)),
            "expected": "Epithelial|Myeloid|T",
            "status": "PASS" if labels == LINEAGES else "FAIL",
        },
    ]
    if any(row["status"] != "PASS" for row in rows):
        raise RuntimeError(f"SATURN output validation failed for {final_path}")
    weights_path = final_path.with_name(
        final_path.stem + "_genes_to_macrogenes.pkl"
    )
    with weights_path.open("rb") as handle:
        weights = pickle.load(handle)
    dimensions = {np.asarray(value).shape for value in weights.values()}
    finite = all(np.isfinite(np.asarray(value)).all() for value in weights.values())
    rows.extend(
        [
            {
                "check": "gene_weight_gene_count",
                "observed": len(weights),
                "expected": 16000,
                "status": "PASS" if len(weights) == 16000 else "FAIL",
            },
            {
                "check": "gene_weight_dimensions",
                "observed": "|".join(map(str, sorted(dimensions))),
                "expected": "(2000,)",
                "status": "PASS" if dimensions == {(2000,)} else "FAIL",
            },
            {
                "check": "gene_weights_finite",
                "observed": finite,
                "expected": True,
                "status": "PASS" if finite else "FAIL",
            },
        ]
    )
    if any(row["status"] != "PASS" for row in rows):
        raise RuntimeError(f"SATURN gene-weight validation failed for {final_path}")
    return rows


def append_hyperparameter_status(
    root: Path, fallback: str, actual_batch_size: int
) -> None:
    path = root / "00_metadata" / "PHASE12_SATURN_HYPERPARAMETERS.csv"
    frame = pd.read_csv(path)
    frame.loc[
        frame["parameter"].eq("resource_fallback"),
        "planned_or_actual_value",
    ] = fallback
    frame.loc[
        frame["parameter"].eq("resource_fallback"),
        "source_or_rationale",
    ] = "actual run outcome; preregistered batch-only fallback if used"
    batch_mask = frame["parameter"].isin(["pretrain_batch_size", "batch_size"])
    frame.loc[batch_mask, "planned_or_actual_value"] = actual_batch_size
    frame.loc[batch_mask, "source_or_rationale"] = (
        "common actual value across all three formal seeds"
    )
    frame.to_csv(path, index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path, default=(Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"]) / "12_saturn")
    )
    parser.add_argument(
        "--mode",
        choices=("pilot", "formal", "all"),
        default="all",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    logs = root.parent / "logs"
    logs.mkdir(parents=True, exist_ok=True)

    for species in ("dog", "human"):
        embedding = root / "protein_embeddings" / f"{species}_gene_embeddings_ESM1b.pt"
        if not embedding.is_file():
            raise RuntimeError(f"Missing validated gene embedding: {embedding}")

    if args.mode in {"pilot", "all"}:
        pilot = root / "method_gate" / "project_pilot"
        pilot_run = pilot / "run.csv"
        pilot_work = pilot / "run"
        pilot_log = logs / "phase12_method_gate_project_pilot.log"
        if pilot_work.exists() and any((pilot_work / "saturn_results").glob("*_seed_20260827.h5ad")):
            final = discover_final(pilot_work, 20260827)
        else:
            code = run_once(
                root, pilot_work, pilot_run, 20260827, 1, 1, 4096,
                "Phase12ProjectPilot", pilot_log
            )
            if code != 0:
                raise RuntimeError(f"Project pilot failed with exit code {code}")
            final = discover_final(pilot_work, 20260827)
        rows = validate_output(final, PILOT_CELLS, 2000)
        for row in rows:
            row.update({"gate": "project_pilot", "seed": 20260827, "result_file": str(final)})
        pd.DataFrame(rows).to_csv(
            root / "method_gate" / "SATURN_PROJECT_PILOT_VALIDATION.csv",
            index=False,
        )

    if args.mode in {"formal", "all"}:
        primary_oom_quarantine = (
            root / "training" / f"seed{SEEDS[0]}" / "oom_attempt_batch4096"
        )
        fallback_used = primary_oom_quarantine.exists()
        actual_batch_size = 1024 if fallback_used else 4096
        all_validation: list[dict[str, object]] = []
        for seed in SEEDS:
            seed_dir = root / "training" / f"seed{seed}"
            run_csv = seed_dir / "run.csv"
            run_work = seed_dir / "run"
            log_path = logs / f"phase12_training_seed{seed}.log"
            completed = False
            if run_work.exists():
                try:
                    final = discover_final(run_work, seed)
                    completed = True
                except RuntimeError:
                    completed = False
            if not completed:
                code = run_once(
                    root, run_work, run_csv, seed, 200, 50,
                    actual_batch_size, "Phase12Primary", log_path
                )
                if code != 0:
                    log_text = log_path.read_text(
                        encoding="utf-8", errors="replace"
                    )
                    if "out of memory" not in log_text.lower():
                        raise RuntimeError(
                            f"Seed {seed} failed for a non-OOM reason; exit={code}"
                        )
                    if actual_batch_size == 1024:
                        raise RuntimeError(
                            "Fallback1 common batch_size=1024 also OOM; stop before "
                            "preregistered Fallback2 cell resampling"
                        )
                    if seed != SEEDS[0]:
                        raise RuntimeError(
                            "A later seed OOMed at batch_size=4096 after an earlier "
                            "seed completed. Stop for a common three-seed rerun at "
                            "Fallback1; never mix batch hyperparameters across seeds."
                        )
                    fallback_used = True
                    actual_batch_size = 1024
                    quarantine = seed_dir / "oom_attempt_batch4096"
                    if quarantine.exists():
                        raise RuntimeError(
                            f"Existing OOM quarantine blocks safe move: {quarantine}"
                        )
                    run_work.rename(quarantine)
                    code = run_once(
                        root, run_work, run_csv, seed, 200, 50,
                        actual_batch_size, "Phase12Primary", log_path
                    )
                    if code != 0:
                        raise RuntimeError(
                            "Fallback1 common batch_size=1024 failed; stop before "
                            "preregistered Fallback2 cell resampling"
                        )
                final = discover_final(run_work, seed)
            rows = validate_output(final, FORMAL_CELLS, 2000)
            for row in rows:
                row.update(
                    {
                        "gate": "formal_training",
                        "seed": seed,
                        "actual_batch_size": actual_batch_size,
                        "result_file": str(final),
                    }
                )
            all_validation.extend(rows)
        pd.DataFrame(all_validation).to_csv(
            root / "audit" / "SATURN_FORMAL_OUTPUT_VALIDATION.csv", index=False
        )
        append_hyperparameter_status(
            root,
            "FALLBACK1_BATCH1024_CELLS_UNCHANGED" if fallback_used else "NONE",
            actual_batch_size,
        )
    print("PASS: requested SATURN pilot/formal workflow")


if __name__ == "__main__":
    main()

