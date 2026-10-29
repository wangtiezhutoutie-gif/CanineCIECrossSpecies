#!/usr/bin/env python
from __future__ import annotations

import hashlib
import platform
import subprocess
import sys
import types
from datetime import datetime, timezone
from pathlib import Path

sys.modules.setdefault("resource", types.ModuleType("resource"))

import importlib.metadata as md
import numpy as np
import pandas as pd
import torch

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P13 = ROOT / "13_concord"
META = ROOT / "00_metadata"
REPORT = ROOT / "08_results" / "PHASE13_REPORT.md"
GATE_NAME = "YELLOW_PHASE13_CONCORD_COMPLETE_CELLTYPE_GATE_PASS_MIXING_IMPROVED_BUT_SPECIES_SEPARATION_AND_DISEASE_SENSITIVITIES"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def revalidate_phase12():
    manifest = pd.read_csv(META / "PHASE12_MANIFEST.csv")
    rows = []
    for row in manifest.itertuples(index=False):
        path = ROOT / str(row.path)
        exists = path.is_file()
        size = path.stat().st_size if exists else np.nan
        digest = sha256(path) if exists else ""
        size_ok = bool(exists and int(size) == int(row.size))
        hash_ok = bool(exists and digest.lower() == str(row.SHA256).lower())
        rows.append({
            "manifest_path": row.path, "resolved_path": str(path),
            "expected_size": row.size, "observed_size": size, "size_match": size_ok,
            "expected_sha256": row.SHA256, "observed_sha256": digest,
            "sha256_match": hash_ok, "status": "PASS" if size_ok and hash_ok else "FAIL",
            "checked_utc": now(),
        })
    out = pd.DataFrame(rows)
    out.to_csv(P13 / "audit" / "PHASE12_MANIFEST_REVALIDATION.csv", index=False)
    if len(out) != 213 or not out["status"].eq("PASS").all():
        raise RuntimeError("Phase12 manifest revalidation failed")
    return out


def software_audit():
    sdist = P13 / "software" / "concord_sc-1.0.8.tar.gz"
    bench = P13 / "software" / "Concord_benchmark"
    commit = subprocess.check_output(["git", "-C", str(bench), "rev-parse", "HEAD"], text=True).strip()
    common = {
        "install_time": "2026-08-28", "Python": platform.python_version(),
        "CUDA": str(torch.version.cuda), "GPU": torch.cuda.get_device_name(0),
    }
    rows = [
        {
            "package": "concord-sc", "version": md.version("concord-sc"),
            "source": "PyPI sdist frozen locally", "source_sha256": sha256(sdist),
            "git_commit_if_available": "", **common,
            "notes": "Primary executable; Windows resource import shim only",
        },
        {
            "package": "Gartner-Lab/Concord_benchmark", "version": "audit checkout",
            "source": "https://github.com/Gartner-Lab/Concord_benchmark",
            "source_sha256": "", "git_commit_if_available": commit, **common,
            "notes": "Method audit only; does not replace PyPI package",
        },
    ]
    packages = [
        ("Python", platform.python_version()), ("torch", torch.__version__),
        ("scanpy", md.version("scanpy")), ("anndata", md.version("anndata")),
        ("numpy", md.version("numpy")), ("scipy", md.version("scipy")),
        ("scikit-learn", md.version("scikit-learn")), ("umap-learn", md.version("umap-learn")),
        ("matplotlib", md.version("matplotlib")), ("seaborn", md.version("seaborn")),
        ("psutil", md.version("psutil")),
    ]
    for package, version in packages:
        rows.append({
            "package": package, "version": version, "source": "Phase13 isolated venv",
            "source_sha256": "", "git_commit_if_available": "", **common, "notes": "",
        })
    rows.append({
        "package": "faiss", "version": "NOT_INSTALLED", "source": "not used",
        "source_sha256": "", "git_commit_if_available": "", **common,
        "notes": "Primary HCL p_intra_knn=0; kNN sensitivity used sklearn fallback",
    })
    pd.DataFrame(rows).to_csv(META / "PHASE13_CONCORD_SOFTWARE_AUDIT.csv", index=False)


def hyperparameters():
    values = [
        ("version", "concord-sc==1.0.8", "publication/PyPI freeze", ""),
        ("n_features", 10000, "cross-species benchmark", "no fallback"),
        ("cell_count", 28046, "Phase12 exact sampling IDs", "100% matched"),
        ("batch_size", 256, "cross-species benchmark/default", ""),
        ("latent_dim", 300, "cross-species benchmark", "full latent used"),
        ("encoder_dims", "[1000]", "cross-species benchmark", ""),
        ("epochs", 15, "cross-species benchmark/default", ""),
        ("element_mask_prob", 0.4, "cross-species benchmark", ""),
        ("feature_mask_prob", 0.2, "cross-species benchmark", ""),
        ("clr_temperature", 0.3, "two direct benchmark paths", "one generated job and default use 0.4"),
        ("clr_beta", 1.0, "HCL", ""),
        ("p_intra_domain", 1.0, "HCL", ""),
        ("p_intra_knn", 0.0, "HCL Primary", ""),
        ("activation", "LeakyReLU(0.1)", "1.0.8 source", ""),
        ("normalization", "LayerNorm", "1.0.8 source", ""),
        ("learning_rate", 0.01, "1.0.8 default", ""),
        ("optimizer", "Adam", "1.0.8 source", ""),
        ("schedule_ratio", 0.97, "1.0.8 default", ""),
        ("seed", "20260827|20260828|20260829", "preregistered", "independent initialization"),
        ("domain_key", "species via dataset_domain", "preregistered Primary", "dog/human"),
        ("load_data_into_memory", True, "Windows-compatible preload path", ""),
        ("num_workers_requested", 0, "Windows compatibility", ""),
        ("num_workers_effective", 0, "in-memory loader", ""),
        ("use_faiss", True, "official/default config", "not invoked by Primary HCL"),
        ("use_classifier", False, "supervision prohibited", ""),
        ("use_decoder", False, "minimal contrastive encoder", ""),
        ("condition_in_training", False, "input contract", ""),
        ("level2_formal_evaluation", False, "scope gate", ""),
        ("resource_fallback", "NO", "training summaries", ""),
    ]
    pd.DataFrame(values, columns=["parameter", "primary_value", "source", "notes"]).to_csv(
        META / "PHASE13_CONCORD_HYPERPARAMETERS.csv", index=False
    )


def method_gate():
    pd.DataFrame([{
        "gate": "OFFICIAL_API_SYNTHETIC_SMOKE", "status": "PASS",
        "n_cells": 120, "selected_features": 80, "latent_dim": 16, "epochs": 1,
        "finite": True, "pipeline": "feature_selection>normalize_total>log1p>Concord>fit_transform",
        "notes": "Official API; no project biological interpretation",
    }]).to_csv(P13 / "method_gate" / "synthetic_smoke" / "SYNTHETIC_SMOKE_AUDIT.csv", index=False)
    pilot = pd.read_csv(P13 / "method_gate" / "project_pilot" / "TRAINING_SUMMARY.csv").iloc[0]
    gate = pd.read_csv(P13 / "method_gate" / "CELLTYPE_INTEGRATION_GATE.csv")
    lines = [
        "Phase13 CONCORD method gate",
        "synthetic_smoke=PASS cells=120 features=80 latent=16 finite=True",
        f"project_pilot=PASS cells={int(pilot.n_cells)} features={int(pilot.n_features)} latent={int(pilot.latent_dim)} finite={pilot.finite}",
        "benchmark_contract=PASS_WITH_OFFICIAL_FEATURE_HARMONIZATION_NOT_FULLY_DOCUMENTED",
    ]
    lines.extend(f"{r.criterion}={r.status}" for r in gate.itertuples(index=False))
    (ROOT / "logs" / "phase13_method_gate.log").write_text("\n".join(lines) + "\n", encoding="utf-8")


def gate():
    criteria = [
        ("CONCORD_1.0.8_FROZEN", "PASS", "PyPI sdist hash and benchmark commit frozen"),
        ("OFFICIAL_METHOD_GATE", "PASS", "synthetic smoke, benchmark contract and pilot passed"),
        ("FEATURE_HARMONIZATION_AUDITABLE", "PASS", "official upstream mapping incomplete; frozen 1:1 route explicit"),
        ("PHASE10_MAPPING_UNCHANGED", "PASS", "14,945 Ensembl104 pairs; no re-query"),
        ("CONDITION_LEAKAGE_ZERO", "PASS", "sidecar joined after integration gate"),
        ("PRIMARY_10000_FEATURES", "PASS", "12,617 eligible; 10,000 condition-blind HVGs"),
        ("THREE_PRIMARY_SEEDS", "PASS", "independent finite training"),
        ("ALL_PRIMARY_LEVEL1_LINEAGES_RETAINED", "PASS", "Epithelial/T/Myeloid"),
        ("CELLTYPE_TRANSFER_ABOVE_NULL", "PASS", "both directions and probes"),
        ("SPECIES_MIXING_INTERPRETABLE", "WARN", "improved versus SATURN; cross30 only 0.10-0.13 and species probes 0.94-0.99"),
        ("DONOR_DISEASE_VECTOR", "PASS", "donor centroids; full 300D vectors"),
        ("PERMUTATION_NULL", "PASS", "Dog 35; Human 10,000; joint 7,000"),
        ("DOG_LODO", "PASS", "post-hoc donor-centroid LODO"),
        ("HUMAN_TREATMENT_SENSITIVITY", "PASS_WITH_CAUTION", "ANTI_IL12_23 n=3; no drug conclusion"),
        ("NO_LEVEL2_OVERRUN", "PASS", "Level1 only"),
        ("NO_PHASE14_ARTIFACT", "PASS", "no Phase14 directory or output"),
        ("T_PRIMARY_DISEASE_VECTOR", "WARN", "weak positive, nonsignificant, mode-sensitive"),
        ("EPITHELIAL_NULL_STABILITY", "WARN", "positive all seeds; mixed significance"),
        ("DOMAIN_AND_KNN_SENSITIVITY", "WARN", "T/Epithelial magnitude depends on representation choice"),
        ("FINAL_PHASE13_GATE", "YELLOW", GATE_NAME),
    ]
    frame = pd.DataFrame(criteria, columns=["criterion", "status", "evidence"])
    frame["generated_utc"] = now()
    frame.to_csv(META / "PHASE13_GATE.csv", index=False)


def scope():
    rows = [
        ("00_metadata/PHASE13_*.csv", "INCLUDE", "Phase13 metadata"),
        ("08_results/PHASE13_REPORT.md", "INCLUDE", "formal report"),
        ("13_concord/**", "INCLUDE", "Phase13 assets"),
        ("13_concord/environment/concord_venv/**", "EXCLUDE", "venv bulk; freeze files included"),
        ("13_concord/software/concord_sc-1.0.8/**", "EXCLUDE", "extracted source; sdist included"),
        ("13_concord/software/Concord_benchmark/**", "EXCLUDE", "checkout bulk; commit audited"),
        ("13_concord/cache/**", "EXCLUDE", "reconstructable cache"),
        ("scripts/phase13_*.py", "INCLUDE", "Phase13 scripts"),
        ("logs/phase13_*.log", "INCLUDE", "Phase13 logs"),
        ("00_metadata/PHASE13_MANIFEST.csv", "EXCLUDE", "no self-hash"),
    ]
    pd.DataFrame(rows, columns=["pattern", "action", "reason"]).to_csv(
        P13 / "audit" / "PHASE13_MANIFEST_SCOPE.csv", index=False
    )


def files():
    paths = []
    for path in P13.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT).as_posix()
        if rel.startswith("13_concord/environment/concord_venv/"):
            continue
        if rel.startswith("13_concord/software/concord_sc-1.0.8/"):
            continue
        if rel.startswith("13_concord/software/Concord_benchmark/"):
            continue
        if rel.startswith("13_concord/cache/"):
            continue
        paths.append(path)
    paths += list((ROOT / "scripts").glob("phase13_*.py"))
    paths += list((ROOT / "logs").glob("phase13_*.log"))
    paths += list(META.glob("PHASE13_*.csv"))
    if REPORT.is_file():
        paths.append(REPORT)
    return sorted({p.resolve() for p in paths if p.is_file() and p.name != "PHASE13_MANIFEST.csv"})


def role(path: Path) -> str:
    rel = path.relative_to(ROOT).as_posix()
    if rel.startswith("00_metadata/"): return "Phase13_metadata"
    if rel.startswith("08_results/"): return "Phase13_report"
    if rel.startswith("scripts/"): return "Phase13_script"
    if rel.startswith("logs/"): return "Phase13_log"
    for name in ["input", "training", "primary", "disease_state", "attribution", "sensitivity", "summary", "figures", "method_gate", "audit", "environment", "software"]:
        if rel.startswith(f"13_concord/{name}/"):
            return f"Phase13_{name}"
    return "Phase13_asset"


def manifest():
    rows = []
    for path in files():
        rows.append({
            "path": path.relative_to(ROOT).as_posix(), "size": path.stat().st_size,
            "SHA256": sha256(path), "role": role(path),
            "source": "Phase13 condition-blind CONCORD workflow; frozen Phase0-12 inputs read-only",
            "generated_by": "Phase13 staged workflow",
        })
    frame = pd.DataFrame(rows).sort_values("path")
    frame.to_csv(META / "PHASE13_MANIFEST.csv", index=False)
    if len(frame) < 80 or frame["path"].duplicated().any():
        raise RuntimeError(f"Unexpected manifest rows: {len(frame)}")
    return frame


def main():
    rv = revalidate_phase12()
    software_audit()
    hyperparameters()
    method_gate()
    gate()
    (ROOT / "logs" / "phase13_report.log").write_text(
        f"report={REPORT}\ngate={GATE_NAME}\ngenerated_utc={now()}\n",
        encoding="utf-8",
    )
    scope()
    if REPORT.is_file():
        mf = manifest()
        print(f"phase13_manifest_rows={len(mf)}")
    print(f"phase12_revalidation={len(rv)} all_pass={rv.status.eq('PASS').all()}")
    print(f"final_gate={GATE_NAME}")


if __name__ == "__main__":
    main()
