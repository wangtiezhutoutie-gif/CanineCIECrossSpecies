#!/usr/bin/env python
"""Finalize method-faithful SATURN ESM1b gene embeddings."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import torch


EXPECTED = {
    "dog": {"unique": 18236, "proteins": 18979, "genes": 8000},
    "human": {"unique": 34791, "proteins": 41431, "genes": 8000},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def load_json(path: Path):
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def run_converter(root: Path, mapping: Path, embedding_dir: Path, output: Path) -> None:
    converter = root / "software" / "SATURN" / "protein_embeddings" / "convert_protein_embeddings_to_gene_embeddings.py"
    env = os.environ.copy()
    env["PYTHONHASHSEED"] = "0"
    env["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
    command = [
        sys.executable, str(converter),
        "--embedding_dir", str(embedding_dir),
        "--gene_symbol_to_protein_ids_path", str(mapping),
        "--embedding_model", "ESM1b",
        "--save_path", str(output),
    ]
    subprocess.run(command, cwd=converter.parent, env=env, check=True)


def ensure_isoform_aliases(
    embedding_dir: Path, protein_to_rep: dict[str, str]
) -> tuple[int, int]:
    created = 0
    preexisting = 0
    for protein_id, representative in sorted(protein_to_rep.items()):
        source = embedding_dir / f"{representative}.pt"
        target = embedding_dir / f"{protein_id}.pt"
        if not source.is_file():
            raise RuntimeError(f"Missing representative tensor: {source}")
        if target.exists():
            preexisting += 1
            continue
        os.link(source, target)
        created += 1
    return created, preexisting


def validate_unique_tensors(
    embedding_dir: Path, representatives: set[str]
) -> tuple[int, int, float, float]:
    dimensions: Counter[int] = Counter()
    finite = 0
    min_value = np.inf
    max_value = -np.inf
    for protein_id in sorted(representatives):
        payload = torch.load(
            embedding_dir / f"{protein_id}.pt",
            map_location="cpu",
            weights_only=False,
        )
        tensor = payload["mean_representations"][33].detach().cpu().float()
        dimensions[int(tensor.numel())] += 1
        array = tensor.numpy()
        if np.isfinite(array).all():
            finite += 1
        min_value = min(min_value, float(array.min()))
        max_value = max(max_value, float(array.max()))
    if dimensions != Counter({1280: len(representatives)}):
        raise RuntimeError(f"Unexpected ESM dimensions: {dimensions}")
    if finite != len(representatives):
        raise RuntimeError("At least one ESM mean representation contains NaN/Inf")
    return len(representatives), finite, min_value, max_value


def unique_symbol_test_map(
    mapping_table: Path, stable_map: dict[str, list[str]], n: int = 100
) -> tuple[dict[str, list[str]], dict[str, list[str]], dict[str, str]]:
    frame = pd.read_csv(mapping_table)
    frame = frame.loc[
        frame["ensembl_gene_id"].isin(stable_map),
        ["ensembl_gene_id", "gene_symbol"],
    ].drop_duplicates()
    frame = frame.dropna(subset=["gene_symbol"])
    frame["gene_symbol"] = frame["gene_symbol"].astype(str).str.strip()
    frame = frame.loc[frame["gene_symbol"].ne("")]
    symbol_counts = frame.groupby("gene_symbol")["ensembl_gene_id"].nunique()
    stable_counts = frame.groupby("ensembl_gene_id")["gene_symbol"].nunique()
    frame = frame.loc[
        frame["gene_symbol"].map(symbol_counts).eq(1)
        & frame["ensembl_gene_id"].map(stable_counts).eq(1)
    ].sort_values(["ensembl_gene_id", "gene_symbol"])
    frame = frame.drop_duplicates("ensembl_gene_id").head(n)
    if len(frame) != n:
        raise RuntimeError(f"Only {len(frame)} symbol-unique adapter genes available")
    stable_test = {
        row.ensembl_gene_id: stable_map[row.ensembl_gene_id]
        for row in frame.itertuples()
    }
    symbol_test = {
        row.gene_symbol: stable_map[row.ensembl_gene_id]
        for row in frame.itertuples()
    }
    symbol_to_stable = dict(zip(frame["gene_symbol"], frame["ensembl_gene_id"]))
    return stable_test, symbol_test, symbol_to_stable


def validate_gene_file(path: Path, expected_genes: int) -> tuple[int, int, float, float]:
    data = torch.load(path, map_location="cpu", weights_only=False)
    if len(data) != expected_genes:
        raise RuntimeError(f"{path}: expected {expected_genes} genes, found {len(data)}")
    dims: Counter[int] = Counter()
    finite = 0
    minimum = np.inf
    maximum = -np.inf
    for tensor in data.values():
        array = tensor.detach().cpu().float().numpy()
        dims[int(array.size)] += 1
        finite += int(np.isfinite(array).all())
        minimum = min(minimum, float(array.min()))
        maximum = max(maximum, float(array.max()))
    if dims != Counter({1280: expected_genes}) or finite != expected_genes:
        raise RuntimeError(f"Invalid gene embedding tensors in {path}: dims={dims}, finite={finite}")
    return len(data), finite, minimum, maximum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path, default=(Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"]) / "12_saturn")
    )
    args = parser.parse_args()
    root = args.root.resolve()
    sequence_dir = root / "protein_sequences"
    output_dir = root / "protein_embeddings"
    gate_dir = root / "method_gate" / "stable_id_adapter"
    gate_dir.mkdir(parents=True, exist_ok=True)
    audit_rows: list[dict[str, object]] = []
    adapter_rows: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []

    for species in ("dog", "human"):
        expected = EXPECTED[species]
        embedding_dir = output_dir / f"{species}_protein_pt"
        stable_json = sequence_dir / f"{species}_stable_id_to_protein_ids.json"
        representative_json = sequence_dir / f"{species}_protein_to_esm_representative.json"
        mapping_table = sequence_dir / f"{species}_stable_id_protein_map.csv.gz"
        stable_map = load_json(stable_json)
        protein_to_rep = load_json(representative_json)
        representatives = set(protein_to_rep.values())

        disk_before = {path.stem for path in embedding_dir.glob("*.pt")}
        if disk_before != representatives:
            missing = sorted(representatives - disk_before)[:5]
            extra = sorted(disk_before - representatives)[:5]
            raise RuntimeError(
                f"{species}: unique extraction mismatch; got {len(disk_before)}, "
                f"expected {len(representatives)}, missing={missing}, extra={extra}"
            )
        if len(representatives) != expected["unique"]:
            raise RuntimeError(f"{species}: preregistered unique count changed")
        if len(protein_to_rep) != expected["proteins"]:
            raise RuntimeError(f"{species}: preregistered protein count changed")
        if len(stable_map) != expected["genes"]:
            raise RuntimeError(f"{species}: preregistered gene count changed")

        tensor_n, finite_n, minimum, maximum = validate_unique_tensors(
            embedding_dir, representatives
        )
        created, preexisting = ensure_isoform_aliases(embedding_dir, protein_to_rep)
        all_expected = set(protein_to_rep)
        all_on_disk = {path.stem for path in embedding_dir.glob("*.pt")}
        if all_expected != all_on_disk:
            raise RuntimeError(
                f"{species}: isoform alias mismatch: expected={len(all_expected)} disk={len(all_on_disk)}"
            )

        output_file = output_dir / f"{species}_gene_embeddings_ESM1b.pt"
        run_converter(root, stable_json, embedding_dir, output_file)
        gene_n, gene_finite, gene_min, gene_max = validate_gene_file(
            output_file, expected["genes"]
        )

        stable_test, symbol_test, symbol_to_stable = unique_symbol_test_map(
            mapping_table, stable_map
        )
        stable_test_path = gate_dir / f"{species}_stable_test_mapping.json"
        symbol_test_path = gate_dir / f"{species}_symbol_test_mapping.json"
        stable_test_path.write_text(
            json.dumps(stable_test, indent=2, sort_keys=True), encoding="utf-8"
        )
        symbol_test_path.write_text(
            json.dumps(symbol_test, indent=2, sort_keys=True), encoding="utf-8"
        )
        stable_output = gate_dir / f"{species}_stable_test_embeddings.pt"
        symbol_output = gate_dir / f"{species}_symbol_test_embeddings.pt"
        run_converter(root, stable_test_path, embedding_dir, stable_output)
        run_converter(root, symbol_test_path, embedding_dir, symbol_output)
        stable_embeddings = torch.load(stable_output, map_location="cpu", weights_only=False)
        symbol_embeddings = torch.load(symbol_output, map_location="cpu", weights_only=False)
        differences = []
        for symbol, stable_id in sorted(symbol_to_stable.items()):
            diff = float(
                torch.max(
                    torch.abs(
                        symbol_embeddings[symbol].float()
                        - stable_embeddings[stable_id].float()
                    )
                ).item()
            )
            differences.append(diff)
            adapter_rows.append(
                {
                    "species": species,
                    "gene_symbol": symbol,
                    "ensembl_gene_id": stable_id,
                    "n_isoforms": len(stable_test[stable_id]),
                    "max_abs_difference": diff,
                    "tolerance": 1e-8,
                    "status": "PASS" if diff <= 1e-8 else "FAIL",
                }
            )
        if max(differences) > 1e-8:
            raise RuntimeError("RED_PROTEIN_EMBEDDING_KEY_ADAPTER_FAILURE")

        audit_rows.extend(
            [
                {
                    "species": species,
                    "check": "unique_official_ESM_outputs",
                    "observed": tensor_n,
                    "expected": expected["unique"],
                    "status": "PASS",
                    "detail": "exact official extract.py outputs before aliasing",
                },
                {
                    "species": species,
                    "check": "unique_tensor_finiteness_and_dimension",
                    "observed": finite_n,
                    "expected": expected["unique"],
                    "status": "PASS",
                    "detail": f"all dim=1280; range=[{minimum:.8g},{maximum:.8g}]",
                },
                {
                    "species": species,
                    "check": "official_isoform_averaging_input_aliases",
                    "observed": len(all_on_disk),
                    "expected": expected["proteins"],
                    "status": "PASS",
                    "detail": f"hardlinks_created={created}; preexisting={preexisting}",
                },
                {
                    "species": species,
                    "check": "official_gene_converter_output",
                    "observed": gene_n,
                    "expected": expected["genes"],
                    "status": "PASS",
                    "detail": (
                        f"finite={gene_finite}; dim=1280; range=[{gene_min:.8g},{gene_max:.8g}]; "
                        f"sha256={sha256(output_file)}"
                    ),
                },
                {
                    "species": species,
                    "check": "stable_ID_adapter_numeric_equivalence",
                    "observed": f"{max(differences):.12g}",
                    "expected": "<=1e-8",
                    "status": "PASS",
                    "detail": "100 unique-symbol genes; frozen official converter on identical protein lists",
                },
            ]
        )

        source_coverage = pd.read_csv(sequence_dir / "PROTEIN_EMBEDDING_COVERAGE.csv")
        row = source_coverage.loc[source_coverage["species"].eq(species)].iloc[0].to_dict()
        row["esm_success_genes"] = gene_n
        row["embedding_status"] = "PASS_OFFICIAL_ESM1B_LAYER33_MEAN_ISOFORM_AVERAGE"
        row["unique_esm_sequences"] = expected["unique"]
        row["protein_isoforms_averaged"] = expected["proteins"]
        row["gene_embedding_dimension"] = 1280
        row["gene_embedding_sha256"] = sha256(output_file)
        coverage_rows.append(row)

    pd.DataFrame(audit_rows).to_csv(
        root / "method_gate" / "PROTEIN_EMBEDDING_VALIDATION.csv", index=False
    )
    pd.DataFrame(adapter_rows).to_csv(
        root / "method_gate" / "STABLE_ID_PROTEIN_EMBEDDING_ADAPTER_VALIDATION.csv",
        index=False,
    )
    pd.DataFrame(coverage_rows).to_csv(
        output_dir / "PROTEIN_EMBEDDING_COVERAGE.csv", index=False
    )
    with (root / "method_gate" / "phase12_gene_embedding_gate.log").open(
        "a", encoding="utf-8", newline=""
    ) as handle:
        csv.writer(handle).writerow(
            [
                datetime.now(timezone.utc).isoformat(),
                "PASS",
                "official ESM unique tensors, isoform aliases, official converter, stable-ID adapter",
            ]
        )
    print("PASS: official ESM1b gene embeddings and stable-ID adapter")


if __name__ == "__main__":
    main()
