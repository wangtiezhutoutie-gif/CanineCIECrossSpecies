#!/usr/bin/env python3
"""Create append-only ontology-annotated Dog and Human h5ad copies."""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P9 = ROOT / "09_cell_ontology"
DOG_SOURCE = ROOT / "05_annotations" / "dog_GSE254005_annotated.h5ad"
HUMAN_SOURCE = ROOT / "07_human_crohn" / "GSE266616" / "h5ad" / "human_GSE266616_TI_annotated.h5ad"
DOG_TARGET = P9 / "dog" / "dog_GSE254005_ontology.h5ad"
HUMAN_TARGET = P9 / "human" / "human_GSE266616_TI_ontology.h5ad"
DOG_MAP = P9 / "dog" / "DOG_CELLTYPE_ONTOLOGY_MAPPING.csv"
HUMAN_MAP = P9 / "human" / "HUMAN_CELLTYPE_ONTOLOGY_MAPPING.csv"
COMPAT = P9 / "harmonized" / "ONTOLOGY_COMPATIBILITY_MATRIX.csv"
NEW_FIELDS = [
    "harmonized_level1",
    "harmonized_level2",
    "harmonized_cl_id",
    "comparison_class",
    "mapping_confidence",
]


def sha256_file(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def write_categorical(group: h5py.Group, key: str, values: pd.Series) -> None:
    if key in group:
        del group[key]
    categorical = pd.Categorical(values.fillna("").astype(str))
    out = group.create_group(key)
    out.attrs["encoding-type"] = "categorical"
    out.attrs["encoding-version"] = "0.2.0"
    out.attrs["ordered"] = bool(categorical.ordered)
    text_dtype = h5py.string_dtype(encoding="utf-8")
    categories = out.create_dataset(
        "categories",
        data=np.asarray(categorical.categories.astype(str), dtype=object),
        dtype=text_dtype,
    )
    categories.attrs["encoding-type"] = "string-array"
    categories.attrs["encoding-version"] = "0.2.0"
    codes_array = np.asarray(categorical.codes)
    max_code = max(0, len(categorical.categories) - 1)
    dtype = np.int8 if max_code <= np.iinfo(np.int8).max else np.int16
    codes = out.create_dataset("codes", data=codes_array.astype(dtype), compression="gzip")
    codes.attrs["encoding-type"] = "array"
    codes.attrs["encoding-version"] = "0.2.0"


def append_obs_fields(path: Path, fields: pd.DataFrame) -> None:
    with h5py.File(path, "r+") as handle:
        obs = handle["obs"]
        existing_order = [
            x.decode() if isinstance(x, bytes) else str(x)
            for x in obs.attrs["column-order"]
        ]
        overlap = set(NEW_FIELDS) & set(existing_order)
        if overlap:
            raise RuntimeError(f"Target unexpectedly already contains Phase9 fields: {sorted(overlap)}")
        if len(fields) != len(obs[obs.attrs["_index"]]):
            raise RuntimeError("Phase9 field length does not match h5ad obs length")
        for key in NEW_FIELDS:
            write_categorical(obs, key, fields[key])
        obs.attrs["column-order"] = np.asarray(existing_order + NEW_FIELDS, dtype=object)
        handle.flush()


def compatibility_lookup(species: str) -> tuple[dict[str, str], dict[str, str]]:
    table = pd.read_csv(COMPAT).fillna("")
    label_col = "dog_native_label" if species == "Dog" else "human_native_label"
    by_label: dict[str, str] = {}
    level1_fallback: dict[str, str] = {}
    for row in table.to_dict("records"):
        value = str(row[label_col])
        if value.startswith("ALL_") and value.endswith("_LEVEL1"):
            level1 = value[4:-7]
            level1_fallback[level1] = str(row["mapping_class"])
            continue
        for label in [x for x in value.split(";") if x]:
            by_label.setdefault(label, str(row["mapping_class"]))
    return by_label, level1_fallback


def dog_fields(obs: pd.DataFrame) -> pd.DataFrame:
    mapping = pd.read_csv(DOG_MAP).fillna("").set_index("dog_annotation_cluster")
    labels, fallback = compatibility_lookup("Dog")
    rows = mapping.loc[obs["annotation_cluster"].astype(str)].reset_index(drop=True)
    harmonized_l1 = np.where(
        rows["ontology_entity_type"].isin(["MIXED", "UNCERTAIN"]),
        "Uncertain",
        rows["dog_level1"].astype(str),
    )
    classes = [
        labels.get(str(l2), fallback.get(str(l1), "UNCERTAIN"))
        for l1, l2 in zip(rows["dog_level1"], rows["dog_level2"])
    ]
    return pd.DataFrame(
        {
            "harmonized_level1": harmonized_l1,
            "harmonized_level2": rows["cl_label"].replace("", "Uncertain").astype(str).to_numpy(),
            "harmonized_cl_id": rows["cl_id"].astype(str).to_numpy(),
            "comparison_class": classes,
            "mapping_confidence": rows["mapping_confidence"].astype(str).to_numpy(),
        },
        index=obs.index,
    )


def human_fields(obs: pd.DataFrame) -> pd.DataFrame:
    mapping = pd.read_csv(HUMAN_MAP).fillna("").set_index(["human_level1", "human_level2"])
    labels, fallback = compatibility_lookup("Human")
    keys = pd.MultiIndex.from_arrays(
        [obs["human_level1"].astype(str), obs["human_level2"].astype(str)]
    )
    rows = mapping.loc[keys].reset_index()
    harmonized_l1 = np.where(
        rows["ontology_entity_type"].isin(["MIXED", "UNCERTAIN"]),
        "Uncertain",
        rows["human_level1"].astype(str),
    )
    classes = []
    for l1, l2 in zip(rows["human_level1"], rows["human_level2"]):
        classes.append(
            fallback.get(str(l1), "UNCERTAIN")
            if str(l2) == "NA_BROAD_ONLY"
            else labels.get(str(l2), fallback.get(str(l1), "UNCERTAIN"))
        )
    return pd.DataFrame(
        {
            "harmonized_level1": harmonized_l1,
            "harmonized_level2": rows["cl_label"].replace("", "Uncertain").astype(str).to_numpy(),
            "harmonized_cl_id": rows["cl_id"].astype(str).to_numpy(),
            "comparison_class": classes,
            "mapping_confidence": rows["mapping_confidence"].astype(str).to_numpy(),
        },
        index=obs.index,
    )


def x_layout_signature(path: Path) -> str:
    """Hash X dataset metadata and edge samples without materializing X."""
    digest = hashlib.sha256()
    with h5py.File(path, "r") as handle:
        x = handle["X"]

        def visit(name, obj):
            digest.update(name.encode())
            if isinstance(obj, h5py.Dataset):
                digest.update(str(obj.shape).encode())
                digest.update(str(obj.dtype).encode())
                size = int(np.prod(obj.shape)) if obj.shape else 1
                if size:
                    if obj.ndim == 1:
                        n = min(4096, obj.shape[0])
                        digest.update(np.asarray(obj[:n]).tobytes())
                        digest.update(np.asarray(obj[-n:]).tobytes())
                    else:
                        n = min(32, obj.shape[0])
                        digest.update(np.asarray(obj[:n]).tobytes())
                        digest.update(np.asarray(obj[-n:]).tobytes())

        if isinstance(x, h5py.Dataset):
            visit("X", x)
        else:
            x.visititems(visit)
    return digest.hexdigest()


def build_one(species: str, source: Path, target: Path) -> dict[str, object]:
    source_adata = ad.read_h5ad(source, backed="r")
    source_obs = source_adata.obs.copy()
    source_shape = source_adata.shape
    source_native_columns = list(source_obs.columns)
    source_adata.file.close()
    fields = dog_fields(source_obs) if species == "Dog" else human_fields(source_obs)
    if fields.isna().any().any():
        raise RuntimeError(f"{species}: Phase9 fields contain missing values")

    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    source_sha = sha256_file(source)
    copied_sha = sha256_file(target)
    if copied_sha != source_sha:
        raise RuntimeError(f"{species}: byte-for-byte source copy verification failed")
    source_x_signature = x_layout_signature(source)
    append_obs_fields(target, fields)

    check = ad.read_h5ad(target, backed="r")
    if check.shape != source_shape:
        raise RuntimeError(f"{species}: shape changed after ontology annotation")
    if not all(field in check.obs.columns for field in NEW_FIELDS):
        raise RuntimeError(f"{species}: one or more Phase9 fields failed readback")
    if not all(column in check.obs.columns for column in source_native_columns):
        raise RuntimeError(f"{species}: native obs field was lost")
    if not check.obs_names.equals(source_obs.index):
        raise RuntimeError(f"{species}: obs index changed")
    check.file.close()
    target_x_signature = x_layout_signature(target)
    if target_x_signature != source_x_signature:
        raise RuntimeError(f"{species}: X layout/sample signature changed")
    return {
        "species": species,
        "source_path": source.relative_to(ROOT).as_posix(),
        "target_path": target.relative_to(ROOT).as_posix(),
        "n_obs": source_shape[0],
        "n_vars": source_shape[1],
        "source_sha256_current": source_sha,
        "pre_annotation_copy_sha256": copied_sha,
        "byte_identical_before_obs_append": "YES",
        "x_layout_sample_signature_source": source_x_signature,
        "x_layout_sample_signature_target": target_x_signature,
        "x_unchanged_verification": "PASS",
        "native_obs_columns_retained": len(source_native_columns),
        "new_obs_columns": ";".join(NEW_FIELDS),
        "readback_status": "PASS",
    }


def main() -> int:
    rows = [
        build_one("Dog", DOG_SOURCE, DOG_TARGET),
        build_one("Human", HUMAN_SOURCE, HUMAN_TARGET),
    ]
    audit = pd.DataFrame(rows)
    audit.to_csv(P9 / "audit" / "ONTOLOGY_OBJECT_INTEGRITY_AUDIT.csv", index=False)
    log_lines = []
    for row in rows:
        log_lines.extend(
            [
                f"{row['species']}_shape={row['n_obs']}x{row['n_vars']}",
                f"{row['species']}_x_unchanged={row['x_unchanged_verification']}",
                f"{row['species']}_readback={row['readback_status']}",
            ]
        )
    log_lines.append("status=PASS")
    (ROOT / "logs" / "phase9_ontology_objects.log").write_text(
        "\n".join(log_lines) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(rows, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
