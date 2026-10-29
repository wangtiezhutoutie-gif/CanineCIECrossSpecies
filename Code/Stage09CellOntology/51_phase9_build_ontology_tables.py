#!/usr/bin/env python3
"""Build Phase9 native-to-CL mappings, compatibility, coverage, and handoffs."""

from __future__ import annotations

import json
import re
import sys
from collections import defaultdict, deque
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P9 = ROOT / "09_cell_ontology"
OBO = P9 / "reference" / "cl-base.obo"
DOG_H5AD = ROOT / "05_annotations" / "dog_GSE254005_annotated.h5ad"
DOG_ANNOT = ROOT / "05_annotations" / "dog_cell_annotation.csv"
HUMAN_H5AD = (
    ROOT
    / "07_human_crohn"
    / "GSE266616"
    / "h5ad"
    / "human_GSE266616_TI_annotated.h5ad"
)
HUMAN_MAP_SOURCE = (
    ROOT
    / "07_human_crohn"
    / "GSE266616"
    / "annotation"
    / "HUMAN_CELLTYPE_MAPPING.csv"
)
REGION_MATRIX = ROOT / "08_region_control" / "summary" / "REGION_INTERPRETATION_MATRIX.csv"
DOG_FORWARD = (
    ROOT
    / "07_human_crohn"
    / "GSE266616"
    / "gsea"
    / "DOG_FORWARD_IN_HUMAN_OBSERVATION.csv"
)
HUMAN_PROGRAMS = (
    ROOT
    / "07_human_crohn"
    / "GSE266616"
    / "gsea"
    / "HUMAN_DISEASE_PROGRAMS.csv"
)
TREATMENT_DIR = ROOT / "07_human_crohn" / "GSE266616" / "gsea" / "treatment_sensitivity"
PRIMARY_GSEA_DIR = ROOT / "07_human_crohn" / "GSE266616" / "gsea" / "full"


LEVEL1_CL = {
    "Epithelial": "CL:0000066",
    "T": "CL:0000084",
    "Myeloid": "CL:0000763",
    "B": "CL:0000236",
    "Plasma": "CL:0000786",
    "Mast": "CL:0000097",
    "Stromal": "CL:0000499",
    "Endothelial": "CL:0000115",
}


def parse_obo(path: Path) -> dict[str, dict[str, object]]:
    terms: dict[str, dict[str, object]] = {}
    current: dict[str, object] | None = None
    for raw in path.read_text(encoding="utf-8").splitlines() + ["[Term]"]:
        if raw == "[Term]":
            if current and current.get("id") and current.get("name"):
                terms[str(current["id"])] = current
            current = {"parents": [], "obsolete": False, "synonyms": []}
            continue
        if current is None:
            continue
        if raw.startswith("id: "):
            current["id"] = raw[4:]
        elif raw.startswith("name: "):
            current["name"] = raw[6:]
        elif raw.startswith("is_a: "):
            parent = raw[6:].split(" ! ")[0].split(" {")[0]
            current["parents"].append(parent)
        elif raw.startswith("is_obsolete: "):
            current["obsolete"] = raw.split(": ", 1)[1].strip().lower() == "true"
        elif raw.startswith("synonym: "):
            match = re.search(r'"([^"]+)"', raw)
            if match:
                current["synonyms"].append(match.group(1))
    return terms


def ancestors(term_id: str, terms: dict[str, dict[str, object]]) -> set[str]:
    seen: set[str] = set()
    queue: deque[str] = deque([term_id])
    while queue:
        node = queue.popleft()
        for parent in terms.get(node, {}).get("parents", []):
            parent = str(parent)
            if parent not in seen:
                seen.add(parent)
                queue.append(parent)
    return seen


def first_cl_parent(term_id: str, terms: dict[str, dict[str, object]]) -> str:
    parents = [str(x) for x in terms[term_id].get("parents", []) if str(x).startswith("CL:")]
    return parents[0] if parents else ""


def rule(
    cl_id: str = "",
    parent_id: str = "",
    entity_type: str = "CELL_TYPE",
    state: str = "",
    confidence: str = "HIGH",
    basis: str = "Native identity and marker evidence support an official CL term",
    comparison_level: str = "FINE_OR_PARENT",
    notes: str = "",
) -> dict[str, str]:
    return {
        "cl_id": cl_id,
        "parent_id": parent_id,
        "entity_type": entity_type,
        "state": state,
        "confidence": confidence,
        "basis": basis,
        "comparison_level": comparison_level,
        "notes": notes,
    }


DOG_LEVEL2_RULES = {
    "B": rule("CL:0000236", "CL:0000945", comparison_level="LEVEL1_BROAD"),
    "Plasma": rule("CL:0000786", "CL:0000946", comparison_level="LEVEL1_BROAD"),
    "Mast": rule("CL:0000097", "CL:0000766", comparison_level="LEVEL1_BROAD"),
    "Stromal": rule("CL:0000499", "CL:0002320", comparison_level="LEVEL1_BROAD"),
    "Enterocyte 1": rule(
        "CL:0000584", "CL:0000066", confidence="MEDIUM",
        basis="Native enterocyte markers support enterocyte; numbered subtype lacks a distinct CL term",
        notes="Project-native numbered subtype retained; official mapping stops at enterocyte",
    ),
    "Enterocyte 2": rule(
        "CL:0000584", "CL:0000066", confidence="MEDIUM",
        basis="Native enterocyte markers support enterocyte; numbered subtype lacks a distinct CL term",
        notes="Project-native numbered subtype retained; official mapping stops at enterocyte",
    ),
    "Enterocyte 3": rule(
        "CL:0000584", "CL:0000066", confidence="MEDIUM",
        basis="Native enterocyte markers support enterocyte; numbered subtype lacks a distinct CL term",
        notes="Project-native numbered subtype retained; official mapping stops at enterocyte",
    ),
    "BEST4+ enterocyte": rule("CL:4030026", "CL:0000584"),
    "Goblet": rule("CL:0019031", "CL:0000066"),
    "Tuft": rule("CL:0019032", "CL:0000066"),
    "Enteroendocrine": rule("CL:1001516", "CL:0000066"),
    "Treg": rule("CL:0000815", "CL:0000084"),
    "Memory T": rule(
        "CL:0000813", "CL:0000084", confidence="MEDIUM",
        basis="Native markers support a memory-like T identity; mapping stops at memory T cell",
    ),
    "Eosinophil": rule("CL:0000771", "CL:0000763"),
    "Macrophage": rule("CL:0000235", "CL:0000763"),
    "Neutrophil": rule("CL:0000775", "CL:0000763"),
    "cDC1": rule(
        "CL:0000990", "CL:0000451", confidence="MEDIUM",
        basis="XCR1/FLT3/CADM1 native evidence supports conventional DC; species-neutral CL mapping stops at parent",
        comparison_level="SHARED_PARENT",
    ),
    "IL22RA2 DC": rule(
        "CL:0000451", "CL:0000113", confidence="MEDIUM",
        basis="Native FCER1A/NAPSA/IL22RA2 evidence supports dendritic cell; no defensible fine CL term",
        comparison_level="SHARED_PARENT",
    ),
    "GZMA-high T": rule(
        "CL:0000084", "CL:0000542", entity_type="CELL_STATE", state="GZMA-high",
        confidence="MEDIUM", basis="T-cell identity is reliable; GZMA-high is retained as a native state",
        comparison_level="PARENT_CELL_TYPE",
    ),
    "IL7R-high T": rule(
        "CL:0000084", "CL:0000542", entity_type="CELL_STATE", state="IL7R-high",
        confidence="MEDIUM", basis="T-cell identity is reliable; IL7R-high is retained as a native state",
        comparison_level="PARENT_CELL_TYPE",
    ),
    "Cycling T": rule(
        "CL:0000084", "CL:0000542", entity_type="CELL_STATE", state="Cycling",
        confidence="MEDIUM", basis="T-cell identity is reliable; cycling is retained as a state",
        comparison_level="PARENT_CELL_TYPE",
    ),
    "Enterocyte 1 (IEG-high)": rule(
        "CL:0000584", "CL:0000066", entity_type="CELL_STATE", state="IEG-high",
        confidence="LOW", basis="Enterocyte parent is supported; IEG-high is a state and cluster is donor-dominated",
        comparison_level="PARENT_CELL_TYPE",
    ),
    "IFN enterocyte": rule(
        "CL:0000584", "CL:0000066", entity_type="CELL_STATE", state="IFN-response",
        confidence="LOW", basis="Enterocyte parent is supported; IFN is a state and clusters are donor-dominated",
        comparison_level="PARENT_CELL_TYPE",
    ),
    "LTF/KRT14 epithelial": rule(
        "CL:0000066", "CL:0000255", entity_type="CELL_STATE", state="LTF/KRT14 injury-like",
        confidence="LOW", basis="Epithelial identity is supported; state lacks a defensible dedicated CL term",
        comparison_level="PARENT_CELL_TYPE",
    ),
    "Cycling macrophage": rule(
        "CL:0000235", "CL:0000763", entity_type="CELL_STATE", state="Cycling",
        confidence="MEDIUM", basis="Macrophage identity is supported; cycling retained as a state",
        comparison_level="PARENT_CELL_TYPE",
    ),
}


HUMAN_FINE_RULES = {
    "BEST4": rule("CL:4030026", "CL:0000584"),
    "EC1-1": rule("CL:0000584", "CL:0000066", "CELL_STATE", "Published EC1-1 transcriptomic subtype", "MEDIUM", "Published reference identity; numbered subtype maps to enterocyte parent"),
    "EC1-2": rule("CL:0000584", "CL:0000066", "CELL_STATE", "Published EC1-2 transcriptomic subtype", "MEDIUM", "Published reference identity; numbered subtype maps to enterocyte parent"),
    "EC2-1": rule("CL:0000584", "CL:0000066", "CELL_STATE", "Published EC2-1 transcriptomic subtype", "MEDIUM", "Published reference identity; numbered subtype maps to enterocyte parent"),
    "EC2-2": rule("CL:0000584", "CL:0000066", "CELL_STATE", "Published EC2-2 transcriptomic subtype", "MEDIUM", "Published reference identity; numbered subtype maps to enterocyte parent"),
    "Goblet1": rule("CL:0019031", "CL:0000066", "CELL_STATE", "Published Goblet1 transcriptomic subtype", "MEDIUM", "Published reference identity; numbered subtype maps to intestine goblet cell"),
    "Goblet2": rule("CL:0019031", "CL:0000066", "CELL_STATE", "Published Goblet2 transcriptomic subtype", "MEDIUM", "Published reference identity; numbered subtype maps to intestine goblet cell"),
    "EEC": rule("CL:1001516", "CL:0000066"),
    "M-like": rule("CL:0000682", "CL:0000066"),
    "Paneth": rule("CL:0000510", "CL:0000066"),
    "Stem": rule("CL:0002250", "CL:0000066", confidence="MEDIUM", basis="Published Stem1/Stem2 references map to intestinal crypt stem cell"),
    "TA": rule("CL:0009012", "CL:0000066", confidence="MEDIUM", basis="Published TA reference maps to transit amplifying cell of small intestine"),
    "Tuft": rule("CL:0019032", "CL:0000066"),
    "CD4 Tcm": rule("CL:0000904", "CL:0000813"),
    "CD103+ CD4 Trm": rule("CL:0000813", "CL:0000084", "CELL_STATE", "CD103-positive CD4 tissue-resident memory phenotype", "MEDIUM", "Published phenotype; no gut-specific Trm CL term, mapping stops at memory T cell"),
    "CD103- CD4 Trm": rule("CL:0000813", "CL:0000084", "CELL_STATE", "CD103-negative CD4 tissue-resident memory phenotype", "MEDIUM", "Published phenotype; no gut-specific Trm CL term, mapping stops at memory T cell"),
    "CD103+ CD8 Trm": rule("CL:0000813", "CL:0000084", "CELL_STATE", "CD103-positive CD8 tissue-resident memory phenotype", "MEDIUM", "Published phenotype; no gut-specific Trm CL term, mapping stops at memory T cell"),
    "KLRG1+ CD8 Trm": rule("CL:0000813", "CL:0000084", "CELL_STATE", "KLRG1-positive CD8 tissue-resident memory phenotype", "MEDIUM", "Published phenotype; no gut-specific Trm CL term, mapping stops at memory T cell"),
    "Treg": rule("CL:0000815", "CL:0000084"),
    "gdT": rule("CL:0000798", "CL:0000084"),
    "MAIT": rule("CL:0000940", "CL:0000084"),
    "NK": rule("CL:0000623", "CL:0000542", confidence="MEDIUM", basis="Published NK label maps to natural killer cell; frozen Human Level1=T is preserved but fine type is excluded from T subtype comparison", comparison_level="SPECIES_SPECIFIC_OR_UNMATCHED"),
    "NK T": rule("CL:0000814", "CL:0000084", confidence="MEDIUM", basis="Published NK T label maps to mature NK T cell", comparison_level="SPECIES_SPECIFIC_OR_UNMATCHED"),
    "ILCs": rule("CL:0001065", "CL:0000542", confidence="MEDIUM", basis="Published ILC label maps to innate lymphoid cell; frozen Human Level1=T is preserved", comparison_level="SPECIES_SPECIFIC_OR_UNMATCHED"),
    "GC B": rule("CL:0000844", "CL:0000236"),
    "Memory B": rule("CL:0000787", "CL:0000236"),
    "Naive B": rule("CL:0000788", "CL:0000236"),
    "IgA plasma": rule("CL:0000987", "CL:0000786"),
    "IgG plasma": rule("CL:0000985", "CL:0000786"),
    "Macrophage": rule("CL:0000235", "CL:0000763"),
    "Monocyte": rule("CL:0000576", "CL:0000763"),
    "Neutrophil": rule("CL:0000775", "CL:0000763"),
    "cDC1": rule("CL:0000990", "CL:0000451", confidence="MEDIUM", basis="Published cDC1 label maps conservatively to conventional dendritic cell parent", comparison_level="SHARED_PARENT"),
    "cDC2": rule("CL:0000990", "CL:0000451", confidence="MEDIUM", basis="Published cDC2 label maps conservatively to conventional dendritic cell parent", comparison_level="SHARED_PARENT"),
    "Lymphoid DC": rule("CL:0000451", "CL:0000113", confidence="LOW", basis="Published label is rare and mapping stops at dendritic cell", comparison_level="SHARED_PARENT"),
    "Fibroblast-1": rule("CL:0000057", "CL:0000499", "CELL_STATE", "Published Fibroblast-1 transcriptomic subtype", "MEDIUM", "Published numbered subtype maps to fibroblast parent"),
    "Fibroblast-2": rule("CL:0000057", "CL:0000499", "CELL_STATE", "Published Fibroblast-2 transcriptomic subtype", "LOW", "Only two TI reference cells; maps to fibroblast parent"),
    "Fibroblast-3": rule("CL:0000057", "CL:0000499", "CELL_STATE", "Published Fibroblast-3 transcriptomic subtype", "MEDIUM", "Published numbered subtype maps to fibroblast parent"),
    "Myofibroblast": rule("CL:0000186", "CL:0002320"),
    "Pericyte": rule("CL:0000669", "CL:0002320"),
    "Contractile pericyte": rule("CL:0000669", "CL:0002320", "CELL_STATE", "Contractile", "MEDIUM", "Published contractile state maps to pericyte parent"),
    "Smooth muscle": rule("CL:0000192", "CL:0008000"),
    "Glial": rule("CL:4040002", "CL:0000125", confidence="MEDIUM", basis="Intestinal published glial cells map to enteroglial cell"),
    "Arterial": rule("CL:1000413", "CL:0000115"),
    "Venous": rule("CL:0002543", "CL:0000115"),
    "Lymphatic endothelium": rule("CL:0002138", "CL:0000115"),
    "Mast": rule("CL:0000097", "CL:0000766"),
}


def uncertain_rule(level1: str, level2: str) -> dict[str, str]:
    if "mixed" in level2.lower():
        return rule(
            entity_type="MIXED", confidence="UNCERTAIN",
            basis="Frozen native annotation explicitly identifies a mixed cluster",
            comparison_level="UNCERTAIN", notes="No formal CL ID assigned",
        )
    return rule(
        entity_type="UNCERTAIN", confidence="UNCERTAIN",
        basis="Frozen native annotation is uncertain or low-confidence unassigned",
        comparison_level="UNCERTAIN", notes="No formal CL ID assigned",
    )


def dog_rule(level1: str, level2: str) -> dict[str, str]:
    if level1 == "Uncertain" or "Uncertain" in level2:
        return uncertain_rule(level1, level2)
    if level2 not in DOG_LEVEL2_RULES:
        raise KeyError(f"Unmapped Dog native label: {level1} / {level2}")
    return DOG_LEVEL2_RULES[level2].copy()


def human_rule(level1: str, level2: str) -> dict[str, str]:
    if level1 == "Uncertain":
        return uncertain_rule(level1, level2)
    if level2 == "NA_BROAD_ONLY":
        cl_id = LEVEL1_CL[level1]
        return rule(
            cl_id,
            parent_id="",
            confidence="HIGH",
            basis="Frozen Human-native broad classifier with marker audit; no fine subtype was imputed",
            comparison_level="LEVEL1_BROAD",
            notes="Broad lineage only",
        )
    if level2 not in HUMAN_FINE_RULES:
        raise KeyError(f"Unmapped Human native label: {level1} / {level2}")
    return HUMAN_FINE_RULES[level2].copy()


def validate_rule(name: str, one: dict[str, str], terms: dict[str, dict[str, object]]) -> None:
    cl_id = one["cl_id"]
    if not cl_id:
        return
    if cl_id not in terms:
        raise KeyError(f"{name}: CL ID absent from frozen ontology: {cl_id}")
    if terms[cl_id]["obsolete"]:
        raise ValueError(f"{name}: CL ID is obsolete: {cl_id}")
    parent = one["parent_id"]
    if not parent:
        parent = first_cl_parent(cl_id, terms)
        one["parent_id"] = parent
    if parent:
        if parent not in terms:
            raise KeyError(f"{name}: parent absent from frozen ontology: {parent}")
        if terms[parent]["obsolete"]:
            raise ValueError(f"{name}: parent is obsolete: {parent}")
        if parent != cl_id and parent not in ancestors(cl_id, terms):
            raise ValueError(f"{name}: {parent} is not an ancestor of {cl_id}")


def label_of(cl_id: str, terms: dict[str, dict[str, object]]) -> str:
    return str(terms[cl_id]["name"]) if cl_id else ""


def join_unique(values) -> str:
    return ";".join(sorted({str(x) for x in values if str(x) not in {"", "nan", "NA"}}))


def mapping_confidence(native: str, ontology: str, entity: str) -> str:
    native = native.upper()
    if entity in {"MIXED", "UNCERTAIN"} or native in {"LOW", "UNCERTAIN"}:
        return "UNCERTAIN"
    rank = {"HIGH": 3, "PUBLISHED": 3, "MEDIUM": 2, "LOW": 1, "UNCERTAIN": 0}
    ontology_rank = rank.get(ontology.upper(), 2)
    native_rank = rank.get(native, 2)
    combined = min(native_rank, ontology_rank)
    return {3: "HIGH", 2: "MEDIUM", 1: "LOW", 0: "UNCERTAIN"}[combined]


def make_dog_mapping(terms) -> tuple[pd.DataFrame, pd.DataFrame]:
    source = pd.read_csv(DOG_ANNOT)
    rows = []
    audit_rows = []
    for item in source.to_dict("records"):
        one = dog_rule(str(item["level1"]), str(item["level2"]))
        validate_rule(f"Dog {item['cluster_id']}", one, terms)
        confidence = mapping_confidence(
            str(item["annotation_confidence"]), one["confidence"], one["entity_type"]
        )
        marker = str(item["positive_markers"])
        rows.append(
            {
                "dog_level1": item["level1"],
                "dog_level2": item["level2"],
                "dog_annotation_cluster": item["cluster_id"],
                "dog_annotation_confidence": item["annotation_confidence"],
                "n_cells": int(item["n_cells"]),
                "n_donors": int(item["n_donors"]),
                "native_marker_support": marker,
                "cl_id": one["cl_id"],
                "cl_label": label_of(one["cl_id"], terms),
                "cl_parent_id": one["parent_id"],
                "cl_parent_label": label_of(one["parent_id"], terms),
                "mapping_confidence": confidence,
                "mapping_basis": one["basis"],
                "comparison_level": one["comparison_level"],
                "ontology_entity_type": one["entity_type"],
                "native_state": one["state"],
                "dominant_donor": item["dominant_donor"],
                "dominant_donor_fraction": item["dominant_donor_fraction"],
                "notes": "; ".join(x for x in [one["notes"], str(item["notes"])] if x),
            }
        )
        audit_rows.append(
            {
                "species": "Dog",
                "native_label": item["level2"],
                "native_level1": item["level1"],
                "native_level2": item["level2"],
                "entity_type": one["entity_type"],
                "cell_type_parent": label_of(one["cl_id"], terms),
                "state_descriptor": one["state"],
                "cl_id": one["cl_id"],
                "mapping_confidence": confidence,
                "phase10_celltype_use": "NO" if one["entity_type"] in {"MIXED", "UNCERTAIN"} else ("PARENT_ONLY" if one["entity_type"] == "CELL_STATE" else "ONTOLOGY_ELIGIBLE_SUBJECT_TO_COVERAGE"),
                "phase10_state_use": "NO_UNMATCHED_STATE" if one["entity_type"] == "CELL_STATE" else "NOT_APPLICABLE",
                "notes": one["notes"],
            }
        )
    mapping = pd.DataFrame(rows).sort_values(["dog_level1", "dog_level2", "dog_annotation_cluster"])
    audit = pd.DataFrame(audit_rows).drop_duplicates().sort_values(["native_level1", "native_level2", "native_label"])
    return mapping, audit


def make_human_mapping(terms) -> tuple[pd.DataFrame, pd.DataFrame]:
    source_map = pd.read_csv(HUMAN_MAP_SOURCE)
    source_marker = (
        source_map.groupby(["human_level1", "human_level2"], dropna=False)["marker_support"]
        .apply(join_unique)
        .to_dict()
    )
    a = ad.read_h5ad(HUMAN_H5AD, backed="r")
    obs = a.obs.copy()
    a.file.close()
    rows = []
    audit_rows = []
    groups = obs.groupby(["human_level1", "human_level2"], observed=True, sort=True)
    for (level1, level2), frame in groups:
        level1, level2 = str(level1), str(level2)
        one = human_rule(level1, level2)
        validate_rule(f"Human {level1}/{level2}", one, terms)
        primary = frame[frame["primary_contrast_eligible"].astype(str) == "YES_SELECTED"]
        conf_counts = frame["annotation_confidence"].astype(str).value_counts()
        native_conf = "PUBLISHED" if "PUBLISHED" in conf_counts.index else str(conf_counts.index[0])
        confidence = mapping_confidence(native_conf, one["confidence"], one["entity_type"])
        marker_counts = frame["marker_support"].astype(str).value_counts().to_dict()
        marker_summary = ";".join(f"{key}={value}" for key, value in sorted(marker_counts.items()))
        supplied_markers = source_marker.get((level1, level2), "")
        rows.append(
            {
                "human_level1": level1,
                "human_level2": level2,
                "published_or_native_label": level1 if level2 == "NA_BROAD_ONLY" else level2,
                "n_cells": len(frame),
                "n_primary_donors": primary["donor"].astype(str).nunique(),
                "native_marker_support": "; ".join(x for x in [marker_summary, supplied_markers] if x),
                "cl_id": one["cl_id"],
                "cl_label": label_of(one["cl_id"], terms),
                "cl_parent_id": one["parent_id"],
                "cl_parent_label": label_of(one["parent_id"], terms),
                "mapping_confidence": confidence,
                "mapping_basis": one["basis"],
                "comparison_level": one["comparison_level"],
                "ontology_entity_type": one["entity_type"],
                "native_state": one["state"],
                "annotation_sources": join_unique(frame["annotation_source"]),
                "native_confidence_distribution": ";".join(f"{key}={value}" for key, value in conf_counts.items()),
                "notes": one["notes"],
            }
        )
        audit_rows.append(
            {
                "species": "Human",
                "native_label": level1 if level2 == "NA_BROAD_ONLY" else level2,
                "native_level1": level1,
                "native_level2": level2,
                "entity_type": one["entity_type"],
                "cell_type_parent": label_of(one["cl_id"], terms),
                "state_descriptor": one["state"],
                "cl_id": one["cl_id"],
                "mapping_confidence": confidence,
                "phase10_celltype_use": "NO" if one["entity_type"] in {"MIXED", "UNCERTAIN"} else ("PARENT_ONLY" if one["entity_type"] == "CELL_STATE" else "ONTOLOGY_ELIGIBLE_SUBJECT_TO_COVERAGE"),
                "phase10_state_use": "NO_NO_ACTIVE_CROHN_FINE_LABELS" if one["entity_type"] == "CELL_STATE" else "NOT_APPLICABLE",
                "notes": one["notes"],
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(audit_rows)


def compatibility_specs() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    def add(level1, level2, cl_id, dog, human, klass, ancestor, analysis, eligible, reason, notes=""):
        rows.append(
            {
                "harmonized_level1": level1,
                "harmonized_level2": level2,
                "harmonized_cl_id": cl_id,
                "dog_labels": dog,
                "human_labels": human,
                "mapping_class": klass,
                "shared_ancestor_cl_id": ancestor,
                "analysis_level": analysis,
                "eligible": eligible,
                "reason": reason,
                "notes": notes,
            }
        )

    for level1 in ["Epithelial", "T", "Myeloid"]:
        add(level1, "NA_LEVEL1", LEVEL1_CL[level1], ["*LEVEL1*"], ["*LEVEL1*"], "EXACT_ONTOLOGY_MATCH", LEVEL1_CL[level1], "PRIMARY_LEVEL1", "YES", "Pre-registered Primary lineage; both species have reliable broad identity and condition-level donor coverage")
    for level1 in ["B", "Plasma", "Mast"]:
        add(level1, "NA_LEVEL1", LEVEL1_CL[level1], ["*LEVEL1*"], ["*LEVEL1*"], "EXACT_ONTOLOGY_MATCH", LEVEL1_CL[level1], "SECONDARY_LEVEL1", "YES", "Both species have reliable broad identity and adequate condition-level donor coverage; not a pre-registered Primary lineage")
    add("Stromal", "NA_LEVEL1", LEVEL1_CL["Stromal"], ["*LEVEL1*"], ["*LEVEL1*"], "EXACT_ONTOLOGY_MATCH", LEVEL1_CL["Stromal"], "EXPLORATORY_LEVEL1", "EXPLORATORY_ONLY", "Ontology match is reliable, but Dog Stromal is strongly donor-dominated and absent in one donor")
    add("Endothelial", "NA_LEVEL1", LEVEL1_CL["Endothelial"], [], ["*LEVEL1*"], "SPECIES_SPECIFIC_OR_UNMATCHED", LEVEL1_CL["Endothelial"], "NOT_ELIGIBLE", "NO", "No reliable Dog Endothelial class in the frozen annotation")
    add("Uncertain", "NA_LEVEL1", "", ["*LEVEL1*"], ["*LEVEL1*"], "UNCERTAIN", "", "NOT_ELIGIBLE", "NO", "Uncertain cells receive no formal CL identity")

    add("Epithelial", "Enterocyte", "CL:0000584", ["Enterocyte 1", "Enterocyte 2", "Enterocyte 3"], ["EC1-1", "EC1-2", "EC2-1", "EC2-2"], "SHARED_PARENT_MATCH", "CL:0000584", "EXPLORATORY_LEVEL2_ONTOLOGY_ONLY", "EXPLORATORY_ONLY", "Both species map to enterocyte, but numbered native subtypes are not asserted equivalent and Human fine labels lack active-Crohn coverage")
    add("Epithelial", "BEST4+ enterocyte", "CL:4030026", ["BEST4+ enterocyte"], ["BEST4"], "EXACT_ONTOLOGY_MATCH", "CL:4030026", "EXPLORATORY_LEVEL2_ONTOLOGY_ONLY", "EXPLORATORY_ONLY", "Same official CL term, but Human fine labels occur only in Healthy published references")
    add("Epithelial", "Intestine goblet cell", "CL:0019031", ["Goblet"], ["Goblet1", "Goblet2"], "SHARED_PARENT_MATCH", "CL:0019031", "EXPLORATORY_LEVEL2_ONTOLOGY_ONLY", "EXPLORATORY_ONLY", "Dog and Human labels share intestine goblet parent; numbered Human subtypes are not forced 1:1")
    add("Epithelial", "Intestinal tuft cell", "CL:0019032", ["Tuft"], ["Tuft"], "EXACT_ONTOLOGY_MATCH", "CL:0019032", "EXPLORATORY_LEVEL2_ONTOLOGY_ONLY", "EXPLORATORY_ONLY", "Same official CL term, but Human fine labels lack active-Crohn coverage")
    add("Epithelial", "Intestinal enteroendocrine cell", "CL:1001516", ["Enteroendocrine"], ["EEC"], "EXACT_ONTOLOGY_MATCH", "CL:1001516", "EXPLORATORY_LEVEL2_ONTOLOGY_ONLY", "EXPLORATORY_ONLY", "Same official CL term, but Human fine labels lack active-Crohn coverage")
    add("T", "Regulatory T cell", "CL:0000815", ["Treg"], ["Treg"], "EXACT_ONTOLOGY_MATCH", "CL:0000815", "EXPLORATORY_LEVEL2_ONTOLOGY_ONLY", "EXPLORATORY_ONLY", "Same official CL term, but Human fine labels lack active-Crohn coverage")
    add("T", "Memory T cell", "CL:0000813", ["Memory T"], ["CD4 Tcm", "CD103+ CD4 Trm", "CD103- CD4 Trm", "CD103+ CD8 Trm", "KLRG1+ CD8 Trm"], "SHARED_PARENT_MATCH", "CL:0000813", "EXPLORATORY_LEVEL2_ONTOLOGY_ONLY", "EXPLORATORY_ONLY", "Only the memory T parent is defensible; Human fine labels lack active-Crohn coverage")
    add("Myeloid", "Macrophage", "CL:0000235", ["Macrophage"], ["Macrophage"], "EXACT_ONTOLOGY_MATCH", "CL:0000235", "EXPLORATORY_LEVEL2_ONTOLOGY_ONLY", "EXPLORATORY_ONLY", "Same official CL term, but Human fine labels lack active-Crohn coverage")
    add("Myeloid", "Neutrophil", "CL:0000775", ["Neutrophil"], ["Neutrophil"], "EXACT_ONTOLOGY_MATCH", "CL:0000775", "EXPLORATORY_LEVEL2_ONTOLOGY_ONLY", "EXPLORATORY_ONLY", "Same official CL term, but Human TI reference contains only two cells and no active-Crohn fine labels")
    add("Myeloid", "Conventional dendritic cell", "CL:0000990", ["cDC1"], ["cDC1"], "SHARED_PARENT_MATCH", "CL:0000990", "EXPLORATORY_LEVEL2_ONTOLOGY_ONLY", "EXPLORATORY_ONLY", "cDC1 names are aligned only at species-neutral conventional dendritic cell parent")
    add("Myeloid", "Dendritic cell broad", "CL:0000451", ["IL22RA2 DC"], ["cDC2", "Lymphoid DC"], "SHARED_PARENT_MATCH", "CL:0000451", "EXPLORATORY_LEVEL2_ONTOLOGY_ONLY", "EXPLORATORY_ONLY", "Only the dendritic-cell parent is defensible; Human fine labels lack active-Crohn coverage")
    add("B", "B-cell subtypes at B parent", "CL:0000236", ["B"], ["GC B", "Memory B", "Naive B"], "BROAD_LINEAGE_ONLY", "CL:0000236", "EXPLORATORY_LEVEL2_ONTOLOGY_ONLY", "EXPLORATORY_ONLY", "Dog has only broad B; Human fine B labels are Healthy-reference-only")
    add("Plasma", "Plasma-cell subtypes at plasma parent", "CL:0000786", ["Plasma"], ["IgA plasma", "IgG plasma"], "BROAD_LINEAGE_ONLY", "CL:0000786", "EXPLORATORY_LEVEL2_ONTOLOGY_ONLY", "EXPLORATORY_ONLY", "Dog has only broad Plasma; Human isotype labels are Healthy-reference-only")
    add("Stromal", "Stromal subtype context", "CL:0000499", ["Stromal"], ["Fibroblast-1", "Fibroblast-2", "Fibroblast-3", "Myofibroblast", "Pericyte", "Contractile pericyte", "Smooth muscle", "Glial"], "BROAD_LINEAGE_ONLY", "CL:0000499", "NOT_ELIGIBLE_LEVEL2", "NO", "Dog has only broad, donor-dominated Stromal; Human fine labels are heterogeneous and Healthy-reference-only")

    for level1, level2, cl_id, dog, human, reason in [
        ("Myeloid", "Eosinophil", "CL:0000771", ["Eosinophil"], [], "No Human eosinophil native label at current resolution"),
        ("Epithelial", "Intestinal crypt stem cell", "CL:0002250", [], ["Stem"], "Human-only at current resolution"),
        ("Epithelial", "Transit amplifying cell", "CL:0009012", [], ["TA"], "Human-only at current resolution"),
        ("Epithelial", "M cell of gut", "CL:0000682", [], ["M-like"], "Human-only at current resolution"),
        ("Epithelial", "Paneth cell", "CL:0000510", [], ["Paneth"], "Human-only at current resolution"),
        ("T", "Gamma-delta T cell", "CL:0000798", [], ["gdT"], "Human-only at current resolution"),
        ("T", "MAIT cell", "CL:0000940", [], ["MAIT"], "Human-only at current resolution"),
        ("T", "NK/ILC/NKT published labels", "CL:0000542", [], ["NK", "NK T", "ILCs"], "Human-only fine identities; frozen Human Level1=T is retained but fine ontology is non-equivalent"),
        ("Myeloid", "Monocyte", "CL:0000576", [], ["Monocyte"], "Human-only at current resolution"),
        ("Endothelial", "Endothelial subtypes", "CL:0000115", [], ["Arterial", "Venous", "Lymphatic endothelium"], "No Dog Endothelial native class"),
    ]:
        add(level1, level2, cl_id, dog, human, "SPECIES_SPECIFIC_OR_UNMATCHED", cl_id, "NOT_ELIGIBLE", "NO", reason)

    for level1, level2, cl_id, dog_label in [
        ("T", "GZMA-high state", "CL:0000084", "GZMA-high T"),
        ("T", "IL7R-high state", "CL:0000084", "IL7R-high T"),
        ("T", "Cycling state", "CL:0000084", "Cycling T"),
        ("Epithelial", "IEG-high state", "CL:0000584", "Enterocyte 1 (IEG-high)"),
        ("Epithelial", "IFN-response state", "CL:0000584", "IFN enterocyte"),
        ("Epithelial", "LTF/KRT14 injury-like state", "CL:0000066", "LTF/KRT14 epithelial"),
        ("Myeloid", "Cycling macrophage state", "CL:0000235", "Cycling macrophage"),
    ]:
        add(level1, level2, cl_id, [dog_label], [], "SPECIES_SPECIFIC_OR_UNMATCHED", cl_id, "NOT_ELIGIBLE_STATE", "NO", "Dog native state has no independently established Human equivalent in Phase9")
    add("Uncertain", "Mixed or uncertain native labels", "", ["Uncertain lymphoid-mixed", "Uncertain epithelial-mixed", "Uncertain T-mixed", "Uncertain"], ["NA_BROAD_ONLY"], "UNCERTAIN", "", "NOT_ELIGIBLE", "NO", "Mixed/uncertain cells are excluded from formal ontology comparison")
    return rows


def select_mask(obs: pd.DataFrame, level1: str, labels: list[str], species: str) -> pd.Series:
    l1col = "level1" if species == "Dog" else "human_level1"
    l2col = "level2" if species == "Dog" else "human_level2"
    if labels == ["*LEVEL1*"]:
        return obs[l1col].astype(str) == level1
    if not labels:
        return pd.Series(False, index=obs.index)
    if level1 == "Uncertain" and species == "Human":
        return obs[l1col].astype(str) == "Uncertain"
    return obs[l2col].astype(str).isin(labels)


def counts_for(obs: pd.DataFrame, mask: pd.Series, species: str) -> dict[str, object]:
    sub = obs.loc[mask]
    if species == "Dog":
        donor_col, cond_col = "donor_id", "condition"
        healthy = sub.loc[sub[cond_col].astype(str) == "Healthy", donor_col].astype(str).nunique()
        disease = sub.loc[sub[cond_col].astype(str) == "CIE", donor_col].astype(str).nunique()
        eligible = sub
    else:
        donor_col, cond_col = "donor", "state"
        eligible = sub.loc[sub["primary_contrast_eligible"].astype(str) == "YES_SELECTED"]
        healthy = eligible.loc[eligible[cond_col].astype(str) == "Healthy", donor_col].astype(str).nunique()
        disease = eligible.loc[eligible[cond_col].astype(str) == "Active_inflamed", donor_col].astype(str).nunique()
    donor_counts = eligible.groupby(donor_col, observed=True).size()
    return {
        "cells_total": len(sub),
        "primary_cells": len(eligible),
        "donors": eligible[donor_col].astype(str).nunique(),
        "healthy_donors": int(healthy),
        "disease_donors": int(disease),
        "minimum_cells_per_positive_donor": int(donor_counts.min()) if len(donor_counts) else 0,
    }


def confidence_for_labels(mapping: pd.DataFrame, labels: list[str], species: str, level1: str) -> str:
    rank = {"HIGH": 3, "MEDIUM": 2, "LOW": 1, "UNCERTAIN": 0}
    if species == "Dog":
        if labels == ["*LEVEL1*"]:
            subset = mapping[mapping["dog_level1"] == level1]
            subset = subset[~subset["ontology_entity_type"].isin(["MIXED", "UNCERTAIN"])]
        else:
            subset = mapping[mapping["dog_level2"].isin(labels)]
    else:
        if labels == ["*LEVEL1*"]:
            subset = mapping[(mapping["human_level1"] == level1) & (mapping["human_level2"] == "NA_BROAD_ONLY")]
        else:
            subset = mapping[mapping["human_level2"].isin(labels)]
    if subset.empty:
        return "NOT_AVAILABLE"
    minimum = min(rank.get(str(x), 0) for x in subset["mapping_confidence"])
    return {3: "HIGH", 2: "MEDIUM", 1: "LOW", 0: "UNCERTAIN"}[minimum]


def build_master_and_coverage(terms, dog_map, human_map, dog_obs, human_obs):
    master_rows, compat_rows, coverage_rows = [], [], []
    for spec in compatibility_specs():
        dm = select_mask(dog_obs, spec["harmonized_level1"], spec["dog_labels"], "Dog")
        hm = select_mask(human_obs, spec["harmonized_level1"], spec["human_labels"], "Human")
        dc = counts_for(dog_obs, dm, "Dog")
        hc = counts_for(human_obs, hm, "Human")
        dog_conf = confidence_for_labels(dog_map, spec["dog_labels"], "Dog", spec["harmonized_level1"])
        human_conf = confidence_for_labels(human_map, spec["human_labels"], "Human", spec["harmonized_level1"])
        cl_id = spec["harmonized_cl_id"]
        ancestor = spec["shared_ancestor_cl_id"]
        if cl_id:
            validate_rule(f"Master {spec['harmonized_level2']}", rule(cl_id, ancestor if ancestor != cl_id else ""), terms)
        dog_labels = ";".join(spec["dog_labels"]).replace("*LEVEL1*", f"ALL_{spec['harmonized_level1']}_LEVEL1")
        human_labels = ";".join(spec["human_labels"]).replace("*LEVEL1*", f"ALL_{spec['harmonized_level1']}_LEVEL1")
        master_rows.append(
            {
                "harmonized_level1": spec["harmonized_level1"],
                "harmonized_level2": spec["harmonized_level2"],
                "harmonized_cl_id": cl_id,
                "harmonized_cl_label": label_of(cl_id, terms),
                "dog_native_labels": dog_labels,
                "human_native_labels": human_labels,
                "dog_cells": dc["cells_total"],
                "human_cells": hc["cells_total"],
                "dog_donors": dc["donors"],
                "human_primary_donors": hc["donors"],
                "dog_mapping_confidence": dog_conf,
                "human_mapping_confidence": human_conf,
                "comparison_mapping_class": spec["mapping_class"],
                "shared_ancestor_cl_id": ancestor,
                "shared_ancestor_cl_label": label_of(ancestor, terms),
                "phase10_eligible": spec["eligible"],
                "phase10_analysis_level": spec["analysis_level"],
                "reason": spec["reason"],
                "notes": spec["notes"],
            }
        )
        dog_ids = join_unique(dog_map.loc[dog_map["dog_level2"].isin(spec["dog_labels"]), "cl_id"]) if spec["dog_labels"] != ["*LEVEL1*"] else cl_id
        human_ids = join_unique(human_map.loc[human_map["human_level2"].isin(spec["human_labels"]), "cl_id"]) if spec["human_labels"] != ["*LEVEL1*"] else cl_id
        compat_rows.append(
            {
                "dog_native_label": dog_labels,
                "human_native_label": human_labels,
                "dog_cl_id": dog_ids,
                "human_cl_id": human_ids,
                "shared_ancestor_cl_id": ancestor,
                "shared_ancestor_label": label_of(ancestor, terms),
                "mapping_class": spec["mapping_class"],
                "dog_confidence": dog_conf,
                "human_confidence": human_conf,
                "comparison_ready": spec["eligible"],
                "recommended_analysis_level": spec["analysis_level"],
                "notes": spec["reason"],
            }
        )
        combined_min = min(dc["minimum_cells_per_positive_donor"], hc["minimum_cells_per_positive_donor"]) if dc["donors"] and hc["donors"] else 0
        if spec["eligible"] == "YES":
            status = "READY_PRIMARY" if spec["analysis_level"] == "PRIMARY_LEVEL1" else "READY_SECONDARY_LEVEL1"
        elif spec["eligible"] == "EXPLORATORY_ONLY":
            status = "EXPLORATORY_ONLY"
        else:
            status = "NOT_ELIGIBLE"
        coverage_rows.append(
            {
                "harmonized_cell_type": f"{spec['harmonized_level1']}::{spec['harmonized_level2']}",
                "dog_cells": dc["cells_total"],
                "dog_donors_total": dc["donors"],
                "dog_healthy_donors": dc["healthy_donors"],
                "dog_cie_donors": dc["disease_donors"],
                "human_cells": hc["primary_cells"],
                "human_healthy_donors": hc["healthy_donors"],
                "human_active_crohn_donors": hc["disease_donors"],
                "minimum_cells_per_donor": combined_min,
                "coverage_status": status,
                "phase10_eligible": spec["eligible"],
                "notes": "Human cell count is restricted to the frozen Primary contrast; " + spec["reason"],
            }
        )
    return pd.DataFrame(master_rows), pd.DataFrame(compat_rows), pd.DataFrame(coverage_rows)


def treatment_status_for_representatives(region: pd.DataFrame) -> dict[tuple[str, str, str], str]:
    statuses = {}
    for item in region.to_dict("records"):
        lineage = str(item["lineage"])
        collection = str(item["representative_collection"])
        pathway = str(item["representative_pathway"])
        key = (lineage, collection, pathway)
        primary_path = PRIMARY_GSEA_DIR / f"Human_{lineage}_{collection}_GSEA.csv"
        primary = pd.read_csv(primary_path)
        prow = primary.loc[primary["pathway_id"].astype(str) == pathway]
        if prow.empty:
            statuses[key] = "PRIMARY_PATHWAY_NOT_FOUND"
            continue
        expected = np.sign(float(prow.iloc[0]["NES"]))
        sensitivity_rows = []
        pattern = f"Human_{lineage}_*_{collection}_GSEA.csv"
        for path in sorted(TREATMENT_DIR.glob(pattern)):
            table = pd.read_csv(path)
            match = table.loc[table["pathway_id"].astype(str) == pathway]
            if not match.empty:
                sensitivity_rows.append(match.iloc[0])
        if not sensitivity_rows:
            statuses[key] = "NOT_EXECUTED"
            continue
        same = [np.sign(float(row["NES"])) == expected for row in sensitivity_rows]
        sig_same = [ok and float(row["padj"]) < 0.05 for ok, row in zip(same, sensitivity_rows)]
        n = len(sensitivity_rows)
        if all(same) and all(sig_same):
            statuses[key] = f"ALL_{n}_FDR05_SAME_DIRECTION"
        elif all(same):
            statuses[key] = f"ALL_{n}_SAME_DIRECTION_DIRECTIONALLY_STABLE"
        else:
            statuses[key] = f"DIRECTION_CHANGE_IN_{sum(not x for x in same)}_OF_{n}_SENSITIVITIES"
    return statuses


def build_handoffs() -> tuple[pd.DataFrame, pd.DataFrame]:
    region = pd.read_csv(REGION_MATRIX)
    dog_forward = pd.read_csv(DOG_FORWARD)
    treatment = treatment_status_for_representatives(region)
    joined = region.merge(
        dog_forward[["lineage", "program", "dog_phase5_stability", "FRE_summary_status"]],
        on=["lineage", "program"], how="left", validate="one_to_one",
    )
    joined["human_treatment_sensitivity_status"] = [
        treatment[(str(row.lineage), str(row.representative_collection), str(row.representative_pathway))]
        for row in joined.itertuples()
    ]
    joined["Phase8_status"] = "GREEN_PHASE8A_REGION_CONTROL_COMPLETE"
    joined["harmonized_level1"] = joined["lineage"]
    joined["phase9_ontology_ready"] = "YES_PRIMARY_LEVEL1"
    joined["notes"] = "Program-level region/treatment/LODO provenance carried forward unchanged; no cross-species score computed"
    provenance_cols = [
        "lineage", "program", "Dog_CIE_direction", "Human_TI_direction", "Human_AC_direction",
        "region_risk", "Phase10_priority", "Phase8_status", "dog_phase5_stability",
        "FRE_summary_status", "human_treatment_sensitivity_status", "harmonized_level1",
        "phase9_ontology_ready", "notes",
    ]
    provenance = joined[provenance_cols].copy()
    forward = pd.DataFrame(
        {
            "lineage": joined["lineage"],
            "program": joined["program"],
            "dog_direction": joined["Dog_CIE_direction"],
            "dog_lodo_stability": joined["dog_phase5_stability"],
            "fre_summary_status": joined["FRE_summary_status"],
            "human_TI_direction": joined["Human_TI_direction"],
            "human_treatment_sensitivity": joined["human_treatment_sensitivity_status"],
            "human_region_status": joined["TI_AC_consistency"],
            "region_risk": joined["region_risk"],
            "phase10_priority": joined["Phase10_priority"],
            "ontology_level": "LEVEL1",
            "ontology_ready": "YES",
            "phase10_forward_status": np.where(
                joined["region_risk"].astype(str) == "LOW_REGION_RISK",
                "READY_PRIMARY_LEVEL1",
                "READY_PRIMARY_LEVEL1_WITH_REGION_CAUTION",
            ),
            "notes": "Provenance handoff only; Phase9 did not compute conservation or gene-level cross-species similarity",
        }
    )
    return provenance, forward


def build_term_snapshot(terms, frames: list[pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    used: set[str] = set()
    for frame in frames:
        for col in frame.columns:
            # Every formal CL identifier is audited, including mapping parents
            # and shared ancestors rather than only the primary mapped term.
            if col.lower().endswith("_id") or col.lower() == "cl_id":
                for value in frame[col].dropna().astype(str):
                    used.update(x for x in value.split(";") if x.startswith("CL:"))
    validation_rows, snapshot_rows = [], []
    for cl_id in sorted(used):
        exists = cl_id in terms
        obsolete = bool(terms.get(cl_id, {}).get("obsolete", True))
        parents = [str(x) for x in terms.get(cl_id, {}).get("parents", []) if str(x).startswith("CL:")]
        validation_rows.append(
            {
                "cl_id": cl_id,
                "cl_label": label_of(cl_id, terms) if exists else "",
                "exists_in_frozen_cl": "YES" if exists else "NO",
                "obsolete": "YES" if obsolete else "NO",
                "parent_relationships_parseable": "YES" if exists and parents else ("YES_ROOT_OR_NO_CL_PARENT" if exists else "NO"),
                "validation_status": "PASS" if exists and not obsolete else "FAIL",
            }
        )
        snapshot_rows.append(
            {
                "cl_id": cl_id,
                "cl_label": label_of(cl_id, terms),
                "parent_id": ";".join(parents),
                "parent_label": ";".join(label_of(x, terms) for x in parents),
                "obsolete": str(obsolete).lower(),
                "source_version": "v2026-06-08",
            }
        )
    validation = pd.DataFrame(validation_rows)
    if not (validation["validation_status"] == "PASS").all():
        raise RuntimeError("At least one used CL ID failed frozen ontology validation")
    return pd.DataFrame(snapshot_rows), validation


def main() -> int:
    for sub in ["dog", "human", "harmonized", "audit", "summary", "reference"]:
        (P9 / sub).mkdir(parents=True, exist_ok=True)
    terms = parse_obo(OBO)
    if len(terms) < 2_000:
        raise RuntimeError(f"Frozen CL parse unexpectedly small: {len(terms)} terms")

    dog_map, dog_audit = make_dog_mapping(terms)
    human_map, human_audit = make_human_mapping(terms)
    dog_map.to_csv(P9 / "dog" / "DOG_CELLTYPE_ONTOLOGY_MAPPING.csv", index=False)
    human_map.to_csv(P9 / "human" / "HUMAN_CELLTYPE_ONTOLOGY_MAPPING.csv", index=False)

    dog_a = ad.read_h5ad(DOG_H5AD, backed="r")
    dog_obs = dog_a.obs.copy()
    dog_a.file.close()
    human_a = ad.read_h5ad(HUMAN_H5AD, backed="r")
    human_obs = human_a.obs.copy()
    human_a.file.close()
    master, compatibility, coverage = build_master_and_coverage(
        terms, dog_map, human_map, dog_obs, human_obs
    )
    master.to_csv(P9 / "harmonized" / "CROSS_SPECIES_CELL_ONTOLOGY_MASTER.csv", index=False)
    compatibility.to_csv(P9 / "harmonized" / "ONTOLOGY_COMPATIBILITY_MATRIX.csv", index=False)
    coverage.to_csv(P9 / "audit" / "CELLTYPE_CROSS_SPECIES_COVERAGE.csv", index=False)
    pd.concat([dog_audit, human_audit], ignore_index=True).to_csv(
        P9 / "audit" / "CELLTYPE_VS_STATE_AUDIT.csv", index=False
    )

    provenance, program_handoff = build_handoffs()
    provenance.to_csv(P9 / "summary" / "PHASE9_PROGRAM_PROVENANCE_HANDOFF.csv", index=False)
    program_handoff.to_csv(P9 / "summary" / "PHASE10_PROGRAM_HANDOFF.csv", index=False)

    cell_handoff = master.rename(
        columns={
            "harmonized_cl_id": "cl_id",
            "harmonized_cl_label": "cl_label",
            "comparison_mapping_class": "mapping_class",
        }
    ).copy()
    cell_handoff["dog_donor_coverage"] = cell_handoff["dog_donors"].map(lambda x: f"{x}/7 donors")
    cell_handoff["human_donor_coverage"] = cell_handoff["human_primary_donors"].map(lambda x: f"{x}/27 Primary donors")
    cell_handoff["mapping_confidence"] = (
        "Dog=" + cell_handoff["dog_mapping_confidence"].astype(str) +
        ";Human=" + cell_handoff["human_mapping_confidence"].astype(str)
    )
    cell_handoff = cell_handoff[
        [
            "harmonized_level1", "harmonized_level2", "cl_id", "cl_label",
            "dog_native_labels", "human_native_labels", "dog_donor_coverage",
            "human_donor_coverage", "mapping_class", "mapping_confidence",
            "phase10_analysis_level", "phase10_eligible", "reason",
        ]
    ]
    cell_handoff.to_csv(P9 / "summary" / "PHASE10_CELLTYPE_HANDOFF.csv", index=False)

    snapshot, validation = build_term_snapshot(
        terms, [dog_map, human_map, master, compatibility, cell_handoff]
    )
    snapshot.to_csv(P9 / "reference" / "CELL_ONTOLOGY_TERM_SNAPSHOT.csv", index=False)
    validation.to_csv(P9 / "audit" / "CL_ID_VALIDATION.csv", index=False)

    (ROOT / "logs" / "phase9_dog_mapping.log").write_text(
        f"dog_mapping_rows={len(dog_map)}\ndog_native_level2_labels={dog_map['dog_level2'].nunique()}\nstatus=PASS\n",
        encoding="utf-8",
    )
    (ROOT / "logs" / "phase9_human_mapping.log").write_text(
        f"human_mapping_rows={len(human_map)}\nhuman_native_level2_labels={human_map['human_level2'].nunique()}\nstatus=PASS\n",
        encoding="utf-8",
    )
    (ROOT / "logs" / "phase9_harmonization.log").write_text(
        "\n".join(
            [
                f"frozen_cl_terms={len(terms)}",
                f"used_cl_ids={len(validation)}",
                f"compatibility_rows={len(compatibility)}",
                f"mapping_class_counts={json.dumps(compatibility['mapping_class'].value_counts().to_dict())}",
                "primary_level1=Epithelial;T;Myeloid",
                "status=PASS",
            ]
        ) + "\n",
        encoding="utf-8",
    )
    (ROOT / "logs" / "phase9_handoff.log").write_text(
        f"program_handoff_rows={len(program_handoff)}\ncelltype_handoff_rows={len(cell_handoff)}\nregion_risk_missing={program_handoff['region_risk'].isna().sum()}\nphase10_priority_missing={program_handoff['phase10_priority'].isna().sum()}\nstatus=PASS\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "dog_mapping_rows": len(dog_map),
                "dog_native_labels": int(dog_map["dog_level2"].nunique()),
                "human_mapping_rows": len(human_map),
                "human_native_labels": int(human_map["human_level2"].nunique()),
                "master_rows": len(master),
                "compatibility_counts": compatibility["mapping_class"].value_counts().to_dict(),
                "phase10_eligibility": master["phase10_eligible"].value_counts().to_dict(),
                "used_cl_ids": len(validation),
                "program_handoff_rows": len(program_handoff),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
