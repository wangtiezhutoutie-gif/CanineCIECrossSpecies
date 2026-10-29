from __future__ import annotations

import csv
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
BASE = ROOT / "10_cross_species_baseline"
PATHWAY = BASE / "pathway_effects"
SUMMARY = BASE / "summary"
LOG = ROOT / "logs" / "phase10_program_matrix.log"

PRE_BASELINE = PATHWAY / "PATHWAY_CROSS_SPECIES_BASELINE_PRE_PROVENANCE.csv.gz"
PROGRAM_HANDOFF = ROOT / "09_cell_ontology" / "summary" / "PHASE10_PROGRAM_HANDOFF.csv"
REGION_MATRIX = ROOT / "08_region_control" / "summary" / "REGION_INTERPRETATION_MATRIX.csv"
IBDVERSE = (
    ROOT / "07_human_crohn" / "IBDverse" / "summary"
    / "IBDVERSE_HUMAN_PROGRAM_REPLICATION.csv"
)


def log(message: str) -> None:
    stamp = datetime.now(timezone.utc).isoformat()
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8", newline="") as handle:
        handle.write(f"[{stamp}] {message}\n")
    print(message, flush=True)


def atomic_csv(
    frame: pd.DataFrame, path: Path, compression: str | None = None
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        frame.to_csv(
            tmp,
            index=False,
            quoting=csv.QUOTE_MINIMAL,
            compression=compression,
        )
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def baseline_status(row: pd.Series) -> str:
    required = ("dog_NES", "human_NES", "dog_FDR", "human_FDR")
    if any(pd.isna(row.get(name)) for name in required):
        return "NOT_TESTABLE"
    if np.sign(float(row["dog_NES"])) != np.sign(float(row["human_NES"])):
        return "OPPOSITE"
    dog_strong = float(row["dog_FDR"]) < 0.10
    human_strong = float(row["human_FDR"]) < 0.10
    if dog_strong and human_strong:
        return "CONCORDANT_STRONG"
    if human_strong and not dog_strong:
        return "CONCORDANT_DOG_LIMITED"
    if dog_strong and not human_strong:
        return "CONCORDANT_HUMAN_LIMITED"
    return "CONCORDANT_DIRECTION_ONLY"


def external_support(status: str) -> str:
    if status == "SAME_DIRECTION_REPLICATED":
        return "YES"
    if status == "MIXED_CELL_SUBTYPE_REPLICATION":
        return "QUALIFIED_MIXED"
    if status.startswith("NOT_ASSESSABLE"):
        return "NOT_ASSESSABLE"
    if status in {"", "NA_NOT_AVAILABLE", "NOT_AVAILABLE_FOR_REPRESENTATIVE_PATHWAY"}:
        return "NOT_AVAILABLE"
    return "NO"


def main() -> None:
    LOG.unlink(missing_ok=True)
    baseline = pd.read_csv(PRE_BASELINE)
    expected = {
        "lineage", "collection", "pathway_id", "pathway_name",
        "dog_NES", "dog_FDR", "human_NES", "human_FDR",
    }
    if not expected.issubset(baseline.columns) or len(baseline) != 4_406:
        raise RuntimeError("Unbiased paired pathway baseline freeze gate failed")
    log(
        "Unbiased common-one2one pathway baseline frozen before handoff read: "
        f"rows={len(baseline)}"
    )

    handoff = pd.read_csv(PROGRAM_HANDOFF, dtype=str, keep_default_na=False)
    region = pd.read_csv(REGION_MATRIX, dtype=str, keep_default_na=False)
    if len(handoff) != 25 or len(region) != 25:
        raise RuntimeError("Expected 25 frozen program handoff rows")
    if not handoff["ontology_level"].eq("LEVEL1").all():
        raise RuntimeError("Non-Level1 program entered Phase10")
    if set(handoff["lineage"]) != {"Epithelial", "T", "Myeloid"}:
        raise RuntimeError("Primary lineage scope violation")
    registered = handoff.merge(
        region[
            [
                "lineage", "program", "FRE_direction",
                "representative_collection", "representative_pathway",
                "region_risk", "Phase10_priority",
            ]
        ],
        on=["lineage", "program"],
        how="left",
        validate="one_to_one",
        suffixes=("_handoff", "_region"),
    )
    if registered["representative_pathway"].eq("").any():
        raise RuntimeError("Missing representative pathway in frozen region matrix")

    ibd = pd.read_csv(IBDVERSE, dtype=str, keep_default_na=False)
    ibd_keyed = ibd[
        [
            "lineage", "collection", "pathway", "ibdverse_status",
            "ibdverse_mapped", "replicated_same_subtypes",
            "replicated_opposite_subtypes", "ibdverse_supporting_cell_labels",
            "ibdverse_opposite_cell_labels",
        ]
    ].drop_duplicates(["lineage", "collection", "pathway"])
    if ibd_keyed.duplicated(["lineage", "collection", "pathway"]).any():
        raise RuntimeError("IBDverse pathway provenance collision")

    final = baseline.drop(
        columns=[
            "region_risk",
            "phase10_priority",
            "ibdverse_replication_if_available",
            "notes",
        ]
    ).copy()
    program_lookup = registered[
        [
            "lineage", "program", "representative_collection",
            "representative_pathway", "region_risk_handoff",
            "phase10_priority", "human_treatment_sensitivity",
            "dog_lodo_stability", "fre_summary_status", "FRE_direction",
        ]
    ].rename(
        columns={
            "representative_collection": "collection",
            "representative_pathway": "pathway_id",
            "region_risk_handoff": "region_risk",
        }
    )
    final = final.merge(
        program_lookup,
        on=["lineage", "collection", "pathway_id"],
        how="left",
        validate="many_to_one",
    )
    final = final.merge(
        ibd_keyed.rename(columns={"pathway": "pathway_id"}),
        on=["lineage", "collection", "pathway_id"],
        how="left",
        validate="many_to_one",
    )
    final["program"] = final["program"].fillna("NOT_REGISTERED_PROGRAM")
    final["region_risk"] = final["region_risk"].fillna(
        "NOT_APPLICABLE_NOT_REGISTERED_PROGRAM"
    )
    final["phase10_priority"] = final["phase10_priority"].fillna(
        "NOT_APPLICABLE_NOT_REGISTERED_PROGRAM"
    )
    final["ibdverse_status"] = final["ibdverse_status"].fillna(
        "NOT_AVAILABLE_FOR_PATHWAY"
    )
    final["ibdverse_replication_if_available"] = final["ibdverse_status"]
    final["region_caution"] = np.select(
        [
            final["region_risk"].eq("HIGH_REGION_RISK"),
            final["region_risk"].eq("REGION_UNKNOWN"),
            final["region_risk"].eq("MODERATE_DIRECTION_ONLY"),
        ],
        ["YES", "UNKNOWN", "YES"],
        default=np.where(
            final["region_risk"].eq("LOW_REGION_RISK"), "NO", "NOT_APPLICABLE"
        ),
    )
    final["treatment_caution"] = np.where(
        final["human_treatment_sensitivity"]
        .fillna("")
        .str.contains("DIRECTION_CHANGE", regex=False),
        "YES",
        np.where(final["program"].eq("NOT_REGISTERED_PROGRAM"), "NOT_APPLICABLE", "NO"),
    )
    final["notes"] = (
        "Unbiased paired baseline; provenance joined only after all 4,406 "
        "pathway comparisons were frozen. FDR is separate by species."
    )
    order = [
        "lineage", "collection", "pathway_id", "pathway_name",
        "dog_NES", "dog_FDR", "human_NES", "human_FDR",
        "same_direction", "delta_NES", "dog_leading_edge",
        "human_leading_edge", "leading_edge_jaccard", "program",
        "region_risk", "region_caution", "phase10_priority",
        "human_treatment_sensitivity", "treatment_caution",
        "dog_lodo_stability", "fre_summary_status", "FRE_direction",
        "ibdverse_replication_if_available", "notes",
    ]
    remaining = [name for name in final.columns if name not in order]
    final = final[order + remaining]
    atomic_csv(
        final,
        PATHWAY / "PATHWAY_CROSS_SPECIES_BASELINE.csv.gz",
        compression="gzip",
    )

    baseline_for_program = baseline.drop(
        columns=[
            "region_risk",
            "phase10_priority",
            "ibdverse_replication_if_available",
        ]
    )
    program = registered.merge(
        baseline_for_program,
        left_on=["lineage", "representative_collection", "representative_pathway"],
        right_on=["lineage", "collection", "pathway_id"],
        how="left",
        validate="one_to_one",
    )
    program = program.merge(
        ibd_keyed,
        left_on=["lineage", "representative_collection", "representative_pathway"],
        right_on=["lineage", "collection", "pathway"],
        how="left",
        validate="one_to_one",
        suffixes=("", "_ibd"),
    )
    program["ibdverse_status"] = program["ibdverse_status"].fillna(
        "NOT_AVAILABLE_FOR_REPRESENTATIVE_PATHWAY"
    )
    program["human_external_support"] = program["ibdverse_status"].map(
        external_support
    )
    program["baseline_status"] = program.apply(baseline_status, axis=1)
    program["human_direction"] = np.where(
        program["human_NES"].isna(),
        "NOT_TESTABLE",
        np.where(
            program["human_NES"].astype(float) > 0,
            "ActiveCrohn_higher",
            "ActiveCrohn_lower",
        ),
    )
    program["dog_common_direction"] = np.where(
        program["dog_NES"].isna(),
        "NOT_TESTABLE",
        np.where(
            program["dog_NES"].astype(float) > 0, "CIE_higher", "CIE_lower"
        ),
    )
    program["region_caution"] = np.select(
        [
            program["region_risk_handoff"].eq("HIGH_REGION_RISK"),
            program["region_risk_handoff"].eq("REGION_UNKNOWN"),
            program["region_risk_handoff"].eq("MODERATE_DIRECTION_ONLY"),
        ],
        ["YES", "UNKNOWN", "YES"],
        default=np.where(
            program["region_risk_handoff"].eq("LOW_REGION_RISK"), "NO", "UNKNOWN"
        ),
    )
    program["treatment_caution"] = np.where(
        program["human_treatment_sensitivity"].str.contains(
            "DIRECTION_CHANGE", regex=False
        ),
        "YES",
        "NO",
    )
    program["same_direction"] = np.where(
        program["dog_NES"].notna() & program["human_NES"].notna(),
        np.sign(program["dog_NES"].astype(float))
        == np.sign(program["human_NES"].astype(float)),
        pd.NA,
    )
    program["notes"] = (
        "Baseline status is traditional 1:1 ortholog common-universe evidence, "
        "not a final translational fidelity class"
    )
    matrix = pd.DataFrame(
        {
            "lineage": program["lineage"],
            "program": program["program"],
            "representative_collection": program["representative_collection"],
            "representative_pathway": program["representative_pathway"],
            "dog_direction": program["dog_direction"],
            "dog_common_direction": program["dog_common_direction"],
            "dog_lodo_stability": program["dog_lodo_stability"],
            "fre_summary_status": program["fre_summary_status"],
            "FRE_direction": program["FRE_direction"],
            "human_direction": program["human_direction"],
            "human_TI_direction_frozen": program["human_TI_direction"],
            "human_treatment_stability": program[
                "human_treatment_sensitivity"
            ],
            "treatment_caution": program["treatment_caution"],
            "ibdverse_status": program["ibdverse_status"],
            "human_external_support": program["human_external_support"],
            "region_risk": program["region_risk_handoff"],
            "region_caution": program["region_caution"],
            "phase10_priority": program["phase10_priority"],
            "dog_common1to1_NES": program["dog_NES"],
            "human_common1to1_NES": program["human_NES"],
            "same_direction": program["same_direction"],
            "dog_FDR": program["dog_FDR"],
            "human_FDR": program["human_FDR"],
            "leading_edge_jaccard": program["leading_edge_jaccard"],
            "baseline_status": program["baseline_status"],
            "notes": program["notes"],
        }
    )
    atomic_csv(matrix, SUMMARY / "PHASE10_PROGRAM_BASELINE_MATRIX.csv")

    interpretations = []
    for row in matrix.itertuples(index=False):
        parts = [row.baseline_status]
        if row.fre_summary_status == "SUMMARY_DISCORDANT":
            parts.append("FRE_BULK_DISCORDANT")
        elif row.fre_summary_status == "SUMMARY_PARTIAL_SUPPORT":
            parts.append("FRE_BULK_PARTIAL_SUPPORT")
        elif row.fre_summary_status == "SUMMARY_DIRECTIONAL_SUPPORT":
            parts.append("FRE_BULK_DIRECTIONAL_SUPPORT")
        else:
            parts.append(row.fre_summary_status)
        if row.region_caution == "YES":
            parts.append("REGION_CAUTION")
        elif row.region_caution == "UNKNOWN":
            parts.append("REGION_UNKNOWN")
        if row.treatment_caution == "YES":
            parts.append("TREATMENT_CAUTION")
        if row.human_external_support == "YES":
            parts.append("IBDVERSE_SAME_DIRECTION_SUPPORT")
        elif row.human_external_support == "QUALIFIED_MIXED":
            parts.append("IBDVERSE_MIXED_SUPPORT")
        interpretations.append(";".join(parts))
    triangulation = pd.DataFrame(
        {
            "lineage": matrix["lineage"],
            "program": matrix["program"],
            "Dog_CIE_direction": matrix["dog_common_direction"],
            "Dog_LODO": matrix["dog_lodo_stability"],
            "FRE_bulk_direction": matrix["FRE_direction"],
            "FRE_status": matrix["fre_summary_status"],
            "Human_Crohn_direction": matrix["human_direction"],
            "Human_treatment_status": matrix["human_treatment_stability"],
            "IBDverse_status": matrix["ibdverse_status"],
            "region_risk": matrix["region_risk"],
            "interpretation": interpretations,
        }
    )
    atomic_csv(triangulation, SUMMARY / "DOG_FRE_HUMAN_TRIANGULATION.csv")
    counts = matrix["baseline_status"].value_counts(dropna=False).to_dict()
    log(f"Program matrix complete: {counts}")
    log(
        "Region cautions preserved: "
        + str(matrix["region_caution"].value_counts(dropna=False).to_dict())
    )


if __name__ == "__main__":
    main()
