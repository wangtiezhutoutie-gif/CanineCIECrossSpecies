from __future__ import annotations

import csv
import hashlib
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
BASE = ROOT / "07_human_crohn" / "GSE266616"
GSEA = BASE / "gsea"
HUMAN_FROZEN = GSEA / "HUMAN_DISEASE_PROGRAMS.csv"
HUMAN_GATE = GSEA / "metadata" / "PHASE7C_GSEA_GATE.csv"
PHASE5 = ROOT / "06_pseudobulk" / "lodo" / "summary" / "PHASE5_PROGRAM_FORWARD_LIST.csv"
PHASE6 = ROOT / "07_external_validation" / "FRE_2026" / "summary" / "PHASE6_HUMAN_FORWARD_LIST.csv"
OUTPUT = GSEA / "DOG_FORWARD_IN_HUMAN_OBSERVATION.csv"
SUMMARY = GSEA / "DOG_FORWARD_IN_HUMAN_OBSERVATION_SUMMARY.csv"
PROVENANCE = GSEA / "DOG_FORWARD_POSTHOC_PROVENANCE.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    tmp = Path(name)
    try:
        frame.to_csv(tmp, index=False, quoting=csv.QUOTE_MINIMAL)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def assert_human_frozen() -> tuple[str, str]:
    if not HUMAN_FROZEN.is_file() or not HUMAN_GATE.is_file():
        raise RuntimeError("Human GSEA is not frozen; Dog-forward observation is prohibited")
    gate = pd.read_csv(HUMAN_GATE)
    if not gate["status"].eq("PASS").all():
        raise RuntimeError("Human GSEA gate did not pass; Dog-forward observation is prohibited")
    full_files = list((GSEA / "full").glob("Human_*_GSEA.csv"))
    if len(full_files) != 9:
        raise RuntimeError(f"Expected nine frozen Human full GSEA files, got {len(full_files)}")
    return sha256(HUMAN_FROZEN), datetime.fromtimestamp(HUMAN_FROZEN.stat().st_mtime, timezone.utc).isoformat()


def load_human_results() -> pd.DataFrame:
    files = sorted((GSEA / "full").glob("Human_*_GSEA.csv"))
    frames = [pd.read_csv(path) for path in files]
    result = pd.concat(frames, ignore_index=True)
    key = ["lineage", "collection", "pathway_id"]
    if result.duplicated(key).any():
        raise RuntimeError("Duplicated Human lineage/collection/pathway GSEA identity")
    return result


def human_equivalent(collection: str, pathway_id: str) -> tuple[str, str]:
    collection = collection.strip()
    pathway_id = pathway_id.strip()
    if collection == "Reactome":
        pathway_id = pathway_id.replace("R-CFA-", "R-HSA-", 1)
    return collection, pathway_id


def parse_representatives(text: str) -> list[tuple[str, str]]:
    reps: list[tuple[str, str]] = []
    for item in str(text).split(";"):
        if not item.strip() or ":" not in item:
            continue
        collection, pathway_id = item.split(":", 1)
        reps.append(human_equivalent(collection, pathway_id))
    return reps


def classify(rows: pd.DataFrame, expected_sign: int) -> tuple[str, str]:
    if rows.empty:
        return "NOT_ASSESSABLE", "No preregistered representative pathway survived Human GSEA size/mapping filters"
    significant = rows.loc[rows["padj"] < 0.05]
    if not significant.empty:
        same = significant.loc[significant["NES"].map(lambda x: (x > 0) - (x < 0)) == expected_sign]
        opposite = significant.loc[significant["NES"].map(lambda x: (x > 0) - (x < 0)) == -expected_sign]
        if not same.empty and opposite.empty:
            return "SAME_DIRECTION", f"{len(same)} representative pathway(s) FDR<0.05 in the Dog-consistent direction"
        if same.empty and not opposite.empty:
            return "OPPOSITE_DIRECTION", f"{len(opposite)} representative pathway(s) FDR<0.05 opposite to Dog"
        if not same.empty and not opposite.empty:
            return "WEAK", "Mixed significant directions across preregistered representatives; no single-direction program call"
    suggestive = rows.loc[rows["padj"] < 0.10]
    if not suggestive.empty:
        signs = suggestive["NES"].map(lambda x: (x > 0) - (x < 0))
        relation = "same" if (signs == expected_sign).sum() >= (signs == -expected_sign).sum() else "opposite"
        return "WEAK", f"Only suggestive FDR<0.10 evidence; predominant direction={relation}"
    return "NOT_PROMINENT", "Representatives were tested but none reached FDR<0.10 in Human"


def main() -> None:
    # The freeze hash/time is captured before either Dog file is opened.
    frozen_hash, frozen_time = assert_human_frozen()
    human = load_human_results()

    # Dog inputs are intentionally read only after the independent Human freeze checks above.
    dog6 = pd.read_csv(PHASE6)
    dog5 = pd.read_csv(PHASE5)[["lineage", "program", "representative_pathways", "full7_direction"]]
    dog = dog6.merge(dog5, on=["lineage", "program"], how="left", validate="one_to_one")
    if len(dog) != 25 or dog["representative_pathways"].isna().any():
        raise RuntimeError("Dog forward list must resolve exactly 25 preregistered program mappings")

    output_rows: list[dict[str, object]] = []
    for row in dog.itertuples(index=False):
        expected_sign = 1 if row.CIE_direction == "CIE_higher" else -1
        reps = parse_representatives(row.representative_pathways)
        observations: list[pd.DataFrame] = []
        for collection, pathway_id in reps:
            match = human.loc[
                human["lineage"].eq(row.lineage)
                & human["collection"].eq(collection)
                & human["pathway_id"].eq(pathway_id)
            ].copy()
            observations.append(match)
        observed = pd.concat(observations, ignore_index=True) if observations else pd.DataFrame()
        status, interpretation = classify(observed, expected_sign)
        evidence = ""
        if not observed.empty:
            evidence = ";".join(
                f"{x.collection}:{x.pathway_id}|NES={x.NES:.4g}|FDR={x.padj:.4g}|{x.direction}"
                for x in observed.sort_values(["padj", "collection", "pathway_id"]).itertuples(index=False)
            )
        output_rows.append({
            "lineage": row.lineage,
            "program": row.program,
            "dog_CIE_direction": row.CIE_direction,
            "dog_phase5_stability": row.phase5_stability,
            "FRE_summary_status": row.FRE_summary_status,
            "FRE_direction": row.FRE_direction,
            "phase6_priority": row.priority,
            "human_observation": status,
            "human_representatives_assessable": int(len(observed)),
            "human_representative_evidence": evidence,
            "interpretation": interpretation,
            "representative_mapping_rule": "GO/Hallmark IDs unchanged; Reactome R-CFA numeric identifier converted to R-HSA",
            "scope_limit": "Post hoc program-level observation only; not formal cross-species fidelity or gene-effect correlation",
        })
    output = pd.DataFrame(output_rows)
    if len(output) != 25:
        raise RuntimeError("Post hoc observation must contain exactly 25 Dog programs")
    atomic_csv(output, OUTPUT)
    order = ["SAME_DIRECTION", "OPPOSITE_DIRECTION", "WEAK", "NOT_PROMINENT", "NOT_ASSESSABLE"]
    summary = output["human_observation"].value_counts().reindex(order, fill_value=0).rename_axis("human_observation").reset_index(name="program_n")
    summary["total_programs"] = len(output)
    atomic_csv(summary, SUMMARY)
    atomic_csv(pd.DataFrame([{
        "human_freeze_file": str(HUMAN_FROZEN),
        "human_freeze_sha256_before_dog_read": frozen_hash,
        "human_freeze_mtime_utc": frozen_time,
        "phase6_dog_forward_file": str(PHASE6),
        "phase6_sha256": sha256(PHASE6),
        "phase5_mapping_file": str(PHASE5),
        "phase5_sha256": sha256(PHASE5),
        "observation_generated_utc": datetime.now(timezone.utc).isoformat(),
        "sequence_assertion": "HUMAN_FROZEN_BEFORE_DOG_READ",
    }]), PROVENANCE)


if __name__ == "__main__":
    main()
