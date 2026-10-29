#!/usr/bin/env python3
"""Phase5 full seven-dog LODO stability summaries and donor-level figures.

Full7 remains the primary result.  Every comparison uses only features/pathways
tested in the corresponding Full7/LODO common universe.  The classifications
below are the pre-registered project-internal descriptive rules recorded in
PHASE5_STABILITY_PARAMETERS.csv; they are not general inferential thresholds.
"""

from __future__ import annotations

import gzip
import hashlib
import platform
import re
import sys
import warnings
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
PB = ROOT / "06_pseudobulk"
LODO = PB / "lodo"
SUMMARY = LODO / "summary"
FIGURES = LODO / "figures" / "stability"
META = LODO / "metadata"
LOG = ROOT / "logs" / "phase5_stability_summary.log"

LINEAGES = ["Epithelial", "T", "Myeloid"]
RUN_ORDER = ["drop_h1", "drop_h2", "drop_h3", "drop_cie2", "drop_cie3", "drop_cie4", "drop_cie6"]
DISPLAY_RUNS = ["Full7", *RUN_ORDER]
DONOR_FROM_SHORT = {
    "drop_h1": "GSM8031698_h_1",
    "drop_h2": "GSM8031699_h_2",
    "drop_h3": "GSM8031700_h_3",
    "drop_cie2": "GSM8031701_cie_2",
    "drop_cie3": "GSM8031702_cie_3",
    "drop_cie4": "GSM8031703_cie_4",
    "drop_cie6": "GSM8031704_cie_6",
}
SHORT_FROM_DONOR = {v: k for k, v in DONOR_FROM_SHORT.items()}
ANCHORS = {
    "Epithelial": ["KCNK16", "MMP7", "GAL3ST2", "CYP1B1"],
    "T": ["GAL3ST2", "TLR9", "TMEM255A", "PTPRT"],
    "Myeloid": ["SERPINE1", "CCL3", "ADAMDEC1", "MMP9", "SAA1"],
}


# Each requested Phase4 program is pre-specified here with an expected direction.
# Regexes only choose representative pathways from the frozen Phase4 results.
PROGRAM_DEFS = {
    "Epithelial": [
        ("Cell cycle / proliferation", "CIE_higher", r"cell cycle|g2m|e2f|mitotic|dna replication|chromosome segregation"),
        ("Ribosome / RNA biogenesis", "CIE_higher", r"ribosom|rrna|ribonucleoprotein|metabolism of rna|rna processing"),
        ("Stress / UPR / MAPK", "CIE_higher", r"unfolded protein|endoplasmic reticulum stress|heat stress|mapk|kras|response to stress"),
        ("TNF / NF-kappaB", "CIE_higher", r"tnfa.signal|tnf signal|tnfr|nf.?kappa|nfkb"),
        ("Mitochondrial respiration", "CIE_lower", r"aerobic respiration|respiratory electron|electron transport chain|mitochondrial respiration"),
        ("Oxidative phosphorylation", "CIE_lower", r"oxidative phosphorylation"),
        ("Fatty acid metabolism", "CIE_lower", r"fatty acid"),
        ("Xenobiotic metabolism", "CIE_lower", r"xenobiotic|drug catabolic|detoxification"),
    ],
    "T": [
        ("TNF / NF-kappaB", "CIE_higher", r"tnfa.signal|tnf signal|tnfr|nf.?kappa|nfkb"),
        ("Inflammatory response", "CIE_higher", r"inflammatory response|inflammation"),
        ("Cytokine / interleukin", "CIE_higher", r"cytokine|interleukin"),
        ("Stress / UPR / MAPK", "CIE_higher", r"unfolded protein|endoplasmic reticulum stress|heat stress|mapk|kras|response to stress"),
        ("Oxidative phosphorylation", "CIE_lower", r"oxidative phosphorylation"),
        ("Mitochondrial respiration", "CIE_lower", r"aerobic respiration|respiratory electron|electron transport chain|mitochondrial respiration"),
        ("Mitochondrial translation", "CIE_lower", r"mitochondrial translat|mitochondrial ribosom"),
        ("Metabolic programs", "CIE_lower", r"fatty acid|xenobiotic|oxoacid|small molecule metabolic|metabolic process"),
    ],
    "Myeloid": [
        ("TNF / NF-kappaB", "CIE_higher", r"tnfa.signal|tnf signal|tnfr|nf.?kappa|nfkb"),
        ("Inflammatory response", "CIE_higher", r"inflammatory response|inflammation"),
        ("TLR4", "CIE_higher", r"toll like receptor 4|toll-like receptor 4|tlr4"),
        ("Cytokine / interleukin / JAK-STAT", "CIE_higher", r"cytokine|interleukin|jak.?stat|stat[0-9]"),
        ("Stress response", "CIE_higher", r"unfolded protein|endoplasmic reticulum stress|heat stress|mapk|kras|response to stress"),
        ("Apoptosis", "CIE_higher", r"apoptosis|apoptotic"),
        ("Respiratory electron transport", "CIE_lower", r"aerobic respiration|respiratory electron|electron transport chain|mitochondrial respiration"),
        ("Oxidative metabolism", "CIE_lower", r"oxidative phosphorylation|oxidative metabolic|aerobic respiration"),
        ("Translation", "CIE_lower", r"cytoplasmic translation|translation$|ribosom|peptide biosynthetic"),
    ],
}


def write_csv_atomic(frame: pd.DataFrame, path: Path, *, gzip_output: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = Path(str(path) + ".tmp")
    if gzip_output:
        with gzip.open(temp, "wt", encoding="utf-8", newline="") as handle:
            frame.to_csv(handle, index=False)
    else:
        frame.to_csv(temp, index=False)
    temp.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_corr(a: pd.Series, b: pd.Series, method: str) -> float:
    valid = a.notna() & b.notna()
    if valid.sum() < 3 or a[valid].nunique() < 2 or b[valid].nunique() < 2:
        return np.nan
    return a[valid].corr(b[valid], method=method)


def jaccard(a: set[str], b: set[str]) -> float:
    union = a | b
    return len(a & b) / len(union) if union else np.nan


def semicolon(values) -> str:
    return ";".join(str(x) for x in values if pd.notna(x) and str(x))


def gene_class(tested: int, same: int, median_delta: float, max_delta: float) -> str:
    if tested < 4:
        return "NOT_ASSESSABLE"
    if same <= 3:
        return "UNSTABLE"
    if tested == 7 and same == 7 and median_delta <= 0.25 and max_delta <= 0.75:
        return "ROBUST"
    if tested >= 6 and same >= 6 and median_delta <= 0.50 and max_delta <= 1.50:
        return "MODERATELY_STABLE"
    return "DONOR_SENSITIVE"


def pathway_class(tested: int, same: int, median_delta: float, max_delta: float,
                  median_jaccard: float, strong_flip_n: int) -> str:
    if tested < 4:
        return "NOT_ASSESSABLE"
    if same <= 3 or strong_flip_n >= 2:
        return "UNSTABLE"
    if (tested == 7 and same == 7 and median_delta <= 0.30 and max_delta <= 0.75
            and pd.notna(median_jaccard) and median_jaccard >= 0.50):
        return "ROBUST"
    if (tested >= 6 and same >= 6 and median_delta <= 0.50 and strong_flip_n == 0
            and pd.notna(median_jaccard) and median_jaccard >= 0.25):
        return "MODERATELY_STABLE"
    return "DONOR_SENSITIVE"


def load_full_de(lineage: str) -> pd.DataFrame:
    path = PB / "de_results" / f"{lineage}_CIE_vs_Healthy_edger.csv"
    frame = pd.read_csv(path)
    if frame["feature_id"].duplicated().any():
        raise RuntimeError(f"Duplicated Full7 feature IDs: {lineage}")
    return frame


def build_gene_outputs(log) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    status = pd.read_csv(META / "LODO_RUN_STATUS.csv")
    if len(status) != 21 or set(status["status"]) != {"SUCCESS"}:
        raise RuntimeError("Gene LODO is not 21/21 SUCCESS")
    master_rows: list[dict] = []
    effects_long: list[pd.DataFrame] = []
    rank_rows: list[dict] = []

    for lineage in LINEAGES:
        full = load_full_de(lineage).copy()
        full_rank = full.sort_values(["rank_metric", "feature_id"], ascending=[False, True]).reset_index(drop=True)
        full_rank["signed_rank"] = np.arange(1, len(full_rank) + 1)
        full_abs = full.assign(abs_rank_metric=full["rank_metric"].abs()).sort_values(
            ["abs_rank_metric", "feature_id"], ascending=[False, True]
        ).reset_index(drop=True)
        full_abs["abs_rank"] = np.arange(1, len(full_abs) + 1)
        full_ranks = full_rank[["feature_id", "signed_rank"]].merge(
            full_abs[["feature_id", "abs_rank"]], on="feature_id"
        )
        run_frames: dict[str, pd.DataFrame] = {}
        lineage_long_parts: list[pd.DataFrame] = []
        for _, run in status.loc[status["lineage"].eq(lineage)].iterrows():
            lodo = pd.read_csv(Path(run["ranked_file"]))
            if set(lodo["lineage"].dropna()) != {lineage}:
                raise RuntimeError(f"LODO ranked lineage mismatch: {run['ranked_file']}")
            lodo = lodo.copy()
            lodo_signed = lodo.sort_values(["rank_metric", "feature_id"], ascending=[False, True]).reset_index(drop=True)
            lodo_signed["lodo_signed_rank"] = np.arange(1, len(lodo_signed) + 1)
            lodo_abs = lodo.assign(abs_rank_metric=lodo["rank_metric"].abs()).sort_values(
                ["abs_rank_metric", "feature_id"], ascending=[False, True]
            ).reset_index(drop=True)
            lodo_abs["lodo_abs_rank"] = np.arange(1, len(lodo_abs) + 1)
            lodo = lodo.merge(lodo_signed[["feature_id", "lodo_signed_rank"]], on="feature_id")
            lodo = lodo.merge(lodo_abs[["feature_id", "lodo_abs_rank"]], on="feature_id")
            run_frames[run["short_label"]] = lodo

            common = full[["feature_id", "rank_metric"]].merge(
                lodo[["feature_id", "rank_metric"]], on="feature_id", suffixes=("_full", "_lodo")
            )
            top_counts = {}
            for n in (100, 200):
                fset = set(common.nlargest(n, "rank_metric_full")["feature_id"])
                lset = set(common.nlargest(n, "rank_metric_lodo")["feature_id"])
                faset = set(common.assign(v=common["rank_metric_full"].abs()).nlargest(n, "v")["feature_id"])
                laset = set(common.assign(v=common["rank_metric_lodo"].abs()).nlargest(n, "v")["feature_id"])
                top_counts[f"top{n}_overlap_n"] = len(fset & lset)
                top_counts[f"top{n}_jaccard"] = jaccard(fset, lset)
                top_counts[f"top{n}_abs_overlap_n"] = len(faset & laset)
                top_counts[f"top{n}_abs_jaccard"] = jaccard(faset, laset)
            rank_rows.append({
                "lineage": lineage,
                "lodo_run": run["lodo_run"],
                "short_label": run["short_label"],
                "dropped_donor": run["dropped_donor"],
                "common_tested_gene_n": len(common),
                "spearman_rho": safe_corr(common["rank_metric_full"], common["rank_metric_lodo"], "spearman"),
                "pearson_r": safe_corr(common["rank_metric_full"], common["rank_metric_lodo"], "pearson"),
                **top_counts,
                "top_definition": "signed rank_metric; abs columns use abs(rank_metric)",
            })

            long = full[["feature_id", "gene_symbol", "logFC", "rank_metric"]].merge(
                lodo[["feature_id", "logFC", "PValue", "FDR", "rank_metric", "lodo_signed_rank", "lodo_abs_rank"]],
                on="feature_id", suffixes=("_full", "_lodo")
            ).merge(full_ranks, on="feature_id")
            long["lineage"] = lineage
            long["lodo_run"] = run["lodo_run"]
            long["short_label"] = run["short_label"]
            long["dropped_donor"] = run["dropped_donor"]
            long["delta_logFC"] = long["logFC_lodo"] - long["logFC_full"]
            long["abs_delta_logFC"] = long["delta_logFC"].abs()
            long["same_direction"] = np.sign(long["logFC_lodo"]) == np.sign(long["logFC_full"])
            long["meaningful_same_direction"] = (
                long["same_direction"] & long["logFC_lodo"].abs().ge(0.25) & long["logFC_full"].abs().ge(0.25)
            )
            long["sign_flip"] = ~long["same_direction"]
            long["meaningful_sign_flip"] = (
                long["sign_flip"] & long["logFC_lodo"].abs().ge(0.25) & long["logFC_full"].abs().ge(0.25)
            )
            effects_long.append(long)
            lineage_long_parts.append(long)

        # Vectorized aggregation avoids repeatedly scanning seven ~10k-row tables
        # for each feature and preserves NOT_TESTED as an absent observation.
        lineage_effects = pd.concat(lineage_long_parts, ignore_index=True)
        aggregate = lineage_effects.groupby("feature_id", sort=False).agg(
            tested_lodo_n=("logFC_lodo", "count"),
            same_direction_n=("same_direction", "sum"),
            meaningful_same_direction_n=("meaningful_same_direction", "sum"),
            sign_flip_n=("sign_flip", "sum"),
            meaningful_sign_flip_n=("meaningful_sign_flip", "sum"),
            median_lodo_logFC=("logFC_lodo", "median"),
            min_lodo_logFC=("logFC_lodo", "min"),
            max_lodo_logFC=("logFC_lodo", "max"),
            median_abs_delta_logFC=("abs_delta_logFC", "median"),
            max_abs_delta_logFC=("abs_delta_logFC", "max"),
            fdr05_retention_n=("FDR", lambda x: int((x < 0.05).sum())),
            fdr10_retention_n=("FDR", lambda x: int((x < 0.10).sum())),
        ).reset_index()
        influence_indices = lineage_effects.groupby("feature_id", sort=False)["abs_delta_logFC"].idxmax()
        influential = lineage_effects.loc[influence_indices, [
            "feature_id", "dropped_donor", "short_label", "delta_logFC"
        ]].rename(columns={
            "dropped_donor": "most_influential_donor",
            "short_label": "most_influential_short_label",
            "delta_logFC": "most_influential_delta",
        })
        lineage_master = full[["feature_id", "gene_symbol", "logFC", "PValue", "FDR", "rank_metric"]].rename(columns={
            "logFC": "full_logFC", "PValue": "full_PValue", "FDR": "full_FDR", "rank_metric": "full_rank_metric"
        }).merge(aggregate, on="feature_id", how="left").merge(influential, on="feature_id", how="left")
        count_columns = ["tested_lodo_n", "same_direction_n", "meaningful_same_direction_n", "sign_flip_n",
                         "meaningful_sign_flip_n", "fdr05_retention_n", "fdr10_retention_n"]
        lineage_master[count_columns] = lineage_master[count_columns].fillna(0).astype(int)
        lineage_master.insert(0, "lineage", lineage)
        lineage_master["stability_class"] = lineage_master.apply(
            lambda r: gene_class(int(r["tested_lodo_n"]), int(r["same_direction_n"]),
                                 float(r["median_abs_delta_logFC"]), float(r["max_abs_delta_logFC"]))
            if int(r["tested_lodo_n"]) else "NOT_ASSESSABLE", axis=1
        )
        lineage_master["notes"] = "Full7 is primary; LODO FDR retention is auxiliary; untested rounds are not imputed."
        master_rows.extend(lineage_master.to_dict("records"))
        print(f"Gene summary {lineage}: {len(full):,} Full7-tested features", file=log)

    master = pd.DataFrame(master_rows)
    effects = pd.concat(effects_long, ignore_index=True)
    ranks = pd.DataFrame(rank_rows)
    write_csv_atomic(master, SUMMARY / "GENE_LODO_STABILITY.csv.gz", gzip_output=True)
    write_csv_atomic(effects, SUMMARY / "GENE_LODO_EFFECTS_LONG.csv.gz", gzip_output=True)
    write_csv_atomic(ranks, SUMMARY / "GENE_RANK_STABILITY_BY_RUN.csv")

    anchor_long = effects.loc[
        effects.apply(lambda r: r["gene_symbol"] in ANCHORS.get(r["lineage"], []), axis=1)
    ].copy()
    anchor_metrics = master.loc[
        master.apply(lambda r: r["gene_symbol"] in ANCHORS.get(r["lineage"], []), axis=1)
    ].copy()
    anchor_out = anchor_long.merge(
        anchor_metrics[["lineage", "feature_id", "tested_lodo_n", "same_direction_n", "median_lodo_logFC",
                        "min_lodo_logFC", "max_lodo_logFC", "fdr05_retention_n", "most_influential_donor",
                        "most_influential_delta", "stability_class"]],
        on=["lineage", "feature_id"], how="left"
    )
    anchor_out["anchor_pre_registered"] = True
    write_csv_atomic(anchor_out, SUMMARY / "ANCHOR_GENE_LODO_SUMMARY.csv")
    return master, effects, ranks


def load_official_gsea() -> pd.DataFrame:
    specs = {
        "Hallmark": ("hallmark", "Hallmark"),
        "Reactome": ("reactome", "Reactome"),
        "GO_BP": ("go_bp", "GOBP"),
    }
    frames = []
    for lineage in LINEAGES:
        for collection, (folder, tag) in specs.items():
            path = PB / "gsea" / "full7" / folder / f"{lineage}_{tag}_GSEA.csv"
            f = pd.read_csv(path)
            f["lineage"] = lineage
            f["collection"] = collection
            frames.append(f)
    return pd.concat(frames, ignore_index=True)


def build_pathway_outputs(log) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    status = pd.read_csv(META / "LODO_GSEA_RUN_STATUS.csv")
    if len(status) != 63 or set(status["status"]) != {"SUCCESS"}:
        raise RuntimeError("LODO GSEA is not 63/63 SUCCESS")
    pair_frames = []
    for path_text in status["comparison_file"]:
        f = pd.read_csv(Path(path_text))
        pair_frames.append(f)
    pairs = pd.concat(pair_frames, ignore_index=True)
    official = load_official_gsea()
    expected_pairs = {(l, c, s) for l in LINEAGES for c in ["Hallmark", "Reactome", "GO_BP"] for s in RUN_ORDER}
    seen_pairs = set(zip(pairs["lineage"], pairs["collection"], pairs["short_label"]))
    if seen_pairs != expected_pairs:
        raise RuntimeError(f"Paired GSEA coverage mismatch: expected {len(expected_pairs)}, observed {len(seen_pairs)}")

    pair_keys = ["lineage", "collection", "pathway_id"]
    grouped = {key: group for key, group in pairs.groupby(pair_keys, sort=False)}
    rows = []
    for _, off in official.iterrows():
        key = (off["lineage"], off["collection"], off["pathway_id"])
        group = grouped.get(key, pd.DataFrame())
        comparable = group.loc[group.get("lodo_NES", pd.Series(dtype=float)).notna() & group.get("pair_full_NES", pd.Series(dtype=float)).notna()].copy()
        tested = len(comparable)
        if tested:
            official_same = np.sign(comparable["lodo_NES"]) == np.sign(float(off["NES"]))
            pair_same = comparable["same_direction"].astype(str).str.lower().eq("true") if comparable["same_direction"].dtype == object else comparable["same_direction"].astype(bool)
            sign_flip = comparable["sign_flip"].astype(str).str.lower().eq("true") if comparable["sign_flip"].dtype == object else comparable["sign_flip"].astype(bool)
            strong = comparable["strong_effect_flip"].astype(str).str.lower().eq("true") if comparable["strong_effect_flip"].dtype == object else comparable["strong_effect_flip"].astype(bool)
            deltas = comparable["delta_NES"].astype(float)
            influence_index = deltas.abs().idxmax()
            influential = comparable.loc[influence_index]
            le = comparable["leading_edge_jaccard"].astype(float)
            med_delta = float(deltas.abs().median())
            max_delta = float(deltas.abs().max())
            med_le = float(le.median()) if le.notna().any() else np.nan
            min_le = float(le.min()) if le.notna().any() else np.nan
            max_le = float(le.max()) if le.notna().any() else np.nan
            same_n = int(pair_same.sum())
            strong_n = int(strong.sum())
            classification = pathway_class(tested, same_n, med_delta, max_delta, med_le, strong_n)
            lodo_values = comparable["lodo_NES"].astype(float)
        else:
            official_same = pd.Series(dtype=bool)
            pair_same = sign_flip = strong = pd.Series(dtype=bool)
            influential = None
            med_delta = max_delta = med_le = min_le = max_le = np.nan
            same_n = strong_n = 0
            classification = "NOT_ASSESSABLE"
            lodo_values = pd.Series(dtype=float)
        rows.append({
            "lineage": off["lineage"],
            "collection": off["collection"],
            "pathway_id": off["pathway_id"],
            "pathway_name": off["pathway_name"],
            "full_NES": off["NES"],
            "full_padj": off["padj"],
            "full_direction": off["direction"],
            "tested_lodo_n": tested,
            "same_direction_n": same_n,
            "same_official_full7_direction_n": int(official_same.sum()),
            "sign_flip_n": int(sign_flip.sum()) if tested else 0,
            "strong_effect_flip_n": strong_n,
            "median_lodo_NES": lodo_values.median() if tested else np.nan,
            "min_lodo_NES": lodo_values.min() if tested else np.nan,
            "max_lodo_NES": lodo_values.max() if tested else np.nan,
            "median_abs_delta_NES": med_delta,
            "max_abs_delta_NES": max_delta,
            "fdr05_retention_n": int((comparable["lodo_padj"] < 0.05).sum()) if tested else 0,
            "fdr10_retention_n": int((comparable["lodo_padj"] < 0.10).sum()) if tested else 0,
            "median_leading_edge_jaccard": med_le,
            "min_leading_edge_jaccard": min_le,
            "max_leading_edge_jaccard": max_le,
            "most_influential_donor": influential["dropped_donor"] if influential is not None else "",
            "most_influential_short_label": influential["short_label"] if influential is not None else "",
            "most_influential_delta_NES": influential["delta_NES"] if influential is not None else np.nan,
            "stability_class": classification,
            "notes": "Official Phase4 Full7 NES is primary; deltas/classes use paired common-universe Full7 versus LODO; FDR retention is auxiliary.",
        })
    master = pd.DataFrame(rows)
    write_csv_atomic(master, SUMMARY / "PATHWAY_LODO_STABILITY.csv")
    qc = pd.read_csv(META / "LODO_GSEA_PAIR_QC.csv")
    qc_out = qc[["lineage", "collection", "lodo_run", "short_label", "dropped_donor",
                 "pathways_comparable", "top20_jaccard", "nes_spearman", "nes_pearson",
                 "direction_concordance", "max_abs_pathway_size_change", "status"]].copy()
    write_csv_atomic(qc_out, SUMMARY / "PATHWAY_RANK_STABILITY_BY_RUN.csv")
    print(f"Pathway master: {len(master):,} official Phase4 pathways; paired rows: {len(pairs):,}", file=log)

    representatives = select_core_representatives(master)
    programs = aggregate_core_programs(representatives)
    write_csv_atomic(programs, SUMMARY / "CORE_PROGRAM_LODO_SUMMARY.csv")
    write_csv_atomic(representatives, SUMMARY / "CORE_PROGRAM_REPRESENTATIVE_PATHWAYS.csv")

    forward = programs[["lineage", "program", "representative_pathways", "full7_direction", "full7_evidence",
                        "direction_stability", "effect_stability", "leading_edge_stability",
                        "most_influential_donor", "stability_class"]].copy()
    forward["forward_to_phase6"] = np.where(
        forward["stability_class"].eq("ROBUST"), "YES",
        np.where(forward["stability_class"].eq("MODERATELY_STABLE"), "YES_WITH_QUALIFICATION", "NO")
    )
    forward["forward_to_human"] = np.where(
        forward["stability_class"].eq("ROBUST"), "YES",
        np.where(forward["stability_class"].eq("MODERATELY_STABLE"), "EXPLORATORY_ONLY", "NO")
    )
    forward["notes"] = "Program-level descriptive LODO status; external canine validation remains required."
    write_csv_atomic(forward, SUMMARY / "PHASE5_PROGRAM_FORWARD_LIST.csv")
    return master, pairs, representatives, programs


def select_core_representatives(master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for lineage, definitions in PROGRAM_DEFS.items():
        lineage_frame = master.loc[master["lineage"].eq(lineage)]
        for program, expected_direction, pattern in definitions:
            matches = lineage_frame.loc[
                lineage_frame["pathway_name"].astype(str).str.contains(pattern, case=False, regex=True, na=False)
                & lineage_frame["full_direction"].eq(expected_direction)
            ].copy()
            if matches.empty:
                rows.append({
                    "lineage": lineage, "program": program, "expected_direction": expected_direction,
                    "selection_pattern": pattern, "collection": "", "pathway_id": "", "pathway_name": "",
                    "full_NES": np.nan, "full_padj": np.nan, "tested_lodo_n": 0, "same_direction_n": 0,
                    "strong_effect_flip_n": 0, "median_lodo_NES": np.nan, "median_abs_delta_NES": np.nan,
                    "max_abs_delta_NES": np.nan, "median_leading_edge_jaccard": np.nan,
                    "most_influential_donor": "", "most_influential_delta_NES": np.nan,
                    "stability_class": "NOT_ASSESSABLE", "selection_note": "No direction-matched frozen Phase4 pathway matched regex."
                })
                continue
            # One strongest Phase4 representative per collection; no LODO result is used for selection.
            matches = matches.sort_values(["full_padj", "full_NES"], ascending=[True, expected_direction == "CIE_lower"])
            chosen = matches.groupby("collection", sort=False).head(1).copy()
            chosen.insert(1, "program", program)
            chosen.insert(2, "expected_direction", expected_direction)
            chosen.insert(3, "selection_pattern", pattern)
            chosen["selection_note"] = "Best official Phase4 FDR-matched pathway per frozen collection; selected without LODO optimization."
            rows.extend(chosen.to_dict("records"))
    return pd.DataFrame(rows)


def aggregate_core_programs(representatives: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (lineage, program), group in representatives.groupby(["lineage", "program"], sort=False):
        assessed = group.loc[group["tested_lodo_n"].ge(4) & group["pathway_id"].astype(str).ne("")].copy()
        expected = group["expected_direction"].iloc[0]
        if assessed.empty:
            pclass = "NOT_ASSESSABLE"
            direction_stability = effect_stability = le_stability = "NOT_ASSESSABLE"
            most_donor = ""
            evidence_reps = assessed
            min_tested = min_same = np.nan
            median_delta = max_delta = median_le = np.nan
        else:
            # Program-level synthesis is deliberately conservative and does not
            # overwrite any pathway class: use only official Full7 FDR<0.10
            # representatives when available, and downgrade mixed evidence.
            evidence_reps = assessed.loc[assessed["full_padj"].lt(0.10)].copy()
            if evidence_reps.empty:
                evidence_reps = assessed.sort_values("full_padj").head(1).copy()
            min_tested = int(evidence_reps["tested_lodo_n"].min())
            min_same = int(evidence_reps["same_direction_n"].min())
            median_delta = float(evidence_reps["median_abs_delta_NES"].median())
            max_delta = float(evidence_reps["max_abs_delta_NES"].max())
            median_le = float(evidence_reps["median_leading_edge_jaccard"].median())
            evidence_classes = set(evidence_reps["stability_class"])
            if "UNSTABLE" in evidence_classes:
                pclass = "UNSTABLE"
            elif evidence_classes == {"ROBUST"}:
                pclass = "ROBUST"
            elif "ROBUST" in evidence_classes and evidence_classes <= {"ROBUST", "MODERATELY_STABLE", "DONOR_SENSITIVE"}:
                pclass = "MODERATELY_STABLE"
            elif evidence_classes <= {"ROBUST", "MODERATELY_STABLE"}:
                pclass = "MODERATELY_STABLE"
            elif "DONOR_SENSITIVE" in evidence_classes:
                pclass = "DONOR_SENSITIVE"
            else:
                pclass = "NOT_ASSESSABLE"
            direction_stability = f"minimum {min_same}/{min_tested} across {len(evidence_reps)} Full7-supported representatives"
            effect_stability = f"median representative median|deltaNES|={median_delta:.3f}; max|deltaNES|={max_delta:.3f}"
            le_stability = f"median representative LE Jaccard={median_le:.3f}"
            donor_score = evidence_reps.groupby("most_influential_donor")["max_abs_delta_NES"].sum().sort_values(ascending=False)
            most_donor = donor_score.index[0] if len(donor_score) else ""
        evidence = semicolon(
            f"{r.collection}:{r.pathway_name} (NES={r.full_NES:.3f}, FDR={r.full_padj:.3g}, {r.stability_class})"
            for r in evidence_reps.itertuples()
        )
        rows.append({
            "lineage": lineage,
            "program": program,
            "representative_pathways": semicolon(f"{r.collection}:{r.pathway_id}" for r in assessed.itertuples()),
            "representative_pathway_n": len(assessed),
            "full7_supported_representative_n": len(evidence_reps),
            "full7_direction": expected,
            "full7_evidence": evidence,
            "direction_stability": direction_stability,
            "effect_stability": effect_stability,
            "leading_edge_stability": le_stability,
            "minimum_tested_lodo_n": min_tested,
            "minimum_same_direction_n": min_same,
            "median_abs_delta_NES": median_delta,
            "maximum_abs_delta_NES": max_delta,
            "median_leading_edge_jaccard": median_le,
            "most_influential_donor": most_donor,
            "stability_class": pclass,
            "notes": "Program synthesis retains individual pathway classes; it uses official Full7 FDR<0.10 representatives when available and conservatively downgrades mixed evidence. Full7 remains primary.",
        })
    return pd.DataFrame(rows)


def build_special_audits(pathway_master: pd.DataFrame, pairs: pd.DataFrame,
                         representatives: pd.DataFrame, programs: pd.DataFrame, log) -> None:
    metabolic_names = {
        "Mitochondrial respiration", "Oxidative phosphorylation", "Fatty acid metabolism",
        "Metabolic programs", "Mitochondrial translation", "Respiratory electron transport", "Oxidative metabolism"
    }
    metabolic = representatives.loc[representatives["program"].isin(metabolic_names)].copy()
    metabolic["cross_lineage_conclusion"] = np.where(
        metabolic["stability_class"].isin(["ROBUST", "MODERATELY_STABLE"]),
        "FORWARD_CANDIDATE", "NOT_FORMAL_FORWARD_CANDIDATE"
    )
    write_csv_atomic(metabolic, SUMMARY / "CROSS_LINEAGE_METABOLIC_LODO_SUMMARY.csv")

    mapping = pd.read_csv(PB / "gsea" / "mapping" / "phase4_feature_to_geneset_id.csv", low_memory=False)
    target_to_symbol = (
        mapping.dropna(subset=["target_geneset_id"]).assign(target_geneset_id=lambda x: x["target_geneset_id"].astype(str).str.replace(r"\.0$", "", regex=True))
        .drop_duplicates("target_geneset_id").set_index("target_geneset_id")["target_symbol"].astype(str).to_dict()
    )
    ieg_genes = {"FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "DUSP1", "DUSP2", "DUSP5"}
    stress_pattern = r"unfolded protein|endoplasmic reticulum stress|heat stress|mapk|ap.?1|immediate early|response to stress"
    epi_official = pathway_master.loc[
        pathway_master["lineage"].eq("Epithelial")
        & pathway_master["pathway_name"].str.contains(stress_pattern, case=False, regex=True, na=False)
        & pathway_master["full_direction"].eq("CIE_higher")
    ].sort_values("full_padj").head(50)
    stress_keys = set(zip(epi_official["collection"], epi_official["pathway_id"]))
    audit_rows = []
    for _, row in pairs.loc[pairs["lineage"].eq("Epithelial")].iterrows():
        if (row["collection"], row["pathway_id"]) not in stress_keys:
            continue
        full_ids = set(str(row.get("pair_full_leading_edge_target_ids", "")).split(";")) - {"", "nan"}
        lodo_ids = set(str(row.get("lodo_leading_edge_target_ids", "")).split(";")) - {"", "nan"}
        full_symbols = {target_to_symbol.get(x, x) for x in full_ids}
        lodo_symbols = {target_to_symbol.get(x, x) for x in lodo_ids}
        audit_rows.append({
            "lineage": "Epithelial", "collection": row["collection"], "pathway_id": row["pathway_id"],
            "pathway_name": row["pathway_name"], "full7_official_NES": epi_official.loc[
                epi_official["collection"].eq(row["collection"]) & epi_official["pathway_id"].eq(row["pathway_id"]), "full_NES"
            ].iloc[0],
            "lodo_run": row["lodo_run"], "short_label": row["short_label"], "dropped_donor": row["dropped_donor"],
            "pair_full_NES": row["pair_full_NES"], "lodo_NES": row["lodo_NES"], "delta_NES": row["delta_NES"],
            "same_direction": row["same_direction"], "leading_edge_jaccard": row["leading_edge_jaccard"],
            "pair_full_IEG_genes": semicolon(sorted(full_symbols & ieg_genes)),
            "lodo_IEG_genes": semicolon(sorted(lodo_symbols & ieg_genes)),
            "pair_full_IEG_n": len(full_symbols & ieg_genes), "lodo_IEG_n": len(lodo_symbols & ieg_genes),
            "cie2_cie4_focus": row["short_label"] in {"drop_cie2", "drop_cie4"},
            "notes": "IEG leading-edge membership is descriptive; pathway comparison uses paired common universe."
        })
    ieg = pd.DataFrame(audit_rows)
    write_csv_atomic(ieg, SUMMARY / "EPITHELIAL_STRESS_IEG_LODO_AUDIT.csv")

    # Pre-registered therapeutic axes: robustness observation only, never a treatment recommendation.
    axis_specs = {
        ("Epithelial", "TNF/NF-kappaB"): ("Hallmark", "TNFA_SIGNALING_VIA_NFKB"),
        ("Epithelial", "JAK-STAT"): ("Hallmark", "IL2_STAT5_SIGNALING"),
        ("Epithelial", "IL-23/STAT3"): ("Hallmark", "IL6_JAK_STAT3_SIGNALING"),
        ("T", "TNF/NF-kappaB"): ("Hallmark", "TNFA_SIGNALING_VIA_NFKB"),
        ("T", "JAK-STAT"): ("Reactome", "Signaling by Interleukins"),
        ("T", "IL-23/STAT3"): ("Hallmark", "IL6_JAK_STAT3_SIGNALING"),
        ("Myeloid", "TNF/NF-kappaB"): ("Hallmark", "TNFA_SIGNALING_VIA_NFKB"),
        ("Myeloid", "JAK-STAT"): ("Hallmark", "IL2_STAT5_SIGNALING"),
        ("Myeloid", "IL-23/STAT3"): ("Hallmark", "IL6_JAK_STAT3_SIGNALING"),
    }
    axis_rows = []
    for (lineage, axis), (collection, exact_name) in axis_specs.items():
        found = pathway_master.loc[
            pathway_master["lineage"].eq(lineage)
            & pathway_master["collection"].eq(collection)
            & pathway_master["pathway_name"].str.casefold().eq(exact_name.casefold())
        ]
        if len(found):
            best = found.iloc[0].to_dict()
            best.update({"axis": axis, "robustness_observation_only": True,
                         "therapeutic_interpretation": "NOT_PERMITTED_IN_PHASE5",
                         "selection_note": "Exact pre-registered Phase4 representative; not chosen by LODO."})
            axis_rows.append(best)
    write_csv_atomic(pd.DataFrame(axis_rows), SUMMARY / "THERAPEUTIC_AXIS_LODO_ROBUSTNESS.csv")
    print(f"Special audits: {len(metabolic)} metabolic representatives; {len(ieg)} Epithelial stress/IEG rows", file=log)


def build_donor_influence(gene_master: pd.DataFrame, gene_effects: pd.DataFrame,
                          pathway_master: pd.DataFrame, pairs: pd.DataFrame,
                          representatives: pd.DataFrame, log) -> pd.DataFrame:
    rows = []
    rep_keys = set(zip(representatives["lineage"], representatives["collection"], representatives["pathway_id"]))
    pair_core = pairs.loc[pairs.apply(lambda r: (r["lineage"], r["collection"], r["pathway_id"]) in rep_keys, axis=1)].copy()
    for lineage in LINEAGES:
        for short in RUN_ORDER:
            donor = DONOR_FROM_SHORT[short]
            gl = gene_effects.loc[gene_effects["lineage"].eq(lineage) & gene_effects["short_label"].eq(short)]
            pp = pairs.loc[pairs["lineage"].eq(lineage) & pairs["short_label"].eq(short)]
            pc = pair_core.loc[pair_core["lineage"].eq(lineage) & pair_core["short_label"].eq(short)]
            gene_most = gene_master.loc[gene_master["lineage"].eq(lineage) & gene_master["most_influential_donor"].eq(donor)]
            path_most = pathway_master.loc[pathway_master["lineage"].eq(lineage) & pathway_master["most_influential_donor"].eq(donor)]
            affected_programs = representatives.loc[
                representatives["lineage"].eq(lineage) & representatives["most_influential_donor"].eq(donor), "program"
            ].dropna().unique()
            rows.append({
                "donor": donor,
                "short_label": short,
                "dropped_condition": "Healthy" if short.startswith("drop_h") else "CIE",
                "lineage": lineage,
                "n_genes_most_influenced": len(gene_most),
                "n_pathways_most_influenced": len(path_most),
                "median_abs_delta_logFC": gl["abs_delta_logFC"].median(),
                "median_abs_delta_NES": pp["delta_NES"].abs().median(),
                "median_core_program_abs_delta_NES": pc["delta_NES"].abs().median(),
                "max_core_program_abs_delta_NES": pc["delta_NES"].abs().max(),
                "core_programs_affected": semicolon(sorted(affected_programs)),
                "notes": "Sensitivity description only; Healthy-drop (2v4) and CIE-drop (3v3) effects are not a mechanically fair individual ranking."
            })
    influence = pd.DataFrame(rows)
    write_csv_atomic(influence, SUMMARY / "DONOR_INFLUENCE_SUMMARY.csv")
    print("Donor influence summary written", file=log)
    return influence


def savefig(fig: plt.Figure, filename: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIGURES / filename, dpi=180, bbox_inches="tight")
    plt.close(fig)


def draw_heatmap(data: pd.DataFrame, title: str, cmap: str, filename: str,
                 center_zero: bool = False, value_format: str = ".2f") -> None:
    if data.empty:
        return
    values = data.to_numpy(dtype=float)
    if center_zero:
        max_abs = np.nanmax(np.abs(values))
        vmin, vmax = -max_abs, max_abs
    else:
        vmin, vmax = np.nanmin(values), np.nanmax(values)
    fig_h = max(4, 0.38 * len(data) + 2)
    fig, ax = plt.subplots(figsize=(10, fig_h))
    image = ax.imshow(values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(data.columns)), data.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(data.index)), data.index)
    ax.set_title(title)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            if np.isfinite(values[i, j]):
                ax.text(j, i, format(values[i, j], value_format), ha="center", va="center", fontsize=6)
    fig.colorbar(image, ax=ax, shrink=0.75)
    savefig(fig, filename)


def build_figures(gene_master: pd.DataFrame, gene_effects: pd.DataFrame,
                  pathway_master: pd.DataFrame, pairs: pd.DataFrame,
                  representatives: pd.DataFrame, influence: pd.DataFrame, log) -> None:
    colors = {"Full7": "#111111", "drop_h1": "#3977A8", "drop_h2": "#4C91C6", "drop_h3": "#73ADD0",
              "drop_cie2": "#C54F3D", "drop_cie3": "#D77855", "drop_cie4": "#A12A36", "drop_cie6": "#E59B78"}
    for lineage in LINEAGES:
        # Anchor gene logFC forest.
        subset = gene_effects.loc[
            gene_effects["lineage"].eq(lineage) & gene_effects["gene_symbol"].isin(ANCHORS[lineage])
        ].copy()
        fig, axes = plt.subplots(len(ANCHORS[lineage]), 1, figsize=(10, 2.1 * len(ANCHORS[lineage])), sharex=True)
        axes = np.atleast_1d(axes)
        for ax, gene in zip(axes, ANCHORS[lineage]):
            g = subset.loc[subset["gene_symbol"].eq(gene)]
            if len(g):
                full_value = float(g["logFC_full"].iloc[0])
                values = {"Full7": full_value, **dict(zip(g["short_label"], g["logFC_lodo"]))}
                ax.axvline(0, color="#AAAAAA", lw=0.8)
                for j, label in enumerate(DISPLAY_RUNS):
                    if label in values:
                        ax.scatter(values[label], j, s=36, color=colors[label], zorder=3)
                ax.set_yticks(range(len(DISPLAY_RUNS)), DISPLAY_RUNS, fontsize=7)
                ax.set_title(gene, loc="left", fontsize=10)
        axes[-1].set_xlabel("edgeR logFC (CIE vs Healthy)")
        fig.suptitle(f"{lineage}: anchor gene Full7 + seven LODO effects", y=1.005)
        savefig(fig, f"{lineage}_anchor_gene_logFC_LODO_forest.png")

        # Gene direction-stability histogram.
        gm = gene_master.loc[gene_master["lineage"].eq(lineage) & gene_master["tested_lodo_n"].eq(7)]
        counts = gm["same_direction_n"].value_counts().reindex(range(8), fill_value=0)
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.bar(counts.index, counts.values, color="#4C78A8")
        ax.set(xlabel="LODO rounds matching Full7 direction", ylabel="Full7-tested genes", title=f"{lineage}: gene direction stability")
        ax.set_xticks(range(8))
        savefig(fig, f"{lineage}_gene_direction_stability_histogram.png")

        # Pathway direction-stability histogram by collection.
        pm = pathway_master.loc[pathway_master["lineage"].eq(lineage) & pathway_master["tested_lodo_n"].eq(7)]
        fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
        for ax, collection in zip(axes, ["Hallmark", "Reactome", "GO_BP"]):
            counts = pm.loc[pm["collection"].eq(collection), "same_direction_n"].value_counts().reindex(range(8), fill_value=0)
            ax.bar(counts.index, counts.values, color="#59A14F")
            ax.set_title(collection)
            ax.set_xlabel("Matching rounds")
            ax.set_xticks(range(8))
        axes[0].set_ylabel("Pathways")
        fig.suptitle(f"{lineage}: pathway direction stability")
        savefig(fig, f"{lineage}_pathway_direction_stability_histogram.png")

        reps = representatives.loc[
            representatives["lineage"].eq(lineage) & representatives["pathway_id"].astype(str).ne("")
        ].copy()
        reps["row_label"] = reps.apply(lambda r: f"{r['program']} | {r['collection']}", axis=1)
        rep_keys = set(zip(reps["collection"], reps["pathway_id"]))
        pair_rep = pairs.loc[pairs["lineage"].eq(lineage) & pairs.apply(lambda r: (r["collection"], r["pathway_id"]) in rep_keys, axis=1)].copy()

        # Core pathway NES forest (one panel per program, LODO values from paired common universe).
        programs = list(dict.fromkeys(reps["program"]))
        fig, axes = plt.subplots(len(programs), 1, figsize=(11, max(6, 1.75 * len(programs))), sharex=False)
        axes = np.atleast_1d(axes)
        for ax, program in zip(axes, programs):
            chosen = reps.loc[reps["program"].eq(program)].sort_values("full_padj").iloc[0]
            q = pair_rep.loc[pair_rep["collection"].eq(chosen["collection"]) & pair_rep["pathway_id"].eq(chosen["pathway_id"])]
            values = {"Full7": chosen["full_NES"], **dict(zip(q["short_label"], q["lodo_NES"]))}
            ax.axvline(0, color="#AAAAAA", lw=0.8)
            for j, label in enumerate(DISPLAY_RUNS):
                if label in values and pd.notna(values[label]):
                    ax.scatter(values[label], j, s=30, color=colors[label])
            ax.set_yticks(range(len(DISPLAY_RUNS)), DISPLAY_RUNS, fontsize=6)
            ax.set_title(f"{program} — {chosen['collection']}:{chosen['pathway_name']}", loc="left", fontsize=8)
        axes[-1].set_xlabel("NES (Full7 official; LODO paired common universe)")
        fig.suptitle(f"{lineage}: core pathway NES LODO forest", y=1.002)
        savefig(fig, f"{lineage}_core_pathway_NES_LODO_forest.png")

        # NES and leading-edge heatmaps for unique representative pathways.
        unique_reps = reps.drop_duplicates(["collection", "pathway_id"]).head(18).copy()
        matrix_rows, le_rows = [], []
        labels = []
        for _, rep in unique_reps.iterrows():
            q = pair_rep.loc[pair_rep["collection"].eq(rep["collection"]) & pair_rep["pathway_id"].eq(rep["pathway_id"])]
            matrix_rows.append([rep["full_NES"], *[q.loc[q["short_label"].eq(s), "lodo_NES"].iloc[0] if q["short_label"].eq(s).any() else np.nan for s in RUN_ORDER]])
            le_rows.append([np.nan, *[q.loc[q["short_label"].eq(s), "leading_edge_jaccard"].iloc[0] if q["short_label"].eq(s).any() else np.nan for s in RUN_ORDER]])
            labels.append(f"{rep['program']} | {rep['collection']}")
        nes = pd.DataFrame(matrix_rows, index=labels, columns=DISPLAY_RUNS)
        draw_heatmap(nes, f"{lineage}: core program NES stability", "RdBu_r", f"{lineage}_NES_stability_heatmap.png", center_zero=True)
        le = pd.DataFrame(le_rows, index=labels, columns=DISPLAY_RUNS)
        draw_heatmap(le.drop(columns="Full7"), f"{lineage}: paired leading-edge Jaccard", "YlGnBu",
                     f"{lineage}_leading_edge_Jaccard_heatmap.png", center_zero=False)

    # Donor x lineage influence heatmap (all pathways and genes; standardized ranks avoided).
    donor_matrix = influence.pivot(index="lineage", columns="short_label", values="median_core_program_abs_delta_NES").reindex(
        index=LINEAGES, columns=RUN_ORDER
    )
    draw_heatmap(donor_matrix, "Donor influence on core programs: median |delta NES|", "OrRd",
                 "donor_influence_core_program_heatmap.png")
    print(f"Figures written: {len(list(FIGURES.glob('*.png')))} PNG files", file=log)


def write_software_versions() -> None:
    versions = pd.DataFrame([
        {"software": "Python", "version": platform.python_version(), "role": "Phase5 aggregation/report runtime"},
        {"software": "pandas", "version": pd.__version__, "role": "Tabular aggregation"},
        {"software": "numpy", "version": np.__version__, "role": "Numerical summaries"},
        {"software": "matplotlib", "version": matplotlib.__version__, "role": "Donor-level stability figures"},
    ])
    write_csv_atomic(versions, ROOT / "00_metadata" / "PHASE5_SOFTWARE_VERSIONS_PYTHON.csv")


def main() -> None:
    SUMMARY.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("w", encoding="utf-8") as log:
        print("Phase5 stability summary", file=log)
        print(f"Python {platform.python_version()}; pandas {pd.__version__}; numpy {np.__version__}; matplotlib {matplotlib.__version__}", file=log)
        print("Full7 remains primary; classifications follow PHASE5_STABILITY_PARAMETERS.csv", file=log)
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            gene_master, gene_effects, gene_ranks = build_gene_outputs(log)
            pathway_master, pairs, representatives, programs = build_pathway_outputs(log)
            build_special_audits(pathway_master, pairs, representatives, programs, log)
            influence = build_donor_influence(gene_master, gene_effects, pathway_master, pairs, representatives, log)
            build_figures(gene_master, gene_effects, pathway_master, pairs, representatives, influence, log)
            write_software_versions()
            if captured:
                print("Warnings:", file=log)
                for warning in captured:
                    print(f"- {warning.message}", file=log)
        print("Output hashes:", file=log)
        for path in sorted([*SUMMARY.glob("*"), *FIGURES.glob("*.png")]):
            if path.is_file():
                print(f"{path.relative_to(ROOT)}\t{path.stat().st_size}\t{sha256(path)}", file=log)
        print("Phase5 stability summaries complete", file=log)


if __name__ == "__main__":
    main()
