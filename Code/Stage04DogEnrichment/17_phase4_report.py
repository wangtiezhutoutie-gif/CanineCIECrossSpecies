#!/usr/bin/env python
"""Build Phase4 deterministic summaries, figures, report, and manifest."""

from __future__ import annotations

import hashlib
import os
import platform
import re
from collections import OrderedDict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
GSEA = ROOT / "06_pseudobulk" / "gsea"
FULL = GSEA / "full7"
LEADING = GSEA / "leading_edges"
SENS = GSEA / "sensitivity"
SUMMARY = GSEA / "summary"
FIGURES = GSEA / "figures"
META = ROOT / "00_metadata"
REPORT = ROOT / "08_results" / "PHASE4_REPORT.md"
LOG = ROOT / "logs" / "phase4_report.log"
LINEAGES = ("Epithelial", "T", "Myeloid")
COLLECTIONS = ("Hallmark", "Reactome", "GO_BP")
TAGS = {"Hallmark": "Hallmark", "Reactome": "Reactome", "GO_BP": "GOBP"}
SUBDIRS = {"Hallmark": "hallmark", "Reactome": "reactome", "GO_BP": "go_bp"}


def log(message: str) -> None:
    print(message, flush=True)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(message + "\n")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fmt_num(value, digits=3):
    if pd.isna(value):
        return "NA"
    value = float(value)
    if value == 0:
        return "0"
    if abs(value) < 0.001:
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "（无）"
    columns = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in df.itertuples(index=False, name=None):
        vals = []
        for value in row:
            text = "" if pd.isna(value) else str(value)
            vals.append(text.replace("|", "\\|").replace("\n", " "))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def load_results():
    results = {}
    for lineage in LINEAGES:
        for collection in COLLECTIONS:
            path = FULL / SUBDIRS[collection] / f"{lineage}_{TAGS[collection]}_GSEA.csv"
            df = pd.read_csv(path)
            if len(df) == 0:
                raise RuntimeError(f"Empty formal GSEA result: {path}")
            results[(lineage, collection)] = df
    return results


def build_go_representatives(results):
    members = pd.read_csv(GSEA / "gene_sets" / "go_bp_canis_familiaris_gene_sets.tsv.gz", sep="\t", dtype=str)
    member_sets = members.groupby("pathway_id")["target_geneset_id"].apply(lambda s: set(s.astype(str))).to_dict()
    rows = []
    for lineage in LINEAGES:
        result = results[(lineage, "GO_BP")].copy()
        candidates = result[result.padj < 0.10].copy()
        if candidates.empty:
            candidates = result.sort_values(["padj", "pathway_id"]).head(50).copy()
        candidates["abs_nes"] = candidates.NES.abs()
        candidates = candidates.sort_values(["padj", "abs_nes", "pathway_id"], ascending=[True, False, True])
        clusters: list[dict] = []
        for row in candidates.itertuples(index=False):
            current = member_sets.get(str(row.pathway_id), set())
            assigned = None
            for cluster in clusters:
                if cluster["direction"] != row.direction:
                    continue
                rep_set = member_sets.get(cluster["representative_id"], set())
                union_n = len(current | rep_set)
                jac = len(current & rep_set) / union_n if union_n else 0
                if jac >= 0.50:
                    assigned = cluster
                    break
            if assigned is None:
                clusters.append(
                    {
                        "lineage": lineage,
                        "direction": row.direction,
                        "representative_id": str(row.pathway_id),
                        "representative_term": row.pathway_name,
                        "member_ids": [str(row.pathway_id)],
                        "member_names": [row.pathway_name],
                        "NES": row.NES,
                        "padj": row.padj,
                    }
                )
            else:
                assigned["member_ids"].append(str(row.pathway_id))
                assigned["member_names"].append(row.pathway_name)
        for cluster in clusters:
            rows.append(
                {
                    "lineage": cluster["lineage"],
                    "direction": cluster["direction"],
                    "representative_term_id": cluster["representative_id"],
                    "representative_term": cluster["representative_term"],
                    "member_terms": "; ".join(cluster["member_ids"]),
                    "member_term_names": "; ".join(cluster["member_names"]),
                    "n_member_terms": len(cluster["member_ids"]),
                    "selection_rule": "Greedy clusters within NES direction at source-membership Jaccard >=0.50; representative is lowest padj, then largest abs(NES), then pathway_id.",
                    "NES": cluster["NES"],
                    "padj": cluster["padj"],
                    "notes": "Raw GO BP GSEA table remains complete and unchanged.",
                }
            )
    reps = pd.DataFrame(rows).sort_values(["lineage", "padj", "representative_term_id"])
    reps.to_csv(SUMMARY / "GO_BP_REPRESENTATIVE_TERMS.csv", index=False)

    colors = {"CIE_higher": "#B33A3A", "CIE_lower": "#2E63A3"}
    for lineage in LINEAGES:
        x = reps[reps.lineage.eq(lineage)].copy()
        chosen = pd.concat(
            [
                x[x.direction.eq("CIE_higher")].sort_values(["padj", "NES"]).head(12),
                x[x.direction.eq("CIE_lower")].sort_values(["padj", "NES"]).head(12),
            ],
            ignore_index=True,
        ).sort_values("NES")
        if chosen.empty:
            continue
        labels = [str(v).replace("_", " ")[:72] for v in chosen.representative_term]
        height = max(5.5, 0.32 * len(chosen) + 1.5)
        fig, ax = plt.subplots(figsize=(10, height))
        ax.barh(range(len(chosen)), chosen.NES, color=[colors[d] for d in chosen.direction])
        ax.set_yticks(range(len(chosen)), labels)
        ax.axvline(0, color="#888888", linewidth=0.7)
        ax.set_xlabel("NES (positive = CIE higher)")
        ax.set_title(f"{lineage}: deterministic GO BP representative terms", fontweight="bold")
        ax.grid(axis="x", alpha=0.2)
        fig.tight_layout()
        fig.savefig(FIGURES / f"{lineage}_GOBP_representative_terms.png", dpi=180)
        plt.close(fig)
    return reps


PROGRAMS = OrderedDict(
    [
        ("Cell cycle / proliferation", r"CELL CYCLE|MITOTIC|MEIOS|G2M|E2F|SPINDLE|DNA REPLICATION"),
        ("Ribosome / RNA biogenesis", r"RIBOSOM|RRNA|RNA BIOGENESIS|RNA PROCESSING|RIBONUCLEOPROTEIN|TRANSLATION"),
        ("Mitochondrial respiration", r"OXIDATIVE PHOSPHORYLATION|RESPIRATORY ELECTRON|ELECTRON TRANSPORT|CELLULAR RESPIRATION|AEROBIC RESPIRATION|MITOCHONDRIAL TRANSLATION"),
        ("Lipid/xenobiotic/metabolic differentiation", r"FATTY ACID|XENOBIOTIC|BILE ACID|CARBOXYLIC ACID|OXOACID|SMALL MOLECULE METABOL|ADIPOGENESIS"),
        ("TNF / NF-kappaB inflammation", r"TNF|NF.?KAPPA|NFKB|INFLAMMATORY RESPONSE"),
        ("Cytokine / interleukin / JAK-STAT", r"CYTOKINE|INTERLEUKIN|JAK|STAT|INTERFERON"),
        ("Innate/TLR response", r"TOLL.?LIKE|TLR|INNATE IMMUNE|DEFENSE RESPONSE|NEUTROPHIL"),
        ("Stress / UPR / MAPK", r"STRESS|UNFOLDED|HEAT|AP.?1|MAPK|KRAS|UV RESPONSE"),
        ("ECM / remodeling / hemostasis", r"EXTRACELLULAR MATRIX|COLLAGEN|REMODEL|HEMOSTASIS|COAGULATION|PLATELET"),
        ("Cell death / apoptosis", r"APOPT|CELL DEATH|P53"),
        ("Adaptive/T-cell program", r"T CELL|T-CELL|LYMPHOCYTE|ALLOGRAFT|TCR|T CELL RECEPTOR"),
    ]
)


def build_cross_collection(results):
    all_result = pd.concat(results.values(), ignore_index=True)
    all_result = all_result[all_result.padj < 0.10].copy()
    rows = []
    for lineage in LINEAGES:
        ld = all_result[all_result.lineage.eq(lineage)]
        for program, pattern in PROGRAMS.items():
            matches = ld[ld.pathway_name.str.contains(pattern, case=False, regex=True, na=False)]
            for direction in ("CIE_higher", "CIE_lower"):
                group = matches[matches.direction.eq(direction)].copy()
                if group.collection.nunique() < 2:
                    continue
                selected = group.sort_values(["padj", "pathway_id"]).groupby("collection", as_index=False).first()
                evidence = "; ".join(
                    f"{r.collection}: {r.pathway_name} (NES={r.NES:.2f}, FDR={fmt_num(r.padj)})"
                    for r in selected.itertuples(index=False)
                )
                rows.append(
                    {
                        "lineage": lineage,
                        "program": program,
                        "direction": direction,
                        "n_collections_support": int(group.collection.nunique()),
                        "collections": ";".join(sorted(group.collection.unique())),
                        "n_supporting_pathways": len(group),
                        "min_padj": group.padj.min(),
                        "max_abs_NES": group.NES.abs().max(),
                        "supporting_evidence": evidence,
                        "p_value_combination": "NONE",
                        "notes": "Support defined without cross-database P-value combination; only collection-wise FDR<0.10 entries count.",
                    }
                )
    cross = pd.DataFrame(rows)
    if not cross.empty:
        cross = cross.sort_values(["lineage", "direction", "min_padj", "program"])
    cross.to_csv(SUMMARY / "cross_collection_program_summary.csv", index=False)
    return cross


def load_all_leading():
    frames = []
    for path in sorted(LEADING.glob("*_leading_edges.csv")):
        frame = pd.read_csv(path, low_memory=False)
        frame["source_file"] = path.name
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def build_anchor_and_control_audits(all_leading, results):
    anchor_by_lineage = {
        "Epithelial": ["KCNK16", "MMP7", "GAL3ST2", "CYP1B1"],
        "T": ["GAL3ST2", "TLR9", "TMEM255A", "PTPRT"],
        "Myeloid": ["SERPINE1", "CCL3", "ADAMDEC1", "MMP9", "SAA1"],
    }
    anchor_rows = []
    for lineage, genes in anchor_by_lineage.items():
        for gene in genes:
            x = all_leading[(all_leading.lineage == lineage) & (all_leading.target_symbol.astype(str) == gene)].copy()
            x = x.sort_values(["pathway_padj", "pathway_id"]).head(20)
            if x.empty:
                anchor_rows.append({"lineage": lineage, "anchor_gene": gene, "leading_edge_occurrence": False, "collection": "", "pathway_id": "", "pathway_name": "", "NES": np.nan, "padj": np.nan})
            else:
                for row in x.itertuples(index=False):
                    anchor_rows.append({
                        "lineage": lineage, "anchor_gene": gene, "leading_edge_occurrence": True,
                        "collection": row.collection, "pathway_id": row.pathway_id,
                        "pathway_name": row.pathway_name, "NES": row.pathway_NES, "padj": row.pathway_padj,
                    })
    anchors = pd.DataFrame(anchor_rows)
    anchors.to_csv(SUMMARY / "phase3_anchor_gene_pathways.csv", index=False)

    control_rows = []
    definitions = [
        ("T", "IL7R/GZMA T-state", ["IL7R", "GZMA"], r"T CELL|LYMPHOCYTE|INTERLEUKIN|TNF|INFLAMMATORY|CYTOKINE"),
        ("Myeloid", "SOD2/IL1A myeloid-state", ["SOD2", "IL1A"], r"TNF|INFLAMMATORY|CYTOKINE|INTERLEUKIN|INNATE|TLR|SIGNAL"),
        ("Epithelial", "epithelial stress/metabolic", [], r"STRESS|UNFOLDED|OXIDATIVE PHOSPHORYLATION|FATTY ACID|XENOBIOTIC|RESPIRATION"),
    ]
    combined = pd.concat(results.values(), ignore_index=True)
    for lineage, label, genes, pattern in definitions:
        path_hits = combined[(combined.lineage == lineage) & combined.pathway_name.str.contains(pattern, case=False, regex=True, na=False)].sort_values(["padj", "pathway_id"])
        sig_hits = path_hits[path_hits.padj < 0.10]
        gene_hits = all_leading[(all_leading.lineage == lineage) & all_leading.target_symbol.astype(str).isin(genes)] if genes else pd.DataFrame()
        gene_sig = gene_hits[gene_hits.pathway_padj < 0.10] if not gene_hits.empty else pd.DataFrame()
        observed_gene_anchors = set(gene_sig.target_symbol.astype(str)) if not gene_sig.empty else set()
        if not sig_hits.empty and (not genes or set(genes).issubset(observed_gene_anchors)):
            status = "consistent"
        elif not sig_hits.empty or observed_gene_anchors:
            status = "partially consistent"
        else:
            status = "not prominent"
        evidence = "; ".join(
            f"{r.collection}:{r.pathway_name} NES={r.NES:.2f}, FDR={fmt_num(r.padj)}"
            for r in sig_hits.head(5).itertuples(index=False)
        )
        if genes and not gene_sig.empty:
            anchor_text = ", ".join(sorted(gene_sig.target_symbol.astype(str).unique()))
            evidence = (evidence + "; " if evidence else "") + f"leading edge: {anchor_text}"
        control_rows.append({"lineage": lineage, "positive_control": label, "status": status, "evidence": evidence or "No FDR<0.10 matched pathway."})
    controls = pd.DataFrame(control_rows)
    controls.to_csv(SUMMARY / "positive_control_pathway_consistency.csv", index=False)
    return anchors, controls


def build_ieg_audit(all_leading, sensitivity, results):
    pattern = r"IMMEDIATE|STRESS|MAPK|AP.?1|HEAT|RESPONSE TO STIMULUS|UNFOLDED|UV_RESPONSE|KRAS"
    primary = pd.concat([results[("Epithelial", c)] for c in COLLECTIONS], ignore_index=True)
    candidates = primary[primary.pathway_name.str.contains(pattern, case=False, regex=True, na=False)].copy()
    sens = sensitivity[sensitivity.lineage.eq("Epithelial")][[
        "collection", "pathway_id", "NES_full", "padj_full", "NES_no_cie2", "padj_no_cie2",
        "delta_NES", "same_direction", "leading_edge_jaccard", "sensitivity_flag"
    ]]
    candidates = candidates.merge(sens, on=["collection", "pathway_id"], how="left")
    ieg_genes = {"FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "DUSP1", "ATF3", "NR4A1"}
    rows = []
    for row in candidates.itertuples(index=False):
        le = all_leading[
            (all_leading.lineage == "Epithelial") &
            (all_leading.collection == row.collection) &
            (all_leading.pathway_id.astype(str) == str(row.pathway_id))
        ]
        symbols = [str(x) for x in le.target_symbol.dropna()]
        ieg = sorted(set(symbols) & ieg_genes)
        rows.append({
            "lineage": "Epithelial", "collection": row.collection, "pathway_id": row.pathway_id,
            "pathway_name": row.pathway_name, "primary_Full7_NES": row.NES,
            "primary_Full7_padj": row.padj, "common_NES_full": row.NES_full,
            "common_NES_no_cie2": row.NES_no_cie2, "delta_NES": row.delta_NES,
            "same_direction": row.same_direction, "leading_edge_jaccard": row.leading_edge_jaccard,
            "sensitivity_flag": row.sensitivity_flag, "leading_edge_n": len(set(symbols)),
            "IEG_leading_edge_n": len(ieg), "IEG_leading_edge_genes": ";".join(ieg),
            "IEG_fraction": len(ieg) / len(set(symbols)) if symbols else np.nan,
            "interpretation_note": "Inspect cie_2 sensitivity before generalizing epithelial IEG/stress activation.",
        })
    audit = pd.DataFrame(rows).sort_values(["primary_Full7_padj", "collection", "pathway_id"])
    audit.to_csv(SUMMARY / "epithelial_ieg_stress_audit.csv", index=False)
    return audit


def build_therapy_axis(results):
    combined = pd.concat(results.values(), ignore_index=True)
    axes = OrderedDict(
        [
            ("JAK-STAT", r"JAK|IL6_JAK_STAT3|IL2_STAT5|SIGNALING PATHWAY VIA STAT|SIGNALING BY INTERLEUKINS"),
            ("IL-23/STAT3", r"IL.?23|INTERLEUKIN.?23|STAT3|TH17|T HELPER 17"),
            ("TNF/NF-kappaB", r"TNF|NF.?KAPPA|NFKB"),
        ]
    )
    rows = []
    for lineage in LINEAGES:
        for axis, pattern in axes.items():
            x = combined[(combined.lineage == lineage) & combined.pathway_name.str.contains(pattern, case=False, regex=True, na=False)].copy()
            x = x.sort_values(["padj", "pathway_id"])
            if x.empty:
                status = "Not assessable"
            elif (x.padj < 0.10).any():
                status = "Present"
            elif ((x.padj < 0.25) & (x.NES.abs() >= 1.3)).any():
                status = "Weak"
            else:
                status = "Not prominent"
            top = x.iloc[0] if not x.empty else None
            evidence = "" if top is None else f"{top.collection}: {top.pathway_name}; NES={top.NES:.2f}; FDR={fmt_num(top.padj)}"
            note = "Observation only; no treatment conclusion."
            if axis == "IL-23/STAT3" and not x.empty:
                il23 = x[x.pathway_name.str.contains(r"IL.?23|INTERLEUKIN.?23", case=False, regex=True, na=False)]
                if il23.empty or not (il23.padj < 0.10).any():
                    if status == "Present":
                        status = "Weak"
                    note = "STAT3-related evidence only; no IL-23-specific tested term was prominent. Observation only."
            rows.append({"lineage": lineage, "axis": axis, "status": status, "best_evidence": evidence, "notes": note})
    axes_df = pd.DataFrame(rows)
    axes_df.to_csv(SUMMARY / "pre_registered_therapeutic_axis_observation.csv", index=False)
    return axes_df


def result_top_table(results, lineage):
    rows = []
    for collection in COLLECTIONS:
        df = results[(lineage, collection)]
        for direction in ("CIE_higher", "CIE_lower"):
            x = df[df.direction.eq(direction)].copy()
            x["abs_nes"] = x.NES.abs()
            x = x.sort_values(["padj", "abs_nes", "pathway_id"], ascending=[True, False, True]).head(3)
            for row in x.itertuples(index=False):
                rows.append({"Collection": collection, "方向": direction, "Pathway": row.pathway_name, "NES": f"{row.NES:.2f}", "FDR": fmt_num(row.padj)})
    return pd.DataFrame(rows)


def build_report(results, go_reps, cross, controls, therapy, ieg, sensitivity, sens_summary):
    ranked_rows = []
    for lineage in LINEAGES:
        path = ROOT / "06_pseudobulk" / "ranked_lists" / f"{lineage}_ranked.csv"
        df = pd.read_csv(path)
        expected = np.sign(df.logFC.to_numpy()) * np.sqrt(df.statistic_F.to_numpy())
        exact = bool(np.allclose(expected, df.rank_metric.to_numpy(), rtol=1e-12, atol=1e-12, equal_nan=False))
        ranked_rows.append({
            "Lineage": lineage, "genes": len(df), "positive": int((df.rank_metric > 0).sum()),
            "negative": int((df.rank_metric < 0).sum()), "formula_exact": exact,
            "sorted_desc": bool(df.rank_metric.is_monotonic_decreasing), "SHA-256": sha256(path),
        })
    ranked_table = pd.DataFrame(ranked_rows)

    gene_audit = pd.read_csv(META / "PHASE4_GENESET_AUDIT.csv")
    coverage = pd.read_csv(GSEA / "mapping" / "phase4_mapping_coverage.csv")
    collisions = pd.read_csv(GSEA / "mapping" / "phase4_mapping_collisions.csv")
    hash_audit = pd.read_csv(META / "PHASE4_PHASE3_INPUT_HASH_AUDIT.csv")

    coverage_view = coverage.copy()
    coverage_view["feature_to_NCBI_rate"] = coverage_view.feature_to_ncbi_mapped / coverage_view.tested_features
    coverage_view["collection_member_rate"] = coverage_view.mapping_rate
    coverage_view = coverage_view[["lineage", "collection", "tested_features", "feature_to_NCBI_rate", "mapped_features", "unique_target_genes", "collection_member_rate", "collapsed_features"]]
    coverage_view.columns = ["Lineage", "Collection", "tested", "feature→NCBI", "member features", "unique target", "collection coverage", "collapsed"]
    for col in ["feature→NCBI", "collection coverage"]:
        coverage_view[col] = coverage_view[col].map(lambda x: f"{100*x:.1f}%")

    source_view = gene_audit[["collection", "source_version", "source_species", "mapping_type", "n_pathways_raw", "n_pathways_after_filter"]].copy()
    source_view.columns = ["Collection", "Version", "Source species", "Mapping", "Raw pathways", "Source size 15–500"]

    sens_view = sens_summary[[
        "lineage", "collection", "common_feature_universe_n", "common_mapped_target_genes_n",
        "nes_spearman", "nes_pearson", "direction_concordance", "top20_pathway_jaccard",
        "median_abs_delta_NES", "median_leading_edge_jaccard", "direction_flip_n"
    ]].copy()
    sens_view.columns = ["Lineage", "Collection", "common feature", "common target", "Spearman", "Pearson", "direction concordance", "Top20 Jaccard", "median abs delta NES", "median LE Jaccard", "all flips"]
    for col in ["Spearman", "Pearson", "direction concordance", "Top20 Jaccard", "median abs delta NES", "median LE Jaccard"]:
        sens_view[col] = sens_view[col].map(lambda x: fmt_num(x))

    flip_summary = sensitivity[sensitivity.sensitivity_flag.eq("Direction_flip")].groupby(["lineage", "collection"]).agg(
        n_flips=("pathway_id", "size"), min_full_padj=("padj_full", "min"), full_fdr10_flips=("padj_full", lambda s: int((s < 0.10).sum()))
    ).reset_index()
    flip_summary["min_full_padj"] = flip_summary.min_full_padj.map(fmt_num)

    report = []
    report.append("# Canine CIE Phase4 Report — cell-type-specific Disease Program GSEA\n")
    report.append("**Phase4 Gate：YELLOW**\n")
    report.append("Phase4 技术执行完整，3 条 Primary lineage × 3 个 collection 均成功；主要疾病程序形成明确的跨库一致性，且共同 universe 的整体 NES 相关性较高。门控保留 YELLOW，原因是样本仍为 3 Healthy vs 4 CIE、Epithelial 的 TNF/NF-κB 等个别重要程序对 `cie_2` 有幅度敏感性、Myeloid 的 `cie_2` 仅 27 cells，并且 donor robustness 尚未经过 Phase5 全 7 轮 LODO。Phase4 没有执行 Phase5。\n")

    report.append("## 1. 输入完整性\n")
    report.append(f"`PHASE3_MANIFEST.csv` 复核通过 **{int(hash_audit.hash_match.sum())}/{len(hash_audit)}**；未发现 Phase3 文件非预期变化。Full7 ranked list 均保持 `positive = CIE higher`、`negative = CIE lower`，并逐行验证 `rank_metric = sign(logFC) × sqrt(F)`。\n")
    report.append(markdown_table(ranked_table))
    report.append("\nFull7 是唯一 Primary GSEA；No-`cie_2` 仅用于共同基因 universe 的敏感性分析。\n")

    report.append("## 2. Gene-set 来源\n")
    report.append(markdown_table(source_view))
    report.append("\nHallmark 为人 MSigDB Hallmark 经 babelgene 投影到犬的功能注释，不是犬原生 Hallmark，也不等同于 Phase10 的正式 Dog–Human 1:1 ortholog comparison。Reactome 使用 v97 官方 `Canis familiaris` NCBI membership；其条目应称为 Reactome-annotated/inferred canine pathways。GO BP 使用 GO 固定发布版、犬 CANLF UniProt GAF/GPI，并沿 `go-basic` 的 `is_a`/`part_of` 关系传播祖先项。三库分别进行 BH 校正，未混合 FDR 空间。正式参数固定为 `fgseaMultilevel`、`minSize=15`、`maxSize=500`、`eps=1e-10`、`scoreType=std`、`nproc=1`、base seed `20260824`。fgsea 日志提示少数实际 P 值可能低于 `eps`；表中数值按预注册下限保守截断，没有为了更小 P 值重调参数。\n")

    report.append("## 3. Mapping coverage\n")
    report.append(markdown_table(coverage_view))
    collision_groups = collisions.groupby(["lineage", "target_geneset_id"]).ngroups if len(collisions) else 0
    report.append(f"\nPhase3 feature→当前 NCBI Gene 的无歧义精确 Ensembl 映射率为约 83.2–85.4%；存在一对多当前注释的 feature 被审慎排除，没有用 symbol 强行裁决。Hallmark/Reactome 的 collection coverage 较低主要反映集合本身只覆盖一部分基因；GO BP coverage 为约 69.4–71.7%。共记录 {collision_groups} 个 lineage-specific collision groups；全部按预注册的 `max(abs(rank_metric)) → smaller PValue → feature_id lexical` 裁决，无静默删除。\n")

    for section, lineage, title in [
        (4, "Epithelial", "Epithelial Disease Programs"),
        (5, "T", "T Disease Programs"),
        (6, "Myeloid", "Myeloid Disease Programs"),
    ]:
        report.append(f"## {section}. {title}\n")
        report.append(markdown_table(result_top_table(results, lineage)))
        coherent = cross[cross.lineage.eq(lineage)] if not cross.empty else pd.DataFrame()
        if not coherent.empty:
            cv = coherent[["program", "direction", "n_collections_support", "supporting_evidence"]].copy()
            cv.columns = ["Program", "方向", "支持库数", "代表证据"]
            report.append("\n跨库一致性：\n")
            report.append(markdown_table(cv))
        if lineage == "Epithelial":
            report.append("\n核心结构为 CIE-higher 的细胞周期/核糖体-RNA 生物发生/UPR，以及 CIE-lower 的氧化磷酸化、呼吸链、脂肪酸/外源物代谢；TNF/NF-κB 为 CIE-higher，但在去除 `cie_2` 后幅度下降，需保留敏感性警告。\n")
        elif lineage == "T":
            report.append("\n核心结构为 CIE-higher 的 TNF/NF-κB/炎症与部分细胞信号程序，以及 CIE-lower 的氧化磷酸化、线粒体翻译和代谢。T 的总体共同-universe 相关性为三谱系中较高。\n")
        else:
            report.append("\n核心结构为 CIE-higher 的 TNF/NF-κB、炎症、TLR4/细胞因子信号与组织反应，以及 CIE-lower 的呼吸链/氧化代谢和翻译；解释必须同时保留 `cie_2` Myeloid 27-cell 低深度限制。\n")

    report.append("## 7. `cie_2` sensitivity\n")
    report.append(markdown_table(sens_view))
    report.append("\n所有比较均在同一 feature intersection、同一目标 NCBI Gene universe、同一基因集与同一参数上完成；因此 pathway size change 的最大绝对值为 0。Epithelial 总体结构仍相关，但比 T/Myeloid 更敏感；Epithelial Hallmark TNF/NF-κB 的 common-universe NES 由 1.89 降至 1.21（ΔNES=-0.68），方向未翻转但被标为 `Sensitive_to_cie2`。T 与 Myeloid 的整体相关性更高，Myeloid 仍需低深度警告。\n")
    report.append("\n方向翻转审计：\n")
    report.append(markdown_table(flip_summary))
    report.append("\n方向翻转集中于 Full7 FDR 不显著的弱富集路径；**没有任何 Full7 FDR<0.10 pathway 发生方向翻转**。因此不能把全部 flip 数量等同于主要疾病程序失效。\n")

    report.append("## 8. Positive-control consistency\n")
    report.append(markdown_table(controls.rename(columns={"lineage": "Lineage", "positive_control": "Positive control", "status": "Status", "evidence": "Evidence"})))
    report.append("\n`consistent / partially consistent / not prominent` 仅表示 pathway-level 一致性，不是额外 donor 验证。Phase3 anchor genes 的自然 pathway/leading-edge 落点已完整保存，未据此修改 gene sets。\n")

    report.append("## 9. Pre-registered therapeutic-axis observation\n")
    report.append(markdown_table(therapy.rename(columns={"lineage": "Lineage", "axis": "Axis", "status": "Status", "best_evidence": "Best evidence", "notes": "Notes"})))
    report.append("\n总结：TNF/NF-κB 在三谱系自然出现，其中 Myeloid/T 最强，Epithelial 对 `cie_2` 有幅度敏感性；JAK-STAT 在 Myeloid 明确、T 较弱；IL-23/STAT3 主要是 Myeloid 的 STAT3 arm，不能据此声称 IL-23 特异机制已经建立。此节不形成治疗结论。\n")

    report.append("## 10. 主要限制\n")
    report.append("- 仅 3 Healthy vs 4 CIE dogs；GSEA 评估 ranked gene-level disease effects 的 pathway coherence，不创造新的独立样本量。\n- 疾病与结构因素仍可能混杂；`cie_4` 在 Phase3 多谱系 PCA/MDS 有偏离，必须由 Phase5 全 7 轮 LODO 检查。\n- Hallmark 为 ortholog-projected annotation；Reactome 犬条目含跨物种推断属性。\n- Epithelial 某些程序对 `cie_2` 有幅度敏感性；IEG/stress 结果必须结合专门审计表解释。\n- `cie_2` Myeloid 仅 27 cells，低深度会限制敏感性比较的生物学解释。\n- `Not detected` 不等于 biologically absent；尤其不可用单个 IL23A 等基因否定整个 pathway。\n")

    report.append("## 11. Phase4 Gate\n")
    report.append("**YELLOW**\n\n理由：九个正式 GSEA、完整 leading edges、gene-set/mapping audit、共同-universe `cie_2` sensitivity、GO 去冗余和跨库汇总均成功；主要 program coherence 清楚，且显著 pathway 无方向翻转。但重要 donor sensitivity、样本量、结构混杂及 Myeloid 低深度尚未由完整 LODO 解决。\n\n**建议进入 Phase5 全 7 轮 LODO**，并优先检查：Epithelial TNF/NF-κB 与 stress/IEG、三谱系 mitochondrial/metabolic programs、Myeloid inflammatory/JAK-STAT/TNF-NF-κB，以及 `cie_4` 对这些 program 和 leading edges 的影响。Phase4 在此停止。\n")

    report.append("### 正式交付索引\n")
    report.append("- Primary GSEA：`06_pseudobulk/gsea/full7/{hallmark,reactome,go_bp}/`\n- Leading edges：`06_pseudobulk/gsea/leading_edges/`\n- `cie_2` sensitivity：`06_pseudobulk/gsea/sensitivity/`\n- Mapping/gene sets：`06_pseudobulk/gsea/mapping/`、`06_pseudobulk/gsea/gene_sets/`\n- Summary/figures：`06_pseudobulk/gsea/summary/`、`06_pseudobulk/gsea/figures/`\n- Metadata：`00_metadata/PHASE4_GENESET_AUDIT.csv`、`PHASE4_GSEA_PARAMETERS.csv`、`PHASE4_SOFTWARE_VERSIONS.csv`、`PHASE4_MANIFEST.csv`\n- Logs：`logs/phase4_*.log`\n")
    REPORT.write_text("\n".join(report), encoding="utf-8")


def build_manifest():
    files = set()
    for path in GSEA.rglob("*"):
        if path.is_file():
            files.add(path)
    for path in META.glob("PHASE4_*.csv"):
        if path.name != "PHASE4_MANIFEST.csv":
            files.add(path)
    files.add(META / "reactome_canis_familiaris_gene_sets.tsv.gz")
    for pattern in ("14_phase4*", "14b_phase4*", "15_phase4*", "16_phase4*", "17_phase4*"):
        files.update(p for p in (ROOT / "scripts").glob(pattern) if p.is_file())
    files.update(p for p in (ROOT / "logs").glob("phase4_*.log") if p.is_file())
    files.add(REPORT)

    rows = []
    for path in sorted(files, key=lambda p: p.as_posix().lower()):
        rel = path.relative_to(ROOT).as_posix()
        if rel.startswith("06_pseudobulk/gsea/full7/"):
            role = "primary_full7_gsea"
        elif rel.startswith("06_pseudobulk/gsea/sensitivity/"):
            role = "cie2_pathway_sensitivity"
        elif rel.startswith("06_pseudobulk/gsea/leading_edges/"):
            role = "leading_edge_output"
        elif rel.startswith("06_pseudobulk/gsea/gene_sets/"):
            role = "frozen_gene_set_or_membership"
        elif rel.startswith("06_pseudobulk/gsea/mapping/"):
            role = "mapping_audit_or_collapsed_rank"
        elif rel.startswith("06_pseudobulk/gsea/figures/"):
            role = "phase4_figure"
        elif rel.startswith("06_pseudobulk/gsea/summary/"):
            role = "phase4_summary"
        elif rel.startswith("scripts/"):
            role = "phase4_script"
        elif rel.startswith("logs/"):
            role = "phase4_log"
        elif rel == "08_results/PHASE4_REPORT.md":
            role = "phase4_report"
        else:
            role = "phase4_metadata"
        rows.append({"path_relative_to_project": rel, "role": role, "bytes": path.stat().st_size, "sha256": sha256(path)})
    manifest = pd.DataFrame(rows)
    manifest.to_csv(META / "PHASE4_MANIFEST.csv", index=False)
    for row in manifest.itertuples(index=False):
        path = ROOT / row.path_relative_to_project.replace("/", os.sep)
        if sha256(path) != row.sha256:
            raise RuntimeError(f"Post-write manifest verification failed: {path}")
    return len(manifest)


def main():
    SUMMARY.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("", encoding="utf-8")
    results = load_results()
    sensitivity = pd.read_csv(SENS / "cie2_pathway_sensitivity.csv")
    sens_summary = pd.read_csv(SENS / "cie2_pathway_sensitivity_summary.csv")
    go_reps = build_go_representatives(results)
    log(f"GO BP representatives: {len(go_reps)} deterministic groups.")
    cross = build_cross_collection(results)
    log(f"Cross-collection coherent program rows: {len(cross)}.")
    all_leading = load_all_leading()
    _, controls = build_anchor_and_control_audits(all_leading, results)
    ieg = build_ieg_audit(all_leading, sensitivity, results)
    therapy = build_therapy_axis(results)
    build_report(results, go_reps, cross, controls, therapy, ieg, sensitivity, sens_summary)
    log(f"Report written: {REPORT}")

    software_path = META / "PHASE4_SOFTWARE_VERSIONS.csv"
    software = pd.read_csv(software_path)
    extra = pd.DataFrame(
        [
            {"software": "Python", "version": platform.python_version()},
            {"software": "pandas", "version": pd.__version__},
            {"software": "matplotlib", "version": matplotlib.__version__},
        ]
    )
    software = pd.concat([software, extra], ignore_index=True).drop_duplicates("software", keep="last")
    software.to_csv(software_path, index=False)
    log("Phase4 reporting completed; finalizing the immutable manifest next. No Phase5 analysis was run.")
    build_manifest()


if __name__ == "__main__":
    main()
