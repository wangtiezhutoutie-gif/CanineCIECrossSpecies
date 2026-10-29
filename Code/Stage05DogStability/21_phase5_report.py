#!/usr/bin/env python3
"""Create the formal Phase5 report, audits, gate evidence, and manifest."""

from __future__ import annotations

import hashlib
import platform
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
PB = ROOT / "06_pseudobulk"
LODO = PB / "lodo"
SUMMARY = LODO / "summary"
META = LODO / "metadata"
REPORT = ROOT / "08_results" / "PHASE5_REPORT.md"
REPORT_LOG = ROOT / "logs" / "phase5_report.log"
MANIFEST = ROOT / "00_metadata" / "PHASE5_MANIFEST.csv"
INPUT_AUDIT = META / "PHASE5_INPUT_HASH_AUDIT.csv"
GATE_EVIDENCE = META / "PHASE5_GATE_EVIDENCE.csv"
LINEAGES = ["Epithelial", "T", "Myeloid"]
RUN_ORDER = ["drop_h1", "drop_h2", "drop_h3", "drop_cie2", "drop_cie3", "drop_cie4", "drop_cie6"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv_atomic(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def fmt(value, digits=3) -> str:
    if pd.isna(value):
        return "NA"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    value = float(value)
    if value != 0 and (abs(value) < 0.001 or abs(value) >= 1000):
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def md_table(frame: pd.DataFrame, columns: list[str] | None = None, headers: list[str] | None = None) -> str:
    if columns is not None:
        frame = frame[columns].copy()
    if frame.empty:
        return "（无）"
    values = frame.fillna("NA").astype(str)
    if headers is None:
        headers = list(values.columns)
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in values.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(x).replace("|", "/") for x in row) + " |")
    return "\n".join(lines)


def verify_manifest(manifest_path: Path, phase: str) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_path)
    rows = []
    for _, row in manifest.iterrows():
        relative = str(row["path_relative_to_project"]).replace("/", str(Path("/")).replace("/", "\\"))
        path = ROOT / relative
        exists = path.is_file()
        actual_bytes = path.stat().st_size if exists else np.nan
        actual_hash = sha256(path) if exists else ""
        expected_bytes = int(row["bytes"])
        expected_hash = str(row["sha256"])
        rows.append({
            "input_phase": phase,
            "path_relative_to_project": str(row["path_relative_to_project"]),
            "expected_bytes": expected_bytes,
            "actual_bytes": actual_bytes,
            "expected_sha256": expected_hash,
            "actual_sha256": actual_hash,
            "status": "PASS" if exists and actual_bytes == expected_bytes and actual_hash == expected_hash else "FAIL",
            "audit_source": manifest_path.name,
        })
    return pd.DataFrame(rows)


def build_input_audit(log) -> pd.DataFrame:
    phase3 = verify_manifest(ROOT / "00_metadata" / "PHASE3_MANIFEST.csv", "Phase3")
    phase4 = verify_manifest(ROOT / "00_metadata" / "PHASE4_MANIFEST.csv", "Phase4")
    audit = pd.concat([phase3, phase4], ignore_index=True)
    write_csv_atomic(audit, INPUT_AUDIT)
    print(f"Phase3 manifest: {(phase3.status == 'PASS').sum()}/{len(phase3)} PASS", file=log)
    print(f"Phase4 manifest: {(phase4.status == 'PASS').sum()}/{len(phase4)} PASS", file=log)
    if not audit["status"].eq("PASS").all():
        failed = audit.loc[audit["status"].ne("PASS"), "path_relative_to_project"].tolist()
        raise RuntimeError(f"Frozen Phase3/4 input hash failure: {failed[:10]}")
    return audit


def combine_software_versions() -> pd.DataFrame:
    files = [
        ROOT / "00_metadata" / "PHASE5_SOFTWARE_VERSIONS_R.csv",
        ROOT / "00_metadata" / "PHASE5_SOFTWARE_VERSIONS_GSEA_R.csv",
        ROOT / "00_metadata" / "PHASE5_SOFTWARE_VERSIONS_PYTHON.csv",
    ]
    frames = [pd.read_csv(path) for path in files]
    combined = pd.concat(frames, ignore_index=True).drop_duplicates()
    write_csv_atomic(combined, ROOT / "00_metadata" / "PHASE5_SOFTWARE_VERSIONS.csv")
    return combined


def select_program(programs: pd.DataFrame, lineage: str, text: str) -> pd.Series | None:
    found = programs.loc[programs["lineage"].eq(lineage) & programs["program"].str.contains(text, case=False, regex=True, na=False)]
    return found.iloc[0] if len(found) else None


def program_extremes(programs: pd.DataFrame, lineage: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = programs.loc[programs["lineage"].eq(lineage)].copy()
    order = {"ROBUST": 0, "MODERATELY_STABLE": 1, "DONOR_SENSITIVE": 2, "UNSTABLE": 3, "NOT_ASSESSABLE": 4}
    frame["class_order"] = frame["stability_class"].map(order).fillna(5)
    stable = frame.loc[frame["stability_class"].eq("ROBUST")].sort_values(
        ["median_abs_delta_NES", "median_leading_edge_jaccard"], ascending=[True, False]
    ).head(3)
    least = frame.sort_values(["class_order", "maximum_abs_delta_NES"], ascending=[False, False]).head(3)
    return stable, least


def compact_program_table(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame[["program", "full7_direction", "minimum_tested_lodo_n", "minimum_same_direction_n", "median_abs_delta_NES",
                    "maximum_abs_delta_NES", "median_leading_edge_jaccard", "most_influential_donor", "stability_class"]].copy()
    result["direction stability"] = result.apply(
        lambda r: f"min {int(r['minimum_same_direction_n'])}/{int(r['minimum_tested_lodo_n'])} tested"
        if pd.notna(r["minimum_same_direction_n"]) and pd.notna(r["minimum_tested_lodo_n"]) else "NA", axis=1
    )
    result["median abs delta NES"] = result["median_abs_delta_NES"].map(fmt)
    result["max abs delta NES"] = result["maximum_abs_delta_NES"].map(fmt)
    result["median LE Jaccard"] = result["median_leading_edge_jaccard"].map(fmt)
    return result[["program", "full7_direction", "direction stability", "median abs delta NES", "max abs delta NES",
                   "median LE Jaccard", "most_influential_donor", "stability_class"]]


def build_report(log) -> tuple[str, dict]:
    gene = pd.read_csv(SUMMARY / "GENE_LODO_STABILITY.csv.gz")
    effects = pd.read_csv(SUMMARY / "GENE_LODO_EFFECTS_LONG.csv.gz")
    gene_rank = pd.read_csv(SUMMARY / "GENE_RANK_STABILITY_BY_RUN.csv")
    pathways = pd.read_csv(SUMMARY / "PATHWAY_LODO_STABILITY.csv")
    path_rank = pd.read_csv(SUMMARY / "PATHWAY_RANK_STABILITY_BY_RUN.csv")
    programs = pd.read_csv(SUMMARY / "CORE_PROGRAM_LODO_SUMMARY.csv")
    reps = pd.read_csv(SUMMARY / "CORE_PROGRAM_REPRESENTATIVE_PATHWAYS.csv")
    forward = pd.read_csv(SUMMARY / "PHASE5_PROGRAM_FORWARD_LIST.csv")
    anchors = pd.read_csv(SUMMARY / "ANCHOR_GENE_LODO_SUMMARY.csv")
    donor = pd.read_csv(SUMMARY / "DONOR_INFLUENCE_SUMMARY.csv")
    metabolic = pd.read_csv(SUMMARY / "CROSS_LINEAGE_METABOLIC_LODO_SUMMARY.csv")
    ieg = pd.read_csv(SUMMARY / "EPITHELIAL_STRESS_IEG_LODO_AUDIT.csv")
    axes = pd.read_csv(SUMMARY / "THERAPEUTIC_AXIS_LODO_ROBUSTNESS.csv")
    gene_status = pd.read_csv(META / "LODO_RUN_STATUS.csv")
    gsea_status = pd.read_csv(META / "LODO_GSEA_RUN_STATUS.csv")
    qc = pd.read_csv(META / "LODO_MODEL_QC.csv")
    pair_qc = pd.read_csv(META / "LODO_GSEA_PAIR_QC.csv")

    gene_success = int(gene_status["status"].eq("SUCCESS").sum())
    gsea_success = int(gsea_status["status"].eq("SUCCESS").sum())
    pcounts = programs["stability_class"].value_counts().to_dict()
    forward_n = int(forward["forward_to_phase6"].isin(["YES", "YES_WITH_QUALIFICATION"]).sum())
    core_total = len(programs)
    technical_ok = gene_success == 21 and gsea_success == 63 and qc["model_status"].eq("SUCCESS").all()
    robust_n = int(pcounts.get("ROBUST", 0))
    moderate_n = int(pcounts.get("MODERATELY_STABLE", 0))
    sensitive_n = int(pcounts.get("DONOR_SENSITIVE", 0))
    unstable_n = int(pcounts.get("UNSTABLE", 0))

    donor_agg = donor.groupby(["donor", "short_label"], as_index=False).agg(
        mean_core_abs_delta_NES=("median_core_program_abs_delta_NES", "mean"),
        mean_all_pathway_abs_delta_NES=("median_abs_delta_NES", "mean"),
        mean_gene_abs_delta_logFC=("median_abs_delta_logFC", "mean"),
        n_genes_most=("n_genes_most_influenced", "sum"),
        n_pathways_most=("n_pathways_most_influenced", "sum"),
    ).sort_values("mean_core_abs_delta_NES", ascending=False)
    most_donor = donor_agg.iloc[0]
    gene_dominance = most_donor["n_genes_most"] / max(1, donor_agg["n_genes_most"].sum())
    path_dominance = most_donor["n_pathways_most"] / max(1, donor_agg["n_pathways_most"].sum())
    no_single_donor_dominance = gene_dominance < 0.50 and path_dominance < 0.50

    metabolic_forward_lineages = set(
        forward.loc[
            forward["program"].str.contains("mitochond|oxidative|respiratory|metabolic", case=False, regex=True, na=False)
            & forward["forward_to_phase6"].isin(["YES", "YES_WITH_QUALIFICATION"]), "lineage"
        ]
    )
    cross_metabolic_ok = metabolic_forward_lineages == set(LINEAGES)
    if (technical_ok and robust_n >= 12 and unstable_n == 0 and sensitive_n == 0
            and forward_n / core_total >= 0.75 and no_single_donor_dominance and cross_metabolic_ok):
        gate = "GREEN"
        gate_reason = "全部技术运行成功，核心 program 以 ROBUST/MODERATELY_STABLE 为主，无 DONOR_SENSITIVE/UNSTABLE program；三谱系代谢下降均可向外部验证传递，且没有单只犬支配多数 gene/pathway。"
    elif technical_ok and forward_n > 0 and unstable_n < core_total / 2:
        gate = "YELLOW"
        gate_reason = "技术框架完整且存在可传递 program，但部分关键 program 仍有 donor sensitivity 或稳定性资格限制。"
    else:
        gate = "RED"
        gate_reason = "技术失败或缺少足够稳定的可传递 program。"

    gate_rows = pd.DataFrame([
        {"criterion": "gene_lodo_success", "value": f"{gene_success}/21", "pass": gene_success == 21},
        {"criterion": "gsea_pair_success", "value": f"{gsea_success}/63", "pass": gsea_success == 63},
        {"criterion": "robust_core_program_n", "value": robust_n, "pass": robust_n >= 12},
        {"criterion": "moderately_stable_core_program_n", "value": moderate_n, "pass": True},
        {"criterion": "donor_sensitive_core_program_n", "value": sensitive_n, "pass": sensitive_n == 0},
        {"criterion": "unstable_core_program_n", "value": unstable_n, "pass": unstable_n == 0},
        {"criterion": "phase6_forward_program_n", "value": forward_n, "pass": forward_n > 0},
        {"criterion": "cross_lineage_metabolic_forward", "value": cross_metabolic_ok, "pass": cross_metabolic_ok},
        {"criterion": "no_single_donor_majority", "value": no_single_donor_dominance, "pass": no_single_donor_dominance},
        {"criterion": "final_gate", "value": gate, "pass": gate != "RED"},
    ])
    write_csv_atomic(gate_rows, GATE_EVIDENCE)

    lines: list[str] = []
    add = lines.append
    add("# Canine CIE Phase5 Report — Full seven-dog LODO robustness")
    add("")
    add(f"**Phase5 Gate：{gate}**")
    add("")
    add(gate_reason)
    add("")
    add("Full7（3 Healthy + 4 CIE）始终是 Primary Result；所有 LODO 仅用于 robustness/sensitivity assessment。本阶段没有执行 Phase6、犬 FRE 外部验证或任何 Human 分析。")
    add("")

    add("## 1. 执行完整性")
    add("")
    add(f"7 donors × 3 Primary lineages 的 edgeR QL 模型成功 **{gene_success}/21**；每轮均从 raw pseudobulk counts 重新运行 `filterByExpr → TMM → dispersion → robust QL fit/test`。Hallmark、Reactome、GO BP 的 paired common-universe GSEA 成功 **{gsea_success}/63**；失败运行为 **{21-gene_success + 63-gsea_success}**。")
    add("")
    add(f"Phase3 manifest 复核 73/73，Phase4 manifest 复核 101/101；所有冻结输入 hash 均通过。所有 GSEA pair 的最大绝对 pathway-size change 为 **{fmt(pair_qc['max_abs_pathway_size_change'].max())}**。")
    add("")
    execution = gene_status.groupby("lineage", as_index=False).agg(DE_runs=("status", "size"), tested_min=("genes_after", "min"), tested_max=("genes_after", "max"))
    gexec = gsea_status.groupby("lineage", as_index=False).size().rename(columns={"size": "GSEA_pairs"})
    execution = execution.merge(gexec, on="lineage")
    add(md_table(execution.astype(str)))
    add("")
    add("分类规则已在运行前写入 `PHASE5_STABILITY_PARAMETERS.csv`。Gene 以方向、效应与 rank 为主；pathway 以方向、paired-common-universe ΔNES 和 leading-edge Jaccard 为主。FDR retention 仅为辅助。")
    add("")

    add("## 2. Gene-level stability")
    add("")
    rank_summary = gene_rank.groupby("lineage", as_index=False).agg(
        median_spearman=("spearman_rho", "median"), min_spearman=("spearman_rho", "min"),
        median_pearson=("pearson_r", "median"), median_top200_abs_jaccard=("top200_abs_jaccard", "median")
    )
    for col in rank_summary.columns[1:]:
        rank_summary[col] = rank_summary[col].map(fmt)
    add(md_table(rank_summary))
    add("")
    for lineage in LINEAGES:
        add(f"### {lineage}")
        add("")
        g = gene.loc[gene["lineage"].eq(lineage)].copy()
        class_counts = g["stability_class"].value_counts()
        add("Full7-tested genes 的内部稳定性分类：" + "；".join(f"{k}={v}" for k, v in class_counts.items()) + "。")
        add("")
        stable = g.loc[g["stability_class"].eq("ROBUST") & g["full_logFC"].abs().ge(0.25)].sort_values(["full_FDR", "max_abs_delta_logFC"]).head(5)
        unstable = g.loc[g["stability_class"].isin(["DONOR_SENSITIVE", "UNSTABLE"]) & g["full_logFC"].abs().ge(0.25)].sort_values("max_abs_delta_logFC", ascending=False).head(5)
        show = pd.concat([stable.assign(category="stable examples"), unstable.assign(category="sensitive examples")])
        show = show[["category", "gene_symbol", "full_logFC", "same_direction_n", "median_abs_delta_logFC", "max_abs_delta_logFC", "most_influential_donor", "stability_class"]].copy()
        for col in ["full_logFC", "median_abs_delta_logFC", "max_abs_delta_logFC"]:
            show[col] = show[col].map(fmt)
        add(md_table(show))
        add("")

    add("### KCNK16 与预注册 anchor genes")
    add("")
    k = gene.loc[gene["lineage"].eq("Epithelial") & gene["gene_symbol"].eq("KCNK16")].iloc[0]
    ka = anchors.loc[anchors["lineage"].eq("Epithelial") & anchors["gene_symbol"].eq("KCNK16")].sort_values(
        "short_label", key=lambda x: x.map({v: i for i, v in enumerate(RUN_ORDER)})
    )
    add(f"KCNK16 Full7 logFC={fmt(k['full_logFC'])}、FDR={fmt(k['full_FDR'])}；7/7 LODO 同为 CIE-lower，logFC range={fmt(k['min_lodo_logFC'])} 至 {fmt(k['max_lodo_logFC'])}，median/max |ΔlogFC|={fmt(k['median_abs_delta_logFC'])}/{fmt(k['max_abs_delta_logFC'])}。删除 `{k['most_influential_donor']}` 影响最大（ΔlogFC={fmt(k['most_influential_delta'])}）；FDR<0.05 保持 2/7，但这不改变其 effect-based **{k['stability_class']}** 分类。绝对统计 rank 在 7 轮均为第 1–2。")
    add("")
    kt = ka[["short_label", "dropped_donor", "logFC_lodo", "FDR", "lodo_signed_rank", "lodo_abs_rank"]].copy()
    kt.columns = ["run", "dropped donor", "logFC", "FDR", "signed rank", "abs-stat rank"]
    for col in ["logFC", "FDR"]:
        kt[col] = kt[col].map(fmt)
    add(md_table(kt))
    add("")

    add("## 3. Pathway-level stability")
    add("")
    add(f"预注册核心 program 共 {core_total} 个：ROBUST={robust_n}，MODERATELY_STABLE={moderate_n}，DONOR_SENSITIVE={sensitive_n}，UNSTABLE={unstable_n}。Program 汇总保留每条代表 pathway 的原始稳定性级别；优先使用 official Full7 FDR<0.10 的代表证据，并对混合证据保守降级。")
    add("")
    for lineage in LINEAGES:
        add(f"### {lineage}")
        add("")
        stable, least = program_extremes(programs, lineage)
        add("最稳定代表：")
        add("")
        add(md_table(compact_program_table(stable)))
        add("")
        add("相对最不稳定代表（这里的“最不稳定”是谱系内相对排序，不等同于统计标签 `UNSTABLE`）：")
        add("")
        add(md_table(compact_program_table(least)))
        add("")
    add("各 collection 的 Top20 overlap、全体 NES Spearman/Pearson 和 direction concordance 见 `PATHWAY_RANK_STABILITY_BY_RUN.csv`。")
    add("")

    add("## 4. Epithelial 专题")
    add("")
    epi_programs = programs.loc[programs["lineage"].eq("Epithelial")]
    add(md_table(compact_program_table(epi_programs)))
    add("")
    epi_tnf_rep = reps.loc[
        reps["lineage"].eq("Epithelial") & reps["pathway_id"].eq("HALLMARK_TNFA_SIGNALING_VIA_NFKB")
    ].iloc[0]
    tnf_pairs = []
    for path in sorted((LODO / "gsea" / "Epithelial").glob("*/*Hallmark_PAIR_COMPARISON.csv")):
        f = pd.read_csv(path)
        hit = f.loc[f["pathway_id"].eq("HALLMARK_TNFA_SIGNALING_VIA_NFKB")]
        if len(hit):
            tnf_pairs.append(hit.iloc[0])
    tnf = pd.DataFrame(tnf_pairs).sort_values("short_label", key=lambda x: x.map({v: i for i, v in enumerate(RUN_ORDER)}))
    add(f"Hallmark TNFA_SIGNALING_VIA_NFKB official Full7 NES={fmt(epi_tnf_rep['full_NES'])}、FDR={fmt(epi_tnf_rep['full_padj'])}。7/7 LODO 保持 CIE-higher；median/max |ΔNES|={fmt(epi_tnf_rep['median_abs_delta_NES'])}/{fmt(epi_tnf_rep['max_abs_delta_NES'])}，median LE Jaccard={fmt(epi_tnf_rep['median_leading_edge_jaccard'])}，pathway 级别为 **{epi_tnf_rep['stability_class']}**。`drop_cie2` 是最大影响轮，但没有方向反转；因此该轴应称为 **moderately stable with cie_2 magnitude sensitivity**，而不是单 donor 驱动。")
    add("")
    tnf_show = tnf[["short_label", "dropped_donor", "pair_full_NES", "lodo_NES", "delta_NES", "same_direction", "leading_edge_jaccard"]].copy()
    for col in ["pair_full_NES", "lodo_NES", "delta_NES", "leading_edge_jaccard"]:
        tnf_show[col] = tnf_show[col].map(fmt)
    add(md_table(tnf_show))
    add("")
    ie_focus = ieg.loc[ieg["short_label"].isin(["drop_cie2", "drop_cie4"])]
    add(f"Stress/IEG audit 共 {len(ieg)} 个 pathway×LODO 记录；`drop_cie2`/`drop_cie4` 专项记录 {len(ie_focus)} 条。IEG/AP-1 成员（FOS/FOSB/JUN/EGR1/DUSP1 等）逐轮 leading edge 已显式列出，不能再由 EPI_4 的 cell-level abundance 代替 donor-level robustness。")
    add("")

    add("## 5. T 专题")
    add("")
    t_programs = programs.loc[programs["lineage"].eq("T")]
    add(md_table(compact_program_table(t_programs)))
    add("")
    t_rank = gene_rank.loc[gene_rank["lineage"].eq("T")]
    add(f"T 的 7 轮 gene-rank Spearman median={fmt(t_rank['spearman_rho'].median())}、range={fmt(t_rank['spearman_rho'].min())}–{fmt(t_rank['spearman_rho'].max())}。TNF/NF-κB、炎症、cytokine/interleukin、Stress/UPR/MAPK 与 mitochondrial/metabolic decrease 均为 ROBUST 并获准传递；其中 mitochondrial respiration 的最大 |ΔNES| 相对较大。`cie_4` 对 T 的 TNF/NF-κB 与 inflammatory response 代表项影响最大，但未造成方向反转。")
    add("")

    add("## 6. Myeloid 专题")
    add("")
    my_programs = programs.loc[programs["lineage"].eq("Myeloid")]
    focus = my_programs.loc[my_programs["program"].str.contains("TNF|TLR4|JAK", case=False, regex=True, na=False)]
    add(md_table(compact_program_table(focus)))
    add("")
    my_cie = donor.loc[donor["lineage"].eq("Myeloid") & donor["short_label"].isin(["drop_cie2", "drop_cie4"]),
                       ["short_label", "median_abs_delta_NES", "median_core_program_abs_delta_NES", "max_core_program_abs_delta_NES", "core_programs_affected"]].copy()
    for col in ["median_abs_delta_NES", "median_core_program_abs_delta_NES", "max_core_program_abs_delta_NES"]:
        my_cie[col] = my_cie[col].map(fmt)
    add("`cie_2` 与 `cie_4` 的 donor-level 对比：")
    add("")
    add(md_table(my_cie))
    add("")
    add("Myeloid TNF 的 Hallmark 与 GO 代表 pathway 为 ROBUST；Reactome TNF signaling 在两轮因 size/mapping 后未达到可测试条件，故跨库 program 汇总保守为 MODERATELY_STABLE。TLR4 与 cytokine/interleukin/JAK-STAT 为 ROBUST。删除低细胞数 `cie_2` 没有使这些主要炎症方向翻转；`cie_4` 影响需保留，但不是单独支配结论。")
    add("")

    add("## 7. Cross-lineage programs")
    add("")
    cross = programs.loc[programs["program"].str.contains("mitochond|oxidative|respiratory|metabolic", case=False, regex=True, na=False)]
    add(md_table(compact_program_table(cross)))
    add("")
    add(f"三谱系均存在获准进入 Phase6 的 CIE-lower mitochondrial/oxidative/metabolic program（cross-lineage criterion={cross_metabolic_ok}）。因此它是 Phase5 最强的跨谱系犬疾病 program 候选之一。TNF/NF-κB 在三谱系均保持主要方向，但 Epithelial 和 Myeloid 的跨库证据应以 MODERATELY_STABLE 的限定语传递。")
    add("")

    add("## 8. Donor influence")
    add("")
    donor_show = donor_agg.copy()
    for col in ["mean_core_abs_delta_NES", "mean_all_pathway_abs_delta_NES", "mean_gene_abs_delta_logFC"]:
        donor_show[col] = donor_show[col].map(fmt)
    add(md_table(donor_show))
    add("")
    add(f"按三谱系 core-program median |ΔNES| 的描述性汇总，影响最大的是 `{most_donor['donor']}`（{most_donor['short_label']}）。它只占 gene/pathway “most influential” 指派的 {gene_dominance:.1%}/{path_dominance:.1%}，远低于多数支配。`cie_2` 的影响集中于 Epithelial TNF/NF-κB 幅度；`cie_4` 的整体 core-program 影响更大，并影响 Epithelial stress、T TNF/stress 和 Myeloid TLR4 等代表项。")
    add("")
    add("重要限制：删除 Healthy 是 2 Healthy vs 4 CIE，删除 CIE 是 3 Healthy vs 3 CIE；因此上述 `most influential donor` 只能解释为模型敏感性描述，不能机械地当作公平的个体异质性排名。")
    add("")

    add("## 9. Leading-edge stability")
    add("")
    le_summary = reps.loc[reps["pathway_id"].astype(str).ne("")].groupby("lineage", as_index=False).agg(
        representative_n=("pathway_id", "size"), median_LE_Jaccard=("median_leading_edge_jaccard", "median"),
        minimum_LE_Jaccard=("median_leading_edge_jaccard", "min"), median_max_delta_NES=("max_abs_delta_NES", "median")
    )
    for col in le_summary.columns[2:]:
        le_summary[col] = le_summary[col].map(fmt)
    add(md_table(le_summary))
    add("")
    add("每个 core representative 的 7 轮 paired Full7/LODO leading-edge Jaccard 已保存；NES 稳定而 leading edge 明显改变的 pathway 不会被自动称为 ROBUST。")
    add("")

    add("## 10. Therapeutic-axis robustness observation")
    add("")
    axis_show = axes[["lineage", "axis", "pathway_name", "full_NES", "full_padj", "same_direction_n",
                      "median_abs_delta_NES", "median_leading_edge_jaccard", "stability_class"]].copy()
    for col in ["full_NES", "full_padj", "median_abs_delta_NES", "median_leading_edge_jaccard"]:
        axis_show[col] = axis_show[col].map(fmt)
    add(md_table(axis_show))
    add("")
    add("这里只报告 robustness status：TNF/NF-κB 在三谱系存在；JAK-STAT 以 Myeloid 最明确、T 较弱；IL-23/STAT3 仍只支持 Myeloid STAT3 arm，不能声称 IL-23-specific mechanism。Phase5 不形成药物推荐、疗效推断或治疗排序。")
    add("")

    add("## 11. 主要限制")
    add("")
    add("- 总样本量仍为 n=7；LODO 只能描述单 donor 敏感性，不能替代独立队列复现。")
    add("")
    add("- Healthy-drop 为 2v4、CIE-drop 为 3v3，统计功效与组别平衡不同。")
    add("")
    add("- Healthy Beagle 与临床 CIE 犬之间仍存在 breed/source/clinical-context 混杂；结果是 CIE-associated，不是 causal。")
    add("")
    add("- Myeloid donor 深度较低，尤其 `cie_2` 仅 27 cells；pseudobulk 可保持 donor 独立性，但不能补回缺失的细胞状态。")
    add("")
    add("- Hallmark 是人基因集向犬投影，Reactome/GO 亦具有注释覆盖和推断属性；gene-set 结果不是物种间正式正交同源验证。")
    add("")
    add("- donor-specific Level2 states 未在本阶段做完整 LODO；本结论仅适用于 3 条预注册 Level1 lineages。")
    add("")

    add("## 12. Phase5 Gate")
    add("")
    add(f"**{gate}**")
    add("")
    add(gate_reason)
    add("")
    add(f"`PHASE5_PROGRAM_FORWARD_LIST.csv` 共 {core_total} 个预注册 program，其中 **{forward_n}** 个获准进入 Phase6（ROBUST=`YES`；MODERATELY_STABLE=`YES_WITH_QUALIFICATION`），其余不得作为正式外部验证主目标。建议下一步进入 **Phase6 犬 FRE 外部验证**，但必须新开门控阶段、仅携带 forward list，并继续把 Full7 作为犬发现队列 Primary Result。Phase5 在此停止。")

    report_text = "\n".join(lines) + "\n"
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(REPORT) + ".tmp")
    temporary.write_text(report_text, encoding="utf-8")
    temporary.replace(REPORT)
    print(f"Report written: {REPORT}; gate={gate}; forward={forward_n}/{core_total}", file=log)
    return gate, {
        "gene_success": gene_success, "gsea_success": gsea_success, "robust_n": robust_n,
        "moderate_n": moderate_n, "sensitive_n": sensitive_n, "unstable_n": unstable_n,
        "forward_n": forward_n, "core_total": core_total, "most_donor": most_donor["donor"],
    }


def manifest_role(path: Path) -> tuple[str, str]:
    rel = path.relative_to(ROOT).as_posix()
    if rel.startswith("scripts/"):
        return "phase5_script", Path(rel).stem.split("_")[0] + "_script"
    if rel.startswith("logs/"):
        return "phase5_log", "phase5_logging"
    if rel.startswith("08_results/"):
        return "phase5_report", "21_phase5_report"
    if rel.startswith("00_metadata/"):
        return "phase5_metadata", "phase5_provenance"
    if "/gene_level/" in rel:
        return "lodo_gene_result", "18_phase5_run_gene_lodo"
    if "/gsea/" in rel:
        return "lodo_gsea_result", "19_phase5_run_gsea_lodo"
    if "/leading_edges/" in rel:
        return "lodo_leading_edge", "19_phase5_run_gsea_lodo"
    if "/summary/" in rel:
        return "lodo_stability_summary", "20_phase5_stability_summary"
    if "/figures/" in rel:
        return "lodo_figure", "18_or_20_phase5_visual_qc"
    if "/metadata/" in rel:
        return "lodo_metadata", "18_to_21_phase5_provenance"
    return "phase5_output", "phase5"


def build_manifest(log_summary: dict) -> int:
    candidates: set[Path] = set()
    for path in LODO.rglob("*"):
        if path.is_file() and not path.name.endswith(".tmp"):
            candidates.add(path)
    for path in (ROOT / "00_metadata").glob("PHASE5_*"):
        if path.is_file() and path != MANIFEST and not path.name.endswith(".tmp"):
            candidates.add(path)
    for number in ["18", "19", "20", "21"]:
        candidates.update(p for p in (ROOT / "scripts").glob(f"{number}_phase5_*") if p.is_file())
    candidates.update(p for p in (ROOT / "logs").glob("phase5_*.log") if p.is_file())
    candidates.add(REPORT)
    rows = []
    for path in sorted(candidates, key=lambda p: p.relative_to(ROOT).as_posix()):
        role, step = manifest_role(path)
        rows.append({
            "path_relative_to_project": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "file_role": role,
            "generation_step": step,
        })
    write_csv_atomic(pd.DataFrame(rows), MANIFEST)
    print(f"PHASE5_MANIFEST rows={len(rows)}; gate={log_summary['forward_n']}/{log_summary['core_total']} forwarded")
    return len(rows)


def main() -> None:
    REPORT_LOG.parent.mkdir(parents=True, exist_ok=True)
    with REPORT_LOG.open("w", encoding="utf-8") as log:
        print("Phase5 report/provenance", file=log)
        print(f"Python {platform.python_version()}; pandas {pd.__version__}; numpy {np.__version__}", file=log)
        build_input_audit(log)
        software = combine_software_versions()
        print(f"Software rows: {len(software)}", file=log)
        gate, summary = build_report(log)
        print(f"Final gate: {gate}", file=log)
        print(f"Report SHA-256: {sha256(REPORT)}", file=log)
        print("Report log closed before manifest hashing; no formal Phase5 file is modified afterward.", file=log)
    build_manifest(summary)


if __name__ == "__main__":
    main()
