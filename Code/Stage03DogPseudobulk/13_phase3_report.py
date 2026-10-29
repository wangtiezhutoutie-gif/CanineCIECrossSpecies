"""Create the Phase 3 summaries, report, and immutable-output manifest."""

from __future__ import annotations

import hashlib
import importlib.metadata
import logging
import math
import os
import platform
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


CODE_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
PB_ROOT = PROJECT_ROOT / "06_pseudobulk"
METADATA_DIR = PB_ROOT / "metadata"
DE_DIR = PB_ROOT / "de_results"
SENSITIVITY_DIR = PB_ROOT / "sensitivity"
RANKED_DIR = PB_ROOT / "ranked_lists"
FIGURE_DIR = PB_ROOT / "figures"
RESULTS_DIR = PROJECT_ROOT / "08_results"
LOG_PATH = PROJECT_ROOT / "logs" / "phase3_report.log"
REPORT_PATH = RESULTS_DIR / "PHASE3_REPORT.md"
MANIFEST_PATH = PROJECT_ROOT / "00_metadata" / "PHASE3_MANIFEST.csv"
ANNOTATED_H5AD = PROJECT_ROOT / "05_annotations" / "dog_GSE254005_annotated.h5ad"
GENE_MAPPING = PROJECT_ROOT / "00_metadata" / "dog_gene_mapping_master.csv"
PHASE2_MANIFEST = PROJECT_ROOT / "00_metadata" / "PHASE2_MANIFEST.csv"
PHASE2_REPORT = RESULTS_DIR / "PHASE2_REPORT.md"
LINEAGES = ("Epithelial", "T", "Myeloid")
SEED = 20240824
MEANINGFUL_EFFECT = 0.25


def configure_logger() -> logging.Logger:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("phase3_report")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    file_handler = logging.FileHandler(LOG_PATH, mode="w", encoding="utf-8")
    stream_handler = logging.StreamHandler(sys.stdout)
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)


def atomic_text(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().eq("true")


def markdown_table(frame: pd.DataFrame, columns: list[str] | None = None) -> str:
    view = frame if columns is None else frame[columns]
    if view.empty:
        return "（无）"
    escaped_headers = [str(column).replace("|", "\\|") for column in view.columns]
    header = "| " + " | ".join(escaped_headers) + " |"
    divider = "| " + " | ".join(["---"] * len(view.columns)) + " |"
    rows = []
    for row in view.itertuples(index=False, name=None):
        values = [str(value).replace("|", "\\|") for value in row]
        rows.append("| " + " | ".join(values) + " |")
    return "\n".join([header, divider, *rows])


def fmt_float(value: float, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}f}"


def validate_phase2(logger: logging.Logger) -> dict[str, object]:
    for path in (ANNOTATED_H5AD, GENE_MAPPING, PHASE2_MANIFEST, PHASE2_REPORT):
        if not path.exists():
            raise FileNotFoundError(path)

    manifest = pd.read_csv(PHASE2_MANIFEST)
    missing: list[str] = []
    mismatched: list[str] = []
    for row in manifest.itertuples(index=False):
        path = PROJECT_ROOT / str(row.path_relative_to_project)
        if not path.exists():
            missing.append(str(path))
        elif sha256(path).lower() != str(row.sha256).lower():
            mismatched.append(str(path))
    if missing or mismatched:
        raise RuntimeError(f"Phase2 manifest failed: missing={missing}; mismatched={mismatched}")
    logger.info("Phase2 manifest verified: %d rows, no missing/hash mismatch", len(manifest))

    adata = ad.read_h5ad(ANNOTATED_H5AD)
    if adata.shape != (38_186, 20_257):
        raise RuntimeError(f"Unexpected h5ad shape: {adata.shape}")
    required_obs = {"donor_id", "sample_id", "condition", "level1", "level2"}
    if required_obs - set(adata.obs.columns):
        raise RuntimeError(f"Missing obs columns: {required_obs - set(adata.obs.columns)}")
    if "counts" not in adata.layers:
        raise RuntimeError("layers['counts'] is absent")
    counts = adata.layers["counts"]
    if not sparse.issparse(counts):
        raise RuntimeError("layers['counts'] is not sparse")
    if counts.data.size and (np.min(counts.data) < 0 or not np.allclose(counts.data, np.rint(counts.data))):
        raise RuntimeError("layers['counts'] is not non-negative integer-valued")
    donors = adata.obs[["donor_id", "condition"]].drop_duplicates()
    if donors["donor_id"].nunique() != 7:
        raise RuntimeError("Expected 7 donors")
    condition_counts = donors.groupby("condition", observed=True)["donor_id"].nunique().to_dict()
    if condition_counts != {"CIE": 4, "Healthy": 3}:
        raise RuntimeError(f"Unexpected donor condition counts: {condition_counts}")
    if "feature_id" not in adata.var or "gene_symbol" not in adata.var:
        raise RuntimeError("feature_id/original gene symbol field missing")
    if not adata.var_names.is_unique or adata.var["feature_id"].duplicated().any():
        raise RuntimeError("Feature identity is not unique")
    mapping = pd.read_csv(GENE_MAPPING)
    required_mapping = {"feature_id", "original_gene_symbol", "unique_var_name"}
    if required_mapping - set(mapping.columns):
        raise RuntimeError("Gene mapping is missing identity columns")
    if len(mapping) != 20_257 or mapping["feature_id"].duplicated().any() or mapping["unique_var_name"].duplicated().any():
        raise RuntimeError("Gene mapping identity validation failed")
    if mapping["original_gene_symbol"].duplicated(keep=False).sum() != 34:
        raise RuntimeError("Expected 34 feature rows affected by duplicated symbols")
    logger.info("Annotated input verified: shape=%s, donors=%s, sparse integer counts", adata.shape, condition_counts)
    del counts, adata
    return {"phase2_manifest_rows": len(manifest), "condition_counts": condition_counts}


def create_runtime_provenance() -> None:
    candidates = [
        (PROJECT_ROOT / "work" / "R-4.4.3-win.exe", "official CRAN R 4.4.3 Windows installer archive"),
        (PROJECT_ROOT / "work" / "innoextract670.zip", "innoextract archive used to unpack installer without system installation"),
    ]
    candidates.extend(
        (path, "downloaded R package binary archive")
        for path in sorted((PROJECT_ROOT / "work" / "phase3_r_packages").glob("*.zip"))
    )
    rows = []
    for path, role in candidates:
        if path.exists():
            rows.append(
                {
                    "file": path.relative_to(PROJECT_ROOT).as_posix(),
                    "role": role,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    atomic_csv(pd.DataFrame(rows), PROJECT_ROOT / "00_metadata" / "PHASE3_RUNTIME_PROVENANCE.csv")


def assemble_results(logger: logging.Logger) -> dict[str, object]:
    sample_meta = pd.read_csv(METADATA_DIR / "pseudobulk_sample_metadata.csv")
    filters = pd.read_csv(METADATA_DIR / "gene_filter_summary.csv")
    normalization = pd.read_csv(METADATA_DIR / "normalization_factors.csv")
    qc = pd.read_csv(METADATA_DIR / "pseudobulk_sample_qc.csv")
    sensitivity = pd.read_csv(SENSITIVITY_DIR / "cie2_sensitivity_summary.csv")
    gene_sensitivity = pd.read_csv(SENSITIVITY_DIR / "cie2_sensitivity_gene_level.csv")
    level2 = pd.read_csv(METADATA_DIR / "level2_pseudobulk_eligibility.csv")
    secondary = pd.read_csv(METADATA_DIR / "secondary_level1_pseudobulk_eligibility.csv")
    level2_eligible = bool_series(level2["eligible"])
    secondary_eligible = bool_series(secondary["eligible"])
    qc_flag = bool_series(qc["descriptive_outlier_flag"])

    summaries: list[dict[str, object]] = []
    top_frames: list[pd.DataFrame] = []
    result_tables: dict[str, pd.DataFrame] = {}
    for lineage in LINEAGES:
        result = pd.read_csv(DE_DIR / f"{lineage}_CIE_vs_Healthy_edger.csv")
        result_tables[lineage] = result
        meta = sample_meta[sample_meta["lineage"].eq(lineage)]
        sens = sensitivity[sensitivity["lineage"].eq(lineage)].iloc[0]
        fdr05 = result["FDR"].lt(0.05)
        fdr10 = result["FDR"].lt(0.10)
        notes = {
            "Epithelial": "Partial cie_2 sensitivity; interpret gene-level extremes cautiously",
            "T": "No FDR<0.10 genes; ranked effects and cie_2 sensitivity are stable",
            "Myeloid": "cie_2 has only 27 cells/low library; Full7 robust QL retained; 3 genes at FDR<0.10",
        }[lineage]
        summaries.append(
            {
                "lineage": lineage,
                "n_cells": int(meta["n_cells"].sum()),
                "healthy_donors": int(meta.loc[meta["condition"].eq("Healthy"), "donor_id"].nunique()),
                "cie_donors": int(meta.loc[meta["condition"].eq("CIE"), "donor_id"].nunique()),
                "tested_genes": len(result),
                "FDR_lt_0_05": int(fdr05.sum()),
                "FDR_lt_0_10": int(fdr10.sum()),
                "up_FDR_lt_0_05": int((fdr05 & result["logFC"].gt(0)).sum()),
                "down_FDR_lt_0_05": int((fdr05 & result["logFC"].lt(0)).sum()),
                "median_abs_logFC": float(result["logFC"].abs().median()),
                "max_abs_logFC": float(result["logFC"].abs().max()),
                "cie2_logFC_spearman": float(sens["logFC_spearman"]),
                "cie2_direction_concordance": float(sens["direction_concordance"]),
                "notes": notes,
            }
        )

        scored = result.copy()
        scored["effect_stat_score"] = scored["logFC"].abs() * np.sqrt(scored["F"].clip(lower=0))
        for direction_label, mask in (
            ("Positive_CIE_higher", scored["logFC"].gt(0)),
            ("Negative_CIE_lower", scored["logFC"].lt(0)),
        ):
            selected = scored.loc[mask].nlargest(10, "effect_stat_score").copy()
            selected.insert(0, "effect_class", direction_label)
            selected.insert(0, "rank_within_class", range(1, len(selected) + 1))
            top_frames.append(
                selected[
                    [
                        "lineage",
                        "effect_class",
                        "rank_within_class",
                        "feature_id",
                        "gene_symbol",
                        "unique_var_name",
                        "logFC",
                        "F",
                        "PValue",
                        "FDR",
                        "effect_stat_score",
                    ]
                ]
            )

    summary = pd.DataFrame(summaries)
    top_effects = pd.concat(top_frames, ignore_index=True)
    atomic_csv(summary, DE_DIR / "PRIMARY_LINEAGE_DE_SUMMARY.csv")
    atomic_csv(top_effects, DE_DIR / "PRIMARY_TOP_CIE_ASSOCIATED_EFFECTS.csv")
    logger.info("Primary summaries and top-effect table written")

    return {
        "sample_meta": sample_meta,
        "filters": filters,
        "normalization": normalization,
        "qc": qc,
        "qc_flag": qc_flag,
        "sensitivity": sensitivity,
        "gene_sensitivity": gene_sensitivity,
        "level2": level2,
        "level2_eligible": level2_eligible,
        "secondary": secondary,
        "secondary_eligible": secondary_eligible,
        "summary": summary,
        "top_effects": top_effects,
        "result_tables": result_tables,
    }


def top_effect_markdown(top_effects: pd.DataFrame, lineage: str, effect_class: str, n: int = 5) -> str:
    view = top_effects[
        top_effects["lineage"].eq(lineage) & top_effects["effect_class"].eq(effect_class)
    ].head(n).copy()
    view["gene"] = view["gene_symbol"].fillna("").astype(str)
    view.loc[view["gene"].eq(""), "gene"] = view.loc[view["gene"].eq(""), "feature_id"]
    view["logFC"] = view["logFC"].map(lambda x: fmt_float(x, 3))
    view["FDR"] = view["FDR"].map(lambda x: f"{x:.3g}")
    view["effect+stat score"] = view["effect_stat_score"].map(lambda x: fmt_float(x, 2))
    return markdown_table(view.rename(columns={"gene": "gene symbol"}), ["gene symbol", "feature_id", "logFC", "FDR", "effect+stat score"])


def build_report(validation: dict[str, object], data: dict[str, object], logger: logging.Logger) -> None:
    sample_meta: pd.DataFrame = data["sample_meta"]
    qc: pd.DataFrame = data["qc"]
    sensitivity: pd.DataFrame = data["sensitivity"]
    gene_sensitivity: pd.DataFrame = data["gene_sensitivity"]
    level2: pd.DataFrame = data["level2"]
    secondary: pd.DataFrame = data["secondary"]
    summary: pd.DataFrame = data["summary"]
    top_effects: pd.DataFrame = data["top_effects"]
    qc_flag: pd.Series = data["qc_flag"]
    level2_eligible: pd.Series = data["level2_eligible"]
    secondary_eligible: pd.Series = data["secondary_eligible"]

    primary_display = summary.copy()
    for column in ("median_abs_logFC", "max_abs_logFC", "cie2_logFC_spearman", "cie2_direction_concordance"):
        primary_display[column] = primary_display[column].map(lambda x: fmt_float(x, 3))

    pb_rows = []
    for lineage in LINEAGES:
        meta = sample_meta[sample_meta["lineage"].eq(lineage)]
        filter_row = data["filters"][data["filters"]["lineage"].eq(lineage)].iloc[0]
        pb_rows.append(
            {
                "lineage": lineage,
                "cells": f"{int(meta['n_cells'].sum()):,}",
                "donors (H/CIE)": "3/4",
                "raw library range": f"{int(meta['library_size'].min()):,}–{int(meta['library_size'].max()):,}",
                "tested genes": f"{int(filter_row['genes_after']):,}",
            }
        )
    pb_display = pd.DataFrame(pb_rows)

    sens_display = sensitivity.copy()
    sens_display = sens_display[
        [
            "lineage",
            "common_tested_genes",
            "logFC_spearman",
            "logFC_pearson",
            "direction_concordance",
            "meaningful_direction_concordance",
            "top_effect_overlap_n",
            "large_effect_sign_flips",
            "median_abs_delta_logFC",
        ]
    ]
    for column in ("logFC_spearman", "logFC_pearson", "direction_concordance", "meaningful_direction_concordance", "median_abs_delta_logFC"):
        sens_display[column] = sens_display[column].map(lambda x: fmt_float(x, 3))
    sens_display["top_effect_overlap_n"] = sens_display["top_effect_overlap_n"].map(lambda x: f"{int(x)}/100")
    sens_display = sens_display.rename(
        columns={
            "common_tested_genes": "common genes",
            "logFC_spearman": "Spearman rho",
            "logFC_pearson": "Pearson r",
            "direction_concordance": "all-gene direction",
            "meaningful_direction_concordance": "meaningful direction",
            "top_effect_overlap_n": "top-effect overlap",
            "large_effect_sign_flips": "large-effect flips",
            "median_abs_delta_logFC": "median |ΔlogFC|",
        }
    )

    qc_details = []
    reason_columns = {
        "low_library_flag": "robust-low library",
        "low_library_ratio_flag": "library <25% lineage median",
        "low_correlation_flag": "low sample correlation",
        "pca_distance_flag": "PCA-distance",
        "mds_distance_flag": "MDS-distance",
        "top_de_contribution_flag": "top-100 DE within-group deviation",
    }
    for row in qc.loc[qc_flag].itertuples(index=False):
        reasons = [label for column, label in reason_columns.items() if str(getattr(row, column)).lower() == "true"]
        qc_details.append(
            {
                "lineage": row.lineage,
                "donor": row.donor_id,
                "cells": int(row.n_cells),
                "raw library": f"{int(row.raw_library_size):,}",
                "min correlation": fmt_float(row.min_sample_correlation, 3),
                "descriptive flags": ", ".join(reasons),
            }
        )
    qc_display = pd.DataFrame(qc_details)

    eligible_level2_count = int(level2_eligible.sum())
    excluded_level2 = level2.loc[~level2_eligible]
    donor_confounded = excluded_level2["reason"].str.contains("DONOR_CONFOUNDED", na=False).sum()
    healthy_coverage = excluded_level2["reason"].str.contains("Healthy donors >=20", regex=False, na=False).sum()
    cie_coverage = excluded_level2["reason"].str.contains("CIE donors >=20", regex=False, na=False).sum()
    eligible_level2_names = ", ".join(level2.loc[level2_eligible, "level2"].astype(str))
    eligible_secondary_names = ", ".join(secondary.loc[secondary_eligible, "level1"].astype(str)) or "无"

    large_flips = gene_sensitivity[bool_series(gene_sensitivity["large_effect_sign_flip"])].copy()
    flip_counts = large_flips.groupby("lineage").size().reindex(LINEAGES, fill_value=0)
    notable_flips = large_flips.assign(abs_delta=large_flips["delta_logFC_no_cie2_minus_full7"].abs()).nlargest(8, "abs_delta")
    if notable_flips.empty:
        flip_markdown = "（无）"
    else:
        notable_flips["gene"] = notable_flips["gene_symbol"].fillna("").astype(str)
        notable_flips.loc[notable_flips["gene"].eq(""), "gene"] = notable_flips.loc[notable_flips["gene"].eq(""), "feature_id"]
        notable_flips["Full7 logFC"] = notable_flips["logFC_full7"].map(lambda x: fmt_float(x, 3))
        notable_flips["No-cie2 logFC"] = notable_flips["logFC_no_cie2"].map(lambda x: fmt_float(x, 3))
        notable_flips["|ΔlogFC|"] = notable_flips["abs_delta"].map(lambda x: fmt_float(x, 3))
        flip_markdown = markdown_table(
            notable_flips.rename(columns={"lineage": "谱系"}),
            ["谱系", "gene", "feature_id", "Full7 logFC", "No-cie2 logFC", "|ΔlogFC|"],
        )

    lines = [
        "# Canine CIE 项目 Phase3 报告",
        "",
        "**Phase3 gate：YELLOW**",
        "",
        "三条预注册 Primary lineages 均已成功建立 7-donor raw-count pseudobulk，condition 设计可估计，edgeR quasi-likelihood 模型与完整 ranked statistics 均已生成。可以进入 Phase4 GSEA，但 Epithelial 对 `cie_2` 的部分基因级极端效应较敏感，T/Myeloid 也存在需要后续 LODO 验证的 donor 结构；因此不授予 GREEN。Phase3 到此停止，本阶段未运行任何 GSEA、完整 LODO、Human、FRE、ortholog 或治疗靶点分析。",
        "",
        "## 1. 输入完整性",
        "",
        f"- `dog_GSE254005_annotated.h5ad` 回读为 38,186 cells × 20,257 features；7 dogs（Healthy=3、CIE=4）。",
        "- `layers['counts']` 为稀疏、非负、整数 raw counts；pseudobulk 未使用 `X` 的 normalized values、PCA、scaled 或 imputed data。",
        "- `donor_id`、`sample_id`、`condition`、Level1、Level2、`feature_id` 均存在。original symbol 使用 `var['gene_symbol']`；稳定 unique display identity 由 `var_names`/mapping 的 `unique_var_name` 保存。",
        "- 20,257 个 feature rows 以唯一 `feature_id` 为统计 row identity。17 个重复 symbols 影响的 34 rows 保持独立，未按 symbol 静默合并。",
        f"- Phase2 manifest 共 {validation['phase2_manifest_rows']} 项全部重新核验，无缺失或 hash mismatch；Phase0–2 正式文件未被覆盖。",
        "",
        "## 2. Pseudobulk 结构",
        "",
        markdown_table(pb_display),
        "",
        "每一谱系单独执行 `filterByExpr` 和 TMM。模型固定为 `~ condition`，Healthy 为 reference，因此 **logFC > 0 表示 CIE higher，logFC < 0 表示 CIE lower**。每只犬是一个生物学重复；没有使用 cell-level P 值。",
        "",
        "## 3. Primary CIE-associated differential expression",
        "",
        markdown_table(
            primary_display,
            [
                "lineage", "n_cells", "healthy_donors", "cie_donors", "tested_genes", "FDR_lt_0_05",
                "FDR_lt_0_10", "up_FDR_lt_0_05", "down_FDR_lt_0_05", "median_abs_logFC",
                "max_abs_logFC", "notes",
            ],
        ),
        "",
        "FDR 数量仅作描述，不代表某谱系更重要。完整测试表和 `sign(logFC) × sqrt(F)` ranked lists 均已保存；后者仅为 Phase4 准备，本阶段未做 enrichment。Top effects 以下按 `abs(logFC) × sqrt(F)` 联合 effect size 与统计证据排序，而不是仅按 P 值。",
    ]

    for lineage in LINEAGES:
        lines.extend(
            [
                "",
                f"### {lineage}",
                "",
                "Top positive CIE-associated effects：",
                "",
                top_effect_markdown(top_effects, lineage, "Positive_CIE_higher"),
                "",
                "Top negative CIE-associated effects：",
                "",
                top_effect_markdown(top_effects, lineage, "Negative_CIE_lower"),
            ]
        )

    lines.extend(
        [
            "",
            "IEG/stress genes（如 FOS/FOSB/JUN/EGR1/DUSP1）若出现在结果中，只能称为 CIE-associated effect；在 pathway coherence、LODO 与外部验证前，不能据此断言疾病机制。",
            "",
            "## 4. Pseudobulk QC",
            "",
            markdown_table(qc_display),
            "",
            "- 全部 7 donors 均保留在 Full7 primary models；上述阈值是描述性 QC 标记，不是自动排除规则。",
            "- `cie_2`：Epithelial 在 PCA/MDS 和相关性上偏离；T 只有 318 cells 且 library 为该谱系中最低；Myeloid 仅 27 cells、raw library 64,176（过滤后 61,680），并在 PCA/MDS 上偏离。",
            "- `cie_4`：三条谱系均有不同程度 PCA/MDS 偏离；T 的 top-100 DE within-group deviation 也被标记，提示部分 top effects 可能有 donor contribution。",
            "- Epithelial `cie_3` 仅触发相对严格的 robust-min-correlation flag；其最低相关系数仍为 0.863，未见样本标识或计数错误。",
            f"- 各谱系 logCPM 最低 pairwise correlation：Epithelial {qc.loc[qc['lineage'].eq('Epithelial'), 'min_sample_correlation'].min():.3f}；T {qc.loc[qc['lineage'].eq('T'), 'min_sample_correlation'].min():.3f}；Myeloid {qc.loc[qc['lineage'].eq('Myeloid'), 'min_sample_correlation'].min():.3f}。",
            "- BCV、QL dispersion、MA、volcano、PCA、MDS、sample-correlation 及 library/detected-gene plots 均已输出。edgeR QL 本流程没有正式 Cook's distance；使用 `estimateDisp(..., robust=TRUE)` 与 `glmQLFit(..., robust=TRUE)` 抗高变异 genes，并另存 PCA/MDS、相关性和 top-DE contribution 描述指标。",
            "- 没有发现需要停止模型或更正样本身份的技术错误，但 n=7 下任何 outlier 判断均不稳定，留待 Phase5 完整 LODO 定量。",
            "",
            "## 5. `cie_2` targeted sensitivity",
            "",
            markdown_table(sens_display),
            "",
            f"这里的 |logFC|≥{MEANINGFUL_EFFECT:.2f} 是预先固定的项目内部敏感性描述阈值，不是显著性阈值。all-gene direction 会被接近 0 的符号改变拉低，因此需与 meaningful-direction 一起读。",
            "",
            "- **Epithelial：部分敏感。** Spearman ρ=0.879，median |ΔlogFC|=0.136，top-100 overlap=78%；双模型均达到 |logFC|≥0.25 的 genes 中方向一致率 98.92%，但仍有 46 个 large-effect sign flips，说明若干极端上/下调效应受 `cie_2` 明显影响。",
            "- **T：稳定。** Spearman ρ=0.977，median |ΔlogFC|=0.045，top-100 overlap=92%；meaningful-direction concordance=100%，large-effect sign flips=0。",
            "- **Myeloid：总体稳定但需低深度警告。** Spearman ρ=0.963，median |ΔlogFC|=0.069，top-100 overlap=86%；meaningful-direction concordance=99.95%，large-effect sign flip=1。由于 `cie_2` Myeloid 只有 27 cells，Full7 和 No-cie2 的 `filterByExpr` tested sets 不同，比较限定于 4,672 个共同 tested genes。",
            "",
            "按 |ΔlogFC| 排序的最明显 large-effect sign flips：",
            "",
            flip_markdown,
            "",
            "Full7 始终是正式 primary result；No-cie2 只作为预设敏感性。未删除其他 donor，也未执行 7 轮 LODO。",
            "",
            "## 6. Level2 eligibility",
            "",
            f"- 31 个 Level2 states 中 **{eligible_level2_count} 个**通过保守 exploratory pseudobulk eligibility screen：{eligible_level2_names}。本阶段只建立 eligibility，不运行这些 state 的 DE，因此没有把 exploratory state 与 Primary Level1 结果等同解释。",
            f"- **{len(excluded_level2)} 个**未通过。非互斥原因计数：Healthy 端 ≥20 cells 的 donor 少于 3（{healthy_coverage}）；CIE 端少于 3（{cie_coverage}）；包含 Phase2 已识别的 donor-confounded annotation component（{donor_confounded}）。",
            "- `Enterocyte 1 (IEG-high)` 明确标记 `DONOR_CONFOUNDED`（含 EPI_4），不得包装为稳定的一般 CIE state。",
            f"- Secondary Level1 中仅 **{eligible_secondary_names}**满足 7-donor/细胞阈值；B、Mast、Stromal 不满足。Secondary DE 非本次核心交付，未强行运行。",
            "",
            "## 7. 主要限制",
            "",
            "- 真正样本量仅 3 Healthy vs 4 CIE；统计功效与 outlier 识别均有限。",
            "- Healthy Beagle 与 clinical CIE 的来源、品种、去势/临床变量结构性混杂；结果只能称为 CIE-associated expression effect，不能称为 causal effect。",
            "- donor 属于唯一 condition，故不能同时拟合 `~ donor + condition`；breed/age/medication 等也不能在 n=7 下可靠校正。",
            "- fine states 广泛存在 donor dominance，Level2 结论边界必须保持 exploratory/descriptive。",
            "- `cie_2` 具有 Phase2 epithelial mt anomaly，且 T/Myeloid pseudobulk 深度偏低；Epithelial gene-level extremes 的敏感性最明显。",
            "- FDR 不显著不等于 effect 无价值；同样，显著但由单 donor 贡献的 gene 也不能自动成为核心机制。",
            "",
            "## 8. Phase3 gate",
            "",
            "**YELLOW — 建议有条件进入 Phase4 GSEA。**",
            "",
            "理由：3 条 Primary lineages 的 Full7 donor-level pseudobulk、QL models、QC、完整 DE tables 与 ranked statistics 均可靠生成，设计可估计，没有关键样本身份/计数失败；但 Epithelial targeted sensitivity 存在 46 个预定义 large-effect sign flips，`cie_2` 与 `cie_4` 在多项描述性 QC 中偏离，且多数 Level2 states 受 coverage/donor-confounding 限制。Phase4 应优先检验 pathway-level coherence，并在 Phase5 用完整 LODO 区分稳健 program 与 donor-driven extremes。",
            "",
            "### 本阶段明确未执行",
            "",
            "Hallmark/Reactome/GO GSEA、完整 7 轮 LODO、FRE validation、Human 数据、ortholog mapping、singIST、SATURN、CONCORD、therapeutic prioritization。",
            "",
            "### 核心交付定位",
            "",
            "- `06_pseudobulk/matrices/`：3 个 raw-count pseudobulk matrices。",
            "- `06_pseudobulk/metadata/`：sample metadata、filter/TMM/QC、Secondary/Level2 eligibility 与 gene identity。",
            "- `06_pseudobulk/de_results/`：完整 Full7 DE、FDR 辅助表、Primary summary 与 top-effect table。",
            "- `06_pseudobulk/sensitivity/`：No-cie2 完整结果、gene-level comparison 与 summary。",
            "- `06_pseudobulk/ranked_lists/`：3 个 Phase4-ready signed-statistic ranked lists。",
            "- `06_pseudobulk/figures/`：各谱系 QC/DE/sensitivity figures。",
            "- `00_metadata/PHASE3_MANIFEST.csv`：Phase3 交付物与核心输入 SHA-256 清单。",
            "",
        ]
    )
    atomic_text("\n".join(lines), REPORT_PATH)
    logger.info("Wrote report: %s", REPORT_PATH)


def build_manifest() -> None:
    explicit_inputs = [ANNOTATED_H5AD, GENE_MAPPING, PHASE2_MANIFEST, PHASE2_REPORT]
    phase3_files: set[Path] = set()
    for directory in (PB_ROOT,):
        phase3_files.update(path for path in directory.rglob("*") if path.is_file())
    for pattern in ("phase3_*.log",):
        phase3_files.update(path for path in (PROJECT_ROOT / "logs").glob(pattern) if path.is_file())
    phase3_files.update(
        path
        for path in (
            CODE_ROOT / "config" / "phase3_config.py",
            CODE_ROOT / "Stage03DogPseudobulk" / "10_phase3_build_pseudobulk.py",
            CODE_ROOT / "Stage03DogPseudobulk" / "11_phase3_edger.R",
            CODE_ROOT / "Stage03DogPseudobulk" / "12_phase3_cie2_sensitivity.R",
            CODE_ROOT / "Stage03DogPseudobulk" / "13_phase3_report.py",
            PROJECT_ROOT / "Environments" / "phase3_environment.txt",
            PROJECT_ROOT / "00_metadata" / "PHASE3_SOFTWARE_VERSIONS.csv",
            PROJECT_ROOT / "00_metadata" / "PHASE3_RUNTIME_PROVENANCE.csv",
            REPORT_PATH,
        )
        if path.exists()
    )
    rows = []
    for path in sorted(set(explicit_inputs) | phase3_files):
        rows.append(
            {
                "path_relative_to_project": path.relative_to(PROJECT_ROOT).as_posix(),
                "role": "immutable_phase2_input" if path in explicit_inputs else "phase3_deliverable",
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    atomic_csv(pd.DataFrame(rows), MANIFEST_PATH)


def main() -> None:
    start = time.time()
    logger = configure_logger()
    logger.info("Phase3 report start")
    logger.info("Project root: %s", PROJECT_ROOT)
    logger.info("Inputs: %s; %s; %s", ANNOTATED_H5AD, GENE_MAPPING, PHASE2_MANIFEST)
    logger.info("Random seed: %d (no stochastic report step)", SEED)
    logger.info(
        "Python=%s; pandas=%s; anndata=%s; scipy=%s; platform=%s",
        platform.python_version(),
        importlib.metadata.version("pandas"),
        importlib.metadata.version("anndata"),
        importlib.metadata.version("scipy"),
        platform.platform(),
    )
    validation = validate_phase2(logger)
    create_runtime_provenance()
    data = assemble_results(logger)
    build_report(validation, data, logger)
    logger.info("Gate assigned: YELLOW")
    logger.info("Phase3 report completed in %.2f seconds", time.time() - start)
    for handler in list(logger.handlers):
        handler.flush()
        handler.close()
        logger.removeHandler(handler)
    logging.shutdown()
    build_manifest()


if __name__ == "__main__":
    main()
