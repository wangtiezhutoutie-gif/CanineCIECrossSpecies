#!/usr/bin/env python3
"""Finalize Phase6A report, provenance, gate, and manifest."""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import platform
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
METADATA = ROOT / "00_metadata"
SUMMARY = ROOT / "07_external_validation" / "FRE_2026" / "summary"
RAW = ROOT / "01_raw_processed" / "dog_FRE_2026"
REPORT = ROOT / "08_results" / "PHASE6_REPORT.md"
LOG = ROOT / "logs" / "phase6_final_report.log"
MANIFEST = METADATA / "PHASE6_MANIFEST.csv"


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def log(message: str) -> None:
    line = f"[{now_utc()}] {message}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(line + "\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_phase5_frozen() -> tuple[int, int]:
    rows = read_csv(METADATA / "PHASE5_MANIFEST.csv")
    failures = []
    for row in rows:
        path = ROOT / Path(row["path_relative_to_project"])
        if (
            not path.is_file()
            or path.stat().st_size != int(row["bytes"])
            or sha256(path).lower() != row["sha256"].lower()
        ):
            failures.append(row["path_relative_to_project"])
    if failures:
        raise RuntimeError(f"Phase5 frozen-file verification failed: {failures}")
    return len(rows), len(failures)


def update_methods_and_gate() -> None:
    methods = [
        {
            "method_item": "authorization",
            "frozen_value": "Option A: published summary-level external validation",
            "interpretation": "No FASTQ acquisition or independent whole-transcriptome reanalysis",
        },
        {
            "method_item": "target_universe",
            "frozen_value": "25 rows from PHASE5_PROGRAM_FORWARD_LIST.csv",
            "interpretation": "No FRE-driven target selection or deletion",
        },
        {
            "method_item": "external_pathway_evidence",
            "frozen_value": "Published Supplementary Table S4: 25 KEGG GSEA rows",
            "interpretation": "NES/padj used only as published summary evidence",
        },
        {
            "method_item": "external_gene_evidence",
            "frozen_value": "12 hypothesis-selected bile-acid genes from Main Table 1 / Supplementary Table S6",
            "interpretation": "Contextual only; not a complete gene-level DE table",
        },
        {
            "method_item": "relationship_types",
            "frozen_value": "EXACT_OR_CLOSE_PATHWAY; BIOLOGICALLY_RELATED; NO_PUBLIC_MATCH",
            "interpretation": "General immune or ribosome terms cannot substitute for TLR4, JAK-STAT, or mitochondrial translation",
        },
        {
            "method_item": "summary_support_classes",
            "frozen_value": "SUMMARY_DIRECTIONAL_SUPPORT; SUMMARY_PARTIAL_SUPPORT; SUMMARY_DISCORDANT; NOT_ASSESSABLE",
            "interpretation": "Every category has evidence_level=PUBLISHED_SUMMARY_ONLY",
        },
        {
            "method_item": "not_performed",
            "frozen_value": "DESeq2; limma; full-rank GSEA; Hallmark/Reactome/GO replication; leading-edge analysis",
            "interpretation": "Unavailable processed data prohibit these analyses",
        },
        {
            "method_item": "bulk_resolution_rule",
            "frozen_value": "Whole-duodenal-mucosa bulk cannot assign signal to lineage or cell state",
            "interpretation": "Discordance is bulk-level and cannot falsify lineage-specific discovery alone",
        },
        {
            "method_item": "human_forward_rule",
            "frozen_value": "ROBUST always YES; MODERATELY_STABLE discordant stays EXPLORATORY_OR_QUALIFIED; other moderate stays YES_WITH_QUALIFICATION",
            "interpretation": "No Phase5 target is erased by FRE summary discordance",
        },
    ]
    write_csv(
        METADATA / "PHASE6_OPTION_A_METHODS.csv",
        ["method_item", "frozen_value", "interpretation"],
        methods,
    )

    gate_rows = [
        {
            "gate_item": "Phase5 frozen-file integrity",
            "status": "PASS",
            "evidence": "389/389 manifest entries match size and SHA-256",
            "consequence": "Discovery evidence remains unchanged",
        },
        {
            "gate_item": "Repository and cohort identity",
            "status": "PASS",
            "evidence": "Zenodo 20330471; Europe PMC PMC13293934; 8 FRE + 4 Healthy independent dogs",
            "consequence": "Published external summary has reliable provenance",
        },
        {
            "gate_item": "Interim data-availability gate",
            "status": "YELLOW_ROUTE_AUTHORIZATION_REQUIRED",
            "evidence": "No whole-transcriptome matrix or complete DE table",
            "consequence": "Phase6 paused before formal validation until route selection",
        },
        {
            "gate_item": "Option A authorization",
            "status": "PASS_AUTHORIZED",
            "evidence": "User explicitly selected published summary-level external validation",
            "consequence": "Summary-only program assessment permitted; FASTQ remains prohibited",
        },
        {
            "gate_item": "25-program classification completeness",
            "status": "PASS",
            "evidence": "25/25 frozen targets classified with relationship type, evidence level, and warnings",
            "consequence": "Summary-level comparison is reproducible",
        },
        {
            "gate_item": "Resolution and subtype limitations",
            "status": "LIMITATION_RETAINED",
            "evidence": "FRE whole-mucosa bulk differs from CIE lineage-specific scRNA pseudobulk",
            "consequence": "Discordance cannot be interpreted as direct falsification",
        },
        {
            "gate_item": "FASTQ / Human access",
            "status": "NOT_PERFORMED",
            "evidence": "No FASTQ payload and no GSE266616/Human data accessed",
            "consequence": "Authorized scope respected",
        },
        {
            "gate_item": "Final Phase6 gate",
            "status": "YELLOW_SUMMARY_LEVEL_EXTERNAL_EVIDENCE_ONLY",
            "evidence": "Authentic reproducible summary assessment, but no independent whole-transcriptome reanalysis",
            "consequence": "Phase7 is scientifically justified with explicit limitations and requires separate authorization",
        },
    ]
    write_csv(
        METADATA / "PHASE6_DATA_AVAILABILITY_GATE.csv",
        ["gate_item", "status", "evidence", "consequence"],
        gate_rows,
    )
    write_csv(
        METADATA / "PHASE6_FINAL_GATE.csv",
        [
            "phase",
            "final_gate",
            "analysis_route",
            "evidence_level",
            "phase7_recommendation",
            "phase7_authorized",
            "notes",
        ],
        [
            {
                "phase": "Phase6A",
                "final_gate": "YELLOW_SUMMARY_LEVEL_EXTERNAL_EVIDENCE_ONLY",
                "analysis_route": "PUBLISHED_SUMMARY_LEVEL_EXTERNAL_SUPPORT_ASSESSMENT",
                "evidence_level": "PUBLISHED_SUMMARY_ONLY",
                "phase7_recommendation": "YES_WITH_EXPLICIT_LIMITATIONS",
                "phase7_authorized": "NO",
                "notes": "Stop before Human data access; use PHASE6_HUMAN_FORWARD_LIST.csv after authorization.",
            }
        ],
    )

    software = [
        {"component": "Python", "version": platform.python_version(), "role": "Phase6 audit, ingest, summary classification, reporting"},
        {"component": "python-docx", "version": importlib.metadata.version("python-docx"), "role": "DOCX table extraction"},
        {"component": "pypdf", "version": importlib.metadata.version("pypdf"), "role": "PDF text/page audit"},
        {"component": "platform", "version": platform.platform(), "role": "execution environment"},
        {"component": "analysis_route", "version": "PUBLISHED_SUMMARY_LEVEL_OPTION_A", "role": "No independent DE/GSEA; evidence_level=PUBLISHED_SUMMARY_ONLY"},
    ]
    write_csv(METADATA / "PHASE6_SOFTWARE_VERSIONS.csv", ["component", "version", "role"], software)


def create_report(phase5_n: int) -> None:
    matrix = read_csv(SUMMARY / "PHASE6_SUMMARY_LEVEL_VALIDATION_MATRIX.csv")
    human = read_csv(SUMMARY / "PHASE6_HUMAN_FORWARD_LIST.csv")
    contextual = read_csv(SUMMARY / "FRE_CONTEXTUAL_FINDINGS.csv")
    samples = read_csv(METADATA / "FRE2026_SAMPLE_METADATA.csv")
    counts = Counter(row["summary_support_class"] for row in matrix)
    forward_counts = Counter(row["forward_to_human"] for row in human)
    priority_counts = Counter(row["priority"] for row in human)

    class_order = [
        "SUMMARY_DIRECTIONAL_SUPPORT",
        "SUMMARY_PARTIAL_SUPPORT",
        "SUMMARY_DISCORDANT",
        "NOT_ASSESSABLE",
    ]
    class_lines = [
        f"| `{name}` | {counts[name]} | {100 * counts[name] / 25:.1f}% |"
        for name in class_order
    ]
    matrix_lines = []
    for row in matrix:
        matrix_lines.append(
            "| {lineage_origin} | {program} | {phase5_direction} | {published_direction} | "
            "{summary_support_class} | {forward_to_human} |".format(**row)
        )

    report = f"""# Phase6 Final Report：独立犬 FRE published summary-level external support

## 最终门控

**YELLOW — SUMMARY-LEVEL EXTERNAL EVIDENCE ONLY**

Phase6A 已按获授权的 Option A 正式完成。结果是可复现的 published summary-level external support assessment，而不是 independent whole-transcriptome reanalysis。25 个 Phase5 frozen programs 均已逐项评估；FASTQ、Human 数据和 Phase7 均未访问。

Phase7 建议：**YES_WITH_EXPLICIT_LIMITATIONS**。该建议不等于 Phase7 已获执行授权。

## 1. Data-availability gate 历史

Phase6 初始审计在 `YELLOW — ROUTE AUTHORIZATION REQUIRED` 暂停：Zenodo 只提供 targeted `Table 6.docx`；Europe PMC supplements 没有 whole-transcriptome sample-level matrix、完整 DE statistics 或 reusable gene rank。随后用户明确授权 Option A，禁止 FASTQ 和 Human 数据访问。

最终使用的公开证据：

- Supplementary Table S4：25 条 KEGG GSEA summary（pathway、NES、p-value、padj）。
- Main Table 1 / Supplementary Table S6：12 个预选胆汁酸基因的 published statistics 和组均值。
- 12 犬样本表：{sum(row['condition'] == 'FRE' for row in samples)} FRE + {sum(row['condition'] == 'Healthy' for row in samples)} Healthy，均为独立犬、duodenal mucosa bulk RNA-seq。

仍然不可执行：DESeq2、limma、完整 rank-based GSEA、Phase4 Hallmark/Reactome/GO replication、leading-edge concordance、全基因 KCNK16 精确检查。

## 2. Option A 方法

- Target universe 严格冻结为 `PHASE5_PROGRAM_FORWARD_LIST.csv` 的 25 行；Phase5 stability、direction 和 forward status 未修改。
- 只允许 `EXACT_OR_CLOSE_PATHWAY`、`BIOLOGICALLY_RELATED`、`NO_PUBLIC_MATCH` 三种关系。
- 不使用泛免疫通路替代 TLR4/JAK-STAT，不使用普通 ribosome 替代 mitochondrial translation，不从 volcano 图片估算基因统计量。
- 每项证据统一标记 `PUBLISHED_SUMMARY_ONLY`。
- `SUMMARY_DISCORDANT` 表示 FRE whole-mucosa bulk 与 CIE lineage-specific discovery 方向不同，不表示 discovery 被否定。

Phase5 frozen-file复核：{phase5_n}/{phase5_n} 项 SHA-256 通过。

## 3. 25-program Summary Support 结果

| Class | n | Percent |
|---|---:|---:|
{chr(10).join(class_lines)}

`CONTEXTUAL_ONLY` 不计入这 25 个 program。

### 逐项结果

| Lineage | Program | CIE direction | Published FRE direction | Summary class | Human status |
|---|---|---|---|---|---|
{chr(10).join(matrix_lines)}

完整 pathway、NES、padj、relationship type 和限制语见 `PHASE6_SUMMARY_LEVEL_VALIDATION_MATRIX.csv`。

## 4. Mitochondrial / metabolic 轴

Phase5 在 Epithelial、T、Myeloid 均显示稳定的 CIE-lower mitochondrial/oxidative/metabolic programs。Published FRE bulk 则显示：

- TCA cycle：NES 2.182，padj 0.005，FRE-higher。
- Carbon metabolism：NES 2.019，padj 0.005，FRE-higher。
- AMPK：NES 2.107，padj 0.005，FRE-higher。
- Mitophagy：NES 1.990，padj 0.005，FRE-higher。
- Peroxisome：NES 1.978，padj 0.005，FRE-higher。

对应 frozen programs 被标记为 `SUMMARY_DISCORDANT`；没有将这些结果包装为支持。最合理的解释是 **CIE lineage-specific metabolic suppression 与 FRE whole-mucosa metabolic activation/enrichment 的异质性**。候选原因包括 FRE subtype、治疗反应表型、whole-mucosa 细胞组成、分辨率、cohort/source/breed 和 pathway/rank 差异。该结果不能单独否定 Phase5 discovery，并被列为 Human 阶段 HIGH priority。

## 5. TNF/NF-kappaB 与 cytokine

- Published `NF-kappa B signaling pathway`：NES -1.886，padj 0.005，FRE-lower。
- Published `Cytokine-cytokine receptor interaction`：NES -2.372，padj 0.005，FRE-lower。

它们分别与 Phase5 的 T/Myeloid CIE-higher，以及 Epithelial moderately stable CIE-higher TNF/NF-kappaB 方向相反，因此标为 `SUMMARY_DISCORDANT`。这是 **bulk-level discordance**：不能写成 “FRE T cell/Myeloid NF-kappaB 下降”，也不能据此否定 discovery lineage signal。

## 6. TLR4、JAK-STAT、IL-23/STAT3 与 KCNK16

- Myeloid TLR4：`NOT_ASSESSABLE`；没有 TLR4-specific published pathway。
- Myeloid JAK-STAT：general cytokine component 为 discordant，但 JAK-STAT 本体 `NOT_ASSESSABLE`；不能用 cytokine-receptor term 等价替代。
- IL-23/STAT3：继续保持 `Weak / STAT3-arm only`；没有 IL-23-specific FRE summary，不能升级或称其 absent。
- Epithelial KCNK16：`NOT_ASSESSABLE`；不在 selected-gene 表中，volcano 没有可复用精确值，未从图片估算 logFC/P 值。

## 7. 其他可解释结果

- 三个 Stress/UPR 相关 lineage-program 通过 FRE-higher `Protein processing in endoplasmic reticulum` 获得 `SUMMARY_PARTIAL_SUPPORT`；该 KEGG term 只覆盖 UPR/ER-processing 部分，不能代表完整 MAPK/stress program。
- Myeloid Translation 通过 FRE-lower `Ribosome, eukaryotes` 获得唯一的 `SUMMARY_DIRECTIONAL_SUPPORT`。
- 同一个 FRE-lower ribosome term 与 Epithelial CIE-higher Ribosome/RNA biogenesis 相反，后者为 `SUMMARY_DISCORDANT`，体现 lineage-dependent interpretation。

## 8. Bile-acid contextual findings

`Bile secretion` 为 FRE-higher（NES 1.935，padj 0.005），并有 12 个 hypothesis-selected bile-acid genes 的公开结果。共 {len(contextual)} 条 pathway/gene contextual records 已保存为 `CONTEXTUAL_ONLY`。这些结果只能作为 FRE subtype-specific transcriptional context，不计入 25-program support counts，也没有回填 Phase5 forward list。

## 9. Discovery-vs-FRE 核心解释

Phase5 CIE programs 经过 LODO 后仍稳定，但 FRE bulk summary 没有统一复制其方向。最突出的是：

1. CIE 三谱系 metabolic/mitochondrial decrease vs FRE bulk metabolic/mitochondrial enrichment。
2. CIE lineage TNF/NF-kappaB increase vs FRE bulk NF-kappaB decrease。
3. CIE T/Myeloid cytokine increase vs FRE bulk cytokine-receptor decrease。

这些结果被解释为 canine disease-subtype / tissue-resolution heterogeneity，而不是技术失败或对 discovery 的直接证伪。Human Crohn 阶段应独立计算 Human effects，再比较 Dog CIE discovery、Dog FRE bulk summary 与 Human Crohn，不能设计成单向寻找“证明相似”的证据。

## 10. Human Forward List

25/25 programs 保留在 Human forward list：

- `YES`：{forward_counts['YES']}。
- `YES_WITH_QUALIFICATION`：{forward_counts['YES_WITH_QUALIFICATION']}。
- `EXPLORATORY_OR_QUALIFIED`：{forward_counts['EXPLORATORY_OR_QUALIFIED']}。
- Priority：HIGH {priority_counts['HIGH']}，MEDIUM {priority_counts['MEDIUM']}。

ROBUST program 即使 FRE summary discordant 仍进入 Human；MODERATELY_STABLE discordant programs 保留为 exploratory/qualified。正式输入为 `PHASE6_HUMAN_FORWARD_LIST.csv`，不得在 Phase7 凭印象重选。

## 11. Final Gate

**YELLOW — SUMMARY-LEVEL EXTERNAL EVIDENCE ONLY**

选择 YELLOW 而非 RED 的理由：repository、样本和 published summary 均真实可核实，25-program 分类可复现，且方向不一致本身是可信的生物学结果。选择 YELLOW 而非 GREEN 的理由：缺少 processed whole-transcriptome data，无法完成独立 DE、冻结 gene-set GSEA 和 leading-edge 验证。

建议进入 Phase7 Human Crohn，但必须携带以下限定：

- Full7 继续是 Primary Canine CIE Discovery Result。
- FRE 只提供 `PUBLISHED_SUMMARY_ONLY` 的 canine external context。
- FRE discordance 不得删除 Phase5 ROBUST programs。
- Phase7 需单独授权；本轮已在 Human 数据下载前停止。

## 12. 正式交付物

- `07_external_validation/FRE_2026/summary/PHASE6_SUMMARY_LEVEL_VALIDATION_MATRIX.csv`
- `07_external_validation/FRE_2026/summary/PHASE6_HUMAN_FORWARD_LIST.csv`
- `07_external_validation/FRE_2026/summary/FRE_CONTEXTUAL_FINDINGS.csv`
- `07_external_validation/FRE_2026/summary/DISCOVERY_VS_FRE_INTERPRETATION_MATRIX.csv`
- `07_external_validation/FRE_2026/summary/PHASE6_SUMMARY_SUPPORT_COUNTS.csv`
- `00_metadata/PHASE6_OPTION_A_METHODS.csv`
- `00_metadata/PHASE6_FINAL_GATE.csv`
- `00_metadata/PHASE6_MANIFEST.csv`

Generated UTC: {now_utc()}
"""
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(report, encoding="utf-8", newline="\n")
    log(f"Wrote final Phase6 report: {REPORT}")


def file_role(relative: str) -> tuple[str, str]:
    if relative.startswith("scripts/"):
        return "phase6_script", "phase6_reproducible_code"
    if relative.startswith("logs/"):
        return "phase6_log", "phase6_execution_log"
    if relative == "08_results/PHASE6_REPORT.md":
        return "phase6_final_report", "phase6_option_a_reporting"
    if relative.startswith("00_metadata/"):
        return "phase6_metadata", "phase6_audit_gate_and_provenance"
    if relative.startswith("07_external_validation/FRE_2026/summary/"):
        return "phase6_summary_level_result", "phase6_option_a_summary_assessment"
    if "/downloaded/" in relative:
        return "external_source_download", "phase6_minimal_public_download"
    if "/extracted/" in relative:
        return "external_source_original", "phase6_supplement_extraction"
    if "/parsed_tables/" in relative:
        return "external_source_parsed_table", "phase6_supplement_ingest"
    if "/reference/" in relative:
        return "external_source_reference", "phase6_repository_and_article_audit"
    return "phase6_output", "phase6_option_a"


def create_manifest() -> int:
    candidates: set[Path] = set()
    metadata_names = [
        "PHASE6_FRE_DATA_AUDIT.csv",
        "FRE2026_DOWNLOAD_MANIFEST.csv",
        "FRE2026_SAMPLE_METADATA.csv",
        "PHASE6_PHASE5_INPUT_HASH_AUDIT.csv",
        "PHASE6_SOFTWARE_VERSIONS.csv",
        "PHASE6_DATA_AVAILABILITY_GATE.csv",
        "PHASE6_OPTION_A_METHODS.csv",
        "PHASE6_FINAL_GATE.csv",
    ]
    candidates.update(METADATA / name for name in metadata_names)
    candidates.update(path for path in RAW.rglob("*") if path.is_file())
    candidates.update(path for path in SUMMARY.rglob("*") if path.is_file())
    candidates.update(path for path in (ROOT / "scripts").glob("*phase6*.py") if path.is_file())
    candidates.update(path for path in (ROOT / "logs").glob("phase6*.log") if path.is_file())
    candidates.add(REPORT)
    missing = sorted(str(path) for path in candidates if not path.is_file())
    if missing:
        raise RuntimeError(f"Missing manifest candidates: {missing}")
    if MANIFEST in candidates:
        candidates.remove(MANIFEST)

    rows = []
    for path in sorted(candidates, key=lambda p: str(p).lower()):
        relative = path.relative_to(ROOT).as_posix()
        role, step = file_role(relative)
        rows.append(
            {
                "path_relative_to_project": relative,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "file_role": role,
                "generation_step": step,
            }
        )
    write_csv(
        MANIFEST,
        ["path_relative_to_project", "bytes", "sha256", "file_role", "generation_step"],
        rows,
    )
    print(f"[{now_utc()}] Wrote final Phase6 manifest with {len(rows)} files", flush=True)
    return len(rows)


def verify_manifest() -> int:
    rows = read_csv(MANIFEST)
    failures = []
    for row in rows:
        path = ROOT / Path(row["path_relative_to_project"])
        if (
            not path.is_file()
            or path.stat().st_size != int(row["bytes"])
            or sha256(path).lower() != row["sha256"].lower()
        ):
            failures.append(row["path_relative_to_project"])
    if failures:
        raise RuntimeError(f"Final Phase6 manifest verification failed: {failures}")
    print(f"[{now_utc()}] Final Phase6 manifest self-check PASS: {len(rows)}/{len(rows)}", flush=True)
    return len(rows)


def main() -> int:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("", encoding="utf-8")
    phase5_n, failures = verify_phase5_frozen()
    if failures or phase5_n != 389:
        raise RuntimeError(f"Expected 389 unchanged Phase5 files; got {phase5_n} with {failures} failures")
    log("Phase5 frozen-file verification PASS: 389/389")
    update_methods_and_gate()
    create_report(phase5_n)
    log("Final gate: YELLOW - SUMMARY-LEVEL EXTERNAL EVIDENCE ONLY")
    log("Phase7 recommended with explicit limitations but not authorized; stopping before Human data access")
    log("Finalizing formal-file manifest; no further formal files will be mutated")
    manifest_n = create_manifest()
    verified_n = verify_manifest()
    if manifest_n != verified_n:
        raise RuntimeError("Manifest row count changed during verification")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        log(f"FATAL: {type(exc).__name__}: {exc}")
        raise
