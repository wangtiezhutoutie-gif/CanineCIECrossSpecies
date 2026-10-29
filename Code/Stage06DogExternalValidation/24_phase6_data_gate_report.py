#!/usr/bin/env python3
"""Create the Phase6 interim data-gate report and formal manifest."""

from __future__ import annotations

import csv
import hashlib
import platform
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
METADATA = ROOT / "00_metadata"
RAW = ROOT / "01_raw_processed" / "dog_FRE_2026"
PARSED = RAW / "reference" / "parsed_tables"
REPORT = ROOT / "08_results" / "PHASE6_REPORT.md"
LOG = ROOT / "logs" / "phase6_data_gate_report.log"
MANIFEST = METADATA / "PHASE6_MANIFEST.csv"


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def log(message: str) -> None:
    line = f"[{now_utc()}] {message}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(line + "\n")


def create_gate_table() -> None:
    rows = [
        {
            "gate_item": "Phase5 frozen-file integrity",
            "status": "PASS",
            "evidence": "389/389 manifest entries match size and SHA-256",
            "consequence": "25-target discovery input can be trusted",
        },
        {
            "gate_item": "Phase5 forward target freeze",
            "status": "PASS",
            "evidence": "25 rows: 20 YES; 5 YES_WITH_QUALIFICATION",
            "consequence": "No FRE-driven target reselection permitted",
        },
        {
            "gate_item": "Repository identity",
            "status": "PASS",
            "evidence": "Zenodo 20330471 and official Europe PMC PMC13293934 supplements",
            "consequence": "Published external cohort provenance is recoverable",
        },
        {
            "gate_item": "Dog/sample identity",
            "status": "PASS",
            "evidence": "8 independent FRE dogs and 4 independent Healthy dogs; duodenal mucosa",
            "consequence": "Sample table and published QC can be frozen",
        },
        {
            "gate_item": "Processed raw gene-count matrix",
            "status": "ABSENT",
            "evidence": "No whole-transcriptome sample-by-gene count matrix in any public file",
            "consequence": "Independent DESeq2 cannot be run",
        },
        {
            "gate_item": "Sample-level normalized expression matrix",
            "status": "ABSENT",
            "evidence": "Only group means for 12 selected bile-acid genes are available",
            "consequence": "Independent limma sensitivity analysis cannot be run",
        },
        {
            "gate_item": "Complete published gene-level DE statistics",
            "status": "ABSENT",
            "evidence": "Only a volcano figure and 12 selected-gene DE rows are public",
            "consequence": "Frozen whole-transcriptome rank/GSEA cannot be reconstructed",
        },
        {
            "gate_item": "Published pathway summary",
            "status": "AVAILABLE_LIMITED",
            "evidence": "25 KEGG GSEA rows with NES, p-value and padj",
            "consequence": "Option A summary-level validation is possible only after authorization",
        },
        {
            "gate_item": "FASTQ payload",
            "status": "NOT_DOWNLOADED",
            "evidence": "No FASTQ is present in Zenodo/Europe PMC; no external FASTQ retrieval attempted",
            "consequence": "Option B requires a separate authorization gate",
        },
        {
            "gate_item": "Phase6 interim gate",
            "status": "YELLOW_ROUTE_AUTHORIZATION_REQUIRED",
            "evidence": "Authentic summary data are usable, but the preregistered independent reanalysis is not possible",
            "consequence": "Stop before formal 25-program classification and before Human Phase7",
        },
    ]
    write_csv(
        METADATA / "PHASE6_DATA_AVAILABILITY_GATE.csv",
        ["gate_item", "status", "evidence", "consequence"],
        rows,
    )


def create_report() -> None:
    phase5_audit = read_csv(METADATA / "PHASE6_PHASE5_INPUT_HASH_AUDIT.csv")
    forward = read_csv(
        ROOT / "06_pseudobulk" / "lodo" / "summary" / "PHASE5_PROGRAM_FORWARD_LIST.csv"
    )
    samples = read_csv(METADATA / "FRE2026_SAMPLE_METADATA.csv")
    qc = read_csv(PARSED / "Supplementary_Table_S2_sample_QC.csv")
    pathways = read_csv(PARSED / "Supplementary_Table_S4_KEGG_GSEA.csv")
    targeted_de = read_csv(PARSED / "Main_Table1_targeted_bile_acid_DE.csv")
    inventory = read_csv(PARSED / "SUPPLEMENTARY_CONTENT_INVENTORY.csv")

    phase5_pass = sum(row["status"] == "PASS" for row in phase5_audit)
    yes_n = sum(row["forward_to_phase6"] == "YES" for row in forward)
    qualified_n = sum(row["forward_to_phase6"] == "YES_WITH_QUALIFICATION" for row in forward)
    fre_n = sum(row["condition"] == "FRE" for row in samples)
    healthy_n = sum(row["condition"] == "Healthy" for row in samples)
    mapping_rates = [float(row["mapping_rate_percent"]) for row in qc]
    proper_rates = [float(row["properly_matched_percent"]) for row in qc]
    gc_values = [float(row["gc_percent"]) for row in qc]

    pathway_by_name = {row["Pathway"]: row for row in pathways}
    key_pathways = [
        "Citrate cycle (TCA cycle)",
        "AMPK signaling pathway",
        "Carbon metabolism",
        "Mitophagy - animal",
        "Peroxisome",
        "NF-kappa B signaling pathway",
        "Cytokine-cytokine receptor interaction",
        "Bile secretion",
    ]
    pathway_lines = []
    for name in key_pathways:
        row = pathway_by_name[name]
        pathway_lines.append(
            f"| {name} | {float(row['NES']):.3f} | {float(row['padj']):.3g} | {row['direction']} |"
        )

    de_by_gene = {row["Gene"]: row for row in targeted_de}
    bile_gene_names = ["NR1H4", "RXRA", "NR0B2", "FGF19", "SLC51A", "ABCC3", "ABCC4", "HSD3B7"]
    bile_lines = []
    for gene in bile_gene_names:
        row = de_by_gene[gene]
        bile_lines.append(
            f"| {gene} | {float(row['log2FC']):.2f} | {float(row['padj']):.4g} | {row['Direction']} |"
        )

    report_text = f"""# Phase6：独立犬 FRE bulk RNA-seq 外部验证

## 当前结论

**Phase6 interim data gate：YELLOW — ROUTE AUTHORIZATION REQUIRED**

正式 Phase6 外部验证尚未完成。实时 repository 审计确认了队列、样本和已发表汇总结果的真实性，但没有发现可支持预注册独立重分析的 whole-transcriptome processed gene-level matrix，也没有完整的 gene-level DE statistics。依照 Phase6 指令第十三条，本阶段在数据门控处停止，等待 Option A 或 Option B 的明确授权；未下载 FASTQ，未进入 Human Phase7。

## 1. Discovery 输入冻结

- `PHASE5_MANIFEST.csv`：{phase5_pass}/{len(phase5_audit)} 项大小和 SHA-256 全部通过。
- `PHASE5_PROGRAM_FORWARD_LIST.csv`：25 个目标保持冻结，其中 {yes_n} 个 `YES`、{qualified_n} 个 `YES_WITH_QUALIFICATION`。
- Full7 仍是 Primary discovery result；FRE 仅可提供 independent canine FRE transcriptomic support。

## 2. Repository 实际提供的文件

### Zenodo DOI 10.5281/zenodo.20330471

实时记录只公开 1 个文件：

- `Table 6.docx`，15,204 bytes；repository MD5 `732b93b8e9d54d9b7e84d5f6bac99854`，本地 MD5 已核验。
- 内容是 12 个预选胆汁酸相关基因的 CTRL/FRE 组均值 raw counts 与 TPM；不是 sample-level matrix，也不是全转录组矩阵。

### Europe PMC PMCID PMC13293934

官方 supplementary ZIP 包含 {len(inventory)} 个成员：4 个补充 PDF 图、6 个 DOCX 表和 4 个主文图像预览。逐项内容如下：

- S1：犬只人口学资料与 FRE 临床摘要。
- S2：12 个样本的测序 QC；表格为 DOCX 内嵌位图。
- S3：sample ID、人口学资料和 PCA1/PCA2 坐标。
- S4：25 条 published KEGG GSEA 汇总，含 NES、p-value、padj。
- S5：subsampling sensitivity 汇总；表格为 DOCX 内嵌位图。
- S6：12 个预选胆汁酸基因的组均值 counts/TPM；与 Zenodo 文件逐字节相同。
- Figures S1-S4：PCA、全局 volcano 图和两张 KEGG pathway map；volcano 没有配套的完整 gene table。

没有发现：processed raw counts、whole-transcriptome normalized expression matrix、sample-level TPM matrix、完整 DE statistics 或 FASTQ payload。

## 3. 12 犬样本身份与已发表 QC

- 总计 12 只独立犬：{fre_n} FRE + {healthy_n} Healthy。
- FRE：8 只 client-owned 犬，经至少 1 个月独占低敏饮食后临床缓解而确诊；取材时均为十二指肠黏膜活检。
- Healthy：4 只 purpose-bred research Beagle，无胃肠病史或临床证据。
- 组织：proximal duodenal mucosa；bulk RNA-seq。
- Published sample QC：mapping rate {min(mapping_rates):.1f}-{max(mapping_rates):.1f}%，properly matched {min(proper_rates):.1f}-{max(proper_rates):.1f}%，GC {min(gc_values):.0f}-{max(gc_values):.0f}%；没有明显 sequencing-QC outlier 的公开迹象。
- 年龄、品种、来源和环境与组别结构性混杂；样本量仅 12，因此即使后续重分析也不应加入复杂协变量。

## 4. 数据能力判定

| 候选统计路线 | 公开数据是否满足 | 当前结论 |
|---|---:|---|
| raw-count DESeq2 `~condition` | 否 | 无 sample-by-gene raw-count matrix |
| normalized-expression limma | 否 | 无 sample-level whole-transcriptome normalized matrix |
| published full-rank replication | 否 | 无完整 gene/log2FC/statistic/p-value 表 |
| published summary-level support | 是，但有限 | 只有 25 条 KEGG GSEA + 12 个 selected-gene DE/组均值；需 Option A 授权 |
| FASTQ reprocessing | 未执行 | 需 Option B 单独授权 |

因此目前无法创建以下正式结果：

- `FRE_to_discovery_gene_mapping.csv`
- `FRE_vs_Healthy_ranked.csv`
- 使用 Phase4 冻结 Hallmark/Reactome/GO BP gene sets 的 independent external GSEA
- stable leading-edge direction concordance
- `PHASE6_EXTERNAL_VALIDATION_MATRIX.csv`

25/25 forward programs 尚未分配 `STRONG_SUPPORT`、`MODERATE_SUPPORT`、`DIRECTIONAL_SUPPORT`、`DISCORDANT`、`NOT_REPLICATED` 或 `NOT_ASSESSABLE` 正式分类。现在强行分类会把 published KEGG summary 冒充冻结 pipeline 的独立重分析。

## 5. Published-summary contextual observations（非正式 Phase6 validation）

这些结果仅用于说明 Option A 可利用的公开证据，不能替代正式 external-validation matrix。

| Published KEGG pathway | NES | padj | FRE direction |
|---|---:|---:|---|
{chr(10).join(pathway_lines)}

### 跨谱系 mitochondrial/metabolic decrease

Phase5 的跨谱系主发现是 CIE-lower mitochondrial/oxidative/metabolic programs；published FRE bulk S4 中 TCA、carbon metabolism、AMPK、mitophagy 和 peroxisome 均为 FRE-higher。现阶段这是一个清楚但尚未正式分级的 discovery-vs-published-bulk direction discordance。它可能来自 FRE subtype、bulk 细胞组成、数据库/rank 定义和 cohort 结构差异，不能据此推翻 Full7 discovery。

### TNF/NF-kappaB

Published FRE bulk 的 `NF-kappa B signaling pathway` 为 FRE-lower（NES=-1.886，padj=0.005），`Cytokine-cytokine receptor interaction` 也为 FRE-lower（NES=-2.372，padj=0.005）。这与 Phase5 多个 lineage 的 CIE-higher inflammatory signal 相反；若 Option A 获授权，应标注为 **bulk-level discordance**，不能写成 T-cell 或 Myeloid 机制被否定。

### Myeloid TLR4/JAK-STAT/cytokine programs

S4 没有提供可直接对应的 TLR4 或 JAK-STAT 行，也没有全基因 rank/leading edge。因此 Myeloid TLR4/JAK-STAT 在当前公开数据层面不可正式评估。`IL6_JAK_STAT3_SIGNALING` 仍只能解释为 STAT3 arm，不得升级为 IL-23-specific mechanism。

### KCNK16

KCNK16 不在公开的 12-gene targeted DE/TPM 表中，volcano 图也不提供可复用的精确数值，因此方向、FDR 和表达水平目前均不可评估。

### Bile-acid axis：仅作 FRE contextual observation

Published `Bile secretion` 为 FRE-higher（NES=1.935，padj=0.005）。主文公开的 selected-gene 结果为：

| Gene | log2FC | padj | Published direction |
|---|---:|---:|---|
{chr(10).join(bile_lines)}

这些是 FRE 队列自身的 hypothesis-driven contextual findings，不能回填或重塑 Phase5 frozen forward list。

## 6. Discovery 与 external cohort 的边界

- Discovery：GSE254005，3 Healthy vs 4 CIE，single-cell RNA-seq，lineage-specific pseudobulk。
- External：4 Healthy vs 8 FRE，whole-duodenal-mucosa bulk RNA-seq。
- FRE 是治疗反应亚型，不等同于整个 CIE。
- Bulk 可以提供 gene/pathway/program 的整体方向支持，不能验证细胞比例、cell-state abundance 或信号来源的具体 lineage。
- 最终措辞必须保持为 `independent canine FRE transcriptomic support for CIE-associated programs`，不能称为 definitive replication。

## 7. 后续二选一（需用户授权）

### Option A：published summary-level validation

使用 S4 的 25 条 KEGG GSEA、主文 12 个 selected-gene DE statistics 和 S6 group means，对 Phase5 forward list 做明确降级的 summary-level support。优点是无需新增大规模下载；限制是不能运行冻结 Hallmark/Reactome/GO pipeline、不能计算 stable leading-edge concordance，预计多个 program 将为 `NOT_ASSESSABLE`。

### Option B：FASTQ acquisition and reprocessing

另行核实原始读段的公开位置和体量，在新门控授权后，按 CanFam3.1/Salmon/gene-level counts/DESeq2 路线重处理，再运行冻结 external GSEA。优点是可完成正式 gene-level rank 和 program matrix；代价是下载量、计算量、reference/release reconstruction 和批次/流程偏差明显增加。

## 8. 当前门控与停止点

**Interim gate：YELLOW — ROUTE AUTHORIZATION REQUIRED**

理由：数据来源、样本身份、published QC 和 summary-level evidence 均可靠，但 processed gene-level data 不足以执行预注册的独立重分析。依据硬门控，本次停止于此：

- 没有下载 FASTQ。
- 没有对 25 个 program 做正式外部支持分类。
- 没有修改 Phase5 forward list。
- 没有下载 GSE266616，也没有进入 Human Crohn Phase7。

在选择 Option A 或 Option B 之前，不建议进入 Phase7。

## 9. 关键交付物

- `00_metadata/PHASE6_FRE_DATA_AUDIT.csv`
- `00_metadata/FRE2026_DOWNLOAD_MANIFEST.csv`
- `00_metadata/FRE2026_SAMPLE_METADATA.csv`
- `00_metadata/PHASE6_PHASE5_INPUT_HASH_AUDIT.csv`
- `00_metadata/PHASE6_DATA_AVAILABILITY_GATE.csv`
- `01_raw_processed/dog_FRE_2026/reference/parsed_tables/`
- `00_metadata/PHASE6_MANIFEST.csv`

Generated UTC: {now_utc()}
"""
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(report_text, encoding="utf-8", newline="\n")
    log(f"Wrote {REPORT}")


def file_role(relative: str) -> tuple[str, str]:
    if relative.startswith("scripts/"):
        return "phase6_script", "phase6_reproducible_code"
    if relative.startswith("logs/"):
        return "phase6_log", "phase6_execution_log"
    if relative == "08_results/PHASE6_REPORT.md":
        return "phase6_interim_report", "phase6_data_gate_reporting"
    if relative.startswith("00_metadata/"):
        return "phase6_metadata", "phase6_audit_and_provenance"
    if "/downloaded/" in relative:
        return "external_source_download", "phase6_minimal_public_download"
    if "/extracted/" in relative:
        return "external_source_original", "phase6_supplement_extraction"
    if "/parsed_tables/" in relative:
        return "external_source_parsed_table", "phase6_supplement_ingest"
    if "/reference/" in relative:
        return "external_source_reference", "phase6_repository_and_article_audit"
    return "phase6_output", "phase6_data_gate"


def create_manifest() -> None:
    candidates: set[Path] = set()
    explicit_metadata = [
        "PHASE6_FRE_DATA_AUDIT.csv",
        "FRE2026_DOWNLOAD_MANIFEST.csv",
        "FRE2026_SAMPLE_METADATA.csv",
        "PHASE6_PHASE5_INPUT_HASH_AUDIT.csv",
        "PHASE6_SOFTWARE_VERSIONS.csv",
        "PHASE6_DATA_AVAILABILITY_GATE.csv",
    ]
    for name in explicit_metadata:
        candidates.add(METADATA / name)
    candidates.update(path for path in RAW.rglob("*") if path.is_file())
    candidates.update(
        ROOT / "scripts" / name
        for name in (
            "22_phase6_fre_download_audit.py",
            "23_phase6_fre_ingest.py",
            "24_phase6_data_gate_report.py",
        )
    )
    candidates.update(
        ROOT / "logs" / name
        for name in (
            "phase6_fre_download_audit.log",
            "phase6_fre_ingest.log",
            "phase6_data_gate_report.log",
        )
    )
    candidates.add(REPORT)
    missing = [str(path) for path in candidates if not path.is_file()]
    if missing:
        raise RuntimeError(f"Manifest candidates missing: {missing}")

    rows: list[dict[str, object]] = []
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
    print(f"[{now_utc()}] Wrote manifest with {len(rows)} formal Phase6 files: {MANIFEST}", flush=True)


def verify_manifest() -> None:
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
        raise RuntimeError(f"Phase6 manifest verification failed: {failures}")
    print(f"[{now_utc()}] Phase6 manifest self-check PASS: {len(rows)}/{len(rows)} files", flush=True)


def main() -> int:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("", encoding="utf-8")
    create_gate_table()
    create_report()
    log("Phase6 stopped at YELLOW data gate; no FASTQ and no Human data were accessed")
    log("Finalizing Phase6 formal-file manifest; no further formal files will be mutated")
    create_manifest()
    verify_manifest()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        log(f"FATAL: {type(exc).__name__}: {exc}")
        raise
