"""Validate Phase 1 deliverables, freeze the environment, and write the gate report.

This script does not perform annotation or any downstream biological analysis.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import logging
import os
import platform
import shutil
import subprocess
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import psutil
import scanpy as sc
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config.paths import (  # noqa: E402
    DOG_RAW_DIR,
    H5AD_CLEAN_DIR,
    H5AD_RAW_DIR,
    LOG_DIR,
    METADATA_DIR,
    PROJECT_ROOT,
    QC_DIR,
    RESULTS_DIR,
    ensure_phase1_directories,
)


SCRIPT_VERSION = "1.0.0"
EXPECTED_PHASE0_HASHES = {
    "MASTER_DATA_AUDIT.csv": "f2f7dd9aa290008adf0fefe0546e8efa11363a280d4a83ea41fabdb35734c99f",
    "PHASE0_AUDIT_REPORT.md": "18d71ddf573988ab2e4b3895e288db89c3e14698fe8c51199c693ad6b399389d",
}
EXPECTED_ARCHIVE_BYTES = 180_408_320
EXPECTED_ARCHIVE_SHA256 = "ebe414a8f87275fe35ddcb2075d93cf594036f663bd1878e883f1cbc87b24cd2"
EXPECTED_DONORS = [
    "GSM8031698_h_1",
    "GSM8031699_h_2",
    "GSM8031700_h_3",
    "GSM8031701_cie_2",
    "GSM8031702_cie_3",
    "GSM8031703_cie_4",
    "GSM8031704_cie_6",
]
EXPECTED_CONDITIONS = {
    "GSM8031698_h_1": "Healthy",
    "GSM8031699_h_2": "Healthy",
    "GSM8031700_h_3": "Healthy",
    "GSM8031701_cie_2": "CIE",
    "GSM8031702_cie_3": "CIE",
    "GSM8031703_cie_4": "CIE",
    "GSM8031704_cie_6": "CIE",
}
EXPECTED_MITO_SYMBOLS = {
    "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6",
    "MT-CO3", "MT-ND3", "ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB",
}
EXPECTED_STEPS = [
    "raw",
    "count filter",
    "feature filter",
    "broad mt filter",
    "doublet",
    "secondary non-epithelial mt filter",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def configure_logging() -> logging.Logger:
    logger = logging.getLogger("phase1_report")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter(
        fmt="%(asctime)sZ %(levelname)s %(message)s", datefmt="%Y-%m-%dT%H:%M:%S"
    )
    handler = logging.FileHandler(LOG_DIR / "phase1_report.log", mode="w", encoding="utf-8")
    handler.setFormatter(formatter)
    handler.formatter.converter = __import__("time").gmtime
    logger.addHandler(handler)
    stream = logging.StreamHandler(sys.stdout)
    stream.setFormatter(formatter)
    stream.formatter.converter = __import__("time").gmtime
    logger.addHandler(stream)
    return logger


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "NOT INSTALLED"


def write_environment_files() -> None:
    required = [
        ("python", platform.python_version(), "CPython executable in project-local environment"),
        ("scanpy", package_version("scanpy"), "imported at report time"),
        ("anndata", package_version("anndata"), "imported at report time"),
        ("numpy", np.__version__, "imported at report time"),
        ("pandas", pd.__version__, "imported at report time"),
        ("scipy", package_version("scipy"), "installed distribution"),
        ("h5py", package_version("h5py"), "installed distribution"),
        ("scrublet", package_version("scrublet"), "exact sklearn neighbors used; annoy not used"),
        ("matplotlib", package_version("matplotlib"), "installed distribution"),
        ("psutil", psutil.__version__, "imported at report time"),
        ("scikit-learn", package_version("scikit-learn"), "Scrublet exact-neighbor backend"),
        ("scikit-image", package_version("scikit-image"), "Scrublet dependency"),
    ]
    pd.DataFrame(required, columns=["software", "version", "notes"]).to_csv(
        METADATA_DIR / "PHASE1_SOFTWARE_VERSIONS.csv", index=False
    )

    freeze = subprocess.run(
        [sys.executable, "-m", "pip", "freeze", "--all"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    lines = [
        "name: canine_cie_phase1",
        "channels:",
        "  - conda-forge",
        "dependencies:",
        f"  - python={platform.python_version()}",
        "  - pip",
        "  - pip:",
    ]
    lines.extend(f"      - {item}" for item in freeze if item.strip())
    (METADATA_DIR / "phase1_environment.yml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def validate_phase0() -> dict[str, str]:
    actual = {}
    for filename, expected_hash in EXPECTED_PHASE0_HASHES.items():
        path = METADATA_DIR / filename
        if not path.is_file():
            raise RuntimeError(f"Missing Phase 0 file: {path}")
        digest = sha256_file(path)
        if digest != expected_hash:
            raise RuntimeError(f"Phase 0 file hash changed: {path}")
        actual[filename] = digest
    return actual


def validate_archive() -> tuple[list[str], str]:
    archive = DOG_RAW_DIR / "GSE254005_RAW.tar"
    if archive.stat().st_size != EXPECTED_ARCHIVE_BYTES:
        raise RuntimeError("Trusted archive byte count changed")
    digest = sha256_file(archive)
    if digest != EXPECTED_ARCHIVE_SHA256:
        raise RuntimeError("Trusted archive SHA-256 changed")

    member_names: list[str] = []
    sample_parts: dict[str, set[str]] = {}
    with tarfile.open(archive, "r") as tar:
        members = tar.getmembers()
        if len(members) != 21:
            raise RuntimeError(f"Expected 21 archive members; found {len(members)}")
        for member in members:
            if not member.isfile() or member.size <= 0:
                raise RuntimeError(f"Invalid archive member: {member.name}")
            if Path(member.name).name != member.name:
                raise RuntimeError(f"Unsafe or nested archive member: {member.name}")
            member_names.append(member.name)
            suffix = next(
                (part for part in ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
                 if member.name.endswith("_" + part)),
                None,
            )
            if suffix is None:
                raise RuntimeError(f"Unexpected archive member: {member.name}")
            prefix = member.name[: -(len(suffix) + 1)]
            sample_parts.setdefault(prefix, set()).add(suffix)
            extracted = tar.extractfile(member)
            if extracted is None:
                raise RuntimeError(f"Cannot read archive member: {member.name}")
            with gzip.GzipFile(fileobj=extracted, mode="rb") as gz:
                while gz.read(8 * 1024 * 1024):
                    pass
    required = {"barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"}
    if set(sample_parts) != set(EXPECTED_DONORS):
        raise RuntimeError(f"Archive donor set mismatch: {sorted(sample_parts)}")
    if any(parts != required for parts in sample_parts.values()):
        raise RuntimeError("One or more donors do not have exactly the 10x triplet")
    return member_names, digest


def validate_extracted(member_names: list[str]) -> None:
    extracted_dir = DOG_RAW_DIR / "extracted"
    for name in member_names:
        path = extracted_dir / name
        if not path.is_file() or path.stat().st_size <= 0:
            raise RuntimeError(f"Missing or empty extracted member: {path}")
        with gzip.open(path, "rb") as handle:
            while handle.read(8 * 1024 * 1024):
                pass


def validate_h5ad() -> dict[str, object]:
    raw_path = H5AD_RAW_DIR / "dog_GSE254005_raw.h5ad"
    clean_path = H5AD_CLEAN_DIR / "dog_GSE254005_clean.h5ad"
    raw = ad.read_h5ad(raw_path, backed="r")
    try:
        if raw.shape != (49_005, 20_257):
            raise RuntimeError(f"Unexpected raw h5ad shape: {raw.shape}")
        raw_donors = set(raw.obs["donor_id"].astype(str))
        if raw_donors != set(EXPECTED_DONORS):
            raise RuntimeError("Raw h5ad donor set mismatch")
    finally:
        raw.file.close()

    clean = ad.read_h5ad(clean_path)
    if clean.shape != (38_186, 20_257):
        raise RuntimeError(f"Unexpected clean h5ad shape: {clean.shape}")
    required_obs = {
        "sample_id", "donor_id", "condition", "dataset", "species", "tissue",
        "total_counts", "n_genes_by_counts", "pct_counts_mt", "scrublet_score",
        "predicted_doublet", "scrublet_evaluated", "coarse_lineage_qc",
    }
    required_var = {"feature_id", "gene_symbol", "feature_type"}
    if not required_obs.issubset(clean.obs.columns):
        raise RuntimeError("Clean h5ad is missing required obs fields")
    if not required_var.issubset(clean.var.columns):
        raise RuntimeError("Clean h5ad is missing required feature identifier fields")
    if "counts" not in clean.layers or not sparse.issparse(clean.layers["counts"]):
        raise RuntimeError("Clean h5ad does not contain sparse raw counts layer")
    if not sparse.issparse(clean.X):
        raise RuntimeError("Clean h5ad X is unexpectedly dense")
    if clean.layers["counts"].data.size and (
        clean.layers["counts"].data.min() < 0
        or not np.allclose(clean.layers["counts"].data, np.round(clean.layers["counts"].data))
    ):
        raise RuntimeError("Clean counts layer is not non-negative integer-valued")
    for donor, condition in EXPECTED_CONDITIONS.items():
        observed = set(clean.obs.loc[clean.obs["donor_id"].astype(str) == donor, "condition"].astype(str))
        if observed != {condition}:
            raise RuntimeError(f"Condition mismatch for {donor}: {observed}")
    return {
        "raw_shape": (49_005, 20_257),
        "clean_shape": clean.shape,
        "healthy_cells": int((clean.obs["condition"].astype(str) == "Healthy").sum()),
        "cie_cells": int((clean.obs["condition"].astype(str) == "CIE").sum()),
        "lineage_counts": clean.obs["coarse_lineage_qc"].astype(str).value_counts().to_dict(),
    }


def validate_tables() -> dict[str, pd.DataFrame]:
    tables = {
        "download": pd.read_csv(METADATA_DIR / "GSE254005_DOWNLOAD_MANIFEST.csv"),
        "sample": pd.read_csv(METADATA_DIR / "GSE254005_SAMPLE_MANIFEST.csv"),
        "ingest": pd.read_csv(METADATA_DIR / "GSE254005_INGEST_SUMMARY.csv"),
        "pre": pd.read_csv(QC_DIR / "dog_qc_pre_filter_per_donor.csv"),
        "flow": pd.read_csv(QC_DIR / "dog_qc_filter_flow.csv"),
        "scrublet": pd.read_csv(QC_DIR / "dog_scrublet_summary.csv"),
        "mito": pd.read_csv(QC_DIR / "dog_mitochondrial_gene_list.csv"),
        "final": pd.read_csv(QC_DIR / "dog_qc_final_per_donor.csv"),
    }
    download = tables["download"]
    if len(download) != 1:
        raise RuntimeError("Download manifest must contain exactly one record")
    download_row = download.iloc[0]
    if (
        int(download_row["actual_bytes"]) != EXPECTED_ARCHIVE_BYTES
        or int(download_row["content_length"]) != EXPECTED_ARCHIVE_BYTES
        or str(download_row["sha256"]).lower() != EXPECTED_ARCHIVE_SHA256
        or str(download_row["archive_readable"]).lower() != "true"
        or str(download_row["status"]) not in {"DOWNLOADED_VALID", "EXISTING_VALID"}
    ):
        raise RuntimeError("Download manifest does not describe the trusted validated archive")
    if len(tables["sample"]) != 7 or len(tables["ingest"]) != 7:
        raise RuntimeError("Sample/ingest manifest does not contain exactly seven donors")
    if set(tables["ingest"]["sample_id"]) != set(EXPECTED_DONORS):
        raise RuntimeError("Ingest manifest donor set mismatch")
    if dict(zip(tables["ingest"]["sample_id"], tables["ingest"]["condition"])) != EXPECTED_CONDITIONS:
        raise RuntimeError("Ingest condition labels mismatch")
    if set(tables["mito"]["gene_symbol"]) != EXPECTED_MITO_SYMBOLS:
        raise RuntimeError("Mitochondrial gene audit does not match the expected 13-gene set")
    if len(tables["flow"]) != 42:
        raise RuntimeError("Filter flow does not contain 6 steps x 7 donors")
    for donor in EXPECTED_DONORS:
        steps = tables["flow"].loc[tables["flow"]["donor"] == donor, "step"].tolist()
        if steps != EXPECTED_STEPS:
            raise RuntimeError(f"Filter flow step mismatch for {donor}")
    if len(tables["scrublet"]) != 7 or tables["scrublet"]["score_threshold"].isna().any():
        raise RuntimeError("Scrublet audit is incomplete")
    if tables["scrublet"]["use_approx_neighbors"].astype(str).str.lower().ne("false").any():
        raise RuntimeError("Scrublet exact-neighbor declaration mismatch")
    retention = tables["final"]["post_qc_cells"] / tables["final"]["raw_cells"]
    if (retention < 0.20).any():
        raise RuntimeError("At least one donor lost nearly all cells")
    return tables


def find_cleanup_targets() -> list[Path]:
    targets: list[Path] = []
    for path in DOG_RAW_DIR.iterdir():
        name = path.name
        if path.is_dir() and name == "GSE254005_RAW.tar.https_segments":
            targets.append(path)
        elif path.is_file() and (
            name.startswith("GSE254005_RAW.tar.corrupt.")
            or name.startswith("GSE254005_RAW.tar.part.corrupt.")
            or name.startswith("GSE254005_RAW.tar.part.untrusted_ftp.")
            or (name.startswith("range_") and name.endswith(".part"))
        ):
            targets.append(path)
    extracted = DOG_RAW_DIR / "extracted"
    if extracted.is_dir():
        targets.extend(path for path in extracted.iterdir() if path.is_file() and ".corrupt." in path.name)
    return sorted(targets)


def cleanup_intermediates(targets: list[Path], do_cleanup: bool) -> pd.DataFrame:
    diagnostics_path = METADATA_DIR / "GSE254005_DOWNLOAD_DIAGNOSTICS.csv"
    if not targets and diagnostics_path.is_file():
        existing = pd.read_csv(diagnostics_path)
        if len(existing) and set(existing["disposition"]) == {
            "removed_after_trusted_inputs_and_h5ad_validated"
        }:
            return existing
    raw_root = DOG_RAW_DIR.resolve()
    final_archive = (DOG_RAW_DIR / "GSE254005_RAW.tar").resolve()
    records: list[dict[str, object]] = []
    for target in targets:
        resolved = target.resolve()
        if not resolved.is_relative_to(raw_root) or resolved == raw_root or resolved == final_archive:
            raise RuntimeError(f"Refusing unsafe cleanup target: {resolved}")
        is_dir_target = resolved.is_dir()
        if is_dir_target:
            files = [item for item in resolved.rglob("*") if item.is_file()]
            byte_count = sum(item.stat().st_size for item in files)
            digest = "NA (directory aggregate)"
            item_count = len(files)
        else:
            byte_count = resolved.stat().st_size
            digest = sha256_file(resolved)
            item_count = 1
        disposition = "recorded_not_removed"
        if do_cleanup:
            if is_dir_target:
                shutil.rmtree(resolved)
            else:
                resolved.unlink()
            disposition = "removed_after_trusted_inputs_and_h5ad_validated"
        records.append(
            {
                "path": str(resolved),
                "type": "directory" if is_dir_target else "file",
                "item_count": item_count,
                "bytes": byte_count,
                "sha256": digest,
                "disposition": disposition,
                "reason": "failed/untrusted transfer, diagnostic range, reusable HTTPS segment, or quarantined corrupt extraction",
            }
        )
    result = pd.DataFrame(
        records,
        columns=["path", "type", "item_count", "bytes", "sha256", "disposition", "reason"],
    )
    result.to_csv(diagnostics_path, index=False)
    return result


def fmt(value: float, digits: int = 1) -> str:
    return f"{value:,.{digits}f}"


def build_report(
    tables: dict[str, pd.DataFrame],
    archive_hash: str,
    h5: dict[str, object],
    cleanup: pd.DataFrame,
) -> str:
    ingest = tables["ingest"].set_index("sample_id").loc[EXPECTED_DONORS]
    pre = tables["pre"].set_index("donor").loc[EXPECTED_DONORS]
    final = tables["final"].set_index("donor").loc[EXPECTED_DONORS]
    scrub = tables["scrublet"].set_index("donor").loc[EXPECTED_DONORS]
    flow = tables["flow"]
    raw_total = int(ingest["cells"].sum())
    final_total = int(final["post_qc_cells"].sum())
    removed_total = raw_total - final_total
    paper_reference = 35_668
    difference = final_total - paper_reference
    difference_pct = 100 * difference / paper_reference
    broad_total = int(
        flow.loc[flow["step"] == "broad mt filter", "cells_after"].sum()
    )
    doublets = int(scrub["predicted_doublets"].sum())
    secondary_removed = int(
        flow.loc[flow["step"] == "secondary non-epithelial mt filter", "cells_removed"].sum()
    )
    worst = final.assign(retention=final["post_qc_cells"] / final["raw_cells"]).sort_values("retention").iloc[0]
    cie2_mt = float(final.loc["GSM8031701_cie_2", "median_mt_pct"])
    gate = "YELLOW"

    lines = [
        "# GSE254005 犬数据 Phase 1 报告",
        "",
        f"生成时间（UTC）：{utc_now()}",
        "",
        "## 执行结论",
        "",
        f"**Phase 1 门控：{gate}。** 7 只犬均成功读取，raw/clean h5ad 均通过独立回读，流程可追溯；但 `GSM8031701_cie_2` 存在显著 donor 级 mt% 异常和较高细胞损失，进入 Phase 2 前必须保留为重点审计对象，不宜将其与其他 donor 等同处理。",
        "",
        "本阶段未进行正式细胞类型注释、Leiden 亚群命名、DEG、pseudobulk、GSEA、LODO、Human 数据处理、singIST、SATURN 或 CONCORD。",
        "",
        "## 1. 下载与完整性",
        "",
        f"- 官方 `GSE254005_RAW.tar` 下载成功：{EXPECTED_ARCHIVE_BYTES:,} bytes（GEO 页面约 172.1 MB）。",
        f"- 冻结 SHA-256：`{archive_hash}`。GEO 未提供已核实 checksum，未伪造 MD5。",
        "- TAR 可读；恰有 21 个非空 gzip member；7 个 donor 各自均恰有 barcodes/features/matrix 三件套；全部 gzip 流通过完整解压/CRC 检查；没有异常额外矩阵。",
        "- 3 Healthy + 4 CIE 的 sample prefix 与 GEO 结构一致；没有把 `cie_2/cie_3/cie_4/cie_6` 擅自映射为论文 CIE_1–CIE_4。",
        "",
        "下载阶段曾因网络/FTP 分段内容不稳定产生不可信或损坏的临时副本；可信 HTTPS 固定范围分段重下并全量 CRC 验证后才原子落盘。所有失败中间件已先登记于 `00_metadata/GSE254005_DOWNLOAD_DIAGNOSTICS.csv`，随后清理；可信 TAR 与验证后的 extracted 三件套均保留。",
        "",
        "## 2. 原始数据结构",
        "",
        "| donor | condition | raw cells | features | matrix shape (cells × features) | nnz |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for donor in EXPECTED_DONORS:
        row = ingest.loc[donor]
        lines.append(
            f"| {donor} | {row['condition']} | {int(row['cells']):,} | {int(row['genes_features']):,} | "
            f"{int(row['matrix_rows_cells']):,} × {int(row['matrix_columns_features']):,} | {int(row['nnz']):,} |"
        )
    lines.extend(
        [
            "",
            f"合并 raw 对象：**{raw_total:,} cells × 20,257 features**。每 donor 有 17 个重复 gene symbol；未删除这些 feature，使用可复现的唯一化 var name，并保留 `feature_id` 与原始 `gene_symbol`。",
            "",
            "## 3. QC 与 Scrublet",
            "",
            "预注册过滤严格按顺序执行：`100 < nCount_RNA < 25000`、`200 < nFeature_RNA < 3000`、`percent.mt < 80%`。未为匹配论文细胞数调参。Scrublet 在每个 donor 内独立运行（expected doublet rate 0.06，seed 20240824，exact sklearn nearest neighbors；未使用 annoy）。",
            "",
            "| donor | raw | final | removed | retention | pre median count | pre median genes | pre median mt% | final median count | final median genes | final median mt% | doublets (%) |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for donor in EXPECTED_DONORS:
        p = pre.loc[donor]
        f = final.loc[donor]
        s = scrub.loc[donor]
        removed = int(f["raw_cells"] - f["post_qc_cells"])
        retention = 100 * f["post_qc_cells"] / f["raw_cells"]
        lines.append(
            f"| {donor} | {int(f['raw_cells']):,} | {int(f['post_qc_cells']):,} | {removed:,} | {retention:.1f}% | "
            f"{fmt(p['median_counts'])} | {fmt(p['median_genes'])} | {fmt(p['median_mt_pct'], 2)} | "
            f"{fmt(f['median_counts'])} | {fmt(f['median_genes'])} | {fmt(f['median_mt_pct'], 2)} | "
            f"{int(s['predicted_doublets'])} ({float(s['predicted_doublet_pct']):.2f}%) |"
        )
    lines.extend(
        [
            "",
            f"宽松三项过滤后 {broad_total:,} cells；Scrublet 共标记 {doublets:,} cells；仅对明确 non-epithelial 的二级 `mt < 12.5%` 又移除 {secondary_removed:,} cells；最终共移除 {removed_total:,}/{raw_total:,}（{100*removed_total/raw_total:.1f}%）。",
            "",
            "犬线粒体基因通过当前 CanFam3.1/Ensembl r104 feature symbol 审计：规则为 symbol 以 `MT-` 开头，另纳入未带前缀的 `ND4L`；共得到可靠的 13 个蛋白编码 mt genes，完整列表已保存。",
            "",
            "QC 分层只使用 marker 检测规则粗分 epithelial/non-epithelial/uncertain，用于决定二级 mt 阈值；这不是正式 annotation，uncertain cells 被保留。规则及 marker 已保存，正式细胞注释留待 Phase 2。",
            "",
            "## 4. 最终数据对象",
            "",
            f"- final cells total：**{final_total:,}**。",
            f"- Healthy：**{int(h5['healthy_cells']):,}**；CIE：**{int(h5['cie_cells']):,}**。",
            f"- 最低保留率 donor：`{worst.name}`，保留 {int(worst['post_qc_cells']):,}/{int(worst['raw_cells']):,}（{100*worst['retention']:.1f}%）。",
            "- clean h5ad 的 `layers['counts']` 保存非负整数 raw counts；`X` 为稀疏的 library-size 10,000 normalization + log1p 矩阵。donor/sample/condition、QC metrics、Scrublet、原始 feature ID 均保留；raw h5ad 未覆盖。",
            "",
            "每 donor 最终细胞数：",
            "",
        ]
    )
    for donor in EXPECTED_DONORS:
        lines.append(f"- `{donor}`：{int(final.loc[donor, 'post_qc_cells']):,}")
    lines.extend(
        [
            "",
            "## 5. 与论文约 35,668 cells 的比较",
            "",
            f"本流程得到 {final_total:,} cells，比论文正参考 35,668 多 {difference:,}（{difference_pct:.2f}%）。该差异未用于反向调节阈值。可能来源包括作者未完全公开的额外 QC/去双细胞细节、对象构建版本差异、以及本流程对 epithelial/uncertain cells 采取的保守 mt 策略。",
            "",
            "## 6. 异常与 Phase 2 前置事项",
            "",
            f"- **主要异常：`GSM8031701_cie_2`。** 最终保留 2,944/5,783（50.9%）；宽松 `mt < 80%` 一步移除 1,257 cells，最终 mt% 中位数仍为 {cie2_mt:.2f}%。其最终细胞以 QC 粗分 epithelial 为主，因此未机械套用 non-epithelial 的 12.5% 阈值。Phase 2 应结合正式细胞类型、复杂度和 donor-level sensitivity analysis 复核，不能为降 mt% 随意删除潜在真实上皮细胞。",
            "- **Doublet：** 各 donor 预测率 0.08%–2.32%，未见极端高 doublet donor；阈值和每 donor 结果均已记录。",
            "- **Feature ID：** 每 donor 有 17 个重复 gene symbol；原始 Ensembl feature ID 与 symbol 已保留，Phase 2 建立 `dog_gene_mapping_master.csv` 时必须显式处理，不能按 symbol 静默丢弃。",
            "- **Metadata：** breed、age、medication、paper dog ID 仍未安全映射，保持 NA；疾病 GEO suffix 未强行对应论文 dog number。",
            "- **下载异常：** 失败/不可信的传输副本只作为中间诊断文件，已登记后清理，不参与任何分析。",
            "",
            "## 7. Phase 1 门控",
            "",
            "### YELLOW",
            "",
            "7/7 donor、Healthy/CIE 标签、mt gene 审计、完整 filter flow、per-donor Scrublet、raw/clean h5ad 回读和 donor 可追溯性均通过。没有 donor 因技术问题接近全部丢失，也没有为复现 35,668 而调参。",
            "",
            "总体数据允许在明确保留上述警告的前提下进入 Phase 2；但 `GSM8031701_cie_2` 的高 mt%/低保留率必须作为正式注释与后续 donor-level sensitivity audit 的优先检查项，因此不标记 GREEN。",
            "",
            "Phase 1 在此停止。",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cleanup-intermediates",
        action="store_true",
        help="After all trusted inputs/outputs validate, remove only recorded failed-transfer intermediates.",
    )
    args = parser.parse_args()
    ensure_phase1_directories()
    logger = configure_logging()
    process = psutil.Process(os.getpid())
    start = datetime.now(timezone.utc)
    logger.info(
        "START script_version=%s python=%s scanpy=%s anndata=%s cleanup_intermediates=%s",
        SCRIPT_VERSION,
        platform.python_version(),
        package_version("scanpy"),
        package_version("anndata"),
        args.cleanup_intermediates,
    )
    try:
        phase0_hashes = validate_phase0()
        logger.info("PHASE0_UNCHANGED hashes=%s", phase0_hashes)
        members, archive_hash = validate_archive()
        validate_extracted(members)
        logger.info(
            "ARCHIVE_VALID bytes=%d sha256=%s members=%d donors=%d gzip_crc=pass",
            EXPECTED_ARCHIVE_BYTES,
            archive_hash,
            len(members),
            len(EXPECTED_DONORS),
        )
        tables = validate_tables()
        h5 = validate_h5ad()
        logger.info(
            "H5AD_VALID raw_shape=%s clean_shape=%s healthy=%d cie=%d raw_counts_layer=pass",
            h5["raw_shape"], h5["clean_shape"], h5["healthy_cells"], h5["cie_cells"],
        )
        write_environment_files()
        targets = find_cleanup_targets()
        logger.info("CLEANUP_TARGETS count=%d targets=%s", len(targets), [str(p) for p in targets])
        cleanup = cleanup_intermediates(targets, args.cleanup_intermediates)
        if args.cleanup_intermediates:
            logger.warning(
                "REMOVED_INTERMEDIATES records=%d bytes=%d; trusted archive/extracted triplets preserved",
                len(cleanup),
                int(cleanup["bytes"].sum()) if len(cleanup) else 0,
            )
        report = build_report(tables, archive_hash, h5, cleanup)
        report_path = RESULTS_DIR / "PHASE1_REPORT.md"
        report_path.write_text(report, encoding="utf-8")
        elapsed = (datetime.now(timezone.utc) - start).total_seconds()
        logger.warning(
            "GATE YELLOW reason=GSM8031701_cie_2 high_mt_and_50.9_percent_retention; Phase2 requires donor-specific audit"
        )
        logger.info(
            "OUTPUT report=%s environment=%s software=%s diagnostics=%s",
            report_path,
            METADATA_DIR / "phase1_environment.yml",
            METADATA_DIR / "PHASE1_SOFTWARE_VERSIONS.csv",
            METADATA_DIR / "GSE254005_DOWNLOAD_DIAGNOSTICS.csv",
        )
        logger.info(
            "END status=success elapsed_seconds=%.3f peak_ram_mb=%.1f",
            elapsed,
            process.memory_info().peak_wset / 1024**2 if hasattr(process.memory_info(), "peak_wset") else process.memory_info().rss / 1024**2,
        )
    except Exception:
        logger.exception("END status=error")
        raise


if __name__ == "__main__":
    main()
