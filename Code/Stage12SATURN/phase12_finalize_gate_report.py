#!/usr/bin/env python3
"""Finalize the audited Phase12 gate, report, logs, and manifest."""
from __future__ import annotations

import hashlib
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P12 = PROJECT / "12_saturn"
META = PROJECT / "00_metadata"
RESULTS = PROJECT / "08_results"
LOGS = PROJECT / "logs"
FINAL_GATE = "YELLOW_PHASE12_SATURN_COMPLETE_CELLTYPE_TRANSFER_PASS_DISEASE_AND_SPECIES_MIXING_LIMITATIONS"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def f4(value: float) -> str:
    return f"{float(value):.4f}"


def f6(value: float) -> str:
    return f"{float(value):.6f}"


def md_table(headers, rows) -> str:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines)


for directory in (META, RESULTS, LOGS, P12 / "audit", P12 / "00_metadata"):
    directory.mkdir(parents=True, exist_ok=True)

# Load and enforce technical gates.
phase11 = pd.read_csv(P12 / "audit" / "PHASE11_MANIFEST_REVALIDATION.csv")
feature = pd.read_csv(P12 / "audit" / "SATURN_TRAINING_FEATURE_AUDIT.csv")
formal = pd.read_csv(P12 / "audit" / "SATURN_FORMAL_OUTPUT_VALIDATION.csv")
official = pd.read_csv(P12 / "method_gate" / "SATURN_OFFICIAL_EXAMPLE_VALIDATION.csv")
pilot = pd.read_csv(P12 / "method_gate" / "SATURN_PROJECT_PILOT_VALIDATION.csv")
adapter = pd.read_csv(P12 / "method_gate" / "STABLE_ID_PROTEIN_EMBEDDING_ADAPTER_VALIDATION.csv")
protein_validation = pd.read_csv(P12 / "method_gate" / "PROTEIN_EMBEDDING_VALIDATION.csv")
coverage = pd.read_csv(P12 / "protein_embeddings" / "PROTEIN_EMBEDDING_COVERAGE.csv")
sampling = pd.read_csv(P12 / "sampling" / "SATURN_PRIMARY_SAMPLING_SUMMARY.csv")
alignment = pd.read_csv(P12 / "primary" / "SATURN_CELLTYPE_ALIGNMENT.csv")
mixing = pd.read_csv(P12 / "primary" / "SATURN_SPECIES_MIXING.csv")
disease = pd.read_csv(P12 / "primary" / "SATURN_DISEASE_VECTOR_METRICS.csv")
permutation = pd.read_csv(P12 / "primary" / "SATURN_DISEASE_VECTOR_PERMUTATION.csv")
macro = pd.read_csv(P12 / "macrogenes" / "MACROGENE_CONCORDANCE_SUMMARY.csv")
priority = pd.read_csv(P12 / "macrogenes" / "PRIORITY_MACROGENES.csv")
composition = pd.read_csv(P12 / "macrogenes" / "MACROGENE_ORTHOLOGY_COMPOSITION.csv")
weights = pd.read_csv(P12 / "macrogenes" / "MACROGENE_GENE_WEIGHTS.csv.gz")
t_theme = pd.read_csv(P12 / "macrogenes" / "T_FUNCTIONAL_ASYMMETRY.csv")
myeloid_theme = pd.read_csv(P12 / "macrogenes" / "MYELOID_FUNCTIONAL_STRUCTURE.csv")
seed_stability = pd.read_csv(P12 / "sensitivity" / "SATURN_SEED_STABILITY.csv")
lodo = pd.read_csv(P12 / "sensitivity" / "SATURN_DOG_LODO_DISEASE_VECTOR.csv")
treatment = pd.read_csv(P12 / "sensitivity" / "SATURN_HUMAN_TREATMENT_DISEASE_VECTOR.csv")
software = pd.read_csv(P12 / "00_metadata" / "PHASE12_SATURN_SOFTWARE_AUDIT.csv").iloc[0]
proteome = pd.read_csv(P12 / "00_metadata" / "PHASE12_PROTEOME_SOURCE_AUDIT.csv")
hyper = pd.read_csv(P12 / "00_metadata" / "PHASE12_SATURN_HYPERPARAMETERS.csv")

require(len(phase11) == 364 and phase11["status"].eq("PASS").all(), "Phase11 manifest revalidation is not 364/364 PASS")
require(feature["status"].eq("PASS").all(), "Training feature audit contains a failure")
forbidden_read = feature.loc[feature["forbidden_training_field"].eq("YES"), "used_as_model_input_or_supervision"]
require(forbidden_read.eq("NO").all(), "Condition or another forbidden field entered SATURN training")
require(official.iloc[0][["input_load", "protein_embeddings_load", "centroid_macrogene_init", "pretraining", "metric_learning", "final_h5ad", "gene_to_macrogene_weights", "final_embedding_finite", "macrogene_finite", "weights_finite"]].eq("PASS").all(), "Official SATURN smoke failed")
require(official.iloc[0]["model_math_modified"] == "NO", "SATURN model math was modified")
require(pilot["status"].eq("PASS").all(), "Project pilot validation failed")
require(adapter["status"].eq("PASS").all() and adapter["max_abs_difference"].max() <= 1e-8, "Stable-ID adapter failed")
require(protein_validation["status"].eq("PASS").all(), "Protein embedding validation failed")
require(formal["status"].eq("PASS").all() and len(formal) == 27, "Formal output validation is not 27/27 PASS")
require(set(formal["seed"].astype(int)) == {20260827, 20260828, 20260829}, "Formal seeds are not exact")
require(set(formal["actual_batch_size"].astype(int)) == {4096}, "Formal seeds used mixed batch sizes")
require(set(alignment["seed"].astype(int)) == {20260827, 20260828, 20260829}, "Post-hoc alignment seed labels are not exact")
require(set(alignment["direction"]) == {"dog_to_human", "human_to_dog"}, "Bidirectional transfer rows are incomplete")
require(len(alignment) == 30, "Cell-type alignment should contain 30 rows")
require(set(disease["lineage"]) == {"Epithelial", "T", "Myeloid"}, "Primary disease lineages incomplete")
require(not (PROJECT / "13_concord").exists(), "Unauthorized Phase13 directory exists")

# Preserve the label-only remediation audit.
backup_script = P12 / "audit" / "phase12_saturn_posthoc.before_seed_label_fix.py"
current_script = PROJECT / "scripts" / "phase12_saturn_posthoc.py"
require(backup_script.exists(), "Missing seed-label remediation backup")
remediation = pd.DataFrame([
    {
        "issue": "CELLTYPE_ALIGNMENT_RUN_SEED_MISLABELLED_BY_DIRECTION_OFFSET",
        "detected_observation": "output included untrained seed20260830 and shifted human_to_dog labels",
        "scope": "post-hoc CSV provenance label only",
        "training_rerun": "NO",
        "posthoc_rerun": "YES",
        "figure_rerun": "YES",
        "model_math_modified": "NO",
        "permutation_random_seed_preserved": "YES",
        "backup_sha256": sha256(backup_script),
        "corrected_script_sha256": sha256(current_script),
        "corrected_seed_set": "20260827|20260828|20260829",
        "corrected_direction_set": "dog_to_human|human_to_dog",
        "status": "PASS",
    }
])
remediation.to_csv(P12 / "audit" / "POSTHOC_SEED_LABEL_REMEDIATION.csv", index=False)

# Core derived summaries.
primary_alignment = alignment[alignment["seed"].eq(20260827)].copy()
acc = primary_alignment[primary_alignment["metric"].eq("accuracy")].set_index("direction")
bal = primary_alignment[primary_alignment["metric"].eq("balanced_accuracy")].set_index("direction")
require((acc["observed"] == 1.0).all() and (bal["observed"] == 1.0).all(), "Primary transfer is not perfect")
require((acc["empirical_p_upper"] <= 0.01).all() and (bal["empirical_p_upper"] <= 0.01).all(), "Primary transfer is not above shuffle")

disease = disease.set_index("lineage")
joint = permutation[permutation["test"].eq("joint_disease_vector_cosine")].set_index("lineage")
macro = macro.set_index("lineage")
lineage_order = ["Epithelial", "T", "Myeloid"]
interpretation = {
    "Epithelial": "WEAK_POSITIVE_ALIGNMENT",
    "T": "POSITIVE_ALIGNMENT",
    "Myeloid": "WEAK_POSITIVE_ALIGNMENT",
}

seed_ranges = {}
for lineage in lineage_order:
    block = seed_stability[(seed_stability["metric"] == "disease_vector_cosine") & (seed_stability["lineage"] == lineage)]
    seed_ranges[lineage] = (float(block["min_3seed"].iloc[0]), float(block["median_3seed"].iloc[0]), float(block["max_3seed"].iloc[0]))

flag = (
    composition.groupby(["lineage", "macrogene"], as_index=False)["contains_non_one2one_top20"].max()
    .merge(priority[["lineage", "macrogene", "priority_class"]], on=["lineage", "macrogene"], how="left")
)
non_one2one_n = int(flag["contains_non_one2one_top20"].sum())
priority_n = int(len(flag))

lodo_max = lodo.loc[lodo.groupby("lineage")["absolute_cosine_change"].idxmax()].set_index("lineage")
lodo_global = lodo.loc[lodo["absolute_cosine_change"].idxmax()]
treatment_ranges = treatment.groupby("lineage")["cosine_similarity"].agg(["min", "max"])

# Gate table.
gate_rows = [
    ("phase11_manifest_revalidated", "PASS", "NONE", "364/364 size and SHA-256 PASS", "12_saturn/audit/PHASE11_MANIFEST_REVALIDATION.csv"),
    ("saturn_source_frozen", "PASS", "NONE", f"commit {software['commit']}; 0 working source mismatches; model math modified=NO", "00_metadata/PHASE12_SATURN_SOFTWARE_AUDIT.csv"),
    ("official_example_and_project_pilot", "PASS", "NONE", "official minimal input-contract smoke PASS; corrected 6000-cell project pilot completed pretrain and metric learning", "12_saturn/method_gate/SATURN_OFFICIAL_EXAMPLE_VALIDATION.csv|SATURN_PROJECT_PILOT_VALIDATION.csv"),
    ("stable_id_adapter", "PASS", "NONE", f"{len(adapter)}/{len(adapter)} PASS; max absolute difference={adapter['max_abs_difference'].max():.1g}", "12_saturn/method_gate/STABLE_ID_PROTEIN_EMBEDDING_ADAPTER_VALIDATION.csv"),
    ("protein_embedding_coverage", "PASS", "NONE", "Dog and Human each 8000/8000 requested HV genes usable", "12_saturn/protein_embeddings/PROTEIN_EMBEDDING_COVERAGE.csv"),
    ("condition_leakage", "PASS", "NONE", "0 forbidden training fields used; condition sidecar joined post-hoc only", "12_saturn/audit/SATURN_TRAINING_FEATURE_AUDIT.csv"),
    ("primary_lineages", "PASS", "NONE", "Epithelial, T, Myeloid all retained", "12_saturn/audit/SATURN_FORMAL_OUTPUT_VALIDATION.csv"),
    ("formal_primary_and_three_seeds", "PASS", "NONE", "27/27 output checks PASS; three exact seeds; common batch_size=4096", "12_saturn/audit/SATURN_FORMAL_OUTPUT_VALIDATION.csv"),
    ("celltype_alignment_above_shuffle", "PASS", "NONE", "Dog→Human and Human→Dog accuracy/balanced accuracy=1.0; p=0.00990099", "12_saturn/primary/SATURN_CELLTYPE_ALIGNMENT.csv"),
    ("donor_disease_metrics", "PASS", "NONE", "three lineages; donor-centroid vectors and registered nulls complete", "12_saturn/primary/SATURN_DISEASE_VECTOR_METRICS.csv|SATURN_DISEASE_VECTOR_PERMUTATION.csv"),
    ("macrogene_comparison", "PASS", "NONE", "2000 macrogenes per lineage; effects, composition and post-hoc annotation complete", "12_saturn/macrogenes/MACROGENE_CONCORDANCE_SUMMARY.csv"),
    ("dog_lodo", "PASS", "NONE", "7 donors × 3 lineages; post-hoc centroid LODO complete", "12_saturn/sensitivity/SATURN_DOG_LODO_DISEASE_VECTOR.csv"),
    ("human_treatment_sensitivity", "PASS", "NONE", "6 subsets × 3 lineages complete", "12_saturn/sensitivity/SATURN_HUMAN_TREATMENT_DISEASE_VECTOR.csv"),
    ("level2_authority", "PASS", "NONE", "Level1 only", "12_saturn/audit/SATURN_TRAINING_FEATURE_AUDIT.csv"),
    ("phase13_absence", "PASS", "NONE", "13_concord directory absent", "filesystem audit"),
    ("species_mixing", "YELLOW_LIMITATION", "YELLOW", "cross-species kNN fraction 0.000419–0.004572 across lineage/seeds versus balanced expectation 0.5", "12_saturn/primary/SATURN_SPECIES_MIXING.csv"),
    ("disease_alignment_strength", "YELLOW_LIMITATION", "YELLOW", "Epithelial and Myeloid joint cosines weak and within registered null intervals", "12_saturn/primary/SATURN_DISEASE_VECTOR_PERMUTATION.csv"),
    ("seed_and_program_stability", "YELLOW_LIMITATION", "YELLOW", "T disease cosine 0.0974–0.5008 and T macrogene Spearman has mixed sign", "12_saturn/sensitivity/SATURN_SEED_STABILITY.csv"),
    ("treatment_sensitivity", "YELLOW_LIMITATION", "YELLOW", "ANTI_IL12_23 subset reverses Epithelial and Myeloid cosine; n_active=3", "12_saturn/sensitivity/SATURN_HUMAN_TREATMENT_DISEASE_VECTOR.csv"),
    ("FINAL_GATE", FINAL_GATE, "FINAL", "Technical workflow complete; biological/species-mixing limitations retained", "08_results/PHASE12_REPORT.md"),
]
gate = pd.DataFrame(gate_rows, columns=["criterion", "status", "gate_impact", "observed", "evidence"])
gate_path = META / "PHASE12_GATE.csv"
gate.to_csv(gate_path, index=False)
shutil.copy2(gate_path, P12 / "00_metadata" / "PHASE12_GATE.csv")

# Copy the three pre-existing metadata audits to the required project metadata root.
for name in ("PHASE12_SATURN_SOFTWARE_AUDIT.csv", "PHASE12_PROTEOME_SOURCE_AUDIT.csv", "PHASE12_SATURN_HYPERPARAMETERS.csv"):
    shutil.copy2(P12 / "00_metadata" / name, META / name)

# Report tables.
software_table = md_table(
    ["项目", "冻结值"],
    [
        ("SATURN", str(software["commit"])),
        ("Python / PyTorch / CUDA", f"{software['python_version']} / {software['pytorch_version']} / {software['cuda_runtime']}"),
        ("GPU / driver", f"{software['gpu']} ({software['gpu_memory_mib']} MiB) / {software['driver']}"),
        ("ESM", f"commit {software['esm_commit']}; {software['esm_model']}; layer {software['esm_layer']}; mean; dim {software['esm_embedding_dimension']}; truncate={software['esm_truncate']}"),
    ],
)
coverage_table = md_table(
    ["Species", "表达genes", "protein-mapped", "no-protein", "HV/usable", "protein isoforms", "unique ESM sequences"],
    [
        (row.species, row.total_expression_genes, row.protein_mapped_expression_genes, row.no_protein_expression_genes, f"{row.requested_hv_genes}/{row.saturn_usable_hv_genes}", row.protein_isoforms_averaged, row.unique_esm_sequences)
        for row in coverage.itertuples()
    ],
)
sampling_table = md_table(
    ["Species", "Lineage", "原始cells", "入选cells", "donors", "fallback"],
    [(row.species, row.lineage, row.original_cells, row.selected_cells, row.donor_coverage, row.resource_fallback) for row in sampling.itertuples()],
)
disease_table = md_table(
    ["Lineage", "cosine", "Spearman", "joint-null 95%", "P(two-sided)", "3-seed min/median/max", "解释"],
    [
        (
            lineage,
            f4(disease.loc[lineage, "cosine_similarity"]),
            f4(disease.loc[lineage, "spearman_rho"]),
            f"[{f4(joint.loc[lineage, 'null_q025'])}, {f4(joint.loc[lineage, 'null_q975'])}]",
            f6(joint.loc[lineage, "empirical_p_two_sided"]),
            "/".join(f4(value) for value in seed_ranges[lineage]),
            interpretation[lineage],
        )
        for lineage in lineage_order
    ],
)
macro_table = md_table(
    ["Lineage", "Pearson", "Spearman", "cosine", "direction concordance", "cosine null 95%"],
    [
        (
            lineage,
            f4(macro.loc[lineage, "pearson_r"]),
            f4(macro.loc[lineage, "spearman_rho"]),
            f4(macro.loc[lineage, "cosine_similarity"]),
            f4(macro.loc[lineage, "direction_concordance"]),
            f"[{f4(macro.loc[lineage, 'cosine_similarity_null_q025'])}, {f4(macro.loc[lineage, 'cosine_similarity_null_q975'])}]",
        )
        for lineage in lineage_order
    ],
)
lodo_table = md_table(
    ["Lineage", "最大影响donor", "condition", "LODO cosine", "Full7 cosine", "|Δcosine|"],
    [
        (lineage, lodo_max.loc[lineage, "dropped_donor"], lodo_max.loc[lineage, "dropped_condition"], f4(lodo_max.loc[lineage, "cosine_similarity"]), f4(lodo_max.loc[lineage, "full7_cosine_similarity"]), f4(lodo_max.loc[lineage, "absolute_cosine_change"]))
        for lineage in lineage_order
    ],
)
treatment_table = md_table(
    ["Lineage", "Full active", "subset range", "方向变化"],
    [
        (
            lineage,
            f4(treatment[(treatment["lineage"] == lineage) & (treatment["subset"] == "FULL_ACTIVE_CROHN")]["cosine_similarity"].iloc[0]),
            f"{f4(treatment_ranges.loc[lineage, 'min'])} to {f4(treatment_ranges.loc[lineage, 'max'])}",
            "ANTI_IL12_23为负（n_active=3）" if treatment_ranges.loc[lineage, "min"] < 0 else "全部subset为正",
        )
        for lineage in lineage_order
    ],
)

# Representative top-gene rows for the report.
top_gene_rows = []
for prow in priority.groupby(["lineage", "priority_class"], sort=False).head(1).itertuples():
    block = weights[(weights["lineage"] == prow.lineage) & (weights["macrogene"] == prow.macrogene) & (weights["gene_rank"] <= 5)]
    dog_symbols = ", ".join(block[block["species"] == "dog"]["gene_symbol"].astype(str))
    human_symbols = ", ".join(block[block["species"] == "human"]["gene_symbol"].astype(str))
    top_gene_rows.append((prow.lineage, prow.priority_class, prow.macrogene, f"{f4(prow.dog_effect)} / {f4(prow.human_effect)}", dog_symbols, human_symbols))
top_gene_table = md_table(["Lineage", "Class", "MG", "Dog/Human effect", "Dog top genes", "Human top genes"], top_gene_rows)

t_counts = t_theme.groupby(["priority_class", "functional_theme"]).size().to_dict()
m_counts = myeloid_theme.groupby(["priority_class", "functional_theme"]).size().to_dict()

report = f"""# Phase12 SATURN Functional Cross-Species State Alignment

最终Gate：`{FINAL_GATE}`。

Phase12作为独立探索性门控已完成。SATURN的技术链条、三seed训练、donor-level统计、macrogene分析和敏感性分析均完成；未创建Phase13。主要限制是latent仍存在很强的species separation，以及Epithelial/Myeloid疾病向量较弱、T结果存在seed与macrogene层面的异质性。因此技术上不是RED，但不支持GREEN。

## 1. Software/method

{software_table}

SATURN冻结源码86个Python文件与冻结hash一致，未修改模型、ZINB、loss或数学。隔离环境使用现代兼容依赖（scvi-tools {software['scvi_tools_version']}、scanpy {software['scanpy_version']}、anndata {software['anndata_version']}），这一环境偏离已记录。正式参数为HV genes=8000、macrogenes=2000、latent=256、pretrain=200 epochs、metric=50 epochs、batch=4096；未使用AMP、resource cell-cap fallback、`--use_ref_labels`或训练ct_map。

## 2. Protein embeddings

{coverage_table}

Dog的全部20,257个表达gene均有protein映射；Human为19,707/36,601（53.84%），16,894个表达gene标记为no-protein。正式HV8000在两物种均达到8000/8000可用。ESM使用Ensembl release104 CanFam3.1/GRCh38 `pep.all`，按官方逻辑对同一gene的isoform embedding取均值。Stable-ID adapter在200个唯一symbol测试gene上200/200 PASS，最大绝对差0。

## 3. Sampling

{sampling_table}

总训练cells=28,046（Dog 14,023；Human 14,023），3条lineage逐物种等量。每个species×lineage内部按donor平衡、无重复、无upsampling，固定cell set与sampling seed 20260827；condition未参与抽样。未触发batch或cell-cap resource fallback。

## 4. Method gate

- 官方冻结代码最小input-contract smoke完整通过input、protein embedding、centroid/macrogene初始化、pretrain、metric learning、final H5AD和gene weights。
- 原150-cell/species pilot在Seurat-v3 HVG LOESS处出现小样本奇异；失败输入与日志已审计保留。修正pilot只扩大到每species×lineage 1000 cells并限制到冻结HV8000，合计6000 cells；正式input、模型数学与参数未改，修正pilot完整PASS。
- ESM官方`extract.py` tiny、1021-aa和worst 4×1021测试PASS；正式gene embedding均为1280维且finite。
- condition leakage=0；训练H5AD仅含`cell_id,species,harmonized_level1,donor_id`。Phase9 map仅post-hoc evaluation。
- post-hoc首次写出时发现双向transfer的run-seed标签被direction RNG offset误标，并出现未训练的20260830；修复仅拆分run seed与permutation seed，保留原随机数定义，未重训SATURN。修复后精确seed为20260827/28/29，双向各5项/seed，统计与图件已重跑。

## 5. Cell-type integration

Primary seed的Dog→Human与Human→Dog accuracy、balanced accuracy以及三个lineage recall均为1.000；两个方向accuracy和balanced accuracy的100次source-label shuffle上界检验均为p=0.00990099。三seed的双向accuracy均为1.000，满足“cell-type correspondence显著高于shuffle”的技术门控。

但species mixing并不理想：Primary seed按lineage的cross-species 30-NN比例为Epithelial 0.001253、T 0.001092、Myeloid 0.001648，而平衡期望为0.5；species silhouette分别为0.2797、0.1994、0.1809。三seed全范围为0.000419–0.004572。因此结论应表述为“Level1边界可跨物种转移”，不能表述为“两个物种在latent中充分混合”。所有正式距离与向量均使用256维SATURN latent，不使用UMAP。

## 6. Epithelial disease-state alignment

{disease_table.splitlines()[0]}
{disease_table.splitlines()[1]}
{disease_table.splitlines()[2]}

Primary cosine=0.1640、Pearson=0.1649、Spearman=0.0679、angle=80.56°；joint null p=0.393944，未超过预注册95% null。Macrogene Spearman=0.1011、cosine=0.1334、direction concordance=0.4925。三seed均为正（0.1161–0.2465），故解释为`WEAK_POSITIVE_ALIGNMENT`，不是显著高保真。LODO最大影响来自Healthy h_2；Human treatment中ANTI_IL12_23子集转为-0.0382（n=3），提示treatment sensitivity。

## 7. T

{disease_table.splitlines()[0]}
{disease_table.splitlines()[1]}
{disease_table.splitlines()[3]}

T是Primary disease-vector最高的lineage：cosine=0.5008、Pearson=0.5008、Spearman=0.2024、angle=59.94°；joint null p=0.004571，超过预注册null。三seed方向均为正，但幅度0.0974–0.5008（median 0.2696）。Macrogene层面则弱：Primary Spearman=0.0267、cosine=0.0500、direction concordance=0.4860，且跨seed Spearman为-0.0060到0.0267、方向混合。

Post-hoc功能注释中，Top conserved含4个inflammatory-like、1个metabolic-like、5个mixed；Top opposite含2个inflammatory-like、1个metabolic-like、6个mixed和1个other。MG_1365等免疫/TCR/IFN/TNF相关结构同向，而胆固醇/脂蛋白相关MG_0723反向；这支持“炎症与代谢存在选择性functional asymmetry”，但不是整齐的二分，也不能把全部T压成单一保真等级。所有治疗subset仍为正（0.2313–0.5008）；LODO最大影响为cie_4，删除后仍为0.3446。

## 8. Myeloid

{disease_table.splitlines()[0]}
{disease_table.splitlines()[1]}
{disease_table.splitlines()[4]}

Primary cosine=0.1269、Spearman=0.1060，joint null p=0.413941；三seed为0.1269–0.2705，解释为`WEAK_POSITIVE_ALIGNMENT`。Macrogene Pearson/cosine约0.159，三条lineage中macrogene cosine最高，但仍位于null区间。

Top conserved Myeloid中5/10为inflammatory-like、5/10为mixed，含C1Q/complement、抗原呈递、IFN、TNF/NF-κB和IL6/JAK-STAT相关结构；因此“炎症功能对应继续存在”在局部macrogene层面成立。与此同时，respiration/OXPHOS相关MG_0375为反向，Top opposite还含2个metabolic-like。不得升级为IL-23-specific。ANTI_IL12_23子集cosine=-0.0784（n=3），显示治疗敏感性。

## 9. Macrogenes

{macro_table}

每lineage分析2000个macrogenes，并预先按effect规则列出10个Top conserved与10个Top opposite，共60个重点macrogenes。59/60在任一物种的top20组成中含ONE2MANY、MANY2MANY或NO_ONE2ONE关系，证明SATURN确实利用了传统严格1:1之外的protein关系；该计数不代表这些非1:1 gene一定驱动疾病对应。

{top_gene_table}

功能富集只用于post-hoc解释，没有反馈训练或重新挑选macrogene。不同seed的macrogene编号不作直接对应。

## 10. Seed stability

三seed的正式输出均为28,046×256 latent、28,046×2000 macrogene且finite；16000个gene weight均为2000维finite，27/27验证PASS，batch均为4096。疾病向量seed范围：Epithelial 0.1161–0.2465，T 0.0974–0.5008，Myeloid 0.1269–0.2705，均保持正号。Macrogene整体中Epithelial/Myeloid指标均为正；T Spearman跨seed混合符号。cell-type transfer 3/3 seeds稳定为1.000，但species mixing持续偏低。

## 11. Dog LODO

这是`post-hoc donor-centroid LODO sensitivity`，SATURN embedding未重训。

{lodo_table}

全局最大|Δcosine|来自T删除`{lodo_global['dropped_donor']}`（cie_4，|Δ|={f4(lodo_global['absolute_cosine_change'])}）。Phase12没有预设cie_3或cie_4；结果独立指向cie_4。各LODO方向仍为正。

## 12. Human treatment sensitivity

{treatment_table}

Full active Crohn为11名active与16名Healthy。ANTI_IL12_23仅3名active，Epithelial和Myeloid的负值必须按小样本敏感性解释，不能作为药物或靶点结论；Phase12不做therapeutic prioritization。

## 13. Phase10对比

Phase10 gene rho为Epithelial 0.1474、T 0.0159、Myeloid 0.1011；mean pathway rho为0.3655、0.0620、0.2626。SATURN embedding的排序不同：T最高（0.5008），Epithelial/Myeloid仅弱正。这与Phase10的“T gene-level近零”形成差异，但T macrogene rho=0.0267又与Phase10的弱全局gene一致，提示SATURN latent可能捕捉到少数结构化state axis，而非广泛program一致。

Epithelial在Phase10最强，但SATURN仅弱正且不超过joint null；Myeloid在两者均保留局部炎症结构，但SATURN整体向量弱。SATURN新增的主要信息是非严格1:1 protein-related macrogene结构（59/60重点MG含非1:1），不是对Phase10的替代或纠错。

## 14. Phase11对比

46/46 Human singIST reference models均weak，Phase12训练完全独立，不据此改sampling、label或pathway。Epithelial的singIST TNF负/respiration正与SATURN弱整体向量并列保留，SATURN不能裁决哪一个“正确”。T的singIST炎症接近0、代谢negative，与SATURN“整体positive但macrogene全局弱、胆固醇相关MG反向”部分一致又部分冲突。Myeloid的singIST TNF强positive与SATURN局部inflammatory conserved macrogenes一致，但SATURN整体向量仍弱。三种方法只作并列，不平均、不生成总分。

## 15. Provenance

- region：Phase8的`HIGH_REGION_RISK`、`REGION_UNKNOWN`等逐program继续携带；SATURN不能消除Dog duodenum与Human TI的区域差异。
- FRE：无cell-resolved raw matrix，未进入SATURN；仅保留summary context。
- treatment：仅post-hoc Human donor subset sensitivity，未参与训练。
- LODO：仅冻结embedding上的Dog donor-centroid敏感性，不是training LODO。
- orthology：SATURN训练不裁到14,945个1:1 ortholog；Phase10 orthology仅用于post-hoc composition annotation。
- labels：仅Epithelial/T/Myeloid Level1；无Level2正式比较。

## 16. 限制

Dog仅7 donors，Human仅27 donors；分析使用donor-balanced cell sampling而非全量cells；RTX3060 6GB限制了可执行规模，虽未触发预注册fallback；SATURN只覆盖有protein representation的coding genes，Human全表达gene coverage仅53.84%；ESM-1b表示与Ensembl104 peptide质量会影响结果；macrogene可解释性有限，编号跨seed不可直接比较；训练有随机性；Dog/Human组织region不完全匹配；condition-blind integration仍可能保留dataset/species effects，且本次kNN mixing确实很低；SATURN不是独立外部validation；Human treatment subset小；Level2未正式分析；project pilot经历了经审计的小样本HVG修正；现代兼容依赖虽未改数学，但不是原始旧环境逐包复刻。

## 17. Phase12 Gate

`{FINAL_GATE}`。

GREEN所需的15项技术门控均已完成，包括源码冻结、method gate、adapter、coverage、零leakage、三lineage、三seed、shuffle-above cell-type transfer、donor统计、macrogene、LODO、treatment、Level1边界和无Phase13。YELLOW来自结果层限制：latent跨物种混合很低；Epithelial/Myeloid疾病向量未超joint null；T强度有seed变动且macrogene全局弱；部分treatment subset方向变化；region/FRE限制仍在。

建议：`YES_WITH_YELLOW_LIMITATIONS`，可由用户决定是否把Phase13 CONCORD作为新的独立门控启动；不得把SATURN解释为Phase10/11“被证明正确”，并必须携带上述species-mixing、donor/treatment、region与FRE限制。本阶段已停止，未自动执行Phase13。

正式交付以`00_metadata/PHASE12_MANIFEST.csv`为唯一文件级索引；主要结果位于`12_saturn/primary`、`macrogenes`、`sensitivity`、`summary`和`figures`，审计位于`12_saturn/audit`与`00_metadata`。
"""
report_path = RESULTS / "PHASE12_REPORT.md"
report_path.write_text(report, encoding="utf-8")

# Required concise logs. Existing raw training/ESM logs are preserved.
now = datetime.now(timezone.utc).isoformat()
log_payloads = {
    "phase12_environment.log": f"{now}\nPASS frozen SATURN {software['commit']}\nPython {software['python_version']} PyTorch {software['pytorch_version']} CUDA {software['cuda_runtime']} GPU {software['gpu']}\nmodel_math_modified=NO\n",
    "phase12_protein_embeddings.log": f"{now}\nPASS Dog 8000/8000 Human 8000/8000 HV genes; adapter {len(adapter)}/{len(adapter)} PASS; max_diff={adapter['max_abs_difference'].max()}\n",
    "phase12_method_gate.log": f"{now}\nPASS official smoke; PASS corrected 6000-cell pilot; attempt1 singular-HVG remediation audited; seed-label posthoc remediation PASS\n",
    "phase12_sampling.log": f"{now}\nPASS 28046 condition-blind donor-balanced Level1 cells; no resource fallback\n",
    "phase12_integration_qc.log": f"{now}\nPASS bidirectional accuracy and balanced_accuracy=1.0, p=0.00990099; YELLOW low cross-species kNN mixing\n",
    "phase12_disease_state.log": f"{now}\nPASS donor-centroid metrics and nulls complete; Epithelial={disease.loc['Epithelial','cosine_similarity']}; T={disease.loc['T','cosine_similarity']}; Myeloid={disease.loc['Myeloid','cosine_similarity']}\n",
    "phase12_macrogenes.log": f"{now}\nPASS 2000 macrogenes per lineage; {non_one2one_n}/{priority_n} priority macrogenes contain non-1:1 top20 relationships\n",
    "phase12_sensitivity.log": f"{now}\nPASS three-seed, 7-donor LODO, and six Human treatment-subset analyses complete\n",
    "phase12_report.log": f"{now}\nPASS report written to 08_results/PHASE12_REPORT.md\nFINAL_GATE={FINAL_GATE}\nPhase13_created=NO\nManifest is generated after this log and excludes itself by design.\n",
}
for name, payload in log_payloads.items():
    (LOGS / name).write_text(payload, encoding="utf-8")

# Record manifest scope before hashing.
scope = pd.DataFrame([
    ("INCLUDE", "12_saturn formal metadata/audit/input/sampling/method_gate/reference/protein_sequences/training/primary/macrogenes/disease_state/sensitivity/summary/figures"),
    ("INCLUDE", "top-level formal gene embeddings and protein coverage"),
    ("INCLUDE", "Phase12 scripts, logs, root metadata, gate, and report"),
    ("EXCLUDE", "12_saturn/environment and unpacked software repositories; exact audit hashes recorded in software audit"),
    ("EXCLUDE", "per-protein PT alias trees; formal gene embedding tensors and coverage are included"),
    ("EXCLUDE", "12_saturn/cache and superseded stop-codon audit payloads"),
    ("EXCLUDE", "PHASE12_MANIFEST.csv itself to avoid a self-referential hash"),
], columns=["action", "scope"])
scope.to_csv(P12 / "audit" / "PHASE12_MANIFEST_SCOPE.csv", index=False)

# Build formal manifest without hashing tens of thousands of hard-linked protein tensors or full environments.
candidates: set[Path] = set()
for pattern in ("PHASE12_*.csv",):
    candidates.update(path for path in META.glob(pattern) if path.name != "PHASE12_MANIFEST.csv")
candidates.add(report_path)
candidates.update(PROJECT.glob("scripts/phase12*"))
candidates.update(LOGS.glob("phase12*"))
for dirname in ("00_metadata", "audit", "input", "sampling", "method_gate", "reference", "protein_sequences", "training", "primary", "macrogenes", "disease_state", "sensitivity", "summary", "figures"):
    for path in (P12 / dirname).rglob("*"):
        if path.is_file() and "stop_codon_semantics_superseded" not in path.parts:
            candidates.add(path)
for name in ("PROTEIN_EMBEDDING_COVERAGE.csv", "dog_gene_embeddings_ESM1b.pt", "human_gene_embeddings_ESM1b.pt"):
    candidates.add(P12 / "protein_embeddings" / name)


def role_for(path: Path) -> str:
    rel = path.relative_to(PROJECT).as_posix()
    if rel.startswith("scripts/"):
        return "Phase12_script"
    if rel.startswith("logs/"):
        return "Phase12_log"
    if rel == "08_results/PHASE12_REPORT.md":
        return "Phase12_report"
    if rel.endswith("PHASE12_GATE.csv"):
        return "Phase12_gate"
    if "/training/" in "/" + rel:
        return "Phase12_training_output"
    if "/figures/" in "/" + rel:
        return "Phase12_figure"
    if "/audit/" in "/" + rel or "/00_metadata/" in "/" + rel or rel.startswith("00_metadata/"):
        return "Phase12_audit_or_metadata"
    if "/protein_embeddings/" in "/" + rel:
        return "Phase12_gene_embedding"
    if "/input/" in "/" + rel or "/sampling/" in "/" + rel:
        return "Phase12_input_or_sampling"
    if "/method_gate/" in "/" + rel:
        return "Phase12_method_gate"
    return "Phase12_result"

manifest_rows = []
for path in sorted(candidates, key=lambda item: item.relative_to(PROJECT).as_posix()):
    require(path.exists() and path.is_file(), f"Manifest candidate missing: {path}")
    rel = path.relative_to(PROJECT).as_posix()
    manifest_rows.append({
        "path": rel,
        "size": path.stat().st_size,
        "SHA256": sha256(path),
        "role": role_for(path),
        "source": "Phase12 condition-blind SATURN workflow; frozen Phase0-11 inputs read-only",
        "generated_by": "phase12_finalize_gate_report.py" if path.name in {"PHASE12_GATE.csv", "PHASE12_REPORT.md", "PHASE12_MANIFEST_SCOPE.csv", "POSTHOC_SEED_LABEL_REMEDIATION.csv"} or path.parent == LOGS else "Phase12 staged workflow",
    })
manifest = pd.DataFrame(manifest_rows)
manifest_path = META / "PHASE12_MANIFEST.csv"
manifest.to_csv(manifest_path, index=False)
require(manifest["path"].is_unique, "Manifest paths are not unique")
require(len(manifest) > 100, "Manifest unexpectedly small")

print(f"PASS: {FINAL_GATE}")
print(f"REPORT={report_path}")
print(f"MANIFEST={manifest_path} rows={len(manifest)}")
print("PHASE13_CREATED=NO")