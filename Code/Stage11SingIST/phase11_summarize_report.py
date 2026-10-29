from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
S11 = ROOT / "11_singIST"
PRIMARY = S11 / "primary" / "method_faithful"
SENS = S11 / "donor_aware_sensitivity"
LODO = S11 / "lodo_sensitivity"
SUMMARY = S11 / "summary"
FIG = S11 / "figures"
REPORT = ROOT / "08_results" / "PHASE11_REPORT.md"
LOG_SUM = ROOT / "logs" / "phase11_summary.log"
LOG_REP = ROOT / "logs" / "phase11_report.log"
SUMMARY.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="notebook")
COLORS = {"Epithelial": "#2878B5", "T": "#F28E2B", "Myeloid": "#59A14F"}

super_df = pd.read_csv(PRIMARY / "SUPERPATHWAY_RECAPITULATION.csv")
cell_df = pd.read_csv(PRIMARY / "CELLTYPE_RECAPITULATION.csv")
gene_df = pd.read_csv(PRIMARY / "GENE_CONTRIBUTIONS.csv.gz")
validity = pd.read_csv(PRIMARY / "HUMAN_REFERENCE_MODEL_VALIDITY.csv")
coverage = pd.read_csv(PRIMARY / "ORTHOLOGY_COVERAGE_BY_SUPERPATHWAY.csv")
sensitivity = pd.read_csv(SENS / "SINGIST_FC_SOURCE_SENSITIVITY.csv")
lodo_super = pd.read_csv(LODO / "SINGIST_DOG_LODO_SUPERPATHWAY.csv")
lodo_summary = pd.read_csv(LODO / "SINGIST_DOG_LODO_SUMMARY.csv")
base = pd.read_csv(ROOT / "10_cross_species_baseline" / "summary" /
                   "PHASE10_PROGRAM_BASELINE_MATRIX.csv")
registered = pd.read_csv(S11 / "superpathways" / "REGISTERED_PROGRAM_PANEL.csv")
hallmark_panel = pd.read_csv(S11 / "superpathways" / "HUMAN_DRIVEN_HALLMARK_PANEL.csv")

lineages = ["Epithelial", "T", "Myeloid"]

def short_pathway(x: str) -> str:
    x = str(x)
    for prefix in ("HALLMARK_", "REACTOME_"):
        if x.startswith(prefix):
            x = x[len(prefix):]
    return x.replace("_", " ").title()


# Phase11 vs Phase10 registered-program matrix.
super_cols = super_df[["pathway", "recapitulation", "reference_model_status",
                       "human_reference_PValue", "human_reference_FDR",
                       "reference_shift_status"]].rename(
    columns={"pathway": "representative_pathway",
             "recapitulation": "singIST_superpathway_recap"}
)
cell_pivot = cell_df.pivot(index="pathway", columns="celltype",
                           values="recapitulation").reset_index().rename(
    columns={"pathway": "representative_pathway",
             "Epithelial": "singIST_epithelial_recap",
             "T": "singIST_T_recap",
             "Myeloid": "singIST_myeloid_recap"}
)
compare = base.merge(super_cols, on="representative_pathway", how="left").merge(
    cell_pivot, on="representative_pathway", how="left")
compare = compare.rename(columns={
    "baseline_status": "phase10_baseline_status",
    "dog_common1to1_NES": "phase10_dog_NES",
    "human_common1to1_NES": "phase10_human_NES",
    "dog_lodo_stability": "dog_lodo",
    "human_treatment_stability": "human_treatment",
    "fre_summary_status": "fre_status",
})

def relationship(row: pd.Series) -> str:
    if row.get("reference_model_status") != "REFERENCE_MODEL_VALIDATED":
        return "REFERENCE_MODEL_WEAK_NOT_CORE"
    recap = row.get("singIST_superpathway_recap")
    baseline = str(row.get("phase10_baseline_status"))
    if not np.isfinite(recap):
        return "NOT_INTERPRETABLE"
    if recap >= 0 and baseline.startswith("CONCORDANT"):
        return "FUNCTIONAL_AND_TRADITIONAL_DIRECTIONALLY_ALIGNED"
    if recap >= 0 and baseline == "OPPOSITE":
        return "FUNCTIONAL_RECAPITULATION_DESPITE_TRADITIONAL_OPPOSITION"
    if recap < 0 and baseline.startswith("CONCORDANT"):
        return "FUNCTIONAL_OPPOSITION_DESPITE_TRADITIONAL_CONCORDANCE"
    if recap < 0 and baseline == "OPPOSITE":
        return "FUNCTIONAL_AND_TRADITIONAL_BOTH_OPPOSED"
    return "MIXED_OR_DIRECTION_ONLY_EVIDENCE"

compare["interpretation"] = compare.apply(relationship, axis=1)
wanted = [
    "lineage", "program", "representative_collection", "representative_pathway",
    "phase10_baseline_status", "phase10_dog_NES", "phase10_human_NES",
    "singIST_superpathway_recap", "singIST_epithelial_recap", "singIST_T_recap",
    "singIST_myeloid_recap", "reference_model_status", "human_reference_PValue",
    "human_reference_FDR", "reference_shift_status", "region_risk", "dog_lodo",
    "human_treatment", "fre_status", "phase10_priority", "interpretation"
]
compare[wanted].to_csv(SUMMARY / "PHASE11_VS_PHASE10_BASELINE.csv", index=False)


# Core positive and negative gene-contribution drivers for registered pathways.
registered_paths = set(registered["pathway_id"])
g = gene_df[gene_df["pathway"].isin(registered_paths)].copy()
driver_parts: list[pd.DataFrame] = []
for (_, _), z in g.groupby(["pathway", "celltype"], sort=False):
    pos = z[z["contribution"] > 0].nlargest(5, "contribution").copy()
    neg = z[z["contribution"] < 0].nsmallest(5, "contribution").copy()
    pos["contribution_direction"] = "POSITIVE"
    neg["contribution_direction"] = "NEGATIVE"
    driver_parts.extend([pos, neg])
drivers = pd.concat(driver_parts, ignore_index=True) if driver_parts else g.head(0)
drivers = drivers.rename(columns={"gene": "human_gene"})

dog_de_parts = []
human_de_parts = []
effect_parts = []
for lineage in lineages:
    dog = pd.read_csv(ROOT / "06_pseudobulk" / "de_results" /
                      f"{lineage}_CIE_vs_Healthy_edger.csv")
    dog = dog.sort_values(["FDR", "PValue"], kind="stable").drop_duplicates("unique_var_name")
    dog = dog[["unique_var_name", "logFC", "PValue", "FDR"]].rename(columns={
        "unique_var_name": "dog_ortholog", "logFC": "dog_phase3_logFC",
        "PValue": "dog_phase3_PValue", "FDR": "dog_phase3_FDR"})
    dog["celltype"] = lineage
    dog_de_parts.append(dog)

    human = pd.read_csv(ROOT / "07_human_crohn" / "GSE266616" / "pseudobulk" /
                        "de_results" / f"Human_{lineage}_ActiveCrohn_vs_Healthy.csv")
    human = human.sort_values(["FDR", "PValue"], kind="stable").drop_duplicates("gene_symbol")
    human = human[["gene_symbol", "logFC", "PValue", "FDR"]].rename(columns={
        "gene_symbol": "human_gene", "logFC": "human_phase7_logFC",
        "PValue": "human_PValue", "FDR": "human_FDR"})
    human["celltype"] = lineage
    human_de_parts.append(human)

    tag = lineage.upper()
    eff = pd.read_csv(ROOT / "10_cross_species_baseline" / "gene_effects" /
                      f"{tag}_DOG_HUMAN_GENE_EFFECTS.csv.gz")
    eff = eff.sort_values(["human_FDR", "dog_FDR"], kind="stable").drop_duplicates(
        ["human_symbol", "dog_symbol"])
    eff = eff[["human_symbol", "dog_symbol", "effect_class",
               "dog_lodo_stability_if_available"]].rename(columns={
        "human_symbol": "human_gene", "dog_symbol": "dog_ortholog",
        "effect_class": "phase10_gene_effect_class",
        "dog_lodo_stability_if_available": "dog_lodo_stability"})
    eff["celltype"] = lineage
    effect_parts.append(eff)

dog_de = pd.concat(dog_de_parts, ignore_index=True)
human_de = pd.concat(human_de_parts, ignore_index=True)
effects = pd.concat(effect_parts, ignore_index=True)
drivers = drivers.merge(dog_de, on=["celltype", "dog_ortholog"], how="left")
drivers = drivers.merge(human_de, on=["celltype", "human_gene"], how="left")
drivers = drivers.merge(effects, on=["celltype", "human_gene", "dog_ortholog"], how="left")

reg_meta = registered.rename(columns={"pathway_id": "pathway", "lineage": "celltype"})[
    ["pathway", "celltype", "program", "region_risk", "phase10_priority"]]
drivers = drivers.merge(reg_meta, on=["pathway", "celltype"], how="left")
drivers["notes"] = (
    "singIST signed contribution; not a differential-expression P value; "
    "top five per sign within registered pathway×cell type"
)
driver_cols = [
    "pathway", "collection", "pathway_name", "celltype", "program",
    "human_gene", "dog_ortholog", "contribution", "contribution_direction",
    "dog_phase3_logFC", "dog_phase3_PValue", "dog_phase3_FDR",
    "dog_lodo_stability", "human_phase7_logFC", "human_PValue", "human_FDR",
    "phase10_gene_effect_class", "region_risk", "phase10_priority", "notes"
]
drivers[driver_cols].to_csv(SUMMARY / "SINGIST_CORE_GENE_DRIVERS.csv", index=False)


# Figures.
def savefig(name: str) -> None:
    plt.tight_layout()
    plt.savefig(FIG / name, dpi=220, bbox_inches="tight")
    plt.close()

# 1. Human-driven superpathway overview.
a_ids = set(hallmark_panel["pathway_id"])
z = super_df[super_df["pathway"].isin(a_ids)].copy().sort_values("recapitulation")
fig, ax = plt.subplots(figsize=(11, max(8, 0.32 * len(z))))
colors = np.where(z["recapitulation"] < 0, "#D95F5F", "#4C78A8")
ax.barh([short_pathway(x) for x in z["pathway"]], z["recapitulation"], color=colors)
ax.axvline(0, color="black", lw=1)
ax.axvline(100, color="#555555", lw=1.2, ls="--")
ax.set_xlabel("singIST superpathway recapitulation (%) — raw, not clipped")
ax.set_ylabel("")
ax.set_title("Human-driven Hallmark panel")
savefig("01_SUPERPATHWAY_RECAPITULATION_OVERVIEW.png")

# 2. Cell-type recapitulation heatmap.
piv = cell_df.pivot(index="pathway", columns="celltype", values="recapitulation")
piv = piv.reindex(columns=lineages)
piv.index = [short_pathway(x) for x in piv.index]
v = np.nanmax(np.abs(piv.to_numpy()))
fig, ax = plt.subplots(figsize=(8, max(11, 0.28 * len(piv))))
sns.heatmap(piv, cmap="vlag", norm=TwoSlopeNorm(vmin=-v, vcenter=0, vmax=v),
            ax=ax, cbar_kws={"label": "Recapitulation (%)"})
ax.set_title("Cell-type recapitulation (unclipped scale)")
ax.set_xlabel("")
ax.set_ylabel("")
savefig("02_CELLTYPE_RECAPITULATION_HEATMAP.png")

# 3. Orthology coverage.
fig, ax = plt.subplots(figsize=(9, 6))
sns.stripplot(data=coverage, x="celltype", y="orthology_coverage", hue="celltype",
              palette=COLORS, jitter=0.22, alpha=0.72, ax=ax, legend=False)
sns.boxplot(data=coverage, x="celltype", y="orthology_coverage", color="white",
            width=0.4, fliersize=0, ax=ax)
ax.set_ylabel("Detected Phase10 one-to-one coverage (%)")
ax.set_xlabel("")
ax.set_title("Orthology coverage by superpathway and Level1 block")
savefig("03_ORTHOLOGY_COVERAGE.png")

# 4. Phase10 direction-aligned Dog NES vs singIST.
plot_cmp = compare.copy()
plot_cmp["phase10_direction_aligned_dog_NES"] = (
    np.sign(plot_cmp["phase10_human_NES"]) * plot_cmp["phase10_dog_NES"]
)
fig, ax = plt.subplots(figsize=(9, 7))
sns.scatterplot(data=plot_cmp, x="phase10_direction_aligned_dog_NES",
                y="singIST_superpathway_recap", hue="phase10_baseline_status",
                style="lineage", s=95, ax=ax)
ax.axhline(0, color="black", lw=1)
ax.axhline(100, color="#555555", lw=1, ls="--")
ax.axvline(0, color="black", lw=1)
ax.set_xlabel("Phase10 Dog NES aligned to Human NES direction")
ax.set_ylabel("singIST superpathway recapitulation (%)")
ax.set_title("Traditional pathway baseline vs functional recapitulation")
savefig("04_PHASE10_VS_SINGIST_SCATTER.png")

# 5. T inflammatory vs metabolic contrast.
t_groups = {
    "Inflammatory": ["HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                     "HALLMARK_INFLAMMATORY_RESPONSE",
                     "HALLMARK_IL2_STAT5_SIGNALING", "R-HSA-449147"],
    "Metabolic/mitochondrial": ["HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                                "GO:0009060", "R-HSA-1428517", "R-HSA-72766"]
}
t_rows = []
for family, ids in t_groups.items():
    for pid in ids:
        hit = cell_df[(cell_df["pathway"] == pid) & (cell_df["celltype"] == "T")]
        if not hit.empty:
            t_rows.append({"family": family, "pathway": pid,
                           "recapitulation": hit.iloc[0]["recapitulation"]})
tplot = pd.DataFrame(t_rows)
fig, ax = plt.subplots(figsize=(11, 6))
sns.barplot(data=tplot, x="pathway", y="recapitulation", hue="family", ax=ax,
            palette={"Inflammatory": "#D95F02", "Metabolic/mitochondrial": "#1B9E77"})
ax.axhline(0, color="black", lw=1)
ax.axhline(100, color="#555555", lw=1, ls="--")
tick_labels = [short_pathway(x.get_text()) for x in ax.get_xticklabels()]
ax.set_xticks(ax.get_xticks())
ax.set_xticklabels(tick_labels, rotation=35, ha="right")
ax.set_ylabel("T-cell recapitulation (%)")
ax.set_xlabel("")
ax.set_title("T inflammatory vs metabolic programs")
savefig("05_T_INFLAMMATORY_VS_METABOLIC.png")

# 6. Core contribution heatmap.
drv_plot = drivers.loc[drivers["contribution"].abs().nlargest(min(60, len(drivers))).index].copy()
drv_plot["feature"] = drv_plot["pathway"].map(short_pathway) + " | " + drv_plot["human_gene"]
drv_piv = drv_plot.pivot_table(index="feature", columns="celltype", values="contribution",
                               aggfunc="sum").reindex(columns=lineages)
v = np.nanmax(np.abs(drv_piv.to_numpy())) if len(drv_piv) else 1
fig, ax = plt.subplots(figsize=(9, max(9, 0.28 * len(drv_piv))))
sns.heatmap(drv_piv, cmap="vlag", norm=TwoSlopeNorm(vmin=-v, vcenter=0, vmax=v),
            ax=ax, cbar_kws={"label": "Signed gene contribution"})
ax.set_title("Largest absolute registered-program gene contributions")
ax.set_xlabel("")
ax.set_ylabel("")
savefig("06_CORE_GENE_CONTRIBUTION_HEATMAP.png")

# 7. LODO delta heatmap.
ld = lodo_super[lodo_super["run_id"] != "Full7"].copy()
full_map = lodo_super[lodo_super["run_id"] == "Full7"].set_index("pathway")["recapitulation"]
ld["delta_from_full7"] = ld["recapitulation"] - ld["pathway"].map(full_map)
lp = ld.pivot(index="pathway", columns="drop_donor", values="delta_from_full7").reindex(columns=[
    "h_1", "h_2", "h_3", "cie_2", "cie_3", "cie_4", "cie_6"])
lp.index = [short_pathway(x) for x in lp.index]
v = np.nanmax(np.abs(lp.to_numpy()))
fig, ax = plt.subplots(figsize=(10, max(7, 0.38 * len(lp))))
sns.heatmap(lp, cmap="vlag", norm=TwoSlopeNorm(vmin=-v, vcenter=0, vmax=v),
            ax=ax, cbar_kws={"label": "LODO − Full7 recapitulation (points)"})
ax.set_title("Dog donor influence on registered superpathways")
ax.set_xlabel("Dropped donor")
ax.set_ylabel("")
savefig("07_DOG_LODO_RECAPITULATION_HEATMAP.png")

# 8. Track A vs Track B.
ss = sensitivity[sensitivity["recap_level"] == "SUPERPATHWAY"].copy()
fig, ax = plt.subplots(figsize=(8, 8))
sns.scatterplot(data=ss, x="method_faithful_recap", y="donor_aware_recap",
                hue="collection", style="same_sign", s=85, ax=ax)
lo = np.nanmin(ss[["method_faithful_recap", "donor_aware_recap"]].to_numpy())
hi = np.nanmax(ss[["method_faithful_recap", "donor_aware_recap"]].to_numpy())
ax.plot([lo, hi], [lo, hi], color="#555555", ls="--", lw=1)
ax.axhline(0, color="black", lw=1)
ax.axvline(0, color="black", lw=1)
ax.set_xlabel("Track A method-faithful recapitulation (%)")
ax.set_ylabel("Track B donor-aware recapitulation (%)")
ax.set_title("FC-source sensitivity")
savefig("08_TRACK_A_VS_TRACK_B.png")


# Key metrics and report.
validated_n = int((validity["reference_model_status"] == "REFERENCE_MODEL_VALIDATED").sum())
weak_n = int((validity["reference_model_status"] == "REFERENCE_MODEL_WEAK").sum())
negative = super_df[super_df["recapitulation"] < 0].sort_values("recapitulation")
over100 = super_df[super_df["recapitulation"] > 100].sort_values("recapitulation", ascending=False)
cell_strong = cell_df.loc[cell_df["recapitulation"].idxmax()]
cell_weak = cell_df.loc[cell_df["recapitulation"].idxmin()]
fc_super = sensitivity[sensitivity["recap_level"] == "SUPERPATHWAY"]
cie4_n = int(lodo_summary["cie_4_is_max_influence"].astype(bool).sum())
exact_zero = super_df[super_df["recapitulation"] == 0].sort_values("pathway")
positive = super_df[super_df["recapitulation"] > 0].sort_values(
    "recapitulation", ascending=False)
lodo_drop = lodo_super[lodo_super["run_id"] != "Full7"].copy()
lodo_full_map = lodo_super[lodo_super["run_id"] == "Full7"].set_index(
    "pathway")["recapitulation"]
lodo_drop["absolute_change"] = (
    lodo_drop["recapitulation"] - lodo_drop["pathway"].map(lodo_full_map)
).abs()
lodo_donor_stats = lodo_drop.groupby("drop_donor", as_index=False).agg(
    sum_absolute_change=("absolute_change", "sum"),
    mean_absolute_change=("absolute_change", "mean"),
    max_absolute_change=("absolute_change", "max"),
).sort_values("sum_absolute_change", ascending=False)
largest_lodo_donor = lodo_donor_stats.iloc[0]
env_record = json.loads((ROOT / "11_singIST" / "environment" /
                         "PHASE11_ENVIRONMENT_HASH.json").read_text(encoding="utf-8"))
environment_hash = env_record["environment_hash"]

metrics = {
    "human_driven_hallmark_n": int(len(hallmark_panel)),
    "registered_program_n": int(len(registered)),
    "unique_superpathway_n": int(len(super_df)),
    "validated_reference_model_n": validated_n,
    "weak_reference_model_n": weak_n,
    "negative_superpathway_n": int(len(negative)),
    "over100_superpathway_n": int(len(over100)),
    "exact_zero_superpathway_n": int(len(exact_zero)),
    "superpathway_max": {"pathway": positive.iloc[0]["pathway"],
                         "recapitulation": float(positive.iloc[0]["recapitulation"])},
    "superpathway_min": {"pathway": negative.iloc[0]["pathway"],
                         "recapitulation": float(negative.iloc[0]["recapitulation"])},
    "celltype_max": {"pathway": cell_strong["pathway"], "celltype": cell_strong["celltype"],
                     "recapitulation": float(cell_strong["recapitulation"])},
    "celltype_min": {"pathway": cell_weak["pathway"], "celltype": cell_weak["celltype"],
                     "recapitulation": float(cell_weak["recapitulation"])},
    "orthology_coverage_min": float(coverage["orthology_coverage"].min()),
    "orthology_coverage_median": float(coverage["orthology_coverage"].median()),
    "orthology_coverage_max": float(coverage["orthology_coverage"].max()),
    "track_same_sign_n": int(fc_super["same_sign"].astype(bool).sum()),
    "track_superpathway_n": int(len(fc_super)),
    "track_b_nonzero_recapitulation_n": int(
        (sensitivity["donor_aware_recap"] != 0).sum()),
    "cie4_max_influence_pathway_n": cie4_n,
    "largest_aggregate_lodo_donor": largest_lodo_donor["drop_donor"],
    "environment_hash": environment_hash,
}
(SUMMARY / "PHASE11_KEY_METRICS.json").write_text(
    json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")

def pathway_lines(df: pd.DataFrame, n: int | None = None) -> str:
    if df.empty:
        return "无。"
    use = df if n is None else df.head(n)
    return "；".join(f"{r.pathway} ({r.recapitulation:.2f}%)" for r in use.itertuples())

def driver_lines(direction: str, n: int = 12) -> str:
    z = drivers[drivers["contribution_direction"] == direction].copy()
    z = z.nlargest(n, "contribution") if direction == "POSITIVE" else z.nsmallest(n, "contribution")
    if z.empty:
        return "无。"
    return "；".join(
        f"{short_pathway(r.pathway)}/{r.celltype}/{r.human_gene} "
        f"({r.contribution:.2f})" for r in z.itertuples())

def recap_for(pid: str, celltype: str | None = None) -> str:
    if celltype is None:
        z = super_df[super_df["pathway"] == pid]
    else:
        z = cell_df[(cell_df["pathway"] == pid) & (cell_df["celltype"] == celltype)]
    return "NA" if z.empty else f"{z.iloc[0]['recapitulation']:.2f}%"

gate = "YELLOW" if weak_n > 0 or len(negative) > 0 or len(over100) > 0 else "GREEN"
report = f"""# Phase11 Report — singIST Functional Recapitulation Analysis

## 1. Software

- 冻结版本：singIST 1.0.2，R 4.6.1，Bioconductor 3.23。
- singIST 使用官方 Bioconductor 3.23 source tarball；安装包未修改。
- 正式流程使用官方 `orthology_cache` 与 `FC_list` API；Windows LOOCV 后端为 `BiocParallel::SnowParam(workers=3)`。
- 环境 SHA-256：`{environment_hash}`。软件来源、package/tree hashes、conda explicit lock 与 session info 已冻结。

## 2. Method gate

- 官方 example：PASS，包括构造器、`fitOptimal()`、superpathway/cell-type/gene 三层输出。
- Frozen Ortholog Adapter：**不需要**。singIST 1.0.2 原生支持 `orthology_cache`；Ensembl BioMart 被跳过。
- 官方 wrapper、逐步官方函数、预计算一次官方 `diff_expressed()` 再经 `FC_list` 复用，最大绝对差均为 **0**（门限 ≤1e-8）。
- Primary mapping 仍为 Phase10 Ensembl release 104 / CanFam3.1：14,945 个 stable-ID 1:1 pairs；其中 14,910 个具有 singIST symbol-level 可表达标识，35 个 Human symbol 缺失条目被显式排除，未重映射。
- singIST 1.0.2 的 pathway 构造器只接受旧式 c2/CP 前缀；Hallmark/GO 使用可审计的 `REACTOME_PHASE11_ALIAS_*` 元数据别名通过校验，真实 pathway ID、collection、gene membership 全部外部冻结。别名不进入任何数学计算。

## 3. Human reference

- GSE266616 terminal ileum：16 Healthy + 11 Active Crohn donors；Epithelial、T、Myeloid 各 27 profiles，共 81 donor×cell-type pseudobulks。
- 36,591 个唯一 Human gene symbols；raw donor-celltype sums 后按 singIST SCE 路径执行 library-factor log normalization。
- Human-driven Hallmark：{len(hallmark_panel)}；注册 programs：{len(registered)}；去重后独立模型：{len(super_df)}。
- Global permutation P 值按 46 个模型统一 BH：FDR<0.10 validated={validated_n}，weak={weak_n}。因此下述 signed recapitulation 只能作为完整、可审计的探索性输出，**没有任何 pathway 可升级为核心生物学结论**。

## 4. Superpathway recapitulation

- 最大正值：{positive.iloc[0]['pathway']} ({positive.iloc[0]['recapitulation']:.2f}%)；其后为 {pathway_lines(positive.iloc[1:], 7)}。
- 精确为 0 的 pathway 共 {len(exact_zero)} 个：{pathway_lines(exact_zero)}
- 负值 {len(negative)} 个：{pathway_lines(negative)}
- >100% {len(over100)} 个：{pathway_lines(over100)}
- 所有值均为 singIST 原始百分比，没有将负值置零，也没有把 >100% 截断。

## 5. Cell-type recapitulation

- 最大：{cell_strong['pathway']} / {cell_strong['celltype']} = {cell_strong['recapitulation']:.2f}%。
- 最小：{cell_weak['pathway']} / {cell_weak['celltype']} = {cell_weak['recapitulation']:.2f}%。
- 逐 pathway×Level1 原始结果见 `CELLTYPE_RECAPITULATION.csv`；未对 T 炎症与代谢结果求平均形成“总体保真度”。

## 6. Gene contributions

- 主要 positive drivers：{driver_lines('POSITIVE')}
- 主要 negative drivers：{driver_lines('NEGATIVE')}
- 每个 pathway×cell type 的 gene contribution 加和均以 1e-8 容差核对 cell-type recapitulation。Contribution 不是新 DE，也没有被解释为跨物种显著性。

## 7. Epithelial

- TNF/NF-κB：super={recap_for('HALLMARK_TNFA_SIGNALING_VIA_NFKB')}，Epithelial={recap_for('HALLMARK_TNFA_SIGNALING_VIA_NFKB','Epithelial')}。
- mitochondrial respiration：super={recap_for('R-HSA-1428517')}，Epithelial={recap_for('R-HSA-1428517','Epithelial')}。
- fatty-acid metabolism：super={recap_for('R-HSA-8978868')}，Epithelial={recap_for('R-HSA-8978868','Epithelial')}。
- xenobiotic metabolism：super={recap_for('HALLMARK_XENOBIOTIC_METABOLISM')}，Epithelial={recap_for('HALLMARK_XENOBIOTIC_METABOLISM','Epithelial')}。
- cell cycle/G2M：super={recap_for('HALLMARK_G2M_CHECKPOINT')}，Epithelial={recap_for('HALLMARK_G2M_CHECKPOINT','Epithelial')}；本次并未复现 Phase10 的强 pathway baseline，且 `HIGH_REGION_RISK` / `REGION_CAUTION=YES` 不因 singIST 值而解除。

## 8. T

- 炎症侧：TNF/NF-κB={recap_for('HALLMARK_TNFA_SIGNALING_VIA_NFKB','T')}；Hallmark inflammatory={recap_for('HALLMARK_INFLAMMATORY_RESPONSE','T')}；IL2/STAT5={recap_for('HALLMARK_IL2_STAT5_SIGNALING','T')}；interleukin={recap_for('R-HSA-449147','T')}。
- 代谢/线粒体侧：Hallmark oxidative phosphorylation={recap_for('HALLMARK_OXIDATIVE_PHOSPHORYLATION','T')}；aerobic respiration={recap_for('GO:0009060','T')}；Reactome respiration={recap_for('R-HSA-1428517','T')}；translation={recap_for('R-HSA-72766','T')}。
- 结论边界：代谢/线粒体的负向信号明确强于炎症侧，但炎症主要接近 0、而非预期的强正向。因此支持“program-selective asymmetry”，**不支持把 T 描述为炎症正向保真**，更不能平均成“T 中等保真”。所有对应 Human models 仍为 weak。

## 9. Myeloid

- TNF/NF-κB={recap_for('HALLMARK_TNFA_SIGNALING_VIA_NFKB','Myeloid')}（全表最强 cell-type 值）；Hallmark inflammatory={recap_for('HALLMARK_INFLAMMATORY_RESPONSE','Myeloid')}；cytokine/JAK={recap_for('R-HSA-449147','Myeloid')}；TLR4={recap_for('R-HSA-166016','Myeloid')}；apoptosis={recap_for('HALLMARK_APOPTOSIS','Myeloid')}；respiration={recap_for('R-HSA-1428517','Myeloid')}；translation={recap_for('R-HSA-72766','Myeloid')}。
- IL6/JAK/STAT3 的结果只解释为 STAT3/JAK-related recapitulation，不升级为 IL-23-specific evidence。

## 10. Orthology coverage

- pathway×cell-type detected one-to-one coverage：{coverage['orthology_coverage'].min():.2f}%–{coverage['orthology_coverage'].max():.2f}%，中位数 {coverage['orthology_coverage'].median():.2f}%。
- coverage 分母为 Human observed features，分子为 Phase10 Ensembl104 1:1 且在相应 Dog Level1 中检测到的 Human genes；另保留 object-presence coverage。
- 正式运行的累计 warnings 经代表性 pathway 重放确认，来自 singIST 对缺少 one-to-one ortholog 的 genes 做预期排除（详见 `phase11_recap_warning_audit.log`），不是拟合失败、非有限值或回退到在线 BioMart。

## 11. Phase10 comparison

- 25 个注册 program 的逐项 traditional baseline 与 singIST 对照见 `PHASE11_VS_PHASE10_BASELINE.csv`。
- 一致/新增信息：T 的 respiration 在 Phase10 为 opposite，singIST T 也为负；Myeloid TNF 的 15.41% 揭示了被 superpathway 汇总值 -0.02% 掩盖的 block-specific 正贡献。
- 不一致：Phase10 的 TNF 三线 lineage pathway baseline 均为 concordant strong，但 singIST T=0、Epithelial=-5.49%；Phase10 Epithelial mitochondrial respiration 同向下调，而 singIST Epithelial=+5.12%；Epithelial G2M 的 Phase10 strong baseline 在 singIST 中近 0。
- 这些差异说明 singIST 提供的是 Human-model-weighted functional perturbation，不等同于 NES concordance。Phase10 gene rho 不作为 singIST 成败标准；两者不合成为最终 fidelity score。

## 12. Dog LODO

- 7 轮均只重算 Dog disease-model change；46 个 Human reference models 从未重拟合。
- `cie_3` 是本阶段按 17 pathways 总绝对变化最大的 donor（sum={largest_lodo_donor['sum_absolute_change']:.2f}，mean={largest_lodo_donor['mean_absolute_change']:.2f}，single-pathway max={largest_lodo_donor['max_absolute_change']:.2f}）；因此 `cie_4` **不再是 singIST LODO 的总体最大影响 donor**。
- `cie_4` 仍是 {cie4_n}/{len(lodo_summary)} 个 pathway 的最大影响 donor（GO:0009060、GO:0042254、KRAS signaling DN/UP）；`cie_2` 仅在 R-HSA-1428517 为最大影响。逐 pathway Full7、median/min/max、sign stability、最大绝对变化均已冻结。

## 13. Donor-aware sensitivity

- Track A official internal DE 的 FDR<0.05 gene 数为 Epithelial/T/Myeloid = 128/92/157；Track B 使用 Phase3 Full7 donor-level edgeR 后为 1/0/0。
- Track B 的 46 个 superpathway、138 个 cell-type recapitulation **全部精确为 0**。所谓同号 {int(fc_super['same_sign'].astype(bool).sum())}/{len(fc_super)} 仅对应 Track A 本身也为 0 的 9 个 pathway，不构成方向复制证据。
- 该完全衰减解释为 Dog n=7 donor-level FDR eligibility sparsity；它限制 Track A 的疾病推断可信度，但不能解释为犬完全不能模拟 Human。

## 14. FRE / treatment / region provenance

- FRE 仅作为 context，不参与训练或 score 修正。
- Human treatment structure 沿用 Phase7/10 provenance，本阶段未扩展为 5 套 Human reference 重拟合。
- Human reference 固定 TI；Dog 为 proximal duodenum。所有 `region_risk`、`phase10_priority`、treatment、FRE 字段均保留在 program summary。

## 15. 主要限制

- Dog n=7；singIST Track A 的内部 cell-level DE 机制不是正式犬疾病推断。
- donor-aware Track B 因 Phase3 FDR 稀疏而强烈衰减。
- Human treatment structure、duodenum vs TI、orthology coverage、Human reference model validity/依赖均限制解释。
- recapitulation 百分比可 <0 或 >100；Level2 未分析；singIST 不是独立外部验证。
- singIST 1.0.2 对现代 Hallmark/GO pathway metadata 的构造器限制通过纯元数据别名处理；gene sets 与数学未改。

## 16. Phase11 Gate

**YELLOW_PHASE11_SINGIST_COMPLETE_REFERENCE_MODEL_AND_DONOR_SENSITIVITY_LIMITATIONS**

技术流程、官方示例、Ensembl104 注入、三层输出、加和门控、Track B、7 轮 LODO 与 scope audit 均完成。YELLOW 的主因是 46/46 Human reference models 均为 weak、Track B 全零、signed recapitulation 混合、最低 coverage 35% 与 donor sensitivity；没有通过删模型、放宽 FDR 或改阈值掩盖这些限制。

建议：技术上可以把项目移交给 Phase12 SATURN 作为**独立、单独门控的探索性分析**，但不能把 Phase11 当作支持核心生物学结论的 GREEN 证据，也不能让 SATURN 继承一个预设的“犬已高保真”结论。

Phase11 在此停止；没有生成 Phase12 SATURN、Phase13 CONCORD、joint embedding、Level2 正式跨物种比较或最终 fidelity score。
"""
REPORT.write_text(report, encoding="utf-8", newline="\n")

LOG_SUM.write_text(
    "PHASE11_SUMMARY_COMPLETE\n" + json.dumps(metrics, ensure_ascii=False, indent=2) + "\n",
    encoding="utf-8")
LOG_REP.write_text(f"PHASE11_REPORT_COMPLETE\nreport={REPORT}\ngate={gate}\n",
                   encoding="utf-8")
print(json.dumps(metrics, ensure_ascii=False, indent=2))
print(f"REPORT_COMPLETE gate={gate} path={REPORT}")
