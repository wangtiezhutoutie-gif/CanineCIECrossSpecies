from pathlib import Path
import hashlib, json, re, sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
OUT = ROOT / "WORKING" / "22_DOG_AUTHOR_TCELL_CROSS_SPECIES"
AOUT = OUT / "01_AUTHOR_ANNOTATION_AUDIT"
BOUT = OUT / "02_DOG_SUBTYPE_ANALYSIS"
COUT = OUT / "03_DOG_HUMAN_CORRESPONDENCE"
WORK = OUT / "work"
for p in (AOUT, BOUT, COUT, WORK):
    p.mkdir(parents=True, exist_ok=True)

H5AD = ROOT / "05_annotations" / "dog_GSE254005_annotated.h5ad"
META = WORK / "author_metadata.csv"
ORTHO = ROOT / "10_cross_species_baseline" / "mapping" / "DOG_HUMAN_ORTHOLOG_MASTER.csv"
AUTHOR_SCRIPT = WORK / "author_3_tcell_analysis.R"

def sha256(path, block=1024*1024):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(block), b""):
            h.update(b)
    return h.hexdigest()

def pick(df, names, required=True):
    for x in names:
        if x in df.columns:
            return x
    if required:
        raise RuntimeError("Missing required metadata column: " + "|".join(names))
    return None

def clean(x):
    if pd.isna(x):
        return ""
    return str(x).strip()

def barcode(x):
    m = re.search(r"([ACGT]{16}-1)(?:_\d+)?$", clean(x).upper())
    return m.group(1) if m else ""

A = sc.read_h5ad(H5AD)
O = A.obs.copy()
for c in ["donor_id", "condition", "level1", "level2", "annotation_cluster"]:
    O[c] = O[c].astype(str)
O["frozen_cell_id"] = O.index.astype(str)
O["barcode"] = O["frozen_cell_id"].map(barcode)
O["donor_id"] = O["donor_id"].astype(str)
F = O.loc[O["level1"].eq("T")].copy()

am = pd.read_csv(META, low_memory=False)
id_col = pick(am, ["author_cell_id"])
sample_col = pick(am, ["sample_name", "name2", "orig.ident"])
condition_col = pick(am, ["cellSource", "condition", "disease"])
l1_col = pick(am, ["celltype.l1", "majorID_sub_big"])
l2_col = pick(am, ["celltype.l2", "majorID_sub_noILC"])
cluster_col = pick(am, ["clusterID_sub", "clusterID_final_noILC", "seurat_clusters"], required=False)
orig_col = pick(am, ["orig.ident"], required=False)
am["author_cell_id"] = am[id_col].map(clean)
am["author_sample"] = am[sample_col].map(clean)
am["author_condition"] = am[condition_col].map(clean)
am["author_celltype_l1"] = am[l1_col].map(clean)
am["author_celltype_l2"] = am[l2_col].map(clean)
am["author_cluster"] = am[cluster_col].map(clean) if cluster_col else ""
am["author_orig_ident"] = am[orig_col].map(clean) if orig_col else ""
am["barcode"] = am["author_cell_id"].map(barcode)
if am["barcode"].eq("").any():
    bad = int(am["barcode"].eq("").sum())
    raise RuntimeError(f"Unable to parse 16bp 10x barcode for {bad} author cells")

frozen_sets = {d: set(g["barcode"]) for d, g in O.groupby("donor_id", observed=True)}
map_rows = []
sample_map = {}
for s, g in am.groupby("author_sample", observed=True):
    aset = set(g["barcode"])
    scores = sorted([(d, len(aset & bs)) for d, bs in frozen_sets.items()], key=lambda z: (-z[1], z[0]))
    best_d, best_n = scores[0]
    second_n = scores[1][1] if len(scores) > 1 else 0
    unique = best_n > 0 and best_n > second_n
    sample_map[s] = best_d if unique else ""
    map_rows.append({
        "author_sample": s,
        "author_cells": len(g),
        "inferred_frozen_donor": best_d if unique else "",
        "best_barcode_overlap": best_n,
        "second_best_overlap": second_n,
        "mapping_status": "UNIQUE_BARCODE_OVERLAP" if unique else "AMBIGUOUS_OR_NO_OVERLAP",
        "all_overlap_counts": ";".join(f"{d}:{n}" for d, n in scores),
    })
sample_map_df = pd.DataFrame(map_rows)
sample_map_df.to_csv(WORK / "author_sample_to_frozen_donor.csv", index=False)

am["inferred_frozen_donor"] = am["author_sample"].map(sample_map).fillna("")
am["join_key"] = am["inferred_frozen_donor"] + "|" + am["barcode"]
dup_author_keys = set(am.loc[am["join_key"].duplicated(False), "join_key"])
flookup = {}
dup_frozen_keys = set()
for idx, r in O.iterrows():
    k = f"{r['donor_id']}|{r['barcode']}"
    if k in flookup:
        dup_frozen_keys.add(k)
    else:
        flookup[k] = idx

rows = []
matched_frozen = set()
for r in am.itertuples(index=False):
    k = r.join_key
    base = {
        "author_cell_id": r.author_cell_id,
        "author_sample": r.author_sample,
        "author_condition": r.author_condition,
        "author_celltype_l1": r.author_celltype_l1,
        "author_celltype_l2": r.author_celltype_l2,
        "author_cluster": r.author_cluster,
        "parsed_barcode": r.barcode,
        "inferred_frozen_donor": r.inferred_frozen_donor,
        "frozen_cell_id": "",
        "frozen_donor_id": "",
        "frozen_condition": "",
        "frozen_level1": "",
        "frozen_level2": "",
        "frozen_annotation_cluster": "",
        "match_status": "AUTHOR_ONLY",
        "identity_note": "",
    }
    if not r.inferred_frozen_donor:
        base["match_status"] = "AMBIGUOUS"
        base["identity_note"] = "Author sample could not be uniquely assigned by barcode overlap"
    elif k in dup_author_keys or k in dup_frozen_keys:
        base["match_status"] = "AMBIGUOUS"
        base["identity_note"] = "Non-unique donor-plus-barcode key"
    elif k in flookup:
        fid = flookup[k]
        fr = O.loc[fid]
        base.update({
            "frozen_cell_id": str(fid),
            "frozen_donor_id": str(fr["donor_id"]),
            "frozen_condition": str(fr["condition"]),
            "frozen_level1": str(fr["level1"]),
            "frozen_level2": str(fr["level2"]),
            "frozen_annotation_cluster": str(fr["annotation_cluster"]),
        })
        ac = clean(r.author_condition).lower()
        fc = clean(fr["condition"]).lower()
        cond_ok = (("health" in ac and "health" in fc) or ("cie" in ac and "cie" in fc))
        if not cond_ok:
            base["match_status"] = "IDENTITY_CONFLICT"
            base["identity_note"] = "Author and frozen disease groups disagree"
        elif r.author_cell_id == str(fid):
            base["match_status"] = "EXACT_MATCH"
            matched_frozen.add(str(fid))
        else:
            base["match_status"] = "VERIFIED_TRANSFORMED_MATCH"
            base["identity_note"] = "Unique author-sample to donor mapping plus identical 16bp barcode"
            matched_frozen.add(str(fid))
    rows.append(base)

for fid, fr in F.iterrows():
    if str(fid) not in matched_frozen:
        rows.append({
            "author_cell_id": "", "author_sample": "", "author_condition": "",
            "author_celltype_l1": "", "author_celltype_l2": "", "author_cluster": "",
            "parsed_barcode": str(fr["barcode"]), "inferred_frozen_donor": "",
            "frozen_cell_id": str(fid), "frozen_donor_id": str(fr["donor_id"]),
            "frozen_condition": str(fr["condition"]), "frozen_level1": str(fr["level1"]),
            "frozen_level2": str(fr["level2"]),
            "frozen_annotation_cluster": str(fr["annotation_cluster"]),
            "match_status": "FROZEN_ONLY", "identity_note": "No author-object counterpart",
        })
cw = pd.DataFrame(rows)
cw.to_csv(AOUT / "02_AUTHOR_FROZEN_CELL_CROSSWALK.csv", index=False)
shared_all = cw[cw["match_status"].isin(["EXACT_MATCH", "VERIFIED_TRANSFORMED_MATCH"])].copy()
shared = shared_all[shared_all["frozen_level1"].eq("T") & shared_all["author_celltype_l2"].ne("")].copy()

def crosstab_rows(df, a, b, table_type):
    ct = pd.crosstab(df[a], df[b])
    out = []
    for av in ct.index:
        for bv in ct.columns:
            n = int(ct.loc[av, bv])
            if n == 0:
                continue
            sub = df[df[a].eq(av) & df[b].eq(bv)]
            out.append({
                "table_type": table_type, "author_label": av, "frozen_label": bv,
                "shared_cell_count": n,
                "row_percentage": 100.0 * n / int(ct.loc[av].sum()),
                "column_percentage": 100.0 * n / int(ct[bv].sum()),
                "healthy_cells": int(sub["frozen_condition"].eq("Healthy").sum()),
                "cie_cells": int(sub["frozen_condition"].eq("CIE").sum()),
                "donor_distribution": ";".join(f"{d}:{v}" for d, v in sub["frozen_donor_id"].value_counts().sort_index().items()),
            })
    return out

ctrows = []
ctrows += crosstab_rows(shared_all, "author_celltype_l1", "frozen_level1", "AUTHOR_L1_X_FROZEN_LEVEL1")
ctrows += crosstab_rows(shared, "author_celltype_l2", "frozen_level2", "AUTHOR_L2_X_FROZEN_LEVEL2")
ctrows += crosstab_rows(shared, "author_celltype_l2", "frozen_annotation_cluster", "AUTHOR_L2_X_FROZEN_CLUSTER")
pd.DataFrame(ctrows).to_csv(AOUT / "04_AUTHOR_FROZEN_LABEL_CROSSTAB.csv", index=False)

# Component A donor counts, targeted marker audit, and annotation-only eligibility.
donors = sorted(F["donor_id"].unique())
donor_condition = F.groupby("donor_id", observed=True)["condition"].first().to_dict()
subtypes = sorted(shared["author_celltype_l2"].dropna().unique())

count_rows = []
for st in subtypes:
    ss = shared[shared["author_celltype_l2"].eq(st)]
    for d in donors:
        z = ss[ss["frozen_donor_id"].eq(d)]
        count_rows.append({
            "author_subtype": st, "donor_id": d, "condition": donor_condition[d],
            "cell_count": int(len(z)),
            "t8_cells": int(z["frozen_annotation_cluster"].eq("T_8").sum()),
            "t9_cells": int(z["frozen_annotation_cluster"].eq("T_9").sum()),
            "matched_statuses": ";".join(sorted(z["match_status"].unique())) if len(z) else "",
        })
donor_counts = pd.DataFrame(count_rows)
donor_counts.to_csv(AOUT / "03_AUTHOR_SUBTYPE_DONOR_COUNTS.csv", index=False)

symbols = A.var["gene_symbol"].astype(str).str.upper().to_numpy()
gene_to_idx = {}
for i, g in enumerate(symbols):
    gene_to_idx.setdefault(g, []).append(i)
counts_layer = A.layers["counts"]
log_layer = A.layers["lognorm"] if "lognorm" in A.layers else None
fid_to_pos = {str(x): i for i, x in enumerate(A.obs_names.astype(str))}
shared["frozen_position"] = shared["frozen_cell_id"].map(fid_to_pos)
if shared["frozen_position"].isna().any():
    raise RuntimeError("Matched frozen cell IDs were not found in h5ad")

marker_panels = {
    "CD8_eff (c0)": ["CD8A", "GZMA", "GZMB", "CCL4", "FASLG"],
    "Tnaive (c1)": ["TCF7", "IL7R", "CCR7", "LTB"],
    "CD8_TRM (c2)": ["CD8A", "ITGAE", "CD69", "CXCR6"],
    "gdT_1 (c3)": ["TRDC", "CD3E", "CD247"],
    "NK_T (c4)": ["CD3E", "NKG7", "KLRD1", "NCR3"],
    "CD8_mem (c5)": ["CD8A", "CCL4", "GZMK", "TCF7"],
    "gdT_2 (c6)": ["TRDC", "CD3E", "CD247"],
    "Treg (c7)": ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "ICOS", "IKZF2"],
    "NK (c8)": ["NKG7", "KLRD1", "NCR3", "CD3E"],
    "T_IFN (c9)": ["ISG15", "OAS1", "IFI44", "IFNG"],
}
marker_rows = []
all_t_pos = shared["frozen_position"].astype(int).to_numpy()
for st in subtypes:
    spos = shared.loc[shared["author_celltype_l2"].eq(st), "frozen_position"].astype(int).to_numpy()
    opos = shared.loc[~shared["author_celltype_l2"].eq(st), "frozen_position"].astype(int).to_numpy()
    for gene in marker_panels.get(st, []):
        ix = gene_to_idx.get(gene, [])
        if len(ix) != 1:
            marker_rows.append({
                "species": "dog", "subtype": st, "gene": gene,
                "mapping_status": "NOT_UNIQUELY_MAPPED", "n_subtype_cells": len(spos),
                "subtype_detection_fraction": np.nan, "other_t_detection_fraction": np.nan,
                "subtype_mean_lognorm": np.nan, "other_t_mean_lognorm": np.nan,
                "mean_lognorm_difference": np.nan,
                "evidence_scope": "Targeted predefined annotation marker; same-dataset support only",
            })
            continue
        j = ix[0]
        xv = np.asarray(counts_layer[spos, j].todense()).ravel() if sp.issparse(counts_layer) else np.asarray(counts_layer[spos, j]).ravel()
        yv = np.asarray(counts_layer[opos, j].todense()).ravel() if sp.issparse(counts_layer) else np.asarray(counts_layer[opos, j]).ravel()
        if log_layer is not None:
            xl = np.asarray(log_layer[spos, j].todense()).ravel() if sp.issparse(log_layer) else np.asarray(log_layer[spos, j]).ravel()
            yl = np.asarray(log_layer[opos, j].todense()).ravel() if sp.issparse(log_layer) else np.asarray(log_layer[opos, j]).ravel()
        else:
            xl, yl = np.log1p(xv), np.log1p(yv)
        marker_rows.append({
            "species": "dog", "subtype": st, "gene": gene,
            "mapping_status": "UNIQUE_GENE_SYMBOL", "n_subtype_cells": len(spos),
            "subtype_detection_fraction": float((xv > 0).mean()) if len(xv) else np.nan,
            "other_t_detection_fraction": float((yv > 0).mean()) if len(yv) else np.nan,
            "subtype_mean_lognorm": float(xl.mean()) if len(xl) else np.nan,
            "other_t_mean_lognorm": float(yl.mean()) if len(yl) else np.nan,
            "mean_lognorm_difference": float(xl.mean() - yl.mean()) if len(xl) and len(yl) else np.nan,
            "evidence_scope": "Targeted predefined annotation marker; same-dataset support only",
        })
marker_df = pd.DataFrame(marker_rows)
marker_df.to_csv(WORK / "dog_targeted_marker_audit.csv", index=False)

elig_rows = []
for st in subtypes:
    dc = donor_counts[donor_counts["author_subtype"].eq(st)].copy()
    h_ge = int(((dc["condition"] == "Healthy") & (dc["cell_count"] >= 20)).sum())
    c_ge = int(((dc["condition"] == "CIE") & (dc["cell_count"] >= 20)).sum())
    total = int(dc["cell_count"].sum())
    max_i = dc["cell_count"].idxmax()
    max_donor = dc.loc[max_i, "donor_id"]
    max_frac = float(dc.loc[max_i, "cell_count"] / total) if total else np.nan
    ss = shared[shared["author_celltype_l2"].eq(st)].copy()
    t8_n = int(ss["frozen_annotation_cluster"].eq("T_8").sum())
    t9_n = int(ss["frozen_annotation_cluster"].eq("T_9").sum())
    t8_frac = t8_n / total if total else 0.0
    t9_frac = t9_n / total if total else 0.0
    comp_flags = []
    gate_dependency = []
    for cl, cg in ss.groupby("frozen_annotation_cluster", observed=True):
        vc = cg["frozen_donor_id"].value_counts()
        dom = float(vc.max() / vc.sum())
        if dom > 0.80:
            comp_flags.append(f"{cl}:{vc.index[0]}:{dom:.4f}")
            reduced = ss[~ss["frozen_annotation_cluster"].eq(cl)]
            rdc = reduced.groupby(["frozen_donor_id", "frozen_condition"], observed=True).size().reset_index(name="n")
            rh = int(((rdc["frozen_condition"] == "Healthy") & (rdc["n"] >= 20)).sum())
            rc = int(((rdc["frozen_condition"] == "CIE") & (rdc["n"] >= 20)).sum())
            if h_ge >= 3 and c_ge >= 3 and (rh < 3 or rc < 3):
                gate_dependency.append(str(cl))
    identity_note = "Official author celltype.l2 label recovered from T-cell object"
    if st == "NK (c8)":
        cls = "IDENTITY_CONFLICT"
        identity_note = "Author NK label conflicts with restriction to frozen Level1 T cells"
    elif h_ge < 3 or c_ge < 3:
        cls = "INSUFFICIENT_COVERAGE"
    elif max_frac > 0.50:
        cls = "DONOR_CONFOUNDED"
    elif gate_dependency:
        cls = "DONOR_CONFOUNDED"
        identity_note += "; donor-dominant component is required to pass the donor gate"
    elif st in {"Treg (c7)", "NK_T (c4)", "T_IFN (c9)"}:
        cls = "TRACEABLE_BUT_ANNOTATION_LIMITED"
    else:
        cls = "TRACEABLE_AND_INTERPRETABLE"
    formal = cls in {"TRACEABLE_AND_INTERPRETABLE", "TRACEABLE_BUT_ANNOTATION_LIMITED"} and h_ge >= 3 and c_ge >= 3
    if st == "NK_T (c4)":
        formal = False
        identity_note += "; mixed NK/T identity not approved for T-state pseudobulk"
    elig_rows.append({
        "author_subtype": st, "total_shared_frozen_t_cells": total,
        "healthy_donors_ge20": h_ge, "cie_donors_ge20": c_ge,
        "minimum_cells_across_seven_donors": int(dc["cell_count"].min()),
        "median_cells_across_seven_donors": float(dc["cell_count"].median()),
        "maximum_cells_across_seven_donors": int(dc["cell_count"].max()),
        "maximum_donor": max_donor, "maximum_donor_fraction": max_frac,
        "one_donor_gt50_flag": "YES" if max_frac > 0.50 else "NO",
        "t8_cells": t8_n, "t8_fraction_of_subtype": t8_frac,
        "t9_cells": t9_n, "t9_fraction_of_subtype": t9_frac,
        "component_gt80_flags": ";".join(comp_flags) if comp_flags else "NONE",
        "gate_dependent_donor_dominant_components": ";".join(gate_dependency) if gate_dependency else "NONE",
        "annotation_classification": cls, "annotation_evidence_note": identity_note,
        "formal_pseudobulk_eligible": "YES" if formal else "NO",
        "gate_reason": "PASS_3H_3CIE_TRACEABLE_ESTIMABLE" if formal else cls,
    })
elig = pd.DataFrame(elig_rows)
elig.to_csv(AOUT / "06_AUTHOR_ANNOTATION_ELIGIBILITY.csv", index=False)
elig.to_csv(BOUT / "01_SUBTYPE_ELIGIBILITY.csv", index=False)

# Freeze the verified shared-cell analysis universe and create raw-count pseudobulks.
for r in elig.itertuples(index=False):
    if r.formal_pseudobulk_eligible != "YES":
        continue
    st = r.author_subtype
    slug = re.sub(r"[^A-Za-z0-9]+", "_", st).strip("_").lower()
    ss = shared[shared["author_celltype_l2"].eq(st)]
    sample_rows, cols = [], []
    for d in donors:
        z = ss[ss["frozen_donor_id"].eq(d)]
        if len(z) < 20:
            continue
        pos = z["frozen_position"].astype(int).to_numpy()
        v = np.asarray(counts_layer[pos].sum(axis=0)).ravel().astype(np.int64)
        cols.append(v)
        sample_rows.append({
            "sample_id": d, "donor_id": d, "condition": donor_condition[d],
            "cell_count": int(len(z)), "library_size": int(v.sum()),
            "author_subtype": st, "analysis_universe": "VERIFIED_AUTHOR_FROZEN_SHARED_T",
        })
    mat = np.column_stack(cols)
    tab = pd.DataFrame(mat, columns=[x["sample_id"] for x in sample_rows])
    tab.insert(0, "gene_symbol", A.var["gene_symbol"].astype(str).to_numpy())
    tab.insert(0, "feature_id", A.var["feature_id"].astype(str).to_numpy())
    tab.to_csv(WORK / f"counts_{slug}.tsv", sep="\t", index=False)
    pd.DataFrame(sample_rows).to_csv(WORK / f"samples_{slug}.tsv", sep="\t", index=False)

# Donor-level composition in the author compartment and shared frozen-T universe.
comp_rows = []
for st in subtypes:
    for d in donors:
        a_den = int((shared_all["frozen_donor_id"].eq(d) & shared_all["author_celltype_l2"].ne("")).sum())
        a_num = int((shared_all["frozen_donor_id"].eq(d) & shared_all["author_celltype_l2"].eq(st)).sum())
        t_den = int(shared["frozen_donor_id"].eq(d).sum())
        t_num = int((shared["frozen_donor_id"].eq(d) & shared["author_celltype_l2"].eq(st)).sum())
        frozen_t_n = int(F["donor_id"].eq(d).sum())
        for universe, num, den in [
            ("AUTHOR_ANNOTATED_MATCHED_COMPARTMENT", a_num, a_den),
            ("SHARED_AUTHOR_FROZEN_T_UNIVERSE", t_num, t_den),
        ]:
            comp_rows.append({
                "author_subtype": st, "donor_id": d, "condition": donor_condition[d],
                "universe": universe, "subtype_cell_count": num,
                "denominator_cell_count": den,
                "proportion": num / den if den else np.nan,
                "author_annotation_coverage_of_frozen_t": t_den / frozen_t_n if frozen_t_n else np.nan,
            })
comp = pd.DataFrame(comp_rows)
summ = comp.groupby(["author_subtype", "universe"], observed=True).apply(
    lambda g: pd.Series({
        "healthy_median_proportion": g.loc[g["condition"].eq("Healthy"), "proportion"].median(),
        "cie_median_proportion": g.loc[g["condition"].eq("CIE"), "proportion"].median(),
    }), include_groups=False).reset_index()
summ["cie_minus_healthy_median_difference"] = summ["cie_median_proportion"] - summ["healthy_median_proportion"]
comp = comp.merge(summ, on=["author_subtype", "universe"], how="left")
comp.to_csv(BOUT / "02_DONOR_COMPOSITION.csv", index=False)

lodo_rows = []
for st in subtypes:
    for universe in comp["universe"].unique():
        base = comp[(comp["author_subtype"] == st) & (comp["universe"] == universe)]
        for omit in donors:
            z = base[base["donor_id"] != omit]
            hm = z.loc[z["condition"].eq("Healthy"), "proportion"].median()
            cm = z.loc[z["condition"].eq("CIE"), "proportion"].median()
            diff = cm - hm
            lodo_rows.append({
                "author_subtype": st, "universe": universe, "omitted_donor": omit,
                "omitted_condition": donor_condition[omit],
                "healthy_donors_retained": int(z["condition"].eq("Healthy").sum()),
                "cie_donors_retained": int(z["condition"].eq("CIE").sum()),
                "healthy_median_proportion": hm, "cie_median_proportion": cm,
                "cie_minus_healthy_median_difference": diff,
                "direction": "CIE_HIGHER" if diff > 0 else ("CIE_LOWER" if diff < 0 else "TIE"),
                "reduced_healthy_n2_caution": "YES" if donor_condition[omit] == "Healthy" else "NO",
            })
pd.DataFrame(lodo_rows).to_csv(BOUT / "03_COMPOSITION_LODO.csv", index=False)

# Component C matching decisions are frozen before any new metabolic results are loaded.
human_rel_path = ROOT / "WORKING" / "20_HUMAN_TCELL_LABEL_FREEZE" / "06_SUBTYPE_ANNOTATION_RELIABILITY.csv"
human_audit_path = ROOT / "WORKING" / "21_HUMAN_TCELL_SUBTYPE_METABOLIC_ANALYSIS" / "06_HUMAN_TCELL_ANNOTATION_AND_TREATMENT_AUDIT.csv"
human_rel = pd.read_csv(human_rel_path)
human_audit = pd.read_csv(human_audit_path, low_memory=False)

ortho = pd.read_csv(ORTHO, low_memory=False)
primary = ortho[ortho["primary_one2one"].astype(str).str.lower().eq("true")].copy()
dog_to_human = dict(zip(primary["dog_original_symbol"].astype(str).str.upper(), primary["human_ortholog_symbol"].astype(str).str.upper()))

marker_evidence_rows = []
for r in marker_df.itertuples(index=False):
    human_symbol = dog_to_human.get(str(r.gene).upper(), "")
    marker_evidence_rows.append({
        "evidence_type": "DOG_TARGETED_MARKER",
        "dog_subtype": r.subtype, "human_subtype": "",
        "dog_gene": r.gene, "human_ortholog_symbol": human_symbol,
        "strict_one_to_one_ortholog": "YES" if human_symbol else "NO",
        "dog_detection_fraction": r.subtype_detection_fraction,
        "dog_other_t_detection_fraction": r.other_t_detection_fraction,
        "dog_mean_lognorm_difference": r.mean_lognorm_difference,
        "human_evidence": "",
        "evidence_boundary": r.evidence_scope,
    })
for r in human_rel.itertuples(index=False):
    marker_evidence_rows.append({
        "evidence_type": "HUMAN_AUTHOR_ANNOTATION",
        "dog_subtype": "", "human_subtype": r.author_subtype,
        "dog_gene": "", "human_ortholog_symbol": "",
        "strict_one_to_one_ortholog": "NOT_APPLICABLE",
        "dog_detection_fraction": np.nan, "dog_other_t_detection_fraction": np.nan,
        "dog_mean_lognorm_difference": np.nan,
        "human_evidence": str(r.author_annotation_evidence),
        "evidence_boundary": str(r.identity_ambiguity_or_boundary),
    })
ht = human_audit[human_audit["record_type"].eq("TREG_MARKER")].copy()
for r in ht.itertuples(index=False):
    marker_evidence_rows.append({
        "evidence_type": "HUMAN_TREG_TARGETED_MARKER",
        "dog_subtype": "", "human_subtype": "Treg",
        "dog_gene": "", "human_ortholog_symbol": r.marker,
        "strict_one_to_one_ortholog": "NOT_APPLICABLE",
        "dog_detection_fraction": np.nan, "dog_other_t_detection_fraction": np.nan,
        "dog_mean_lognorm_difference": np.nan,
        "human_evidence": f"detection={r.treg_detection_fraction}; other_T={r.other_T_detection_fraction}; direction_support={r.marker_direction_support}",
        "evidence_boundary": "Predefined same-dataset marker support; not an independent reannotation",
    })
marker_ortholog = pd.DataFrame(marker_evidence_rows)
marker_ortholog.to_csv(COUT / "02_MARKER_AND_ORTHOLOG_EVIDENCE.csv", index=False)

elig_map = elig.set_index("author_subtype").to_dict("index")
m_by = marker_df[marker_df["mapping_status"].eq("UNIQUE_GENE_SYMBOL")].copy()

def dog_marker_support(st, genes, minimum):
    z = m_by[(m_by["subtype"] == st) & (m_by["gene"].isin(genes))]
    good = z[(z["subtype_detection_fraction"] >= 0.005) & (z["mean_lognorm_difference"] > 0)]
    return int(len(good)), sorted(good["gene"].tolist()), int(minimum)

candidates = [
    ("CD8_TRM (c2)", "CD103+ CD8", "resident CD8 state", ["CD8A", "ITGAE", "CD69", "CXCR6"], 2),
    ("Treg (c7)", "Treg", "regulatory-like state", ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "ICOS", "IKZF2"], 3),
    ("CD8_eff (c0)", "KLRG1+ CD8", "effector/cytotoxic analogy", ["CD8A", "GZMA", "GZMB", "CCL4"], 2),
    ("CD8_mem (c5)", "KLRG1+ CD8", "CD8 memory versus KLRG1-positive state", ["CD8A", "CCL4", "GZMK", "TCF7"], 2),
    ("Tnaive (c1)", "CD103- CD4", "non-resident conventional T analogy", ["TCF7", "IL7R", "CCR7"], 2),
    ("gdT_1 (c3)", "", "gamma-delta T state", ["TRDC", "CD3E"], 2),
    ("gdT_2 (c6)", "", "gamma-delta T state", ["TRDC", "CD3E"], 2),
    ("NK_T (c4)", "", "mixed NK/T state", ["CD3E", "NKG7"], 2),
    ("NK (c8)", "", "NK state", ["NKG7", "KLRD1"], 2),
    ("T_IFN (c9)", "", "interferon activation state", ["ISG15", "OAS1", "IFI44"], 2),
]
match_rows = []
freeze_time = datetime.now(timezone.utc).isoformat()
for dog_st, human_st, identity_basis, genes, min_support in candidates:
    if dog_st not in elig_map:
        continue
    e = elig_map[dog_st]
    ngood, goodgenes, needed = dog_marker_support(dog_st, genes, min_support)
    formal_dog = e["formal_pseudobulk_eligible"] == "YES"
    human_ok = False
    human_boundary = "NO_NAMED_HUMAN_COUNTERPART"
    if human_st:
        hz = human_rel[human_rel["author_subtype"].eq(human_st)]
        human_ok = len(hz) == 1 and str(hz.iloc[0]["phase5b_eligibility"]).startswith("ELIGIBLE")
        human_boundary = str(hz.iloc[0]["identity_ambiguity_or_boundary"]) if len(hz) else "HUMAN_STATE_NOT_FOUND"
    if dog_st == "CD8_TRM (c2)" and human_ok and formal_dog and ngood >= needed:
        decision = "SUPPORTED_EXPLORATORY_MATCH"
        rationale = "Independent author labels support resident CD8 identity in both species; multiple predefined canine markers and donor gates support an exploratory match."
    elif dog_st == "Treg (c7)" and human_ok and formal_dog and ngood >= needed:
        decision = "SUPPORTED_EXPLORATORY_MATCH"
        rationale = "Both are original-author Treg labels with multi-marker regulatory-state support and valid donor-level estimates; canonical identity remains annotation-limited."
    elif human_st and ngood >= max(1, needed - 1):
        decision = "PARTIAL_FUNCTIONAL_ANALOGY"
        rationale = "Broad functional identity overlaps, but author terminology, marker definition, residency or donor evidence does not establish the same cellular state."
    elif human_st:
        decision = "NOT_COMPARABLE"
        rationale = "Named states lack sufficient multi-marker and/or donor-supported evidence for cellular equivalence."
    else:
        decision = "NOT_ASSESSABLE"
        rationale = "No directly corresponding frozen human author-defined candidate state is available."
    if not formal_dog and decision == "SUPPORTED_EXPLORATORY_MATCH":
        decision = "PARTIAL_FUNCTIONAL_ANALOGY"
        rationale += " Dog donor-level disease estimation did not pass the frozen gate."
    match_rows.append({
        "dog_author_subtype": dog_st, "human_author_subtype": human_st,
        "matching_classification": decision, "identity_basis": identity_basis,
        "dog_supporting_markers": ";".join(goodgenes) if goodgenes else "NONE",
        "dog_marker_support_count": ngood, "dog_marker_support_required": needed,
        "dog_annotation_classification": e["annotation_classification"],
        "dog_formal_pseudobulk_eligible": e["formal_pseudobulk_eligible"],
        "dog_healthy_donors_ge20": e["healthy_donors_ge20"],
        "dog_cie_donors_ge20": e["cie_donors_ge20"],
        "human_phase5b_eligible": "YES" if human_ok else "NO",
        "human_annotation_boundary": human_boundary,
        "matching_rationale": rationale,
        "tissue_boundary": "Dog duodenum versus human terminal ileum",
        "decision_basis_excludes_metabolic_outcomes": "YES",
        "freeze_timestamp_utc": freeze_time,
    })
matches = pd.DataFrame(match_rows)
match_path = COUT / "01_MATCHING_DECISIONS_FROZEN.csv"
matches.to_csv(match_path, index=False)
freeze_record = {
    "path": str(match_path), "sha256": sha256(match_path),
    "freeze_timestamp_utc": freeze_time,
    "metabolic_outcomes_loaded_before_freeze": False,
    "supported_pairs": int(matches["matching_classification"].eq("SUPPORTED_EXPLORATORY_MATCH").sum()),
}
(WORK / "matching_decision_freeze.json").write_text(json.dumps(freeze_record, indent=2), encoding="utf-8")

# Program compatibility is established without using disease effects.
dog_react = pd.read_csv(ROOT / "06_pseudobulk" / "gsea" / "gene_sets" / "reactome_canis_familiaris_gene_sets.tsv.gz", sep="\t")
dog_hall = pd.read_csv(ROOT / "06_pseudobulk" / "gsea" / "gene_sets" / "hallmark_dog_ortholog_projected_membership.tsv.gz", sep="\t")
dog_members_all = pd.concat([dog_react, dog_hall], ignore_index=True)
feature_map = pd.read_csv(ROOT / "06_pseudobulk" / "gsea" / "mapping" / "phase4_feature_to_geneset_id.csv", low_memory=False)
feature_map["target_geneset_id"] = feature_map["target_geneset_id"].astype(str)
dog_members_all["target_geneset_id"] = dog_members_all["target_geneset_id"].astype(str)
human_members = pd.read_csv(ROOT / "WORKING" / "21_HUMAN_TCELL_SUBTYPE_METABOLIC_ANALYSIS" / "work" / "frozen_target_memberships.tsv", sep="\t")
program_pairs = [
    ("R-CFA-1428517", "R-HSA-1428517", "Aerobic respiration and respiratory electron transport", "Reactome"),
    ("R-CFA-9937383", "R-HSA-9937383", "Mitochondrial ribosome-associated quality control", "Reactome"),
    ("HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_XENOBIOTIC_METABOLISM", "Xenobiotic Metabolism", "Hallmark"),
]
compat_rows, shared_member_rows = [], []
primary2 = primary[["dog_original_feature_id", "human_ortholog_symbol", "mapping_confidence"]].drop_duplicates()
for dog_pid, human_pid, pname, coll in program_pairs:
    dm = dog_members_all[dog_members_all["pathway_id"].astype(str).eq(dog_pid)]
    ids = set(dm["target_geneset_id"].astype(str))
    fm = feature_map[feature_map["target_geneset_id"].isin(ids)][["feature_id"]].drop_duplicates()
    mo = fm.merge(primary2, left_on="feature_id", right_on="dog_original_feature_id", how="inner")
    hm = human_members[human_members["pathway_id"].astype(str).eq(human_pid)]
    hgenes = set(hm["gene_symbol"].astype(str).str.upper())
    mapped_h = set(mo["human_ortholog_symbol"].astype(str).str.upper())
    shared_genes = sorted(mapped_h & hgenes)
    for g in shared_genes:
        shared_member_rows.append({
            "dog_pathway_id": dog_pid, "human_pathway_id": human_pid,
            "pathway_name": pname, "human_gene_symbol": g,
        })
    compat_rows.append({
        "dog_pathway_id": dog_pid, "human_pathway_id": human_pid,
        "pathway_name": pname, "collection": coll,
        "dog_native_members": int(len(ids)),
        "dog_members_with_strict_one2one_ortholog": int(len(mapped_h)),
        "human_native_members": int(len(hgenes)),
        "shared_strict_one2one_pathway_members": int(len(shared_genes)),
        "compatibility_status": "STRICT_ORTHOLOG_MATCHED_ANALYSIS_POSSIBLE" if len(shared_genes) >= 15 else "INSUFFICIENT_SHARED_PATHWAY_MEMBERS",
        "comparison_boundary": "Native species memberships differ; matched analysis uses only shared one-to-one ortholog members.",
    })
compat = pd.DataFrame(compat_rows)
compat.to_csv(COUT / "03_GENESET_COMPATIBILITY.csv", index=False)
pd.DataFrame(shared_member_rows).to_csv(WORK / "shared_program_memberships.tsv", sep="\t", index=False)

preflight = {
    "author_cells": int(len(am)),
    "frozen_t_cells": int(len(F)),
    "shared_frozen_t_cells": int(len(shared)),
    "author_subtypes": subtypes,
    "eligible_subtypes": elig.loc[elig["formal_pseudobulk_eligible"].eq("YES"), "author_subtype"].tolist(),
    "supported_matches": matches.loc[matches["matching_classification"].eq("SUPPORTED_EXPLORATORY_MATCH"), ["dog_author_subtype", "human_author_subtype"]].to_dict("records"),
    "matching_decision_sha256": freeze_record["sha256"],
    "source_checksums": {
        "author_metadata_sha256": sha256(META),
        "frozen_h5ad_sha256": sha256(H5AD),
        "ortholog_master_sha256": sha256(ORTHO),
        "author_script_sha256": sha256(AUTHOR_SCRIPT),
    },
}
(WORK / "phase6_preflight.json").write_text(json.dumps(preflight, indent=2), encoding="utf-8")
print(json.dumps({"status": "PHASE6_PREFLIGHT_COMPLETE", **{k: preflight[k] for k in ["author_cells", "frozen_t_cells", "shared_frozen_t_cells", "eligible_subtypes", "supported_matches"]}}, indent=2))
