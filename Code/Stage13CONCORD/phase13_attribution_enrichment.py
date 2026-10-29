#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P13 = ROOT / "13_concord"
GENESET_DIR = ROOT / "07_human_crohn" / "GSE266616" / "gsea" / "gene_sets"
FEATURES = P13 / "attribution" / "CONCORD_DISEASE_ALIGNED_FEATURES.csv"
MAPPING = P13 / "input" / "CONCORD_SHARED_FEATURE_MAPPING.csv"
OUT = P13 / "attribution" / "CONCORD_ATTRIBUTION_ENRICHMENT.csv"
THEME_OUT = P13 / "attribution" / "CONCORD_ATTRIBUTION_THEME_SUMMARY.csv"
SHARED_OUT = P13 / "attribution" / "CONCORD_ATTRIBUTION_SHARED_PATHWAYS.csv"
SETS = {
    "Hallmark": GENESET_DIR / "HUMAN_HALLMARK_GENE_SETS.tsv.gz",
    "Reactome": GENESET_DIR / "HUMAN_REACTOME_GENE_SETS.tsv.gz",
    "GO_BP": GENESET_DIR / "HUMAN_GO_BP_GENE_SETS.tsv.gz",
}


def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty(len(p), dtype=float)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def theme(name: str) -> str:
    x = name.upper()
    rules = [
        ("IMMUNE_TNF_IFN_JAK", ["TNF", "NFKB", "NF-K", "INTERFERON", "JAK", "STAT", "CYTOKINE", "INTERLEUKIN", "INFLAM", "IMMUNE", "T CELL", "T-CELL"]),
        ("ANTIGEN_COMPLEMENT", ["ANTIGEN", "COMPLEMENT", "C1Q", "MHC CLASS", "HLA-"]),
        ("OXPHOS_METABOLISM", ["OXIDATIVE PHOSPHORYLATION", "RESPIRATORY", "MITOCHON", "METABOL", "FATTY ACID", "GLYCOL", "TRICARBOXYLIC", "TCA CYCLE", "LIPID", "XENOBIOTIC"]),
        ("TRANSLATION", ["RIBOSOM", "TRANSLATION", "PROTEIN SYNTHESIS", "PEPTIDE BIOSYNTH"]),
        ("PROLIFERATION", ["CELL CYCLE", "MITOTIC", "MITOSIS", "DNA REPLICATION", "E2F", "G2M", "PROLIFER"]),
        ("EPITHELIAL_BARRIER_STRUCTURE", ["EPITHEL", "TIGHT JUNCTION", "CELL JUNCTION", "KERATIN", "INTERMEDIATE FILAMENT", "EXTRACELLULAR MATRIX", "ECM", "DIGEST", "ABSORPTION"]),
    ]
    for label, terms in rules:
        if any(t in x for t in terms):
            return label
    return "OTHER"


def load_sets(universe: set[str]):
    result = {}
    for collection, path in SETS.items():
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        grouped = []
        for keys, block in frame.groupby(["pathway_id", "pathway_name", "source_version"], sort=False):
            genes = set(block["gene_symbol"].dropna().astype(str)) & universe
            if 15 <= len(genes) <= 500:
                grouped.append((*keys, genes))
        result[collection] = grouped
        if not grouped:
            raise RuntimeError(f"No eligible {collection} sets")
    return result


def main():
    features = pd.read_csv(FEATURES, low_memory=False)
    mapping = pd.read_csv(MAPPING, dtype=str)
    selected_mapping = mapping.loc[mapping["selected_hvg"].eq("True")]
    universe = set(selected_mapping["human_symbol"].dropna().astype(str))
    universe.discard("")
    if len(selected_mapping) != 10000 or len(universe) != 9999:
        raise RuntimeError(f"Unexpected symbol universe: {len(universe)}")
    sets = load_sets(universe)
    n_universe = len(universe)
    rows = []
    group_cols = ["lineage", "species", "latent_dimension", "latent_index", "attribution_tail"]
    for keys, block in features.groupby(group_cols, sort=False):
        query = set(block["human_symbol"].dropna().astype(str)) & universe
        if len(query) < 40:
            raise RuntimeError(f"Small query {keys}: {len(query)}")
        loading = float(block["dog_disease_loading"].iloc[0] if keys[1] == "dog" else block["human_disease_loading"].iloc[0])
        aligned_tail = "TOP_POSITIVE" if loading >= 0 else "TOP_NEGATIVE"
        tail_relation = "DISEASE_ALIGNED" if keys[4] == aligned_tail else "DISEASE_OPPOSED"
        for collection, pathways in sets.items():
            tested = []
            p_values = []
            for pathway_id, pathway_name, source_version, genes in pathways:
                overlap = query & genes
                pvalue = float(hypergeom.sf(len(overlap) - 1, n_universe, len(genes), len(query)))
                tested.append((pathway_id, pathway_name, source_version, genes, overlap, pvalue))
                p_values.append(pvalue)
            fdr = bh(p_values)
            for item, qvalue in zip(tested, fdr):
                pathway_id, pathway_name, source_version, genes, overlap, pvalue = item
                if len(overlap) < 2:
                    continue
                a = len(overlap)
                b = len(query) - a
                c = len(genes) - a
                d = n_universe - a - b - c
                odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
                rows.append({
                    "lineage": keys[0], "species": keys[1],
                    "latent_dimension": keys[2], "latent_index": keys[3],
                    "attribution_tail": keys[4], "disease_loading": loading,
                    "tail_relation_to_disease": tail_relation,
                    "collection": collection, "pathway_id": pathway_id,
                    "pathway_name": pathway_name, "theme": theme(pathway_name),
                    "source_version": source_version, "universe_size": n_universe,
                    "pathway_size_in_universe": len(genes), "query_size": len(query),
                    "overlap_n": len(overlap), "odds_ratio_haldane": odds,
                    "p_value": pvalue, "fdr_bh_within_context_collection": qvalue,
                    "overlap_genes": "|".join(sorted(overlap)),
                    "posthoc_only": "YES", "feedback_to_training": "NO",
                })
    enrich = pd.DataFrame(rows).sort_values(
        ["lineage", "species", "latent_index", "attribution_tail", "collection", "fdr_bh_within_context_collection", "p_value"]
    )
    enrich.to_csv(OUT, index=False)

    sig = enrich.loc[enrich["fdr_bh_within_context_collection"] <= 0.10].copy()
    theme_rows = []
    for keys, block in sig.groupby(["lineage", "species", "tail_relation_to_disease", "theme"], sort=False):
        theme_rows.append({
            "lineage": keys[0], "species": keys[1],
            "tail_relation_to_disease": keys[2], "theme": keys[3],
            "significant_pathway_dimension_hits": len(block),
            "n_unique_dimensions": block["latent_dimension"].nunique(),
            "n_unique_pathways": block["pathway_id"].nunique(),
            "minimum_fdr": block["fdr_bh_within_context_collection"].min(),
            "top_pathways": "|".join(block.sort_values("fdr_bh_within_context_collection")["pathway_name"].drop_duplicates().head(10)),
        })
    pd.DataFrame(theme_rows).to_csv(THEME_OUT, index=False)

    aligned = sig.loc[sig["tail_relation_to_disease"] == "DISEASE_ALIGNED"].copy()
    dog = aligned.loc[aligned["species"] == "dog"]
    human = aligned.loc[aligned["species"] == "human"]
    shared = dog.merge(
        human,
        on=["lineage", "latent_dimension", "pathway_id"],
        suffixes=("_dog", "_human"),
    )
    keep = [
        "lineage", "latent_dimension", "pathway_id", "pathway_name_dog", "theme_dog",
        "collection_dog", "disease_loading_dog", "disease_loading_human",
        "attribution_tail_dog", "attribution_tail_human",
        "fdr_bh_within_context_collection_dog", "fdr_bh_within_context_collection_human",
        "overlap_genes_dog", "overlap_genes_human",
    ]
    shared[keep].to_csv(SHARED_OUT, index=False)

    annotations = {}
    for row in sig.sort_values("fdr_bh_within_context_collection").itertuples(index=False):
        for gene in str(row.overlap_genes).split("|"):
            key = (row.lineage, row.species, row.latent_dimension, row.attribution_tail, gene)
            annotations.setdefault(key, [])
            label = f"{row.collection}:{row.pathway_name}"
            if label not in annotations[key] and len(annotations[key]) < 5:
                annotations[key].append(label)
    features["pathway_annotation"] = [
        "|".join(annotations.get((r.lineage, r.species, r.latent_dimension, r.attribution_tail, str(r.human_symbol)), []))
        or "NO_FDR_0.10_ENRICHED_PATHWAY_MEMBERSHIP"
        for r in features.itertuples(index=False)
    ]
    features.to_csv(FEATURES, index=False)
    print(f"enrichment_rows={len(enrich)} significant_fdr10={len(sig)} shared_aligned={len(shared)}")


if __name__ == "__main__":
    main()
