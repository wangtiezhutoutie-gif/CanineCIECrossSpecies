#!/usr/bin/env python
"""Freeze canine Phase4 gene sets, feature-ID mapping, and audit tables."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import platform
import re
from collections import defaultdict
from datetime import date
from pathlib import Path

import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
GSEA = ROOT / "06_pseudobulk" / "gsea"
RAW = GSEA / "gene_sets" / "raw"
GENESETS = GSEA / "gene_sets"
MAPPING = GSEA / "mapping"
META = ROOT / "00_metadata"
LOG = ROOT / "logs" / "phase4_prepare_genesets.log"
LINEAGES = ("Epithelial", "T", "Myeloid")
RETRIEVAL_DATE = "2026-08-24"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def log(message: str) -> None:
    print(message, flush=True)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(message + "\n")


def require(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)


def parse_gene_info(path: Path):
    gene_info = pd.read_csv(path, sep="\t", dtype=str, compression="gzip")
    gene_info.columns = [c.lstrip("#") for c in gene_info.columns]
    gene_to_symbol = dict(zip(gene_info["GeneID"], gene_info["Symbol"]))
    ens_to_gene: dict[str, set[str]] = defaultdict(set)
    symbol_to_gene: dict[str, set[str]] = defaultdict(set)
    for row in gene_info.itertuples(index=False):
        gene_id = str(row.GeneID)
        symbol = str(row.Symbol)
        if symbol and symbol != "-":
            symbol_to_gene[symbol.upper()].add(gene_id)
        synonyms = getattr(row, "Synonyms", "-")
        if synonyms and synonyms != "-":
            for syn in str(synonyms).split("|"):
                if syn and syn != "-":
                    symbol_to_gene[syn.upper()].add(gene_id)
        for xref in str(row.dbXrefs).split("|"):
            if xref.startswith("Ensembl:"):
                ens_to_gene[xref.split(":", 1)[1].split(".", 1)[0]].add(gene_id)
    return gene_info, gene_to_symbol, ens_to_gene, symbol_to_gene


def parse_obo(path: Path):
    terms: dict[str, dict] = {}
    current: dict | None = None
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "[Term]":
                if current and current.get("id"):
                    terms[current["id"]] = current
                current = {"parents": set(), "obsolete": False}
            elif line.startswith("["):
                if current and current.get("id"):
                    terms[current["id"]] = current
                current = None
            elif current is not None:
                if line.startswith("id: GO:"):
                    current["id"] = line[4:]
                elif line.startswith("name: "):
                    current["name"] = line[6:]
                elif line.startswith("namespace: "):
                    current["namespace"] = line[11:]
                elif line.startswith("is_a: GO:"):
                    current["parents"].add(line.split()[1])
                elif line.startswith("relationship: part_of GO:"):
                    current["parents"].add(line.split()[2])
                elif line == "is_obsolete: true":
                    current["obsolete"] = True
    if current and current.get("id"):
        terms[current["id"]] = current
    bp = {
        go_id: term
        for go_id, term in terms.items()
        if term.get("namespace") == "biological_process" and not term.get("obsolete")
    }
    return bp


def ancestors_for(go_id: str, terms: dict[str, dict], memo: dict[str, set[str]], active=None):
    if go_id in memo:
        return memo[go_id]
    if active is None:
        active = set()
    if go_id in active:
        return {go_id}
    active = active | {go_id}
    result = {go_id}
    for parent in terms.get(go_id, {}).get("parents", set()):
        if parent in terms:
            result.update(ancestors_for(parent, terms, memo, active))
    memo[go_id] = result
    return result


def parse_gpi_mapping(path: Path, ens_to_gene, symbol_to_gene):
    uniprot_to_gene: dict[str, set[str]] = defaultdict(set)
    method: dict[tuple[str, str], str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            uniprot = fields[0].split(":", 1)[-1]
            symbol = fields[1].upper()
            ensembl_ids = {
                token.split(":", 1)[1].split(".", 1)[0]
                for token in fields[9].split("|")
                if token.startswith("Ensembl:")
            }
            ids = set().union(*(ens_to_gene.get(x, set()) for x in ensembl_ids)) if ensembl_ids else set()
            mapping_method = "GPI_ENSEMBL_TO_NCBI_GENE"
            if not ids and len(symbol_to_gene.get(symbol, set())) == 1:
                ids = set(symbol_to_gene[symbol])
                mapping_method = "GPI_UNIQUE_SYMBOL_TO_NCBI_GENE"
            for gene_id in ids:
                uniprot_to_gene[uniprot].add(gene_id)
                method[(uniprot, gene_id)] = mapping_method
    return uniprot_to_gene, method


def build_go_membership(gaf_path, gpi_path, obo_path, ens_to_gene, symbol_to_gene, gene_to_symbol):
    terms = parse_obo(obo_path)
    memo: dict[str, set[str]] = {}
    uniprot_to_gene, mapping_methods = parse_gpi_mapping(gpi_path, ens_to_gene, symbol_to_gene)
    direct: dict[str, set[str]] = defaultdict(set)
    direct_method: dict[tuple[str, str], set[str]] = defaultdict(set)
    annotated_rows = mapped_rows = 0
    with gzip.open(gaf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 15 or fields[8] != "P" or "NOT" in fields[3].split("|"):
                continue
            go_id = fields[4]
            if go_id not in terms:
                continue
            annotated_rows += 1
            uniprot = fields[1]
            genes = uniprot_to_gene.get(uniprot, set())
            if genes:
                mapped_rows += 1
            for gene_id in genes:
                direct[go_id].add(gene_id)
                direct_method[(go_id, gene_id)].add(mapping_methods.get((uniprot, gene_id), "GPI_MAPPING"))

    propagated: dict[str, set[str]] = defaultdict(set)
    for go_id, genes in direct.items():
        for ancestor in ancestors_for(go_id, terms, memo):
            propagated[ancestor].update(genes)

    rows = []
    for go_id, genes in propagated.items():
        name = terms[go_id].get("name", go_id)
        for gene_id in sorted(genes, key=lambda x: (len(x), x)):
            rows.append(
                {
                    "collection": "GO_BP",
                    "pathway_id": go_id,
                    "pathway_name": name,
                    "target_geneset_id": gene_id,
                    "target_symbol": gene_to_symbol.get(gene_id, ""),
                    "mapping_type": "CANINE_NATIVE_GO_UNIPROT_WITH_GO_ANCESTOR_PROPAGATION",
                    "source_version": "GO release 2026-08-05; CANLF GAF generated 2026-07-28",
                    "source_species": "Canis lupus familiaris",
                    "target_species": "Canis lupus familiaris",
                }
            )
    return pd.DataFrame(rows), {
        "n_bp_terms_ontology": len(terms),
        "n_direct_terms_mapped": len(direct),
        "n_propagated_terms": len(propagated),
        "n_gaf_bp_rows": annotated_rows,
        "n_gaf_bp_rows_mapped": mapped_rows,
        "n_uniprot_mapped": len(uniprot_to_gene),
    }


def build_feature_mapping(ens_to_gene, gene_to_symbol):
    master = pd.read_csv(META / "dog_gene_mapping_master.csv", dtype=str).fillna("")
    ranks = {lineage: pd.read_csv(ROOT / "06_pseudobulk" / "ranked_lists" / f"{lineage}_ranked.csv") for lineage in LINEAGES}
    used_features = set().union(*(set(df["feature_id"]) for df in ranks.values()))
    master = master[master["feature_id"].isin(used_features)].copy()
    rows = []
    for row in master.itertuples(index=False):
        targets = sorted(ens_to_gene.get(row.feature_id, set()), key=lambda x: (len(x), x))
        if len(targets) == 1:
            status = "MAPPED_EXACT_ENSEMBL_TO_NCBI_GENE"
            target = targets[0]
            note = "Phase3 feature_id matched an Ensembl dbXref in frozen canine NCBI gene_info."
        elif len(targets) > 1:
            status = "AMBIGUOUS_ENSEMBL_TO_MULTIPLE_NCBI_GENE"
            target = ""
            note = "Ambiguous current annotation; excluded rather than selecting by symbol."
        else:
            status = "UNMAPPED_CURRENT_NCBI_GENE_INFO"
            target = ""
            note = "No symbol fallback was used for Phase3 feature identity."
        rows.append(
            {
                "feature_id": row.feature_id,
                "original_gene_symbol": row.original_gene_symbol,
                "unique_var_name": row.unique_var_name,
                "source_ensembl_id": row.feature_id,
                "target_geneset_id": target,
                "target_symbol": gene_to_symbol.get(target, ""),
                "mapping_status": status,
                "mapping_source": "NCBI Canis familiaris gene_info dated 2026-08-19",
                "mapping_type": "DOG_ENSEMBL_TO_NCBI_GENE",
                "collision_group": "",
                "used_for_gsea": "TRUE_AS_MAPPING_CANDIDATE" if target else "FALSE",
                "notes": note,
            }
        )
    mapping = pd.DataFrame(rows)
    counts = mapping.loc[mapping.target_geneset_id.ne(""), "target_geneset_id"].value_counts()
    duplicated_targets = set(counts[counts > 1].index)
    mapping.loc[mapping.target_geneset_id.isin(duplicated_targets), "collision_group"] = (
        "NCBI_GENE_" + mapping.loc[mapping.target_geneset_id.isin(duplicated_targets), "target_geneset_id"]
    )
    mapping.to_csv(MAPPING / "phase4_feature_to_geneset_id.csv", index=False)

    collision_rows = []
    collapsed_by_lineage = {}
    for lineage, rank in ranks.items():
        merged = rank.merge(
            mapping[["feature_id", "target_geneset_id", "target_symbol", "mapping_status"]],
            on="feature_id",
            how="left",
            validate="one_to_one",
        )
        mapped = merged[merged["target_geneset_id"].fillna("").ne("")].copy()
        mapped["abs_rank_metric"] = mapped["rank_metric"].abs()
        mapped = mapped.sort_values(
            ["target_geneset_id", "abs_rank_metric", "PValue", "feature_id"],
            ascending=[True, False, True, True],
            kind="mergesort",
        )
        mapped["mapping_rank"] = mapped.groupby("target_geneset_id").cumcount() + 1
        mapped["kept_after_collision"] = mapped.mapping_rank.eq(1)
        for target, group in mapped.groupby("target_geneset_id"):
            if len(group) <= 1:
                continue
            kept = group.iloc[0]["feature_id"]
            for x in group.itertuples(index=False):
                collision_rows.append(
                    {
                        "lineage": lineage,
                        "target_geneset_id": target,
                        "target_symbol": x.target_symbol,
                        "feature_id": x.feature_id,
                        "rank_metric": x.rank_metric,
                        "abs_rank_metric": x.abs_rank_metric,
                        "PValue": x.PValue,
                        "kept_feature_id": kept,
                        "is_kept": bool(x.kept_after_collision),
                        "decision_rule": "max(abs(rank_metric)); then smaller PValue; then feature_id lexical",
                    }
                )
        collapsed = mapped[mapped.kept_after_collision].copy()
        collapsed = collapsed.sort_values(["rank_metric", "feature_id"], ascending=[False, True])
        collapsed.to_csv(MAPPING / f"{lineage}_Full7_collapsed_rank.csv", index=False)
        collapsed_by_lineage[lineage] = collapsed
    pd.DataFrame(collision_rows).to_csv(MAPPING / "phase4_mapping_collisions.csv", index=False)
    return mapping, ranks, collapsed_by_lineage


def main():
    for directory in [GENESETS, MAPPING, META, LOG.parent]:
        directory.mkdir(parents=True, exist_ok=True)
    LOG.write_text("", encoding="utf-8")

    required = [
        RAW / "h.all.v2026.1.Hs.symbols.gmt",
        RAW / "NCBI2Reactome_All_Levels_v97.txt",
        RAW / "go-basic_2026-08-05.obo",
        RAW / "CANLF-uniprot_2026-08-05.gaf.gz",
        RAW / "CANLF-uniprot_2026-08-05.gpi.gz",
        RAW / "Canis_familiaris.gene_info_2026-08-19.gz",
        GENESETS / "hallmark_dog_ortholog_projected.tsv.gz",
    ]
    for path in required:
        require(path)

    manifest = pd.read_csv(META / "PHASE3_MANIFEST.csv", dtype=str)
    manifest_audit = []
    for row in manifest.itertuples(index=False):
        path = ROOT / row.path_relative_to_project.replace("/", os.sep)
        actual_hash = sha256(path) if path.is_file() else "MISSING"
        manifest_audit.append(
            {
                "path_relative_to_project": row.path_relative_to_project,
                "role": row.role,
                "expected_sha256": row.sha256,
                "actual_sha256": actual_hash,
                "hash_match": actual_hash == row.sha256,
            }
        )
    manifest_audit = pd.DataFrame(manifest_audit)
    manifest_audit.to_csv(META / "PHASE4_PHASE3_INPUT_HASH_AUDIT.csv", index=False)
    if not manifest_audit.hash_match.all():
        bad = manifest_audit.loc[~manifest_audit.hash_match, "path_relative_to_project"].tolist()
        raise RuntimeError(f"Phase3 manifest validation failed: {bad}")
    log(f"PHASE3_MANIFEST verified: {len(manifest_audit)}/{len(manifest_audit)} files.")

    gene_info, gene_to_symbol, ens_to_gene, symbol_to_gene = parse_gene_info(required[5])
    log(f"NCBI gene_info: {len(gene_info)} canine genes; {len(ens_to_gene)} Ensembl IDs.")

    hallmark = pd.read_csv(GENESETS / "hallmark_dog_ortholog_projected.tsv.gz", sep="\t", dtype=str)
    hallmark_core = hallmark[[
        "collection", "pathway_id", "pathway_name", "target_geneset_id", "target_symbol",
        "mapping_type", "source_version", "source_species", "target_species"
    ]].drop_duplicates()
    hallmark_core.to_csv(GENESETS / "hallmark_dog_ortholog_projected_membership.tsv.gz", sep="\t", index=False)

    reactome = pd.read_csv(
        RAW / "NCBI2Reactome_All_Levels_v97.txt",
        sep="\t",
        header=None,
        names=["target_geneset_id", "pathway_id", "url", "pathway_name", "evidence", "species"],
        dtype=str,
    )
    reactome = reactome[reactome.species.eq("Canis familiaris")].copy()
    reactome["target_symbol"] = reactome.target_geneset_id.map(gene_to_symbol).fillna("")
    reactome["collection"] = "Reactome"
    reactome["mapping_type"] = "CANINE_NATIVE_REACTOME_ANNOTATED_OR_INFERRED"
    reactome["source_version"] = "Reactome v97"
    reactome["source_species"] = "Canis familiaris"
    reactome["target_species"] = "Canis lupus familiaris"
    reactome_core = reactome[[
        "collection", "pathway_id", "pathway_name", "target_geneset_id", "target_symbol",
        "mapping_type", "source_version", "source_species", "target_species"
    ]].drop_duplicates()
    reactome_core.to_csv(GENESETS / "reactome_canis_familiaris_gene_sets.tsv.gz", sep="\t", index=False)
    reactome.to_csv(META / "reactome_canis_familiaris_gene_sets.tsv.gz", sep="\t", index=False)
    log(f"Reactome v97 dog subset: {len(reactome_core)} memberships; {reactome_core.pathway_id.nunique()} pathways.")

    go_bp, go_stats = build_go_membership(
        RAW / "CANLF-uniprot_2026-08-05.gaf.gz",
        RAW / "CANLF-uniprot_2026-08-05.gpi.gz",
        RAW / "go-basic_2026-08-05.obo",
        ens_to_gene,
        symbol_to_gene,
        gene_to_symbol,
    )
    go_bp = go_bp.drop_duplicates()
    go_bp.to_csv(GENESETS / "go_bp_canis_familiaris_gene_sets.tsv.gz", sep="\t", index=False)
    log(f"GO BP: {len(go_bp)} propagated memberships; {go_bp.pathway_id.nunique()} pathways; {json.dumps(go_stats)}")

    mapping, ranks, collapsed_by_lineage = build_feature_mapping(ens_to_gene, gene_to_symbol)
    log(f"Feature mapping: {(mapping.target_geneset_id != '').sum()}/{len(mapping)} exact Ensembl-to-NCBI mappings.")

    collections = {"Hallmark": hallmark_core, "Reactome": reactome_core, "GO_BP": go_bp}
    coverage_rows = []
    audit_rows = []
    source_meta = {
        "Hallmark": ("MSigDB Hallmark official GMT", "MSigDB 2026.1.Hs; babelgene 22.9", "Homo sapiens", "ORTHOLOG_PROJECTED_HALLMARK_BABELGENE", "NCBI Gene ID"),
        "Reactome": ("Reactome NCBI2Reactome_All_Levels", "Reactome v97", "Canis familiaris", "CANINE_NATIVE_REACTOME_ANNOTATED_OR_INFERRED", "NCBI Gene ID"),
        "GO_BP": ("GO Consortium CANLF UniProt GAF + go-basic", "GO 2026-08-05; CANLF GAF 2026-07-28", "Canis lupus familiaris", "CANINE_NATIVE_GO_UNIPROT_WITH_GO_ANCESTOR_PROPAGATION", "NCBI Gene ID"),
    }
    for collection, members in collections.items():
        path_sizes = members.groupby("pathway_id").target_geneset_id.nunique()
        audit_rows.append(
            {
                "collection": collection,
                "source": source_meta[collection][0],
                "source_version": source_meta[collection][1],
                "source_species": source_meta[collection][2],
                "target_species": "Canis lupus familiaris",
                "mapping_type": source_meta[collection][3],
                "retrieval_date": RETRIEVAL_DATE,
                "n_pathways_raw": int(members.pathway_id.nunique()),
                "n_pathways_after_filter": int(((path_sizes >= 15) & (path_sizes <= 500)).sum()),
                "gene_identifier": source_meta[collection][4],
                "notes": "After-filter count uses source membership size 15-500; fgsea applies the same bounds after intersecting each ranked universe.",
            }
        )
        member_genes = set(members.target_geneset_id.astype(str))
        for lineage, rank in ranks.items():
            line_map = rank[["feature_id"]].merge(
                mapping[["feature_id", "target_geneset_id"]], on="feature_id", how="left"
            )
            in_collection = line_map.target_geneset_id.fillna("").isin(member_genes)
            mapped_features = int(in_collection.sum())
            unique_targets = int(line_map.loc[in_collection, "target_geneset_id"].nunique())
            coverage_rows.append(
                {
                    "lineage": lineage,
                    "collection": collection,
                    "tested_features": len(rank),
                    "mapped_features": mapped_features,
                    "unique_target_genes": unique_targets,
                    "mapping_rate": mapped_features / len(rank),
                    "unmapped_features": len(rank) - mapped_features,
                    "collapsed_features": mapped_features - unique_targets,
                    "feature_to_ncbi_mapped": int(line_map.target_geneset_id.fillna("").ne("").sum()),
                    "notes": "Collection-specific rate counts features mapped to an NCBI Gene ID present in that collection.",
                }
            )
    pd.DataFrame(coverage_rows).to_csv(MAPPING / "phase4_mapping_coverage.csv", index=False)
    pd.DataFrame(audit_rows).to_csv(META / "PHASE4_GENESET_AUDIT.csv", index=False)

    source_files = [
        (RAW / "h.all.v2026.1.Hs.symbols.gmt", "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt", "MSigDB 2026.1.Hs"),
        (RAW / "h.all.v2026.1.Hs.symbols.official.gmt", "Official-copy byte identity check", "MSigDB 2026.1.Hs"),
        (RAW / "NCBI2Reactome_All_Levels_v97.txt", "https://download.reactome.org/97/NCBI2Reactome_All_Levels.txt", "Reactome v97"),
        (RAW / "go-basic_2026-08-05.obo", "https://release.geneontology.org/2026-08-05/ontology/go-basic.obo", "GO 2026-08-05"),
        (RAW / "CANLF-uniprot_2026-08-05.gaf.gz", "https://release.geneontology.org/2026-08-05/annotations/gaf/CANLF-uniprot.gaf.gz", "GO 2026-08-05"),
        (RAW / "CANLF-uniprot_2026-08-05.gpi.gz", "https://release.geneontology.org/2026-08-05/annotations/gpi/CANLF-uniprot.gpi.gz", "GO 2026-08-05"),
        (RAW / "Canis_familiaris.gene_info_2026-08-19.gz", "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Canis_familiaris.gene_info.gz", "NCBI file dated 2026-08-19"),
    ]
    source_rows = []
    for path, url, version in source_files:
        require(path)
        source_rows.append({
            "file_relative_to_project": path.relative_to(ROOT).as_posix(),
            "source_url": url,
            "source_version": version,
            "retrieval_date": RETRIEVAL_DATE,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    pd.DataFrame(source_rows).to_csv(META / "PHASE4_GENESET_SOURCE_FILES.csv", index=False)

    parameters = [
        ("algorithm", "fgseaMultilevel"), ("rank_metric", "sign(logFC) * sqrt(F)"),
        ("positive_direction", "CIE higher"), ("negative_direction", "CIE lower"),
        ("minSize", "15"), ("maxSize", "500"), ("eps", "1e-10"),
        ("scoreType", "std"), ("nproc", "1"), ("seed", "20260824"),
        ("multiple_testing", "Benjamini-Hochberg separately within lineage and collection"),
        ("primary_universe", "Full7 complete Phase3 tested rank after deterministic ID collapse"),
        ("sensitivity_universe", "Full7/No-cie2 common feature universe; same target IDs and gene sets"),
        ("collision_rule", "max(abs(rank_metric)); smaller PValue; feature_id lexical"),
        ("sensitivity_flag_direction_flip", "NES signs differ"),
        ("sensitivity_flag_sensitive", "same direction and (abs(delta_NES)>=0.5 or leading-edge Jaccard<0.25)"),
        ("sensitivity_flag_stable", "same direction and abs(delta_NES)<0.5 and leading-edge Jaccard>=0.25"),
    ]
    pd.DataFrame(parameters, columns=["parameter", "value"]).to_csv(META / "PHASE4_GSEA_PARAMETERS.csv", index=False)

    log("Phase4 gene-set preparation and mapping audit completed.")


if __name__ == "__main__":
    main()
