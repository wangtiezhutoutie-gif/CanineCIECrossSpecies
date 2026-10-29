from __future__ import annotations

import csv
import gzip
import hashlib
import os
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
RAW = ROOT / "06_pseudobulk" / "gsea" / "gene_sets" / "raw"
BASE = ROOT / "07_human_crohn" / "GSE266616" / "gsea"
OUT = BASE / "gene_sets"
META = ROOT / "00_metadata"

HALLMARK = RAW / "h.all.v2026.1.Hs.symbols.gmt"
REACTOME = RAW / "NCBI2Reactome_All_Levels_v97.txt"
GENE_INFO = RAW / "Homo_sapiens.gene_info.gz"
GO_GAF = RAW / "goa_human.gaf.gz"
GO_OBO = RAW / "go-basic_2026-08-05.obo"

OUTPUTS = {
    "Hallmark": OUT / "HUMAN_HALLMARK_GENE_SETS.tsv.gz",
    "Reactome": OUT / "HUMAN_REACTOME_GENE_SETS.tsv.gz",
    "GO_BP": OUT / "HUMAN_GO_BP_GENE_SETS.tsv.gz",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_tsv_gz(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        frame.to_csv(tmp, sep="\t", index=False, compression="gzip", quoting=csv.QUOTE_MINIMAL)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        frame.to_csv(tmp, index=False)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def require_inputs() -> None:
    missing = [str(path) for path in (HALLMARK, REACTOME, GENE_INFO, GO_GAF, GO_OBO) if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing frozen Human gene-set inputs: " + "; ".join(missing))


def load_human_symbols() -> tuple[dict[str, str], set[str]]:
    gene_info = pd.read_csv(GENE_INFO, sep="\t", compression="gzip", dtype=str, low_memory=False)
    gene_info = gene_info.loc[gene_info["#tax_id"].eq("9606")].copy()
    entrez_to_symbol = dict(zip(gene_info["GeneID"], gene_info["Symbol"], strict=False))
    official = set(gene_info["Symbol"].dropna().astype(str))
    if len(entrez_to_symbol) < 60_000:
        raise RuntimeError(f"Unexpectedly small Human gene_info mapping: {len(entrez_to_symbol)}")
    return entrez_to_symbol, official


def hallmark_frame() -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    with HALLMARK.open("rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            pathway_id, description, *symbols = fields
            name = pathway_id.removeprefix("HALLMARK_").replace("_", " ").title()
            for symbol in symbols:
                if symbol:
                    rows.append({
                        "collection": "Hallmark", "pathway_id": pathway_id,
                        "pathway_name": name, "gene_symbol": symbol,
                        "source_version": "MSigDB 2026.1.Hs", "source_species": "Homo sapiens",
                        "source_reference": description,
                    })
    frame = pd.DataFrame(rows).drop_duplicates(["pathway_id", "gene_symbol"])
    if frame["pathway_id"].nunique() != 50:
        raise RuntimeError(f"Expected 50 Hallmark pathways, got {frame['pathway_id'].nunique()}")
    return frame


def reactome_frame(entrez_to_symbol: dict[str, str]) -> pd.DataFrame:
    raw = pd.read_csv(
        REACTOME, sep="\t", header=None,
        names=["entrez", "pathway_id", "url", "pathway_name", "evidence", "species"],
        dtype=str, low_memory=False,
    )
    human = raw.loc[raw["species"].eq("Homo sapiens")].copy()
    human["gene_symbol"] = human["entrez"].map(entrez_to_symbol)
    human = human.loc[human["gene_symbol"].notna() & human["pathway_id"].str.startswith("R-HSA-", na=False)]
    frame = human.assign(
        collection="Reactome", source_version="Reactome v97", source_species="Homo sapiens",
        source_reference=human["url"],
    )[["collection", "pathway_id", "pathway_name", "gene_symbol", "source_version", "source_species", "source_reference"]]
    frame = frame.drop_duplicates(["pathway_id", "gene_symbol"])
    if frame["pathway_id"].nunique() < 2_000:
        raise RuntimeError(f"Unexpectedly few Human Reactome pathways: {frame['pathway_id'].nunique()}")
    return frame


def parse_obo() -> tuple[dict[str, str], dict[str, str]]:
    names: dict[str, str] = {}
    aliases: dict[str, str] = {}
    block: dict[str, list[str]] = defaultdict(list)

    def commit() -> None:
        term_id = block.get("id", [None])[0]
        if not term_id or block.get("namespace", [None])[0] != "biological_process":
            return
        if block.get("is_obsolete", ["false"])[0] == "true":
            return
        names[term_id] = block.get("name", [term_id])[0]
        for alt in block.get("alt_id", []):
            aliases[alt] = term_id

    with GO_OBO.open("rt", encoding="utf-8") as handle:
        in_term = False
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "[Term]":
                if in_term:
                    commit()
                block = defaultdict(list)
                in_term = True
                continue
            if line.startswith("["):
                if in_term:
                    commit()
                block = defaultdict(list)
                in_term = False
                continue
            if in_term and ": " in line:
                key, value = line.split(": ", 1)
                block[key].append(value)
        if in_term:
            commit()
    if len(names) < 20_000:
        raise RuntimeError(f"Unexpectedly few GO BP terms: {len(names)}")
    return names, aliases


def go_bp_frame(official_symbols: set[str]) -> pd.DataFrame:
    names, aliases = parse_obo()
    rows: set[tuple[str, str]] = set()
    with gzip.open(GO_GAF, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 15 or fields[8] != "P" or "NOT" in fields[3].split("|"):
                continue
            symbol = fields[2]
            go_id = aliases.get(fields[4], fields[4])
            if go_id in names and symbol in official_symbols:
                rows.add((go_id, symbol))
    frame = pd.DataFrame(rows, columns=["pathway_id", "gene_symbol"])
    frame["collection"] = "GO_BP"
    frame["pathway_name"] = frame["pathway_id"].map(names)
    frame["source_version"] = "GO basic 2026-08-05; GOA Human current acquisition 2026-08-25"
    frame["source_species"] = "Homo sapiens"
    frame["source_reference"] = "https://current.geneontology.org/annotations/goa_human.gaf.gz"
    frame = frame[["collection", "pathway_id", "pathway_name", "gene_symbol", "source_version", "source_species", "source_reference"]]
    if frame["pathway_id"].nunique() < 10_000:
        raise RuntimeError(f"Unexpectedly few annotated Human GO BP pathways: {frame['pathway_id'].nunique()}")
    return frame


def main() -> None:
    require_inputs()
    entrez_to_symbol, official_symbols = load_human_symbols()
    frames = {
        "Hallmark": hallmark_frame(),
        "Reactome": reactome_frame(entrez_to_symbol),
        "GO_BP": go_bp_frame(official_symbols),
    }
    audit_rows: list[dict[str, object]] = []
    for collection, frame in frames.items():
        frame = frame.sort_values(["pathway_id", "gene_symbol"], kind="stable").reset_index(drop=True)
        atomic_tsv_gz(frame, OUTPUTS[collection])
        sizes = frame.groupby("pathway_id", sort=False).size()
        audit_rows.append({
            "collection": collection,
            "membership_file": str(OUTPUTS[collection]),
            "pathways": int(frame["pathway_id"].nunique()),
            "memberships": int(len(frame)),
            "unique_gene_symbols": int(frame["gene_symbol"].nunique()),
            "pathways_size_15_500": int(sizes.between(15, 500).sum()),
            "sha256": sha256(OUTPUTS[collection]),
            "status": "PASS",
        })
    inputs = pd.DataFrame([
        {"asset": p.name, "path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p), "acquired_or_frozen": "2026-08-25" if p in {GENE_INFO, GO_GAF} else "PREEXISTING_FROZEN"}
        for p in (HALLMARK, REACTOME, GENE_INFO, GO_GAF, GO_OBO)
    ])
    atomic_csv(pd.DataFrame(audit_rows), META / "PHASE7_HUMAN_GENE_SET_AUDIT.csv")
    atomic_csv(inputs, META / "PHASE7_HUMAN_GENE_SET_INPUTS.csv")
    atomic_csv(pd.DataFrame([{
        "gate": "PHASE7_HUMAN_GENE_SETS",
        "status": "GREEN",
        "criteria": "Human-native Hallmark/Reactome/GO BP memberships frozen independently of Dog results",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
    }]), META / "PHASE7_HUMAN_GENE_SET_GATE.csv")


if __name__ == "__main__":
    main()
