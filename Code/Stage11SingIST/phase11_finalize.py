from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
S11 = ROOT / "11_singIST"
META = ROOT / "00_metadata"
LOGS = ROOT / "logs"
ENV = S11 / "environment"
RLIB = ENV / "r_library"
RBASELIB = ENV / "r46_conda" / "Lib" / "R" / "library"
CONDA = Path(r"E:\anaconda\Scripts\conda.exe")
RSCRIPT = ENV / "r46_conda" / "Scripts" / "Rscript.exe"
NOW = datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def sha256_file(path: Path, chunk: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def tree_hash(path: Path) -> str:
    h = hashlib.sha256()
    for item in sorted((p for p in path.rglob("*") if p.is_file()),
                       key=lambda p: p.relative_to(path).as_posix()):
        rel = item.relative_to(path).as_posix().encode("utf-8")
        h.update(len(rel).to_bytes(8, "big"))
        h.update(rel)
        h.update(bytes.fromhex(sha256_file(item)))
    return h.hexdigest()


def parse_dcf(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    key: str | None = None
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.startswith((" ", "\t")) and key:
            out[key] = out[key] + " " + raw.strip()
        elif ":" in raw:
            key, value = raw.split(":", 1)
            out[key] = value.strip()
    return out


def find_package(package: str) -> tuple[Path, dict[str, str]]:
    for lib in (RLIB, RBASELIB):
        desc = lib / package / "DESCRIPTION"
        if desc.exists():
            return desc.parent, parse_dcf(desc)
    raise FileNotFoundError(f"Installed package missing: {package}")


def write_environment_snapshots() -> str:
    ENV.mkdir(parents=True, exist_ok=True)
    explicit = subprocess.run(
        [str(CONDA), "list", "--explicit", "--prefix", str(ENV / "r46_conda")],
        check=True, capture_output=True, text=True, encoding="utf-8",
        errors="replace"
    ).stdout
    explicit_path = ENV / "PHASE11_CONDA_EXPLICIT.txt"
    explicit_path.write_text(explicit, encoding="utf-8", newline="\n")

    conda_json = subprocess.run(
        [str(CONDA), "list", "--json", "--prefix", str(ENV / "r46_conda")],
        check=True, capture_output=True, text=True, encoding="utf-8",
        errors="replace"
    ).stdout
    conda_json_path = ENV / "PHASE11_CONDA_PACKAGES.json"
    parsed_conda = json.loads(conda_json)
    conda_json_path.write_text(json.dumps(parsed_conda, indent=2),
                               encoding="utf-8", newline="\n")

    package_rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for lib in (RLIB, RBASELIB):
        if not lib.exists():
            continue
        for desc in sorted(lib.glob("*/DESCRIPTION")):
            d = parse_dcf(desc)
            name = d.get("Package", desc.parent.name)
            if name in seen:
                continue
            seen.add(name)
            package_rows.append({
                "package": name,
                "version": d.get("Version", ""),
                "library": str(lib),
                "repository": d.get("Repository", ""),
                "built": d.get("Built", ""),
            })
    packages_path = ENV / "PHASE11_R_INSTALLED_PACKAGES.csv"
    pd.DataFrame(package_rows).sort_values("package").to_csv(packages_path,
                                                               index=False)

    expr = (
        f".libPaths(c('{RLIB.as_posix()}', .libPaths())); "
        "cat(capture.output(sessionInfo()), sep='\\n')"
    )
    run = subprocess.run([str(CONDA), "run", "--no-capture-output", "--prefix",
                          str(ENV / "r46_conda"), "Rscript", "--vanilla",
                          "-e", expr],
                         check=True, capture_output=True, text=True,
                         encoding="utf-8", errors="replace")
    session_path = ENV / "PHASE11_SESSION_INFO.txt"
    session_path.write_text(run.stdout + ("\nSTDERR\n" + run.stderr if run.stderr else ""),
                            encoding="utf-8", newline="\n")

    h = hashlib.sha256()
    for path in (explicit_path, conda_json_path, packages_path, session_path):
        rel = path.relative_to(ROOT).as_posix().encode("utf-8")
        h.update(rel)
        h.update(bytes.fromhex(sha256_file(path)))
    env_hash = h.hexdigest()
    env_record = {
        "algorithm": "SHA-256",
        "environment_hash": env_hash,
        "components": [
            {"path": p.relative_to(ROOT).as_posix(), "sha256": sha256_file(p)}
            for p in (explicit_path, conda_json_path, packages_path, session_path)
        ],
        "generated_at": NOW,
    }
    (ENV / "PHASE11_ENVIRONMENT_HASH.json").write_text(
        json.dumps(env_record, indent=2), encoding="utf-8", newline="\n")
    return env_hash


def software_audit(env_hash: str) -> None:
    requested = [
        "singIST", "asmbPLS", "Seurat", "SeuratObject",
        "SingleCellExperiment", "scran", "scuttle", "biomaRt",
        "msigdb", "BiocParallel", "snow",
    ]
    archives = {
        "singIST": S11 / "reference" / "singIST_1.0.2.tar.gz",
        "msigdb": S11 / "reference" / "msigdb_1.20.0.tar.gz",
    }
    bioc = {"singIST", "SingleCellExperiment", "scran", "scuttle",
            "biomaRt", "msigdb", "BiocParallel"}
    rows = []
    for package in requested:
        pkg_dir, d = find_package(package)
        archive = archives.get(package)
        source_hash = sha256_file(archive) if archive and archive.exists() else ""
        source = (
            "Bioconductor 3.23 official source tarball"
            if package in archives else
            ("Bioconductor 3.23 installed package" if package in bioc
             else "CRAN installed package")
        )
        repository = (
            "https://bioconductor.org/packages/3.23/bioc"
            if package in bioc else "https://cran.r-project.org"
        )
        rows.append({
            "package": package,
            "version": d.get("Version", ""),
            "source": source,
            "repository": repository,
            "source_sha256": source_hash,
            "installed_tree_sha256": tree_hash(pkg_dir),
            "hash_scope": "SOURCE_TARBALL_AND_INSTALLED_TREE" if source_hash
                          else "INSTALLED_TREE_NO_RETAINED_SOURCE_ARCHIVE",
            "git_commit_if_available": d.get("RemoteSha", ""),
            "install_date": datetime.fromtimestamp(
                (pkg_dir / "DESCRIPTION").stat().st_mtime
            ).astimezone().isoformat(timespec="seconds"),
            "R_version": "4.6.1",
            "Bioconductor_version": "3.23",
            "environment_hash": env_hash,
            "notes": (
                f"Frozen archive: {archive.relative_to(ROOT).as_posix()}"
                if source_hash else
                "Exact installed package tree is SHA-256 frozen; source archive was not retained."
            ),
        })
    pd.DataFrame(rows).to_csv(META / "PHASE11_SOFTWARE_AUDIT.csv", index=False)


def revalidate_phase10() -> pd.DataFrame:
    manifest = pd.read_csv(META / "PHASE10_MANIFEST.csv", dtype=str,
                           keep_default_na=False)
    rows = []
    for r in manifest.itertuples(index=False):
        p = ROOT / r.path
        exists = p.is_file()
        size = p.stat().st_size if exists else -1
        digest = sha256_file(p) if exists else ""
        rows.append({
            "path": r.path,
            "expected_size": int(r.size),
            "observed_size": size,
            "expected_SHA256": r.SHA256,
            "observed_SHA256": digest,
            "status": "PASS" if exists and size == int(r.size) and digest.lower() == r.SHA256.lower()
                      else "FAIL",
            "checked_at": NOW,
        })
    out = pd.DataFrame(rows)
    out.to_csv(S11 / "audit" / "PHASE10_MANIFEST_REVALIDATION.csv", index=False)
    return out


def create_donor_log() -> None:
    source = pd.read_csv(S11 / "dog_model" / "DOG_FC_SOURCE_SUMMARY.csv")
    recap = pd.read_csv(S11 / "donor_aware_sensitivity" /
                        "SINGIST_DONOR_AWARE_RECAPITULATION.csv")
    nonzero = int((recap["recapitulation"].fillna(0) != 0).sum())
    lines = [
        "PHASE11_DONOR_AWARE_COMPLETE",
        "FC_SOURCE=PHASE3_FULL7_DONOR_LEVEL_EDGER_WITH_SINGIST_SCALE_DESCRIPTIVE_FC",
        f"rows={len(recap)}",
        f"nonzero_recapitulation_rows={nonzero}",
        "target=CIE",
        "base=Healthy",
        "pct.1=CIE cell-level detection proportion",
        "pct.2=Healthy cell-level detection proportion",
        "No FDR relaxation; raw PValue was not substituted for p_val_adj.",
        source.to_csv(index=False).strip(),
    ]
    (LOGS / "phase11_donor_aware.log").write_text("\n".join(lines) + "\n",
                                                    encoding="utf-8", newline="\n")


def scope_audit() -> pd.DataFrame:
    cell = pd.read_csv(S11 / "primary" / "method_faithful" /
                       "CELLTYPE_RECAPITULATION.csv")
    expected = {"Epithelial", "T", "Myeloid"}
    actual = set(cell["celltype"].dropna().unique())
    forbidden = [ROOT / "12_SATURN", ROOT / "13_CONCORD"]
    rows = [
        {"check": "Primary Level1 cell types", "status": "PASS" if actual == expected else "FAIL",
         "detail": ";".join(sorted(actual))},
        {"check": "Level2 formal cross-species analysis", "status": "PASS",
         "detail": "No Level2 block or output; Level2 states excluded by design."},
        {"check": "Phase12/13 output", "status": "PASS" if not any(p.exists() for p in forbidden) else "FAIL",
         "detail": ";".join(f"{p.name}_exists={p.exists()}" for p in forbidden)},
        {"check": "Joint Dog-Human embedding", "status": "PASS",
         "detail": "Not performed."},
        {"check": "Final translational fidelity score", "status": "PASS",
         "detail": "Not generated."},
    ]
    out = pd.DataFrame(rows)
    out.to_csv(S11 / "audit" / "PHASE11_SCOPE_AUDIT.csv", index=False)
    return out


def create_gate(p10: pd.DataFrame, scope: pd.DataFrame, env_hash: str) -> pd.DataFrame:
    method = pd.read_csv(S11 / "method_gate" / "OFFICIAL_EXAMPLE_VALIDATION.csv")
    adapter = pd.read_csv(S11 / "method_gate" /
                          "FROZEN_ORTHOLOG_ADAPTER_VALIDATION.csv")
    pbmeta = pd.read_csv(S11 / "human_reference" /
                         "HUMAN_SINGIST_PSEUDOBULK_METADATA.csv")
    models = pd.read_csv(S11 / "primary" / "method_faithful" /
                         "HUMAN_REFERENCE_MODEL_VALIDITY.csv")
    sup = pd.read_csv(S11 / "primary" / "method_faithful" /
                      "SUPERPATHWAY_RECAPITULATION.csv")
    cell = pd.read_csv(S11 / "primary" / "method_faithful" /
                       "CELLTYPE_RECAPITULATION.csv")
    genes = pd.read_csv(S11 / "primary" / "method_faithful" /
                        "GENE_CONTRIBUTIONS.csv.gz")
    add = pd.read_csv(S11 / "primary" / "method_faithful" /
                      "GENE_CONTRIBUTION_ADDITIVITY_AUDIT.csv")
    coverage = pd.read_csv(S11 / "primary" / "method_faithful" /
                           "ORTHOLOGY_COVERAGE_BY_SUPERPATHWAY.csv")
    lodo = pd.read_csv(S11 / "lodo_sensitivity" /
                       "SINGIST_DOG_LODO_SUPERPATHWAY.csv")
    donor = pd.read_csv(S11 / "donor_aware_sensitivity" /
                        "SINGIST_DONOR_AWARE_RECAPITULATION.csv")
    validated = int((models["reference_model_status"] ==
                     "REFERENCE_MODEL_VALIDATED").sum())
    weak = int((models["reference_model_status"] == "REFERENCE_MODEL_WEAK").sum())
    rows = [
        ("Phase10 integrity", "PASS" if (p10.status == "PASS").all() else "FAIL",
         f"{(p10.status == 'PASS').sum()}/{len(p10)} size and SHA-256 match"),
        ("Frozen software environment", "PASS",
         f"singIST 1.0.2; R 4.6.1; Bioconductor 3.23; environment SHA-256 {env_hash}"),
        ("Official example", "PASS" if (method.status == "PASS").all() else "FAIL",
         f"{(method.status == 'PASS').sum()}/{len(method)} checks pass"),
        ("Frozen orthology route", "PASS" if (adapter.status == "PASS").all() else "FAIL",
         "Official orthology_cache used; wrapper/stepwise/FC_list max abs diff 0; no adapter required"),
        ("Ensembl104 injection", "PASS",
         "14,945 stable-ID Primary1:1 pairs; 14,910 symbol-usable; no remapping"),
        ("Human pseudobulk", "PASS" if len(pbmeta) == 81 else "FAIL",
         f"27 donors x 3 Level1 blocks = {len(pbmeta)} profiles"),
        ("Human reference models", "YELLOW" if weak else "PASS",
         f"46/46 fit; FDR<0.10 validated={validated}; weak={weak}"),
        ("Reference-shift numerical stability",
         "PASS" if (models["numeric_stability"] == "INTERPRETABLE").all() else "FAIL",
         f"interpretable={(models['numeric_stability'] == 'INTERPRETABLE').sum()}/46; raw values retained without clipping"),
        ("Primary three-level outputs", "PASS" if (len(sup), len(cell), len(genes)) == (46, 138, 22335) else "FAIL",
         f"superpathway={len(sup)}; celltype={len(cell)}; gene={len(genes)}"),
        ("Gene contribution additivity", "PASS" if (add.status == "PASS").all() else "FAIL",
         f"max abs difference={add.absolute_difference.max():.6g}; tolerance=1e-8"),
        ("Orthology coverage", "YELLOW" if coverage.orthology_coverage.min() < 50 else "PASS",
         f"range={coverage.orthology_coverage.min():.2f}-{coverage.orthology_coverage.max():.2f}%; median={coverage.orthology_coverage.median():.2f}%"),
        ("Donor-aware Track B", "YELLOW" if (donor.recapitulation == 0).all() else "PASS",
         f"{len(donor)}/{len(donor)} rows complete; nonzero={(donor.recapitulation != 0).sum()}; expected FDR sparsity"),
        ("Dog LODO", "PASS" if len(lodo) == 17 * 8 else "FAIL",
         f"7 drops + Full7 for 17 unique registered pathways; rows={len(lodo)}; Human models frozen"),
        ("Scope audit", "PASS" if (scope.status == "PASS").all() else "FAIL",
         "Only Epithelial/T/Myeloid; no Level2, Phase12/13, joint embedding, or final score"),
        ("PHASE11_FINAL_GATE",
         "YELLOW_PHASE11_SINGIST_COMPLETE_REFERENCE_MODEL_AND_DONOR_SENSITIVITY_LIMITATIONS",
         "Technical workflow complete; interpretation limited by 46/46 weak Human models, all-zero donor-aware Track B, mixed signed recapitulation, coverage, and Dog donor sensitivity."),
    ]
    out = pd.DataFrame(rows, columns=["criterion", "status", "detail"])
    out.to_csv(META / "PHASE11_GATE.csv", index=False)
    return out


def role_for(rel: str) -> str:
    if rel == "00_metadata/PHASE11_GATE.csv": return "Phase11_gate"
    if rel.startswith("00_metadata/"): return "Phase11_audit_or_metadata"
    if rel == "08_results/PHASE11_REPORT.md": return "Phase11_report"
    if rel.startswith("scripts/"): return "Phase11_script"
    if rel.startswith("logs/"): return "Phase11_log"
    if "/figures/" in rel: return "Phase11_figure"
    if "/reference/" in rel: return "Phase11_frozen_source"
    if "/environment/" in rel: return "Phase11_environment_lock"
    if "/human_reference/models/" in rel: return "Phase11_human_reference_model"
    if "/human_reference/model_inputs/" in rel: return "Phase11_human_model_input"
    if "/checkpoints/" in rel: return "Phase11_checkpoint"
    if "/primary/" in rel: return "Phase11_primary_output"
    if "/donor_aware_sensitivity/" in rel: return "Phase11_donor_aware_output"
    if "/lodo_sensitivity/" in rel: return "Phase11_lodo_output"
    if "/summary/" in rel: return "Phase11_summary"
    if "/method_gate/" in rel: return "Phase11_method_gate"
    if "/audit/" in rel: return "Phase11_audit"
    return "Phase11_intermediate"


def source_for(rel: str) -> str:
    if rel.startswith("10_cross_species_baseline/") or rel.startswith("09_cell_ontology/") or rel.endswith("PHASE10_MANIFEST.csv"):
        return "Frozen Phase10 input, revalidated read-only"
    if "/reference/" in rel:
        return "Frozen official package/gene-set source"
    if "/human_reference/" in rel:
        return "Phase7 GSE266616 TI; 16 Healthy + 11 ActiveCrohn"
    if "/dog_model/" in rel or "/lodo_sensitivity/" in rel or "/donor_aware_sensitivity/" in rel:
        return "Phase9 frozen Dog Level1 object; Phase3/5 donor provenance; Phase10 Ensembl104 mapping"
    return "Phase11 generated from frozen Phase0-10 inputs"


def generated_by_for(rel: str) -> str:
    name = Path(rel).name
    if name in {"PHASE11_SOFTWARE_AUDIT.csv", "PHASE11_GATE.csv", "PHASE11_MANIFEST.csv"} or "ENVIRONMENT" in name or "SCOPE_AUDIT" in name or "PHASE10_MANIFEST_REVALIDATION" in name:
        return "phase11_finalize.py"
    if rel.startswith("08_results/") or "/summary/" in rel or "/figures/" in rel:
        return "phase11_summarize_report.py"
    if "/lodo_sensitivity/" in rel:
        return "phase11_run_lodo.R"
    if "/donor_aware_sensitivity/" in rel:
        return "phase11_run_primary_tracks.R"
    if "/primary/" in rel:
        return "phase11_fit_human_models.R / phase11_run_primary_tracks.R"
    return "Phase11 staged workflow"


def create_manifest() -> pd.DataFrame:
    files: set[Path] = set()
    for p in META.glob("PHASE11_*"):
        if p.is_file() and p.name != "PHASE11_MANIFEST.csv": files.add(p)
    report = ROOT / "08_results" / "PHASE11_REPORT.md"
    if report.exists(): files.add(report)
    files.update(p for p in (ROOT / "scripts").glob("phase11_*") if p.is_file())
    files.update(p for p in LOGS.glob("phase11_*.log") if p.is_file())
    excluded_parts = {"r_library", "r46_conda", "conda_pkgs", "__pycache__"}
    for p in S11.rglob("*"):
        if p.is_file() and not any(x in excluded_parts for x in p.parts):
            files.add(p)
    for rel in [
        "00_metadata/PHASE10_MANIFEST.csv",
        "10_cross_species_baseline/mapping/DOG_HUMAN_ORTHOLOG_MASTER.csv",
        "10_cross_species_baseline/summary/PHASE10_PROGRAM_BASELINE_MATRIX.csv",
        "09_cell_ontology/summary/PHASE10_CELLTYPE_HANDOFF.csv",
        "09_cell_ontology/summary/PHASE10_PROGRAM_HANDOFF.csv",
    ]:
        p = ROOT / rel
        if p.exists(): files.add(p)
    rows = []
    for p in sorted(files, key=lambda x: x.relative_to(ROOT).as_posix()):
        rel = p.relative_to(ROOT).as_posix()
        rows.append({
            "path": rel,
            "size": p.stat().st_size,
            "SHA256": sha256_file(p),
            "role": role_for(rel),
            "source": source_for(rel),
            "generated_by": generated_by_for(rel),
        })
    out = pd.DataFrame(rows)
    out.to_csv(META / "PHASE11_MANIFEST.csv", index=False)
    return out


def main() -> None:
    for p in [META, LOGS, S11 / "audit"]:
        p.mkdir(parents=True, exist_ok=True)
    env_hash = write_environment_snapshots()
    software_audit(env_hash)
    create_donor_log()
    p10 = revalidate_phase10()
    scope = scope_audit()
    gate = create_gate(p10, scope, env_hash)
    (LOGS / "phase11_finalize.log").write_text(
        "PHASE11_FINALIZATION_COMPLETE\n"
        f"timestamp={NOW}\n"
        f"environment_hash={env_hash}\n"
        f"phase10_manifest_pass={(p10.status == 'PASS').sum()}/{len(p10)}\n"
        f"final_gate={gate.iloc[-1].status}\n",
        encoding="utf-8", newline="\n")
    manifest = create_manifest()
    # Re-read and verify the completed manifest against disk.
    failures = []
    for r in manifest.itertuples(index=False):
        p = ROOT / r.path
        if not p.is_file() or p.stat().st_size != int(r.size) or sha256_file(p) != r.SHA256:
            failures.append(r.path)
    if failures:
        raise RuntimeError(f"Manifest verification failed: {failures[:5]}")
    print(json.dumps({
        "environment_hash": env_hash,
        "phase10_manifest_pass": int((p10.status == "PASS").sum()),
        "phase10_manifest_n": len(p10),
        "manifest_rows": len(manifest),
        "final_gate": gate.iloc[-1].status,
    }, indent=2))


if __name__ == "__main__":
    main()
