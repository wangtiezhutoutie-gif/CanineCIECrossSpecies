# Canine CIE–human Crohn cross-species single-cell analysis code

This repository is the curated code and environment-metadata package for a cross-species comparison of canine chronic inflammatory enteropathy (CIE) and human Crohn disease. It supports lineage- and program-specific analyses; it does not claim global disease equivalence.

## Repository scope

Included:

- `Code/`: 147 Python, R, PowerShell and JavaScript analysis scripts, plus four shared configuration modules, organized as Stage00–Stage16 plus human-subtype, matched-background and figure-generation extensions.
- `Environments/`: stage-specific environment records retained from the analysis. They are not presented as one universal environment.
- `MethodsCodeMapping.csv`: Methods-to-code, input, output, environment and display-item mapping.
- `FigureSourceMapping.csv`: figure-panel plotting and analysis provenance.
- `ReproductionGuide.md`: required inputs, execution order, expected outputs and verification limits.

Not included:

- the manuscript, internal audits, review drafts, logs, backups or invalidated runs;
- raw or complete expression matrices, H5AD/loom objects, pseudobulk result collections, model weights or full analysis outputs;
- journal Source Data files or frozen figure/table data. Those files remain controlled separately and must be obtained from the journal deposit or the data owners when released.

## Public datasets

The workflow uses public repository inputs, including:

- canine CIE: GEO `GSE254005`;
- human Crohn disease: GEO `GSE266616`;
- healthy anatomical-region evidence: GEO `GSE185224`;
- coeliac-disease context: GEO `GSE277276`;
- disease-context sensitivity: GEO `GSE190564`.

Download data from the authoritative repositories. Do not substitute normalized expression where a script requires raw counts, and preserve downloaded raw data read-only.

## Configuration

Set `CANINE_CIE_PROJECT_ROOT` to a writable analysis root containing the expected input layout before running scripts. The variable replaces the original workstation path without changing analysis algorithms.

Optional script-specific settings are:

- `CANINE_CIE_QA_ROOT`: output directory for the Phase 7 workbook visual QA artifacts;
- `CANINE_CIE_PHASE10_REVALIDATION`: optional cached Phase 9 revalidation CSV;
- `CANINE_CIE_PYTHON`: Python executable used by the PowerShell transfer wrapper.

Example in PowerShell:

```powershell
$env:CANINE_CIE_PROJECT_ROOT = 'D:\path\to\analysis-root'
```

The code expects the historical analysis-root layout described by script inputs and `MethodsCodeMapping.csv`; cloning this repository alone does not provide the scientific inputs.

## Reproduction order

Follow [ReproductionGuide.md](ReproductionGuide.md). In summary: obtain the public inputs, select the stage-specific environment, configure the project root, then run Stage01–Stage16 in order. Stage00 is a documented preflight checkpoint, not a packaged executable. Run the human-subtype and matched-background extensions only where the mapping requires them, and generate figures only from authorized frozen source data.

## Verification scope

The release process performs static, low-cost checks only: Python compilation, R parsing when a compatible R executable is available, PowerShell/JavaScript syntax checks, environment-file parsing, mapping-path checks, secret scanning and Git inclusion checks. It does not rerun the full single-cell workflow, recompute the 25 frozen programs, retrain models or regenerate frozen statistics.

`PASS` values in the mapping and environment records describe the pre-existing project provenance checks named in those files. They do not mean that all 147 scripts were rerun in one environment for this repository release. Missing data, archived environments, GPU requirements and external services remain `NOT VERIFIED` unless explicitly tested.

## Scientific boundaries

- Donor, not cell, is the biological replicate.
- Canine biological replication is limited to 3 healthy and 4 CIE donors.
- Dog samples are primarily duodenum; the main human Crohn reference is terminal ileum.
- singIST human reference models were weak, and SATURN T-cell alignment was seed/representation sensitive.
- Myeloid correspondence can reflect generic mucosal inflammation rather than Crohn specificity.
- TNF is a translational benchmark, not a newly established canine CIE therapy.

## License and release status

No authorized `LICENSE` is currently present. Until the code owners confirm project-code licensing and third-party redistribution rights and add an appropriate license, no reuse permission is granted by this package and public release remains blocked. See [RepositoryReleaseGuide.md](RepositoryReleaseGuide.md).
