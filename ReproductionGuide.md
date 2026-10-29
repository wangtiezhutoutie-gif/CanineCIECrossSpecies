# Reproduction guide

## Scope and verification boundary

This guide describes how to reconstruct the analysis from authorized public inputs. It does not execute the analyses automatically, redistribute the scientific data, or claim that every historical stage can be reproduced in one environment.

The repository release is validated with static and low-cost checks. The complete workflow, GPU models, 25 frozen lineage-by-program results, final statistics and figures are not recomputed during release. If a required dataset, archived environment, GPU, external service or licensed asset is unavailable, record the stage as `NOT VERIFIED`; do not substitute a different input or threshold.

## Required inputs

Obtain the following from their authoritative repositories:

| Role | Accession |
|---|---|
| Canine CIE primary data | `GSE254005` |
| Human Crohn primary reference | `GSE266616` |
| Healthy anatomical-region evidence | `GSE185224` |
| Coeliac-disease context | `GSE277276` |
| Disease-context sensitivity | `GSE190564` |

Several scripts also require frozen intermediate tables, ortholog mappings, gene sets or figure Source Data that are not redistributed in this code repository. Consult `MethodsCodeMapping.csv` and `FigureSourceMapping.csv`. Inputs using `${CANINE_CIE_PROJECT_ROOT}` are logical locations beneath the configured analysis root, not files bundled here.

## Environment setup

1. Review `Environments/EnvironmentManifest.csv` and `Environments/README.md`.
2. Select the record for the stage being reproduced. Do not merge the records into an unverified universal environment.
3. Set `CANINE_CIE_PROJECT_ROOT` to the analysis root. Keep raw inputs read-only and write outputs to the stage locations expected by the scripts.
4. Where applicable, set `CANINE_CIE_QA_ROOT`, `CANINE_CIE_PHASE10_REVALIDATION` or `CANINE_CIE_PYTHON` as described in `README.md`.

The Phase 1 and Phase 2 YAML files are Conda records; Phase 3 and Phase 13 have separate package records; Phase 11 retains its explicit Conda and R package/session information. Historical package prefixes are provenance, not portable executable defaults.

## Preflight (Stage00)

The packaged `Code/Stage00DataAudit/` directory contains documentation only; no Stage00 executable was present in the curated release copy. Before analysis, independently verify accessions, download hashes, raw-versus-normalized matrix roles, donor/sample identifiers and expected file locations. This checkpoint is therefore `NOT VERIFIED` as an executable stage in this repository.

## Execution sequence

1. `Stage01DogQualityControl` — acquire/ingest `GSE254005`, distinguish raw counts from normalized expression and apply the frozen QC workflow. Expected output: audited canine object and QC report.
2. `Stage02DogAnnotation` — reproduce the author/Level1 lineage mapping and annotation audits. Expected output: frozen Epithelial/T/Myeloid annotations.
3. `Stage03DogPseudobulk` — aggregate raw UMI counts by donor and lineage and fit donor-level edgeR models. Expected output: differential-expression tables.
4. `Stage04DogEnrichment` — prepare the defined gene sets and run enrichment. Expected output: lineage/program enrichment tables.
5. `Stage05DogStability` — perform leave-one-dog-out gene/program stability assessment.
6. `Stage06DogExternalValidation` — reproduce the limited external-data availability and summary-level assessment.
7. `Stage07HumanCrohnAnalysis` — analyze `GSE266616` at participant level.
8. `Stage08RegionControl` — assess healthy anatomical-region sensitivity.
9. `Stage09CellOntology` — freeze the ontology crosswalk.
10. `Stage10OrthologComparison` — construct the strict ortholog universe and gene/program comparisons.
11. `Stage11SingIST`, `Stage12SATURN`, `Stage13CONCORD` — run method-specific state analyses using their archived inputs, versions and seeds. GPU/model reruns are not part of release validation.
12. `Stage14CeliacDiseaseControl` — reproduce `GSE277276` disease-context comparisons.
13. `Stage15TargetPrioritization` — use the frozen candidate pools and run005 evidence correction; do not use invalidated run004 as final evidence.
14. `Stage16EvidenceMap` — integrate independent evidence axes without treating confidence as statistical significance.
15. `HumanSubtypeAnalysis` — use the Phase12-corrected human pseudobulk/index chain for affected human epithelial/myeloid subtype results.
16. `MatchedBackgroundCalibration` — use the saved 5,000 K20 draws; do not choose sensitivity exclusions to obtain significance.
17. `FigureGeneration` — render only from authorized frozen source data. The source data are not in this repository.

## Expected checkpoints

The following are frozen expectations to compare against authorized outputs; they were not newly recomputed for repository release:

- 7 canine donors and 38,186 final canine cells;
- 25 frozen lineage×program combinations: Epithelial 8, T 8, Myeloid 9;
- leave-one-dog-out classes: 20 `ROBUST`, 5 `MODERATELY_STABLE`;
- final evidence categories: 10/7/4/3/1 and confidence: 3/10/12;
- Supplementary Figure 1 K20 saved-draw count: 5,000;
- human subtype outputs trace to Phase12-corrected scripts and source tables;
- target prioritization traces to run005.

## Known limitations

Large GPU models should be rerun only when the archived environment, seed and frozen input hashes are available. R package snapshots, external URLs and APIs may no longer resolve identically. Windows-only BITS wrappers require Windows. JavaScript workbook and rendering scripts require their declared packages. Report any unmet prerequisite as `NOT VERIFIED` rather than changing algorithms, thresholds or frozen results.
