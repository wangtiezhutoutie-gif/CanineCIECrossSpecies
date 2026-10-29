# Environment records

These files preserve stage-specific environment evidence; they are not one tested universal environment.

- `phase1_environment.yml` and `phase2_environment.yml`: Conda environment records for the early canine stages.
- `phase3_environment.txt`: Phase 3 environment record.
- `PHASE11_CONDA_EXPLICIT.txt`, `PHASE11_CONDA_PACKAGES.json`, `PHASE11_R_INSTALLED_PACKAGES.csv`, `PHASE11_SESSION_INFO.txt` and `PHASE11_ENVIRONMENT_HASH.json`: the archived Phase 11 Conda/R environment inventory and integrity metadata.
- `pip-freeze-phase13.txt` and `requirements-phase13.txt`: Phase 13 Python package records.
- `gpu-audit.txt`: retained GPU environment evidence for the model-dependent workflow.

Paths embedded in explicit package URLs, R library inventories or build metadata record the historical analysis environment. They are not executable path defaults and do not replace `CANINE_CIE_PROJECT_ROOT`.

The packaging host record is deliberately excluded because it is not an analysis environment. Exact reproduction still depends on the public scientific inputs, frozen intermediate inputs where authorized, compatible R/Bioconductor repositories, external services and—in model stages—the archived GPU/software stack. A present environment record is not evidence that all scripts were rerun during repository release.

## Static dependency inventory

A release-time import scan found the following principal third-party dependencies across the packaged scripts:

- Python: Pillow, anndata, h5py, joblib, matplotlib, NumPy, pandas, psutil, Scanpy, SciPy, Scrublet, seaborn, scikit-learn, PyTorch and umap-learn, plus the CONCORD and SATURN `score_adata` code used in their respective archived environments.
- R: babelgene, BiocParallel, data.table, digest, edgeR, fgsea, ggplot2, limma, Matrix, msigdbr, S4Vectors, scuttle, singIST, SingleCellExperiment and SummarizedExperiment.
- JavaScript: `@oai/artifact-tool` for workbook processing and `sharp` for raster export.
- PowerShell: Windows BITS cmdlets for the staged transfer wrappers.

The stage records provide versions only where captured by the original analysis environment. In particular, JavaScript package versions and the upstream SATURN `score_adata` module are not pinned by this package and remain `NOT VERIFIED` for clean-room execution.
