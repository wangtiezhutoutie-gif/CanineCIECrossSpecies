# Repository release guide

Public release is permitted only after every gate below is documented as satisfied.

1. Confirm ownership of the project code, third-party redistribution rights and the intended license. Add only a license authorized by the rights holders.
2. Confirm that the manuscript, Source Data, internal audits, logs, backups, raw/complete expression matrices, model weights and invalidated outputs are not tracked.
3. Scan both the proposed worktree and the complete Git history for credentials, private keys, personal paths and sensitive data. A `.gitignore` rule alone is not a history cleanup.
4. Verify that executable paths use `CANINE_CIE_PROJECT_ROOT` or documented script-specific configuration. Historical paths in non-executable environment provenance records must be clearly identified as such.
5. Run the static checks described in `README.md`; mark missing software/data/GPU/external-service tests `NOT VERIFIED`.
6. Confirm every path in `MethodsCodeMapping.csv` and `FigureSourceMapping.csv`, including which scientific inputs are not redistributed.
7. Review the exact staged file list, file sizes and SHA256 manifest before committing.
8. Create or select the intended GitHub repository, verify its owner/name, and avoid force-pushing unknown history.
9. After push, query the remote repository to verify `Public` visibility, `main`, the latest commit, required paths, file count and absence of excluded content.

## Current blocker

This release candidate has no authorized `LICENSE`, and the supplied project materials do not establish project-code ownership or third-party redistribution permission. Do not make the repository public until the rights holders resolve this gate. Do not paste credentials into logs or chat; use GitHub's official browser authentication flow when publishing is authorized.
