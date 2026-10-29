param(
    [Parameter(Mandatory = $true)]
    [int]$UpstreamWatcherPid
)
$ErrorActionPreference = 'Stop'
$project = $env:CANINE_CIE_PROJECT_ROOT
if ([string]::IsNullOrWhiteSpace($project)) { throw 'Set CANINE_CIE_PROJECT_ROOT' }
$root = Join-Path $project '12_saturn'
$python = Join-Path $root 'environment\saturn_gpu_conda\python.exe'
$log = Join-Path $project 'logs\phase12_downstream_watcher.log'

"UTC_START=$([DateTime]::UtcNow.ToString('o'))" | Set-Content -LiteralPath $log -Encoding UTF8
"WAITING_FOR_ESM_WATCHER_PID=$upstreamWatcherPid" | Add-Content -LiteralPath $log
while (Get-Process -Id $upstreamWatcherPid -ErrorAction SilentlyContinue) {
    Start-Sleep -Seconds 30
}
"DOG_COUNT=$((Get-ChildItem -LiteralPath (Join-Path $root 'protein_embeddings\dog_protein_pt') -Filter *.pt -File).Count)" | Add-Content -LiteralPath $log
"HUMAN_COUNT=$((Get-ChildItem -LiteralPath (Join-Path $root 'protein_embeddings\human_protein_pt') -Filter *.pt -File).Count)" | Add-Content -LiteralPath $log
$dogCount = (Get-ChildItem -LiteralPath (Join-Path $root 'protein_embeddings\dog_protein_pt') -Filter *.pt -File).Count
$humanCount = (Get-ChildItem -LiteralPath (Join-Path $root 'protein_embeddings\human_protein_pt') -Filter *.pt -File).Count
if ($dogCount -ne 18236 -or $humanCount -ne 34791) {
    "STOP_ESM_COUNT_GATE dog=$dogCount human=$humanCount" | Add-Content -LiteralPath $log
    exit 31
}
if (Test-Path -LiteralPath (Join-Path $project '13_concord')) {
    "STOP_PHASE13_DIRECTORY_PRESENT" | Add-Content -LiteralPath $log
    exit 32
}

"START_GENE_EMBEDDING_FINALIZER=$([DateTime]::UtcNow.ToString('o'))" | Add-Content -LiteralPath $log
& $python (Join-Path $project 'scripts\phase12_finalize_gene_embeddings.py') *>> $log
if ($LASTEXITCODE -ne 0) {
    "STOP_GENE_EMBEDDING_FINALIZER_EXIT=$LASTEXITCODE" | Add-Content -LiteralPath $log
    exit 33
}

"START_SATURN_GATED_RUNNER=$([DateTime]::UtcNow.ToString('o'))" | Add-Content -LiteralPath $log
& $python (Join-Path $project 'scripts\phase12_run_saturn_gated.py') --mode all *>> $log
if ($LASTEXITCODE -ne 0) {
    "STOP_SATURN_GATED_RUNNER_EXIT=$LASTEXITCODE" | Add-Content -LiteralPath $log
    exit 34
}

"START_POSTHOC=$([DateTime]::UtcNow.ToString('o'))" | Add-Content -LiteralPath $log
& $python (Join-Path $project 'scripts\phase12_saturn_posthoc.py') *>> $log
if ($LASTEXITCODE -ne 0) {
    "STOP_POSTHOC_EXIT=$LASTEXITCODE" | Add-Content -LiteralPath $log
    exit 35
}
if (Test-Path -LiteralPath (Join-Path $project '13_concord')) {
    "STOP_PHASE13_DIRECTORY_CREATED_UNEXPECTEDLY" | Add-Content -LiteralPath $log
    exit 36
}
"PHASE12_POSTHOC_COMPLETE=$([DateTime]::UtcNow.ToString('o'))" | Add-Content -LiteralPath $log
exit 0
