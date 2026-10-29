param(
    [Parameter(Mandatory = $true)]
    [int]$DownstreamPid
)
$ErrorActionPreference = 'Stop'
$project = $env:CANINE_CIE_PROJECT_ROOT
if ([string]::IsNullOrWhiteSpace($project)) { throw 'Set CANINE_CIE_PROJECT_ROOT' }
$root = Join-Path $project '12_saturn'
$python = Join-Path $root 'environment\saturn_gpu_conda\python.exe'
$log = Join-Path $project 'logs\phase12_figure_watcher.log'
"UTC_START=$([DateTime]::UtcNow.ToString('o'))" | Set-Content -LiteralPath $log -Encoding UTF8
"WAITING_FOR_DOWNSTREAM_PID=$downstreamPid" | Add-Content -LiteralPath $log
while (Get-Process -Id $downstreamPid -ErrorAction SilentlyContinue) {
    Start-Sleep -Seconds 30
}
$downstreamLog = Join-Path $project 'logs\phase12_downstream_watcher.log'
if (-not (Select-String -LiteralPath $downstreamLog -SimpleMatch 'PHASE12_POSTHOC_COMPLETE=' -Quiet)) {
    "STOP_DOWNSTREAM_GATE_NOT_COMPLETE" | Add-Content -LiteralPath $log
    exit 41
}
if (Test-Path -LiteralPath (Join-Path $project '13_concord')) {
    "STOP_PHASE13_DIRECTORY_PRESENT" | Add-Content -LiteralPath $log
    exit 42
}
"START_FIGURES=$([DateTime]::UtcNow.ToString('o'))" | Add-Content -LiteralPath $log
& $python (Join-Path $project 'scripts\phase12_saturn_figures.py') *>> $log
if ($LASTEXITCODE -ne 0) {
    "STOP_FIGURES_EXIT=$LASTEXITCODE" | Add-Content -LiteralPath $log
    exit 43
}
"PHASE12_FIGURES_COMPLETE=$([DateTime]::UtcNow.ToString('o'))" | Add-Content -LiteralPath $log
exit 0
