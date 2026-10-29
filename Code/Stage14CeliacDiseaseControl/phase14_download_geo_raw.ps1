param(
    [int]$SegmentCount = 8,
    [int]$ThrottleLimit = 3,
    [string]$PinnedIPv4 = '130.14.250.10'
)

$ErrorActionPreference = 'Stop'
$projectRoot = $env:CANINE_CIE_PROJECT_ROOT
if ([string]::IsNullOrWhiteSpace($projectRoot)) { throw 'Set CANINE_CIE_PROJECT_ROOT' }
$phaseRoot = Join-Path $projectRoot '14_ced_specificity'
$rawRoot = Join-Path $phaseRoot 'raw_processed'
$partsRoot = Join-Path $rawRoot 'GSE277276_RAW.parts'
$tarPath = Join-Path $rawRoot 'GSE277276_RAW.tar'
$extractRoot = Join-Path $rawRoot 'GSE277276_cellbender'
$logPath = Join-Path $projectRoot 'logs\phase14_download_geo_raw.log'
$url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277276/suppl/GSE277276_RAW.tar'
$totalBytes = [int64]1688524800

foreach ($directory in @($rawRoot, $partsRoot, $extractRoot, (Split-Path -Parent $logPath))) {
    New-Item -ItemType Directory -Path $directory -Force | Out-Null
}

$segmentSize = [int64][math]::Ceiling($totalBytes / [double]$SegmentCount)
$segments = foreach ($index in 0..($SegmentCount - 1)) {
    $start = [int64]$index * $segmentSize
    $end = [math]::Min($totalBytes - 1, $start + $segmentSize - 1)
    [pscustomobject]@{
        index = $index
        start = [int64]$start
        end = [int64]$end
        expected = [int64]($end - $start + 1)
        path = Join-Path $partsRoot ('part_{0:D2}.bin' -f $index)
    }
}

$header = @(
    "started_utc=$([DateTime]::UtcNow.ToString('o'))",
    "url=$url",
    "pinned_ipv4=$PinnedIPv4",
    "segments=$SegmentCount",
    "throttle_limit=$ThrottleLimit",
    "total_bytes=$totalBytes"
)
$header | Set-Content -LiteralPath $logPath -Encoding utf8

$results = $segments | ForEach-Object -Parallel {
    $segment = $_
    $downloadUrl = $using:url
    $serverIp = $using:PinnedIPv4
    $hostName = 'ftp.ncbi.nlm.nih.gov'
    $partPath = $segment.path
    $expected = [int64]$segment.expected
    $attempt = 0
    while ((Get-Item -LiteralPath $partPath -ErrorAction SilentlyContinue).Length -lt $expected) {
        $existingItem = Get-Item -LiteralPath $partPath -ErrorAction SilentlyContinue
        $existing = if ($null -eq $existingItem) { [int64]0 } else { [int64]$existingItem.Length }
        if ($existing -gt $expected) { throw "Segment $($segment.index) exceeds expected size" }
        $offset = [int64]$segment.start + $existing
        $chunkPath = "$partPath.chunk"
        if (Test-Path -LiteralPath $chunkPath) { Remove-Item -LiteralPath $chunkPath -Force }
        $attempt++
        & curl.exe -4 --resolve "$hostName`:443`:$serverIp" --range "$offset-$($segment.end)" `
            --connect-timeout 30 --max-time 2400 --speed-time 120 --speed-limit 1024 --fail --silent --show-error `
            --output $chunkPath $downloadUrl
        $curlExit = $LASTEXITCODE
        $chunkItem = Get-Item -LiteralPath $chunkPath -ErrorAction SilentlyContinue
        $chunkBytes = if ($null -eq $chunkItem) { [int64]0 } else { [int64]$chunkItem.Length }
        if ($chunkBytes -gt 0) {
            $destination = [IO.File]::Open($partPath, [IO.FileMode]::Append, [IO.FileAccess]::Write, [IO.FileShare]::Read)
            try {
                $source = [IO.File]::OpenRead($chunkPath)
                try { $source.CopyTo($destination) } finally { $source.Dispose() }
            } finally { $destination.Dispose() }
            Remove-Item -LiteralPath $chunkPath -Force
        }
        $now = (Get-Item -LiteralPath $partPath -ErrorAction SilentlyContinue).Length
        "segment=$($segment.index) attempt=$attempt curl_exit=$curlExit appended=$chunkBytes bytes=$now expected=$expected"
        if ($now -lt $expected) {
            if ($chunkBytes -eq 0) { Start-Sleep -Seconds ([math]::Min(30, 2 * $attempt)) }
            else { Start-Sleep -Seconds 2 }
        }
    }
    $finalSize = (Get-Item -LiteralPath $partPath).Length
    if ($finalSize -ne $expected) { throw "Segment $($segment.index) final size mismatch" }
    "segment=$($segment.index) status=PASS bytes=$finalSize"
} -ThrottleLimit $ThrottleLimit

$results | Add-Content -LiteralPath $logPath -Encoding utf8

if (Test-Path -LiteralPath $tarPath) {
    $existingTar = Get-Item -LiteralPath $tarPath
    if ($existingTar.Length -ne $totalBytes) {
        throw "Existing RAW tar has unexpected size $($existingTar.Length); refusing overwrite"
    }
} else {
    $output = [IO.File]::Open($tarPath, [IO.FileMode]::CreateNew, [IO.FileAccess]::Write, [IO.FileShare]::Read)
    try {
        foreach ($segment in $segments | Sort-Object index) {
            $input = [IO.File]::OpenRead($segment.path)
            try { $input.CopyTo($output) } finally { $input.Dispose() }
        }
    } finally { $output.Dispose() }
}

$tarSize = (Get-Item -LiteralPath $tarPath).Length
if ($tarSize -ne $totalBytes) { throw "Combined RAW tar size $tarSize does not equal $totalBytes" }
& tar.exe -xf $tarPath -C $extractRoot
if ($LASTEXITCODE -ne 0) { throw "tar extraction failed with exit code $LASTEXITCODE" }
$files = Get-ChildItem -LiteralPath $extractRoot -File
$summary = [pscustomobject]@{
    completed_utc = [DateTime]::UtcNow.ToString('o')
    tar_bytes = $tarSize
    extracted_files = $files.Count
    extracted_bytes = ($files | Measure-Object Length -Sum).Sum
}
$summary | ConvertTo-Json -Compress | Add-Content -LiteralPath $logPath -Encoding utf8
$summary | ConvertTo-Json -Compress
