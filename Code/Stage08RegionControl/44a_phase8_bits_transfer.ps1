$ErrorActionPreference = 'Stop'

$Root = $env:CANINE_CIE_PROJECT_ROOT
if ([string]::IsNullOrWhiteSpace($Root)) { throw 'Set CANINE_CIE_PROJECT_ROOT' }
$Audit = Join-Path $Root '08_region_control\metadata\HUMAN_AC_REGION_AUDIT.csv'
$RawRoot = Join-Path $Root '08_region_control\human_AC_TI\raw_processed'
$TransferManifest = Join-Path $Root '08_region_control\metadata\HUMAN_AC_BITS_TRANSFER_MANIFEST.csv'
$LogPath = Join-Path $Root 'logs\phase8_region_analysis.log'
$LocalLogPath = Join-Path $Root '08_region_control\logs\phase8_region_analysis.log'
$Forbidden = @('GSE266616_RAW.tar', '.fastq', '.sra', 'scRNA_gene_counts.mtx')

function Write-Phase8Log {
    param([Parameter(Mandatory = $true)][string]$Message)
    $line = '[{0}] {1}' -f ([DateTime]::UtcNow.ToString('o')), $Message
    Write-Output $line
    Add-Content -LiteralPath $LogPath -Value $line -Encoding utf8
    Add-Content -LiteralPath $LocalLogPath -Value $line -Encoding utf8
}

New-Item -ItemType Directory -Force -Path (Split-Path $LogPath -Parent) | Out-Null
New-Item -ItemType Directory -Force -Path (Split-Path $LocalLogPath -Parent) | Out-Null

try {
    $allRows = @(Import-Csv -LiteralPath $Audit)
    $selected = @($allRows | Where-Object { $_.eligible.StartsWith('YES_') })
    if ($selected.Count -ne 22) {
        throw "Expected 22 frozen AC analysis archives, observed $($selected.Count)"
    }
    $totalBytes = ($selected | Measure-Object -Property processed_file_size -Sum).Sum
    if ([int64]$totalBytes -ne 944416204) {
        throw "Frozen AC byte total changed: $totalBytes"
    }
    if (@($selected | Where-Object donor_id -eq 'OR95').Count -ne 0) {
        throw 'OR95 entered the AC transfer scope'
    }

    Write-Phase8Log 'Starting BITS transfer for 22 AC processed archives; OR95/adjacent/noninflamed/RAW/SRA/FASTQ excluded'
    $records = New-Object System.Collections.Generic.List[object]
    $index = 0
    foreach ($row in $selected) {
        $index++
        $sourceUri = [Uri]$row.processed_download_url
        if ($sourceUri.Scheme -ne 'https' -or $sourceUri.Host -ne 'ftp.ncbi.nlm.nih.gov' -or -not $sourceUri.AbsolutePath.Contains('/geo/samples/')) {
            throw "Unexpected source URL: $($row.processed_download_url)"
        }
        $fileName = [Uri]::UnescapeDataString([System.IO.Path]::GetFileName($sourceUri.AbsolutePath))
        if (-not $fileName.EndsWith('.tar.gz', [StringComparison]::OrdinalIgnoreCase)) {
            throw "Refusing non-tar.gz file: $fileName"
        }
        foreach ($token in $Forbidden) {
            if ($fileName.IndexOf($token, [StringComparison]::OrdinalIgnoreCase) -ge 0) {
                throw "Forbidden target encountered: $fileName"
            }
        }
        $cohortDir = if ($row.disease -eq 'Healthy') { 'healthy_AC' } elseif ($row.inflammation_state -eq 'Active_inflamed') { 'crohn_AC_active' } else { throw "Unexpected cohort: $($row.GSM)" }
        $destinationDir = Join-Path $RawRoot $cohortDir
        New-Item -ItemType Directory -Force -Path $destinationDir | Out-Null
        $resolvedDir = (Resolve-Path -LiteralPath $destinationDir).Path
        if (-not $resolvedDir.StartsWith(($RawRoot + '\'), [StringComparison]::OrdinalIgnoreCase)) {
            throw "Destination escaped Phase8 raw root: $resolvedDir"
        }
        $destination = Join-Path $resolvedDir $fileName
        $expected = [int64]$row.processed_file_size
        $status = 'REUSED_SIZE_VERIFIED'
        if (Test-Path -LiteralPath $destination) {
            if ((Get-Item -LiteralPath $destination).Length -ne $expected) {
                throw "Existing final file has unexpected size: $destination"
            }
        }
        else {
            $status = 'BITS_DOWNLOADED'
            $displayName = "CaninePhase8A_$($row.GSM)"
            Write-Phase8Log ("BITS {0}/22 START {1} {2} bytes" -f $index, $row.GSM, $expected)
            $job = Get-BitsTransfer -Name $displayName -ErrorAction SilentlyContinue | Select-Object -First 1
            if ($null -eq $job) {
                $job = Start-BitsTransfer -Source $row.processed_download_url -Destination $destination -DisplayName $displayName -Priority Foreground -Asynchronous
            }
            $deadline = (Get-Date).AddHours(2)
            $nextProgress = Get-Date
            while ($true) {
                $job = Get-BitsTransfer -JobId $job.JobId -ErrorAction Stop
                if ($job.JobState -eq 'Transferred') {
                    $job | Complete-BitsTransfer
                    break
                }
                if ($job.JobState -in @('Error', 'Cancelled', 'Acknowledged')) {
                    throw "BITS permanent state for $($row.GSM): $($job.JobState) $($job.ErrorDescription)"
                }
                if ($job.JobState -in @('Suspended', 'TransientError')) {
                    Write-Phase8Log ("BITS RESUME {0} state={1} bytes={2}/{3} error={4}" -f $row.GSM, $job.JobState, $job.BytesTransferred, $job.BytesTotal, $job.ErrorDescription)
                    $job | Resume-BitsTransfer -Asynchronous -ErrorAction SilentlyContinue
                }
                if ((Get-Date) -ge $nextProgress) {
                    Write-Phase8Log ("BITS PROGRESS {0} state={1} bytes={2}/{3}" -f $row.GSM, $job.JobState, $job.BytesTransferred, $job.BytesTotal)
                    $nextProgress = (Get-Date).AddSeconds(30)
                }
                if ((Get-Date) -gt $deadline) {
                    throw "BITS transfer exceeded two-hour per-file deadline: $($row.GSM)"
                }
                Start-Sleep -Seconds 5
            }
            if (-not (Test-Path -LiteralPath $destination) -or (Get-Item -LiteralPath $destination).Length -ne $expected) {
                throw "BITS destination validation failed: $destination"
            }
            Write-Phase8Log ("BITS {0}/22 COMPLETE {1} {2} bytes" -f $index, $row.GSM, $expected)
        }
        $records.Add([pscustomobject]@{
            GSM = $row.GSM
            donor_id = $row.donor_id
            disease = $row.disease
            state = $row.inflammation_state
            source = $row.processed_download_url
            destination = $destination
            expected_bytes = $expected
            actual_bytes = (Get-Item -LiteralPath $destination).Length
            status = $status
            completed_utc = [DateTime]::UtcNow.ToString('o')
        })
        $records | Export-Csv -LiteralPath $TransferManifest -NoTypeInformation -Encoding utf8
    }
    Write-Phase8Log 'BITS transfer completed for 22/22 AC archives; formal archive/SHA verification remains a separate step'
    exit 0
}
catch {
    Write-Phase8Log ("FATAL: {0}: {1}" -f $_.Exception.GetType().Name, $_.Exception.Message)
    throw
}
