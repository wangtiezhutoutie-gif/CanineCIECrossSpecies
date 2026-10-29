$ErrorActionPreference = 'Stop'

$Root = $env:CANINE_CIE_PROJECT_ROOT
if ([string]::IsNullOrWhiteSpace($Root)) { throw 'Set CANINE_CIE_PROJECT_ROOT' }
$SampleMaster = Join-Path $Root '00_metadata\GSE266616_HUMAN_SAMPLE_MASTER.csv'
$RawRoot = Join-Path $Root '07_human_crohn\GSE266616\raw_processed'
$TransferManifest = Join-Path $Root '00_metadata\GSE266616_BITS_TRANSFER_MANIFEST.csv'
$LogPath = Join-Path $Root 'logs\phase7_BITS_transfer.log'
$Python = if ($env:CANINE_CIE_PYTHON) { $env:CANINE_CIE_PYTHON } else { 'python' }
$Verifier = Join-Path $Root 'scripts\32_phase7_download_TI.py'
$Forbidden = @('GSE266616_RAW.tar', '.fastq', '.sra', 'scRNA_gene_counts.mtx')

function Write-Phase7Log {
    param([Parameter(Mandatory = $true)][string]$Message)
    $line = '[{0}] {1}' -f ([DateTime]::UtcNow.ToString('o')), $Message
    Write-Output $line
    Add-Content -LiteralPath $LogPath -Value $line -Encoding utf8
}

function Get-CohortDirectory {
    param([Parameter(Mandatory = $true)]$Row)
    if ($Row.region -ne 'TI' -or $Row.tissue -ne 'terminal_ileum') {
        throw "Non-TI sample reached BITS transfer: $($Row.GSM)"
    }
    if ($Row.healthy_or_crohn -eq 'Healthy' -and $Row.inflammation_state -eq 'Healthy') {
        return 'healthy_TI'
    }
    if ($Row.healthy_or_crohn -ne 'Crohn') {
        throw "Unexpected TI disease class: $($Row.GSM) $($Row.healthy_or_crohn)"
    }
    switch ($Row.inflammation_state) {
        'Active_inflamed' { return 'crohn_TI_active' }
        'Adjacent' { return 'crohn_TI_adjacent' }
        'Non_inflamed' { return 'crohn_TI_noninflamed' }
        default { throw "Unexpected Crohn TI inflammation state: $($Row.GSM) $($Row.inflammation_state)" }
    }
}

New-Item -ItemType Directory -Force -Path (Split-Path $LogPath -Parent) | Out-Null
Set-Content -LiteralPath $LogPath -Value '' -Encoding utf8

try {
    $allRows = @(Import-Csv -LiteralPath $SampleMaster)
    $selected = @($allRows | Where-Object { $_.region -eq 'TI' })
    if ($selected.Count -ne 48) {
        throw "Expected 48 frozen TI samples, observed $($selected.Count)"
    }
    $totalBytes = ($selected | Measure-Object -Property processed_file_size -Sum).Sum
    if ([int64]$totalBytes -ne 2302117661) {
        throw "Frozen TI byte total changed: $totalBytes"
    }

    Write-Phase7Log 'Starting BITS transfer for 48 per-GSM TI processed archives; AC/RAW/SRA/FASTQ/Zenodo matrix excluded'
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
        $cohortDir = Get-CohortDirectory -Row $row
        $destinationDir = Join-Path $RawRoot $cohortDir
        New-Item -ItemType Directory -Force -Path $destinationDir | Out-Null
        $resolvedDir = (Resolve-Path -LiteralPath $destinationDir).Path
        if (-not $resolvedDir.StartsWith(($RawRoot + '\'), [StringComparison]::OrdinalIgnoreCase)) {
            throw "Destination escaped raw_processed root: $resolvedDir"
        }
        $destination = Join-Path $resolvedDir $fileName
        $expected = [int64]$row.processed_file_size
        $status = 'REUSED_SIZE_VERIFIED'
        if (Test-Path -LiteralPath $destination) {
            $existingSize = (Get-Item -LiteralPath $destination).Length
            if ($existingSize -ne $expected) {
                throw "Existing final file has unexpected size and was not overwritten: $destination ($existingSize != $expected)"
            }
        }
        else {
            $status = 'BITS_DOWNLOADED'
            $transferUrl = $row.processed_download_url
            $displayName = "CaninePhase7_$($row.GSM)"
            Write-Phase7Log ("BITS {0}/48 START {1} {2} bytes" -f $index, $row.GSM, $expected)
            $job = Get-BitsTransfer -Name $displayName -ErrorAction SilentlyContinue | Select-Object -First 1
            if ($null -eq $job) {
                $job = Start-BitsTransfer -Source $transferUrl -Destination $destination -DisplayName $displayName -Priority Foreground -Asynchronous
            }
            $deadline = (Get-Date).AddHours(1)
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
                    Write-Phase7Log ("BITS RESUME {0} state={1} bytes={2}/{3} error={4}" -f $row.GSM, $job.JobState, $job.BytesTransferred, $job.BytesTotal, $job.ErrorDescription)
                    $job | Resume-BitsTransfer -Asynchronous -ErrorAction SilentlyContinue
                }
                if ((Get-Date) -ge $nextProgress) {
                    Write-Phase7Log ("BITS PROGRESS {0} state={1} bytes={2}/{3}" -f $row.GSM, $job.JobState, $job.BytesTransferred, $job.BytesTotal)
                    $nextProgress = (Get-Date).AddSeconds(30)
                }
                if ((Get-Date) -gt $deadline) {
                    throw "BITS transfer exceeded one-hour per-file deadline: $($row.GSM) state=$($job.JobState) bytes=$($job.BytesTransferred)/$($job.BytesTotal)"
                }
                Start-Sleep -Seconds 5
            }
            if (-not (Test-Path -LiteralPath $destination)) {
                throw "BITS reported completion without destination: $destination"
            }
            $actual = (Get-Item -LiteralPath $destination).Length
            if ($actual -ne $expected) {
                throw "BITS byte mismatch for $($row.GSM): $actual != $expected"
            }
            Write-Phase7Log ("BITS {0}/48 COMPLETE {1} {2} bytes" -f $index, $row.GSM, $actual)
        }

        $records.Add([pscustomobject]@{
            GSM = $row.GSM
            sample_title = $row.sample_title
            patient_id = $row.patient_id
            inflammation_state = $row.inflammation_state
            source_metadata_url = $row.processed_download_url
            transfer_url = $row.processed_download_url
            destination = $destination
            expected_bytes = $expected
            actual_bytes = (Get-Item -LiteralPath $destination).Length
            transfer_status = $status
            completed_utc = [DateTime]::UtcNow.ToString('o')
        })
        $records | Export-Csv -LiteralPath $TransferManifest -NoTypeInformation -Encoding utf8
    }

    Write-Phase7Log 'BITS transfer complete for 48/48 TI archives; starting tar/SHA-256 verification'
    & $Python $Verifier
    if ($LASTEXITCODE -ne 0) {
        throw "Formal TI archive verifier failed with exit code $LASTEXITCODE"
    }
    Write-Phase7Log 'Phase7B BITS transfer and formal archive verification complete'
    exit 0
}
catch {
    Write-Phase7Log ("FATAL: {0}: {1}" -f $_.Exception.GetType().Name, $_.Exception.Message)
    throw
}
