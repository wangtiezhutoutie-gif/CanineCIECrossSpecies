param(
    [ValidateSet('Start','Status','Finalize')]
    [string]$Mode='Status',
    [string]$JobName='Phase14_GSE277276_PROCESSED_101'
)
$ErrorActionPreference='Stop'
$root=$env:CANINE_CIE_PROJECT_ROOT
if ([string]::IsNullOrWhiteSpace($root)) { throw 'Set CANINE_CIE_PROJECT_ROOT' }
$phase=Join-Path $root '14_ced_specificity'
$expectPath=Join-Path $phase 'metadata\GSE277276_PROCESSED_FILE_EXPECTATIONS.csv'
$stage=Join-Path $phase 'raw_processed\GSE277276_bits_staging'
$final=Join-Path $phase 'raw_processed\GSE277276_cellbender'
$log=Join-Path $root 'logs\phase14_download.log'
New-Item -ItemType Directory -Path $stage,$final,(Split-Path $log -Parent) -Force|Out-Null
$expect=Import-Csv -LiteralPath $expectPath
if($expect.Count -ne 102){throw "Expected 102 GEO files; found $($expect.Count)"}
if($Mode -eq 'Start'){
    $remaining=@($expect|Where-Object{-not(Test-Path -LiteralPath(Join-Path $stage $_.file_name))})
    if($remaining.Count -eq 0){'All expected files already staged';exit 0}
    $src=@();$dst=@()
    foreach($r in $remaining){
        $src+=('https://www.ncbi.nlm.nih.gov/geo/download/?acc='+$r.gsm+'&format=file&file='+[Uri]::EscapeDataString($r.file_name))
        $dst+=(Join-Path $stage $r.file_name)
    }
    $job=Start-BitsTransfer -Source $src -Destination $dst -Asynchronous -DisplayName $JobName -Description 'Official GEO GSE277276 processed CellBender files'
    Set-BitsTransfer -BitsJob $job -Priority Foreground
    Resume-BitsTransfer -BitsJob $job -Asynchronous|Out-Null
}
if($Mode -eq 'Status' -or $Mode -eq 'Start'){
    $job=Get-BitsTransfer -Name $JobName -ErrorAction SilentlyContinue
    if($null -eq $job){[pscustomobject]@{State='NO_ACTIVE_JOB';StagedFiles=(Get-ChildItem $stage -File).Count}}
    else{[pscustomobject]@{
        JobId=$job.JobId;State=$job.JobState;BytesTransferred=$job.BytesTransferred
        BytesTotal=$job.BytesTotal
        Percent=if($job.BytesTotal -gt 0 -and $job.BytesTotal -lt [uint64]::MaxValue){[math]::Round(100*$job.BytesTransferred/$job.BytesTotal,2)}else{0}
        FilesTransferred=$job.FilesTransferred;QueuedFiles=$job.FileList.Count
        Error=$job.ErrorDescription}}
    exit 0
}
$job=Get-BitsTransfer -Name $JobName -ErrorAction SilentlyContinue
if($null -ne $job){
    if($job.JobState -ne 'Transferred'){throw "BITS job is $($job.JobState), not Transferred"}
    Complete-BitsTransfer -BitsJob $job
}
$staged=Get-ChildItem -LiteralPath $stage -File
if($staged.Count -ne 102){throw "Expected 102 staged files; found $($staged.Count)"}
$expectedNames=@($expect.file_name|Sort-Object)
$observedNames=@($staged.Name|Sort-Object)
if(Compare-Object $expectedNames $observedNames){throw 'Staged names differ from GEO expectation'}
$audit=@()
foreach($r in $expect){
    $p=Join-Path $stage $r.file_name;$item=Get-Item -LiteralPath $p
    if($item.Length -le 0){throw "Zero-byte staged file: $p"}
    $valid=$false;$uncompressed=[int64]0;$stream=[IO.File]::OpenRead($p)
    try{
        $gz=[IO.Compression.GZipStream]::new($stream,[IO.Compression.CompressionMode]::Decompress)
        try{$buffer=New-Object byte[](8MB);while(($n=$gz.Read($buffer,0,$buffer.Length))-gt 0){$uncompressed+=$n};$valid=$true}
        finally{$gz.Dispose()}
    }finally{$stream.Dispose()}
    if(-not $valid){throw "gzip validation failed: $p"}
    $audit+=[pscustomobject]@{gsm=$r.gsm;filename=$r.file_name;compressed_bytes=$item.Length;uncompressed_bytes=$uncompressed;gzip_status='PASS'}
}
$audit|Export-Csv -LiteralPath(Join-Path $phase 'audit\GSE277276_BITS_GZIP_AUDIT.csv') -NoTypeInformation
foreach($r in $expect){Move-Item -LiteralPath(Join-Path $stage $r.file_name) -Destination(Join-Path $final $r.file_name) -Force}
$finalFiles=Get-ChildItem -LiteralPath $final -File|Where-Object{$_.Name -in $expectedNames}
if($finalFiles.Count -ne 102){throw "Final directory lacks expected 102 files"}
$summary=[pscustomobject]@{
    completed_utc=[DateTime]::UtcNow.ToString('o');source='NCBI GEO web endpoint via BITS'
    expected_files=102;validated_gzip_files=102
    compressed_bytes=($finalFiles|Measure-Object Length -Sum).Sum
    status='PASS_ATOMIC_STAGE_THEN_PUBLISH'}
$summary|ConvertTo-Json|Set-Content -LiteralPath $log -Encoding utf8
$summary

