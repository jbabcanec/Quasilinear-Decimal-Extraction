param(
    [long[]]$Sizes=@(100000,1000000),
    [int]$Repeats=3,
    [long[]]$LargeSizes=@(10000000),
    [int]$LargeRepeats=3,
    [string]$OutputPath="$PSScriptRoot/machin_benchmark_2026-09-07.json"
)
$ErrorActionPreference='Stop'
$env:PATH='C:/msys64/ucrt64/bin;'+$env:PATH
$projectRoot=Split-Path $PSScriptRoot -Parent
function Invoke-Run([string]$Executable,[string[]]$Arguments) {
    $info=[Diagnostics.ProcessStartInfo]::new()
    $info.FileName=Join-Path $projectRoot $Executable
    $info.WorkingDirectory=$projectRoot
    $info.UseShellExecute=$false
    $info.CreateNoWindow=$true
    $info.RedirectStandardOutput=$true
    $info.RedirectStandardError=$true
    foreach($argument in $Arguments) {$info.ArgumentList.Add($argument)}
    $process=[Diagnostics.Process]::new()
    $process.StartInfo=$info
    $watch=[Diagnostics.Stopwatch]::StartNew()
    [void]$process.Start()
    $output=$process.StandardOutput.ReadToEndAsync()
    $errors=$process.StandardError.ReadToEndAsync()
    [long]$peak=0
    [double]$cpu=0
    while(-not $process.WaitForExit(25)) {
        try {$process.Refresh();$peak=[Math]::Max($peak,$process.PeakWorkingSet64);$cpu=$process.TotalProcessorTime.TotalSeconds} catch {}
        if($watch.Elapsed.TotalSeconds -gt 600) {$process.Kill($true);$process.WaitForExit();throw 'Own benchmark child exceeded 600 seconds'}
    }
    $watch.Stop()
    try {$peak=[Math]::Max($peak,$process.PeakWorkingSet64);$cpu=$process.TotalProcessorTime.TotalSeconds} catch {}
    $stdout=$output.GetAwaiter().GetResult().Trim()
    $stderr=$errors.GetAwaiter().GetResult().Trim()
    $code=$process.ExitCode
    $process.Dispose()
    if($code -ne 0) {throw "$Executable failed: $stderr"}
    return @{stdout=$stdout;stderr=$stderr;wall_seconds=$watch.Elapsed.TotalSeconds;cpu_seconds=$cpu;peak_working_set_bytes=$peak}
}
function Save-Report {[IO.File]::WriteAllText($OutputPath,($report|ConvertTo-Json -Depth 15)+"`n")}
function Median([double[]]$Values) {
    $sorted=@($Values|Sort-Object)
    if($sorted.Count%2) {return $sorted[[int][Math]::Floor($sorted.Count/2)]}
    return ($sorted[$sorted.Count/2-1]+$sorted[$sorted.Count/2])/2
}
$hashes=@{}
foreach($file in @('decimal_extract_machin.cpp','decimal_extract_shared_residues.cpp','validation/mpfr_window_interval.cpp','validation/benchmark_machin.ps1')) {
    $hashes[$file]=(Get-FileHash (Join-Path $projectRoot $file) -Algorithm SHA256).Hash.ToLowerInvariant()
}
$report=[ordered]@{
    started_utc=[DateTime]::UtcNow.ToString('o')
    cpu=(Get-CimInstance Win32_Processor|Select-Object -First 1 -ExpandProperty Name)
    compiler_flags='-O3 -march=native -std=c++17'
    source_sha256=$hashes
    methodology='Sequential fresh processes, rotating order, no discarded samples. Both new modes use a 128-bit initial guard and the same rigorous interval certificate; the archived Euler control uses its existing fixed 192-bit guard. MPFR references are outside timers. No prefix is computed by either new mode.'
    checks=@{};references=@{};runs=@();summary=@();complete=$false
}
try {
    $report.checks.algebra=(Invoke-Run 'extract_machin.exe' @('--selftest')).stdout|ConvertFrom-Json
    $positions=@(1..100)+@(127,128,129,239,240,241,997,1000,4095,4096,9999,10000)
    foreach($position in $positions) {
        $ref=(Invoke-Run 'validation/mpfr_window_interval.exe' @("$position")).stdout|ConvertFrom-Json
        foreach($mode in @('machin','euler')) {
            $arguments=@("$position",'--timing')
            if($mode -eq 'euler') {$arguments+='--euler'}
            $run=Invoke-Run 'extract_machin.exe' $arguments
            $metrics=$run.stderr|ConvertFrom-Json
            if($run.stdout -cne $ref.digits -or -not $metrics.certified) {throw "Smoke mismatch at $position in $mode"}
        }
    }
    $adaptive=@()
    foreach($position in @(1,17,239,1000,10000)) {
        $ref=(Invoke-Run 'validation/mpfr_window_interval.exe' @("$position")).stdout|ConvertFrom-Json
        foreach($mode in @('machin','euler')) {
            $arguments=@("$position",'--guard=1','--timing')
            if($mode -eq 'euler') {$arguments+='--euler'}
            $run=Invoke-Run 'extract_machin.exe' $arguments
            $metrics=$run.stderr|ConvertFrom-Json
            if($run.stdout -cne $ref.digits -or -not $metrics.certified -or $metrics.passes -le 1) {throw 'Adaptive certificate check failed'}
            $adaptive+=$metrics
        }
    }
    $report.checks.smoke=@{status='PASS';positions=$positions;modes=2}
    $report.checks.adaptive=$adaptive
    Write-Host 'PASS: exact core/residue checks, 224 independent smoke windows, 10 forced guard retries.'
    Save-Report
    $modes=@('machin','euler-certified','euler-archived')
    foreach($position in @($Sizes)+@($LargeSizes)) {
        $report.references["$position"]=(Invoke-Run 'validation/mpfr_window_interval.exe' @("$position")).stdout|ConvertFrom-Json
        $count=if($Sizes -contains $position) {$Repeats} else {$LargeRepeats}
        for($repeat=0;$repeat -lt $count;++$repeat) {
            for($offset=0;$offset -lt $modes.Count;++$offset) {
                $mode=$modes[($repeat+$offset)%$modes.Count]
                $arguments=@("$position",'--timing')
                $executable='extract_machin.exe'
                if($mode -eq 'euler-certified') {$arguments+='--euler'}
                if($mode -eq 'euler-archived') {$executable='extract_shared_residues.exe'}
                $run=Invoke-Run $executable $arguments
                if($run.stdout -cne $report.references["$position"].digits) {throw "Timed mismatch at $position in $mode"}
                $metrics=$run.stderr|ConvertFrom-Json -AsHashtable
                $metrics.mode=$mode
                $metrics.repeat=$repeat
                $metrics.digits=$run.stdout
                foreach($field in @('wall_seconds','cpu_seconds','peak_working_set_bytes')) {$metrics[$field]=$run[$field]}
                $report.runs+=$metrics
                Save-Report
                Write-Host ($metrics|ConvertTo-Json -Compress)
            }
        }
    }
    foreach($position in @($Sizes)+@($LargeSizes)) {foreach($mode in $modes) {
        $runs=@($report.runs|Where-Object { $_.N -eq $position -and $_.mode -eq $mode })
        $report.summary+=@{
            N=$position;mode=$mode;count=$runs.Count
            wall_median=(Median $runs.wall_seconds)
            wall_min=($runs.wall_seconds|Measure-Object -Minimum).Minimum
            wall_max=($runs.wall_seconds|Measure-Object -Maximum).Maximum
            peak_working_set_bytes=($runs.peak_working_set_bytes|Measure-Object -Maximum).Maximum
        }
    }}
    $report.complete=$true
    $report.completed_utc=[DateTime]::UtcNow.ToString('o')
    Save-Report
    Write-Host 'PASS: all timed windows independently verified.'
} catch {$report.failure=$_.Exception.Message;Save-Report;throw}
