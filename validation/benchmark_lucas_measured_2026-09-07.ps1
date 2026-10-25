param(
    [long[]]$Sizes = @(100000,1000000),
    [int]$Repeats = 3,
    [long[]]$LargeSizes = @(10000000),
    [int]$RunLimitSeconds = 900,
    [string]$OutputPath = "$PSScriptRoot/lucas_benchmark_2026-09-08.json"
)
$ErrorActionPreference='Stop'
$env:PATH='C:/msys64/ucrt64/bin;'+$env:PATH
$projectRoot=Split-Path $PSScriptRoot -Parent
function Invoke-Run([string]$Executable,[string[]]$Arguments) {
    $start=[Diagnostics.ProcessStartInfo]::new()
    $start.FileName=Join-Path $projectRoot $Executable
    $start.WorkingDirectory=$projectRoot
    $start.UseShellExecute=$false; $start.CreateNoWindow=$true
    $start.RedirectStandardOutput=$true; $start.RedirectStandardError=$true
    foreach($argument in $Arguments) { $start.ArgumentList.Add($argument) }
    $process=[Diagnostics.Process]::new(); $process.StartInfo=$start
    $timer=[Diagnostics.Stopwatch]::StartNew()
    [void]$process.Start()
    $stdoutTask=$process.StandardOutput.ReadToEndAsync()
    $stderrTask=$process.StandardError.ReadToEndAsync()
    [long]$peak=0; [double]$cpu=0
    while(-not $process.WaitForExit(100)) {
        try { $process.Refresh(); $peak=[Math]::Max($peak,$process.PeakWorkingSet64); $cpu=$process.TotalProcessorTime.TotalSeconds } catch {}
        if($timer.Elapsed.TotalSeconds -gt $RunLimitSeconds) {
            $process.Kill($true); $process.WaitForExit()
            throw "$Executable exceeded the ${RunLimitSeconds}-second per-run limit."
        }
    }
    $timer.Stop()
    try { $peak=[Math]::Max($peak,$process.PeakWorkingSet64); $cpu=$process.TotalProcessorTime.TotalSeconds } catch {}
    $stdout=$stdoutTask.GetAwaiter().GetResult().Trim()
    $stderr=$stderrTask.GetAwaiter().GetResult().Trim()
    $exitCode=$process.ExitCode; $process.Dispose()
    if($exitCode -ne 0) { throw "$Executable failed: $stderr" }
    return @{stdout=$stdout; stderr=$stderr; wall_seconds=$timer.Elapsed.TotalSeconds; cpu_seconds=$cpu; peak_working_set_bytes=$peak}
}
function Save-Report {
    [IO.File]::WriteAllText($OutputPath,($report | ConvertTo-Json -Depth 15)+"`n")
}
function Get-Median([double[]]$Values) {
    $ordered=@($Values | Sort-Object)
    if($ordered.Count%2) { return $ordered[[int][Math]::Floor($ordered.Count/2)] }
    return ($ordered[$ordered.Count/2-1]+$ordered[$ordered.Count/2])/2
}
$hashes=@{}
foreach($source in @('decimal_extract_lucas.cpp','decimal_extract_weighted.cpp','decimal_extract_shared_residues.cpp','validation/mpfr_window_interval.cpp','validation/benchmark_lucas.ps1')) {
    $hashes[$source]=(Get-FileHash (Join-Path $projectRoot $source) -Algorithm SHA256).Hash.ToLowerInvariant()
}
$report=[ordered]@{
    started_utc=[DateTime]::UtcNow.ToString('o')
    cpu=(Get-CimInstance Win32_Processor | Select-Object -First 1 -ExpandProperty Name)
    compiler_flags='-O3 -march=native -std=c++17'
    source_sha256=$hashes
    methodology='Sequential fresh processes with rotating order. MPFR directed-rounding references run outside extraction timing. New recursive and direct-control modes both perform exact word aggregation and the same internal decimal certificate. The existing Euler executable uses its fixed 192-bit guard. Large sizes have one run per mode.'
    checks=@{}; references=@{}; runs=@(); summary=@(); complete=$false
}
try {
    $check=Invoke-Run 'extract_lucas.exe' @('--selftest')
    $report.checks.algebra=$check.stdout | ConvertFrom-Json
    Write-Host $check.stdout
    $smoke=@(1..32)+@(63,64,97,997,1000,10000)
    foreach($position in $smoke) {
        $reference=(Invoke-Run 'validation/mpfr_window_interval.exe' @("$position")).stdout | ConvertFrom-Json
        $result=Invoke-Run 'extract_lucas.exe' @("$position",'--timing')
        if($result.stdout -cne $reference.digits) { throw "Digit mismatch at $position" }
        $metrics=$result.stderr | ConvertFrom-Json
        if(-not $metrics.certified) { throw 'Missing internal certificate' }
    }
    $adaptive=@()
    foreach($position in @(1,17,97,1000)) {
        $reference=(Invoke-Run 'validation/mpfr_window_interval.exe' @("$position")).stdout | ConvertFrom-Json
        $result=Invoke-Run 'extract_lucas.exe' @("$position",'--guard=1','--timing')
        if($result.stdout -cne $reference.digits) { throw "Adaptive certificate mismatch at $position" }
        $adaptive+=($result.stderr | ConvertFrom-Json)
    }
    if(-not ($adaptive | Where-Object certificate_passes -gt 1)) { throw 'Adaptive test did not exercise retries' }
    $report.checks.digits=@{status='PASS'; positions=$smoke; method='MPFR directed-rounding intervals'}
    $report.checks.adaptive_certificate=$adaptive
    Save-Report
    Write-Host 'PASS: smoke digits and forced certificate retries.'
    $targets=@()
    foreach($position in $Sizes) { $targets+=@{N=$position; repetitions=$Repeats; category='repeated'} }
    foreach($position in $LargeSizes) { $targets+=@{N=$position; repetitions=1; category='single-large-run'} }
    $modes=@('recursive','direct-control','existing-euler')
    foreach($target in $targets) {
        $position=$target.N
        $reference=(Invoke-Run 'validation/mpfr_window_interval.exe' @("$position")).stdout | ConvertFrom-Json
        $report.references["$position"]=$reference
        for($rep=0;$rep -lt $target.repetitions;++$rep) {
            for($index=0;$index -lt $modes.Count;++$index) {
                $mode=$modes[($index+$rep)%$modes.Count]
                $executable=if($mode -eq 'existing-euler') {'extract_shared_residues.exe'} else {'extract_lucas.exe'}
                $runArguments=@("$position",'--timing')
                if($mode -eq 'direct-control') { $runArguments+='--direct' }
                $loadBefore=Get-CimInstance Win32_Processor | Select-Object -First 1 -ExpandProperty LoadPercentage
                $result=Invoke-Run $executable $runArguments
                if($result.stdout -cne $reference.digits) { throw "Timed digit mismatch at $position ($mode)" }
                $metrics=$result.stderr | ConvertFrom-Json -AsHashtable
                if($mode -eq 'recursive' -and $position -ge 100000) {
                    foreach($counter in @('levels','prime_queries','power_queries','nonunit_divisions')) {
                        if($metrics[$counter] -le 0) { throw "New algorithm path $counter was not exercised at $position" }
                    }
                }
                $metrics.mode=$mode; $metrics.repeat=$rep; $metrics.category=$target.category
                $metrics.digits=$result.stdout; $metrics.wall_seconds=$result.wall_seconds
                $metrics.cpu_seconds=$result.cpu_seconds; $metrics.peak_working_set_bytes=$result.peak_working_set_bytes
                $metrics.host_cpu_percent_before=$loadBefore
                $report.runs+=$metrics; Save-Report
                Write-Host ($metrics | ConvertTo-Json -Compress)
            }
        }
    }
    foreach($group in ($report.runs | Group-Object N,mode)) {
        $rows=@($group.Group)
        $report.summary+=@{
            N=$rows[0].N; mode=$rows[0].mode; repetitions=$rows.Count
            wall_median=(Get-Median @($rows | ForEach-Object wall_seconds))
            wall_min=($rows | ForEach-Object wall_seconds | Measure-Object -Minimum).Minimum
            wall_max=($rows | ForEach-Object wall_seconds | Measure-Object -Maximum).Maximum
            core_median=(Get-Median @($rows | ForEach-Object core))
            word_median=(Get-Median @($rows | ForEach-Object word))
            peak_bytes=($rows | ForEach-Object peak_working_set_bytes | Measure-Object -Maximum).Maximum
        }
    }
    $report.complete=$true; $report.finished_utc=[DateTime]::UtcNow.ToString('o'); Save-Report
    Write-Host 'PASS: all benchmark outputs match independent intervals.'
} catch {
    $report.error=$_.ToString(); Save-Report; throw
}
