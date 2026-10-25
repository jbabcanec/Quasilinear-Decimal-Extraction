param(
    [long[]]$Sizes = @(100000,1000000,10000000),
    [int]$Repeats = 3,
    [string]$OutputPath = "$PSScriptRoot/shared_residue_benchmark_2026-09-07.json"
)
$ErrorActionPreference = 'Stop'
$env:PATH = 'C:/msys64/ucrt64/bin;' + $env:PATH
$projectRoot = Split-Path $PSScriptRoot -Parent

function Invoke-NativeRun([string]$Executable, [string[]]$RunArguments) {
    $startInfo = [System.Diagnostics.ProcessStartInfo]::new()
    $startInfo.FileName = Join-Path $projectRoot $Executable
    $startInfo.WorkingDirectory = $projectRoot
    $startInfo.UseShellExecute = $false
    $startInfo.CreateNoWindow = $true
    $startInfo.RedirectStandardOutput = $true
    $startInfo.RedirectStandardError = $true
    foreach ($argument in $RunArguments) { $startInfo.ArgumentList.Add($argument) }
    $process = [System.Diagnostics.Process]::new()
    $process.StartInfo = $startInfo
    $timer = [System.Diagnostics.Stopwatch]::StartNew()
    [void]$process.Start()
    $stdoutTask = $process.StandardOutput.ReadToEndAsync()
    $stderrTask = $process.StandardError.ReadToEndAsync()
    $process.WaitForExit()
    $timer.Stop()
    $stdout = $stdoutTask.GetAwaiter().GetResult().Trim()
    $stderr = $stderrTask.GetAwaiter().GetResult().Trim()
    $code = $process.ExitCode
    $process.Dispose()
    if ($code -ne 0) { throw "$Executable exited ${code}: $stderr" }
    return @{ stdout=$stdout; stderr=$stderr; wall_seconds=$timer.Elapsed.TotalSeconds }
}

function Save-Report { [IO.File]::WriteAllText($OutputPath, ($report | ConvertTo-Json -Depth 12)+"`n") }

$sourceFiles = @('decimal_extract_v5.cpp','decimal_extract_shared_residues.cpp',
                 'validation/mpfr_window_interval.cpp','validation/compare_shared_residues.ps1')
$hashes = @{}
foreach ($file in $sourceFiles) { $hashes[$file]=(Get-FileHash (Join-Path $projectRoot $file) -Algorithm SHA256).Hash.ToLowerInvariant() }
$cpu = Get-CimInstance Win32_Processor | Select-Object -First 1 -ExpandProperty Name
$report = [ordered]@{
    date_utc=[DateTime]::UtcNow.ToString('o'); cpu=$cpu;
    method='Sequential fresh processes; rotating order; three repetitions. Wall time includes process launch. MPFR verification runs are outside extraction timing.';
    compiler_flags='-O3 -march=native -std=c++17';
    source_sha256=$hashes; checks=@{}; references=@{}; runs=@(); verified=$false
}
$check = Invoke-NativeRun 'extract_shared_residues.exe' @('--check-shared')
$report.checks.shared = $check.stdout | ConvertFrom-Json
$check = Invoke-NativeRun 'extract_shared_residues.exe' @('--check-core')
$report.checks.core = $check.stdout | ConvertFrom-Json
$smokePositions = @(1..64)+@(97,997,10000)
foreach ($digitPosition in $smokePositions) {
    $reference = (Invoke-NativeRun 'validation/mpfr_window_interval.exe' @("$digitPosition")).stdout | ConvertFrom-Json
    $digits = (Invoke-NativeRun 'extract_shared_residues.exe' @("$digitPosition")).stdout
    if ($digits -cne $reference.digits) { throw "Smoke mismatch at $digitPosition" }
}
$report.checks.digits = @{status='PASS';positions=$smokePositions;method='MPFR directed-rounding enclosure'}
Save-Report
Write-Host ($report.checks | ConvertTo-Json -Depth 4 -Compress)

foreach ($digitPosition in $Sizes) {
    $reference = (Invoke-NativeRun 'validation/mpfr_window_interval.exe' @("$digitPosition")).stdout | ConvertFrom-Json
    $report.references["$digitPosition"] = $reference
    $modes = if ($digitPosition -eq 1000000) { @('baseline','shared-only','joint-only','shared-joint') } else { @('baseline','shared-joint') }
    for ($rep=0; $rep -lt $Repeats; ++$rep) {
        for ($index=0; $index -lt $modes.Count; ++$index) {
            $mode=$modes[($index+$rep)%$modes.Count]
            $executable=if ($mode -eq 'baseline') { 'extract_euler_baseline.exe' } else { 'extract_shared_residues.exe' }
            $runArguments=@("$digitPosition",'--timing')
            if ($mode -in @('shared-only','joint-only')) { $runArguments += "--$mode" }
            $result=Invoke-NativeRun $executable $runArguments
            if ($result.stdout -cne $reference.digits) { throw "Benchmark mismatch at $digitPosition ($mode)" }
            $stats=$result.stderr | ConvertFrom-Json -AsHashtable
            $stats.mode=$mode; $stats.repeat=$rep; $stats.digits=$result.stdout
            $stats.wall_seconds=$result.wall_seconds
            $report.runs += $stats
            Save-Report
            Write-Host ($stats | ConvertTo-Json -Compress)
        }
    }
}
$report.verified=$true
Save-Report
Write-Host "PASS: all outputs agree with independent directed-rounding intervals."
