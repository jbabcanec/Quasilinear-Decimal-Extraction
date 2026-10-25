<#
.SYNOPSIS
Build, check, and compare weighted Machin against both Machin controls; retain the faster measured comparison.
.EXAMPLE
./validation/run_weighted_machin_best.ps1 -OutputPath validation/weighted_machin_run2.json
.NOTES
All engines use the same compiler, flags, runtime, and 192-bit initial guard.
Weighted uses the original Newton formula; both original and half-width controls are measured.
Reports are exclusively created; existing results are never overwritten.
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory=$true)][string]$OutputPath,
    [ValidateRange(1,10000000)][long[]]$Positions=@(100000,1000000,10000000),
    [ValidateRange(1,2147483647)][int]$Repeats=5,
    [string]$CompilerPath='C:/msys64/ucrt64/bin/g++.exe',
    [string]$PythonCommand='python'
)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$reportPath=$ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($OutputPath)
if(Test-Path -LiteralPath $reportPath) { throw "Refusing to overwrite: $reportPath" }
if(-not (Test-Path -LiteralPath (Split-Path -Parent $reportPath) -PathType Container)) {
    throw 'The report directory must already exist.'
}
if(@($Positions | Select-Object -Unique).Count -ne $Positions.Count) { throw 'Positions must be unique.' }
$compiler=(Get-Item -LiteralPath $CompilerPath).FullName
$python=(Get-Command -Name $PythonCommand -CommandType Application | Select-Object -First 1).Source
$runtimeDirectory=Split-Path -Parent $compiler
$buildDirectory=Join-Path $projectRoot 'tmp/weighted_machin_best'
$sourceNames=@(
    'decimal_extract_weighted_machin.cpp',
    'validation/decimal_extract_machin_half_newton_2026-10-24.cpp',
    'validation/decimal_extract_weighted_refined_2026-10-24.cpp',
    'decimal_extract_weighted_words.hpp',
    'validation/weighted_controls_common/decimal_extract_machin.cpp',
    'validation/weighted_controls_common/decimal_extract_shared_residues.cpp',
    'validation/mpfr_window_interval.cpp',
    'validation/benchmark_weighted_optimization.py',
    'validation/run_weighted_machin_best.ps1'
)
$sourcePaths=@($sourceNames | ForEach-Object { Join-Path $projectRoot $_ })
foreach($source in $sourcePaths) {
    if(-not (Test-Path -LiteralPath $source -PathType Leaf)) { throw "Missing source: $source" }
}
function Invoke-Checked {
    param([string]$Executable,[string[]]$Arguments)
    & $Executable @Arguments
    if($LASTEXITCODE -ne 0) { throw "Exit ${LASTEXITCODE}: $Executable" }
}
$builds=@(
    @{Name='machin_previous'; Source=$sourcePaths[4]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='machin_half_newton'; Source=$sourcePaths[1]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='weighted_machin'; Source=$sourcePaths[0]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='mpfr_window_interval'; Source=$sourcePaths[6]; Libraries=@('-lmpfr','-lgmpxx','-lgmp')}
)
$savedPath=$env:PATH
Push-Location -LiteralPath $projectRoot
try {
    $env:PATH=$runtimeDirectory+[IO.Path]::PathSeparator+$savedPath
    New-Item -ItemType Directory -Path $buildDirectory -Force | Out-Null
    Invoke-Checked $compiler @('--version')
    foreach($build in $builds) {
        $build.Executable=Join-Path $buildDirectory ($build.Name+'.exe')
        Write-Host "Building $($build.Name)"
        Invoke-Checked $compiler (@('-O3','-march=native','-std=c++17',$build.Source)+$build.Libraries+@('-o',$build.Executable))
    }
    foreach($build in $builds[0..2]) {
        $lines=@(Invoke-Checked $build.Executable @('--selftest'))
        if(-not $lines.Count) { throw "No selftest output: $($build.Name)" }
        foreach($line in $lines) {
            $record=$line | ConvertFrom-Json
            if($record.selftest -ne 'PASS' -and $record.weighted_selftest -ne 'PASS') {
                throw "Unexpected selftest: $line"
            }
            Write-Host "$($build.Name): $line"
        }
    }
    Invoke-Checked $builds[2].Executable @('--core-check','10000')
    $arguments=@((Join-Path $projectRoot 'validation/benchmark_weighted_optimization.py'))
    foreach($build in $builds[0..2]) {
        $arguments+=@('--engine',"$($build.Name)=$($build.Executable)",'--engine-arg',"$($build.Name)=--guard=192")
    }
    $arguments+=@('--positions')+@($Positions | ForEach-Object { $_.ToString([Globalization.CultureInfo]::InvariantCulture) })
    $arguments+=@('--repeats',"$Repeats",'--output',$reportPath,'--runtime-dir',$runtimeDirectory,
        '--mpfr-exe',$builds[3].Executable,'--source')+$sourcePaths
    Invoke-Checked $python $arguments
    Write-Host "Verified report: $reportPath"
} finally {
    $env:PATH=$savedPath
    Pop-Location
}
