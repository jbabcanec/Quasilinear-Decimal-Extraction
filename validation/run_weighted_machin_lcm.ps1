<#
.SYNOPSIS
Build and verify four matched Machin engines before comparing their timings.
.EXAMPLE
./validation/run_weighted_machin_lcm.ps1 -OutputPath validation/weighted_machin_lcm_run1.json
.NOTES
All engines use the same compiler, flags, runtime, and 192-bit initial guard.
The improved log-cubed control includes the shared inverse and scaling changes.
Reports are exclusively created; existing results are never overwritten.
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory=$true)][string]$OutputPath,
    [ValidateRange(1,10000000)][long[]]$Positions=@(100000,1000000,10000000),
    [ValidateRange(1,2147483647)][int]$Repeats=3,
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
$buildDirectory=Join-Path $projectRoot 'tmp/weighted_machin_lcm_final'
$sourceNames=@(
    'decimal_extract_weighted_machin.cpp',
    'validation/decimal_extract_machin_improved_2026-10-24.cpp',
    'validation/decimal_extract_weighted_machin_pre_lcm_2026-10-24.cpp',
    'decimal_extract_weighted_words.hpp',
    'validation/weighted_controls_common/decimal_extract_machin.cpp',
    'validation/weighted_controls_common/decimal_extract_shared_residues.cpp',
    'validation/mpfr_window_interval.cpp',
    'validation/weighted_machin_arithmetic_checks.cpp',
    'validation/benchmark_weighted_optimization.py',
    'validation/run_weighted_machin_lcm.ps1'
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
function Assert-PassRecords {
    param([string]$Name,[object[]]$Lines)
    if(-not $Lines.Count) { throw "No check output: $Name" }
    foreach($line in $Lines) {
        $record=$line | ConvertFrom-Json
        $passProperties=@($record.PSObject.Properties | Where-Object { $_.Value -ceq 'PASS' })
        if(-not $passProperties.Count) { throw "Unexpected check result: $line" }
        Write-Host "${Name}: $line"
    }
}
$builds=@(
    @{Name='machin_previous'; Source=$sourcePaths[4]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='machin_improved'; Source=$sourcePaths[1]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='weighted_previous'; Source=$sourcePaths[2]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='weighted_lcm'; Source=$sourcePaths[0]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='mpfr_window_interval'; Source=$sourcePaths[6]; Libraries=@('-lmpfr','-lgmpxx','-lgmp')},
    @{Name='weighted_machin_arithmetic_checks'; Source=$sourcePaths[7]; Libraries=@('-lgmpxx','-lgmp')}
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
    foreach($build in $builds[0..3]) {
        $lines=@(Invoke-Checked $build.Executable @('--selftest'))
        Assert-PassRecords $build.Name $lines
    }
    foreach($build in $builds[2..3]) {
        $lines=@(Invoke-Checked $build.Executable @('--core-check','10000'))
        Assert-PassRecords $build.Name $lines
    }
    $arithmeticLines=@(Invoke-Checked $builds[5].Executable @())
    Assert-PassRecords $builds[5].Name $arithmeticLines
    $arithmeticRecord=$arithmeticLines[-1] | ConvertFrom-Json
    if($arithmeticRecord.signed_scale_cases -ne 540) { throw 'Expected 540 exact signed scale checks.' }
    $arguments=@((Join-Path $projectRoot 'validation/benchmark_weighted_optimization.py'))
    foreach($build in $builds[0..3]) {
        $arguments+=@('--engine',"$($build.Name)=$($build.Executable)",'--engine-arg',"$($build.Name)=--guard=192")
    }
    $arguments+=@('--positions')+@($Positions | ForEach-Object { $_.ToString([Globalization.CultureInfo]::InvariantCulture) })
    $arguments+=@('--repeats',"$Repeats",'--output',$reportPath,'--runtime-dir',$runtimeDirectory,
        '--mpfr-exe',$builds[4].Executable,'--source')+$sourcePaths
    Invoke-Checked $python $arguments
    Write-Host "Verified report: $reportPath"
} finally {
    $env:PATH=$savedPath
    Pop-Location
}
