<#
.SYNOPSIS
Build, validate, and compare preserved LCM and final reconstruction extractors.
.EXAMPLE
./validation/run_weighted_machin_reconstruction.ps1 -OutputPath validation/weighted_machin_reconstruction_run1.json
.NOTES
Both engines use identical compiler flags, runtime, and initial guard 192.
All correctness checks run before the sequential fresh-process benchmark.
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
$buildDirectory=Join-Path $projectRoot 'tmp/weighted_machin_reconstruction_final'
$sourceNames=@(
    'decimal_extract_weighted_machin.cpp',
    'validation/decimal_extract_weighted_machin_pre_reconstruction_2026-10-24.cpp',
    'validation/mpfr_window_interval.cpp',
    'validation/weighted_machin_arithmetic_checks.cpp',
    'validation/weighted_machin_position_core_checks.cpp',
    'validation/benchmark_weighted_optimization.py',
    'validation/weighted_machin_transfer_normalization_checks_2026-10-24.py',
    'validation/weighted_machin_reconstruction_derivation_2026-10-24.md',
    'validation/run_weighted_machin_reconstruction.ps1'
)
$sourcePaths=@($sourceNames | ForEach-Object { Join-Path $projectRoot $_ })
foreach($source in $sourcePaths) {
    if(-not (Test-Path -LiteralPath $source -PathType Leaf)) { throw "Missing source: $source" }
}
$expectedBaselineHash='7F84ECC227EFF3CC012C0300FD06045591CB4C097986AD0B4BD071D52902C5B2'
if((Get-FileHash -LiteralPath $sourcePaths[1] -Algorithm SHA256).Hash -ne $expectedBaselineHash) {
    throw 'The preserved LCM baseline source hash differs from the recorded snapshot.'
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
        $passProperties=@($record.PSObject.Properties | Where-Object {
            ($_.Value -is [string]) -and ($_.Value -ceq 'PASS')
        })
        if(-not $passProperties.Count) { throw "Unexpected check result: $line" }
        Write-Host "${Name}: $line"
    }
}
$builds=@(
    @{Name='weighted_lcm'; Source=$sourcePaths[1]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='weighted_reconstruction'; Source=$sourcePaths[0]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='mpfr_window_interval'; Source=$sourcePaths[2]; Libraries=@('-lmpfr','-lgmpxx','-lgmp')},
    @{Name='weighted_machin_arithmetic_checks'; Source=$sourcePaths[3]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='weighted_machin_position_core_checks'; Source=$sourcePaths[4]; Libraries=@('-lgmpxx','-lgmp')}
)
$savedPath=$env:PATH
Push-Location -LiteralPath $projectRoot
try {
    $env:PATH=$runtimeDirectory+[IO.Path]::PathSeparator+$savedPath
    New-Item -ItemType Directory -Path $buildDirectory -Force | Out-Null
    $mathLines=@(Invoke-Checked $python @($sourcePaths[6]))
    $mathRecord=($mathLines -join [Environment]::NewLine) | ConvertFrom-Json
    if($mathRecord.status -cne 'PASS') { throw 'Transfer and normalization algebra checks failed.' }
    Write-Host ('Transfer and normalization: '+($mathRecord | ConvertTo-Json -Compress))
    Invoke-Checked $compiler @('--version')
    foreach($build in $builds) {
        $build.Executable=Join-Path $buildDirectory ($build.Name+'.exe')
        Write-Host "Building $($build.Name)"
        Invoke-Checked $compiler (@('-O3','-march=native','-std=c++17',$build.Source)+$build.Libraries+@('-o',$build.Executable))
    }
    foreach($build in $builds[0..1]) {
        Assert-PassRecords $build.Name @(Invoke-Checked $build.Executable @('--selftest'))
        Assert-PassRecords $build.Name @(Invoke-Checked $build.Executable @('--core-check','10000'))
    }
    $arithmeticLines=@(Invoke-Checked $builds[3].Executable @())
    Assert-PassRecords $builds[3].Name $arithmeticLines
    $arithmeticRecord=$arithmeticLines[-1] | ConvertFrom-Json
    if($arithmeticRecord.signed_scale_cases -ne 540) { throw 'Expected 540 exact signed scale checks.' }
    $coreLines=@(Invoke-Checked $builds[4].Executable @())
    Assert-PassRecords $builds[4].Name $coreLines
    $coreRecord=$coreLines[-1] | ConvertFrom-Json
    if($coreRecord.complete_position_cores -ne 14) { throw 'Expected 14 complete position-core checks.' }
    $arguments=@($sourcePaths[5])
    foreach($build in $builds[0..1]) {
        $arguments+=@('--engine',"$($build.Name)=$($build.Executable)",'--engine-arg',"$($build.Name)=--guard=192")
    }
    $arguments+=@('--positions')+@($Positions | ForEach-Object { $_.ToString([Globalization.CultureInfo]::InvariantCulture) })
    $arguments+=@('--repeats',"$Repeats",'--output',$reportPath,'--runtime-dir',$runtimeDirectory,
        '--mpfr-exe',$builds[2].Executable,'--source')+$sourcePaths
    Invoke-Checked $python $arguments
    Write-Host "Verified report: $reportPath"
} finally {
    $env:PATH=$savedPath
    Pop-Location
}
