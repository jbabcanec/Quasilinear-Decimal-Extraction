<#
.SYNOPSIS
Build preserved and optimized extractors, validate their cores, and benchmark them.
.EXAMPLE
./validation/run_weighted_optimization.ps1 -OutputPath validation/weighted_optimization_run1.json
.EXAMPLE
./validation/run_weighted_optimization.ps1 -OutputPath validation/weighted_smoke_run1.json -Positions 100,997 -Repeats 1
.NOTES
Requires Python 3, a GCC compiler with GMP/MPFR, and the corresponding runtime DLLs.
All builds use -O3 -march=native -std=c++17. Outputs go to tmp/weighted_optimization.
The two certified comparison engines start with guard 192, matching weighted.
Their existing certificate retry behavior remains enabled and is recorded in timing JSON.
An existing report is rejected before any compilation; the Python harness also uses
exclusive report creation. Relative OutputPath values resolve from the caller's directory.
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidateNotNullOrEmpty()]
    [string]$OutputPath,

    [ValidateNotNullOrEmpty()]
    [ValidateRange(1, 10000000)]
    [long[]]$Positions = @(100000, 1000000, 10000000),

    [ValidateRange(1, 2147483647)]
    [int]$Repeats = 3,

    [ValidateNotNullOrEmpty()]
    [string]$CompilerPath = 'C:/msys64/ucrt64/bin/g++.exe',

    [ValidateNotNullOrEmpty()]
    [string]$PythonCommand = 'python',

    [ValidateRange(1, 2147483647)]
    [int]$TimeoutSeconds = 600
)

$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$reportPath = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($OutputPath)
if (Test-Path -LiteralPath $reportPath) {
    throw "Refusing to overwrite an existing report or path: $reportPath"
}
if (-not (Test-Path -LiteralPath (Split-Path -Parent $reportPath) -PathType Container)) {
    throw "The output directory must already exist: $(Split-Path -Parent $reportPath)"
}
if (@($Positions | Select-Object -Unique).Count -ne $Positions.Count) {
    throw 'Positions must be unique.'
}

$compiler = (Get-Item -LiteralPath $CompilerPath).FullName
if (-not (Test-Path -LiteralPath $compiler -PathType Leaf)) {
    throw "Compiler is not a file: $compiler"
}
$python = (Get-Command -Name $PythonCommand -CommandType Application -ErrorAction Stop | Select-Object -First 1).Source
$runtimeDirectory = Split-Path -Parent $compiler
$buildDirectory = Join-Path $projectRoot 'tmp/weighted_optimization'
$sourceNames = @(
    'validation/decimal_extract_weighted_baseline_2026-10-24.cpp',
    'decimal_extract_weighted_optimized.cpp',
    'decimal_extract_weighted_words.hpp',
    'validation/weighted_controls_common/decimal_extract_machin.cpp',
    'validation/weighted_controls_common/decimal_extract_shared_residues.cpp',
    'validation/mpfr_window_interval.cpp',
    'validation/benchmark_weighted_optimization.py',
    'validation/run_weighted_optimization.ps1'
)
$sourcePaths = @($sourceNames | ForEach-Object { Join-Path $projectRoot $_ })
foreach ($source in $sourcePaths) {
    if (-not (Test-Path -LiteralPath $source -PathType Leaf)) {
        throw "Missing source dependency: $source"
    }
}

function Invoke-NativeChecked {
    param([string]$Executable, [string[]]$Arguments)
    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Native command failed with exit code ${LASTEXITCODE}: $Executable"
    }
}

$originalExecutable = Join-Path $buildDirectory 'extract_weighted_original.exe'
$optimizedExecutable = Join-Path $buildDirectory 'extract_weighted_optimized.exe'
$eulerExecutable = Join-Path $buildDirectory 'extract_euler_common.exe'
$machinExecutable = Join-Path $buildDirectory 'extract_machin_common.exe'
$verifierExecutable = Join-Path $buildDirectory 'mpfr_window_interval.exe'
$commonMachinSource = Join-Path $projectRoot 'validation/weighted_controls_common/decimal_extract_machin.cpp'
$compilerFlags = @('-O3', '-march=native', '-std=c++17')
$gmpLibraries = @('-lgmpxx', '-lgmp')
$builds = @(
    @{ Source = $sourcePaths[0]; Output = $originalExecutable; Extra = @(); Libraries = $gmpLibraries },
    @{ Source = $sourcePaths[1]; Output = $optimizedExecutable; Extra = @(); Libraries = $gmpLibraries },
    @{ Source = $commonMachinSource; Output = $eulerExecutable; Extra = @('-DCOMMON_EULER_DEFAULT=1'); Libraries = $gmpLibraries },
    @{ Source = $commonMachinSource; Output = $machinExecutable; Extra = @(); Libraries = $gmpLibraries },
    @{ Source = (Join-Path $projectRoot 'validation/mpfr_window_interval.cpp'); Output = $verifierExecutable; Extra = @(); Libraries = @('-lmpfr', '-lgmpxx', '-lgmp') }
)

$savedPath = $env:PATH
Push-Location -LiteralPath $projectRoot
try {
    $env:PATH = $runtimeDirectory + [IO.Path]::PathSeparator + $savedPath
    New-Item -ItemType Directory -Path $buildDirectory -Force | Out-Null
    Invoke-NativeChecked -Executable $compiler -Arguments @('--version')
    foreach ($build in $builds) {
        Write-Host "Building $(Split-Path -Leaf $build.Output)"
        $buildArguments = $compilerFlags + $build.Extra + @($build.Source) + $build.Libraries + @('-o', $build.Output)
        Invoke-NativeChecked -Executable $compiler -Arguments $buildArguments
    }

    foreach ($executable in @($originalExecutable, $optimizedExecutable)) {
        $selftestOutput = (Invoke-NativeChecked -Executable $executable -Arguments @('--selftest')) -join "`n"
        $selftest = $selftestOutput | ConvertFrom-Json
        if ($selftest.selftest -ne 'PASS') {
            throw "Self-test did not report PASS: $executable"
        }
        Write-Host "$(Split-Path -Leaf $executable): $selftestOutput"
    }

    $benchmarkArguments = @(
        (Join-Path $projectRoot 'validation/benchmark_weighted_optimization.py'),
        '--engine', "original=$originalExecutable",
        '--engine', "weighted_optimized=$optimizedExecutable",
        '--engine', "euler_common=$eulerExecutable",
        '--engine', "machin_common=$machinExecutable",
        '--engine-arg', 'euler_common=--guard=192',
        '--engine-arg', 'machin_common=--guard=192',
        '--positions'
    ) + @($Positions | ForEach-Object { $_.ToString([Globalization.CultureInfo]::InvariantCulture) }) + @(
        '--repeats', "$Repeats",
        '--output', $reportPath,
        '--timeout', "$TimeoutSeconds",
        '--reference-timeout', "$TimeoutSeconds",
        '--runtime-dir', $runtimeDirectory,
        '--mpfr-exe', $verifierExecutable,
        '--check-baseline-core', $originalExecutable,
        '--core-engine', 'weighted_optimized',
        '--source'
    ) + $sourcePaths
    Invoke-NativeChecked -Executable $python -Arguments $benchmarkArguments
    Write-Host "Verified benchmark report: $reportPath"
} finally {
    $env:PATH = $savedPath
    Pop-Location
}
