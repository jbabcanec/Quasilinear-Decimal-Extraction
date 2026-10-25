<#
.SYNOPSIS
Build, check, and compare serial weighted, parallel weighted, improved Machin,
and the experimental direct-H Lucas reconstruction.
.EXAMPLE
pwsh -File validation/run_weighted_machin_parallel.ps1 -OutputPath validation/weighted_machin_parallel_local.json
.NOTES
The default run contains 24 sequential fresh-process samples: four engines,
two positions, and three repetitions. Engine order rotates. All builds use
the same compiler flags and initial guard 192; independent MPFR enclosure
checks run outside the timers. Reports and build manifests are never overwritten.
Parallel total/parallel_wall are wall times; its phase fields sum worker
elapsed times and overlap. Parallelism does not change the log^(5/2) work bound.
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory=$true)][string]$OutputPath,
    [ValidateRange(1,10000000)][long[]]$Positions=@(1000000,10000000),
    [ValidateRange(1,2147483647)][int]$Repeats=3,
    [ValidateRange(1,64)][int]$Threads=8,
    [string]$CompilerPath='C:/msys64/ucrt64/bin/g++.exe',
    [string]$PythonCommand='python'
)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$reportPath=$ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($OutputPath)
$manifestPath=$reportPath+'.build.json'
foreach($path in @($reportPath,$manifestPath)) {
    if(Test-Path -LiteralPath $path) { throw "Refusing to overwrite: $path" }
}
if(-not (Test-Path -LiteralPath (Split-Path -Parent $reportPath) -PathType Container)) {
    throw 'The report directory must already exist.'
}
if(-not $Positions.Count -or @($Positions | Select-Object -Unique).Count -ne $Positions.Count) {
    throw 'Provide a nonempty list of unique positions.'
}
$compiler=(Get-Item -LiteralPath $CompilerPath).FullName
$python=(Get-Command -Name $PythonCommand -CommandType Application | Select-Object -First 1).Source
$runtimeDirectory=Split-Path -Parent $compiler
$buildDirectory=Join-Path $projectRoot ('tmp/weighted_machin_parallel/'+[Guid]::NewGuid().ToString('N'))
$compileFlags=@('-O3','-march=native','-std=c++17','-pthread')
$sourceNames=@(
    'validation/decimal_extract_weighted_machin_pre_reconstruction_2026-10-24.cpp',
    'decimal_extract_weighted_machin_parallel.cpp',
    'validation/decimal_extract_machin_improved_2026-10-24.cpp',
    'decimal_extract_weighted_machin_lucas.cpp',
    'validation/machin_lucas_common.hpp',
    'validation/machin_lucas_coefficients.hpp',
    'validation/machin_lucas_reconstruction.hpp',
    'validation/machin_lucas_adic.hpp',
    'validation/weighted_controls_common/decimal_extract_shared_residues.cpp',
    'decimal_extract_weighted_words.hpp',
    'validation/weighted_machin_parallel_checks.cpp',
    'validation/machin_lucas_checks.cpp',
    'validation/machin_lucas_adic_checks.cpp',
    'validation/mpfr_window_interval.cpp',
    'validation/benchmark_weighted_optimization.py',
    'validation/run_weighted_machin_parallel.ps1'
)
$sourcePaths=@($sourceNames | ForEach-Object { Join-Path $projectRoot $_ })
$sourceHashes=[ordered]@{}
foreach($source in $sourcePaths) {
    if(-not (Test-Path -LiteralPath $source -PathType Leaf)) { throw "Missing source: $source" }
    $sourceHashes[$source]=(Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash
}
$expectedBaselineHash='7F84ECC227EFF3CC012C0300FD06045591CB4C097986AD0B4BD071D52902C5B2'
if($sourceHashes[$sourcePaths[0]] -ne $expectedBaselineHash) {
    throw 'The preserved serial LCM baseline differs from its recorded source snapshot.'
}
function Invoke-Checked {
    param([string]$Executable,[string[]]$Arguments)
    & $Executable @Arguments
    if($LASTEXITCODE -ne 0) { throw "Exit ${LASTEXITCODE}: $Executable" }
}
$testRecords=[System.Collections.Generic.List[object]]::new()
function Add-PassRecords {
    param([string]$Name,[object[]]$Lines)
    if(-not $Lines.Count) { throw "No check output: $Name" }
    foreach($line in $Lines) {
        $record=$line | ConvertFrom-Json
        $passed=@($record.PSObject.Properties | Where-Object {
            ($_.Value -is [string]) -and ($_.Value -ceq 'PASS')
        })
        if(-not $passed.Count) { throw "Unexpected check result: $line" }
        $testRecords.Add([ordered]@{test=$Name; result=$record})
        Write-Host "${Name}: $line"
    }
}
$builds=@(
    @{Name='weighted_serial'; Source=$sourcePaths[0]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='weighted_parallel'; Source=$sourcePaths[1]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='machin_improved'; Source=$sourcePaths[2]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='lucas_direct_H'; Source=$sourcePaths[3]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='parallel_checks'; Source=$sourcePaths[10]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='lucas_checks'; Source=$sourcePaths[11]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='adic_checks'; Source=$sourcePaths[12]; Libraries=@('-lgmpxx','-lgmp')},
    @{Name='mpfr_window_interval'; Source=$sourcePaths[13]; Libraries=@('-lmpfr','-lgmpxx','-lgmp')}
)
$savedPath=$env:PATH
Push-Location -LiteralPath $projectRoot
try {
    $env:PATH=$runtimeDirectory+[IO.Path]::PathSeparator+$savedPath
    New-Item -ItemType Directory -Path $buildDirectory -Force | Out-Null
    $compilerVersion=@(Invoke-Checked $compiler @('--version'))
    $compilerTarget=@(Invoke-Checked $compiler @('-dumpmachine'))
    foreach($build in $builds) {
        $build.Executable=Join-Path $buildDirectory ($build.Name+'.exe')
        $build.Arguments=$compileFlags+@($build.Source)+$build.Libraries+@('-o',$build.Executable)
        Write-Host "Building $($build.Name)"
        Invoke-Checked $compiler $build.Arguments
    }
    foreach($build in $builds[0..3]) {
        Add-PassRecords ($build.Name+' selftest') @(Invoke-Checked $build.Executable @('--selftest'))
    }
    foreach($index in @(0,1,3)) {
        $build=$builds[$index]
        Add-PassRecords ($build.Name+' core-check') @(Invoke-Checked $build.Executable @('--core-check','10000'))
    }
    foreach($build in $builds[4..6]) {
        Add-PassRecords $build.Name @(Invoke-Checked $build.Executable @())
    }
    foreach($source in $sourcePaths) {
        if((Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash -ne $sourceHashes[$source]) {
            throw "Source changed during compilation or validation: $source"
        }
    }
    $runtimeHashes=[ordered]@{}
    foreach($name in @('libgmp-10.dll','libgmpxx-4.dll','libmpfr-6.dll','libwinpthread-1.dll',
                       'libgcc_s_seh-1.dll','libstdc++-6.dll')) {
        $path=Join-Path $runtimeDirectory $name
        if(Test-Path -LiteralPath $path -PathType Leaf) {
            $runtimeHashes[$path]=(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash
        }
    }
    $manifest=[ordered]@{
        schema_version=1
        created_utc=[DateTime]::UtcNow.ToString('o')
        compiler=$compiler
        compiler_sha256=(Get-FileHash -LiteralPath $compiler -Algorithm SHA256).Hash
        compiler_version=$compilerVersion
        compiler_target=$compilerTarget
        common_compile_flags=$compileFlags
        initial_guard=192
        parallel_threads=$Threads
        managed_processor_count=[Environment]::ProcessorCount
        processor_identifier=$env:PROCESSOR_IDENTIFIER
        powershell_version=$PSVersionTable.PSVersion.ToString()
        operating_system=[Environment]::OSVersion.ToString()
        runtime_directory=$runtimeDirectory
        runtime_sha256=$runtimeHashes
        source_sha256=$sourceHashes
        builds=@($builds | ForEach-Object {
            [ordered]@{name=$_.Name; source=$_.Source; command=@($compiler)+$_.Arguments;
                executable=$_.Executable; sha256=(Get-FileHash -LiteralPath $_.Executable -Algorithm SHA256).Hash}
        })
        exact_checks=@($testRecords.ToArray())
        phase_timing_note='Parallel phase fields sum overlapping worker elapsed times. Compare total or process wall time, not phase sums.'
        complexity_note='Fixed-thread parallelism preserves the serial O(N log^(5/2) N) one-pass work and O(N) stored-data bounds; no lower exponent is claimed.'
    }
    $json=$manifest | ConvertTo-Json -Depth 20
    $stream=[IO.File]::Open($manifestPath,[IO.FileMode]::CreateNew,[IO.FileAccess]::Write,[IO.FileShare]::None)
    try {
        $writer=[IO.StreamWriter]::new($stream,[Text.UTF8Encoding]::new($false))
        try { $writer.WriteLine($json) } finally { $writer.Dispose() }
    } finally { $stream.Dispose() }
    $arguments=@($sourcePaths[14])
    foreach($build in $builds[0..3]) {
        $arguments+=@('--engine',"$($build.Name)=$($build.Executable)",
                      '--engine-arg',"$($build.Name)=--guard=192")
    }
    $arguments+=@('--engine-arg',"weighted_parallel=--threads=$Threads")
    $arguments+=@('--positions')+@($Positions | ForEach-Object { $_.ToString([Globalization.CultureInfo]::InvariantCulture) })
    $arguments+=@('--repeats',"$Repeats",'--output',$reportPath,'--runtime-dir',$runtimeDirectory,
        '--mpfr-exe',$builds[7].Executable,'--source')+$sourcePaths+@($manifestPath)
    Invoke-Checked $python $arguments
    Write-Host "Verified report: $reportPath"
    Write-Host "Build, environment, and exact-check manifest: $manifestPath"
} finally {
    $env:PATH=$savedPath
    Pop-Location
}
