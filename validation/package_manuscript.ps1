param(
    [string]$PdfLaTex = 'pdflatex',
    [string]$PdfInfo = 'pdfinfo',
    [string]$PdfToText = 'pdftotext',
    [string]$PdfFonts = 'pdffonts'
)
# Build from a fresh extraction of an explicit source-only archive.
# Benchmark programs and historical datasets are read, never rerun or rewritten.
$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$projectRoot = Split-Path $PSScriptRoot -Parent
$main = 'decimal_extraction_paper.tex'
$pdfName = 'decimal_extraction_paper.pdf'
$culture = [Globalization.CultureInfo]::InvariantCulture
function Require([bool]$Condition, [string]$Message) {
    if (-not $Condition) { throw $Message }
}
function HashFile([string]$Path) {
    (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
}
function ReadJson([string]$Relative) {
    Get-Content -Raw -LiteralPath (Join-Path $projectRoot $Relative) | ConvertFrom-Json -AsHashtable
}
function WriteUtf8([string]$Path, [string]$Content) {
    [IO.File]::WriteAllText($Path, $Content + [char]10, [Text.UTF8Encoding]::new($false))
}
function F([double]$Value, [int]$Places) {
    $Value.ToString("F$Places", $culture)
}
function ResolvePdfTool([string]$Name) {
    $command = Get-Command $Name -ErrorAction SilentlyContinue
    if ($command) { return $command.Source }
    $bundled = Join-Path $env:USERPROFILE ".cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/$Name.exe"
    Require (Test-Path -LiteralPath $bundled) "Tool unavailable: $Name; pass its path as a parameter."
    return $bundled
}
$PdfLaTex = ResolvePdfTool $PdfLaTex
$PdfInfo = ResolvePdfTool $PdfInfo
$PdfToText = ResolvePdfTool $PdfToText
$PdfFonts = ResolvePdfTool $PdfFonts
$source = [IO.File]::ReadAllText((Join-Path $projectRoot $main))
$labels = @([regex]::Matches($source, '\\label\{([^}]+)\}') | ForEach-Object {$_.Groups[1].Value})
$refs = @([regex]::Matches($source, '\\(?:eqref|ref)\{([^}]+)\}') | ForEach-Object {$_.Groups[1].Value})
$bib = @([regex]::Matches($source, '\\bibitem\{([^}]+)\}') | ForEach-Object {$_.Groups[1].Value})
$cites = @([regex]::Matches($source, '\\cite(?:\[[^\]]*\])?\{([^}]+)\}') | ForEach-Object {$_.Groups[1].Value.Split(',')} | Sort-Object -Unique)
Require ($labels.Count -eq @($labels | Sort-Object -Unique).Count) 'Duplicate labels.'
Require (@($refs | Where-Object {$_ -notin $labels}).Count -eq 0) 'Missing references.'
Require (@(Compare-Object ($bib | Sort-Object) $cites).Count -eq 0) 'Bibliography/citation mismatch.'
Require ($source -notmatch 'Python|CPython|gmpy2|\bv[345]\b|\\path\{|README|anc/|\.cpp|\.py|\.json|\bTODO\b|\bFIXME\b|\bTBD\b') 'Development notation or placeholder in manuscript.'
Require ($source -notmatch '\\(?:input|include|write18)\s*\{') 'Unexpected external TeX dependency.'
Require ([regex]::Matches($source, '\\begin\{table\}').Count -eq 4) 'Expected three performance tables and one complexity table.'

$backup = Join-Path $projectRoot 'notes/manuscript_backups/before_consolidation_2026-09-07/decimal_extraction_paper.tex'
$proofPreserved = $null
if (Test-Path -LiteralPath $backup) {
    $old = [IO.File]::ReadAllText($backup).Replace([string][char]13, '')
    $current = $source.Replace([string][char]13, '')
    $proofPattern = '(?s)\\section\{Computational model and notation\}.*?(?=\\section\{Comparison with earlier methods\})'
    $theoremPattern = '(?s)\\begin\{theorem\}\[Certified decimal extraction\].*?\\end\{theorem\}'
    Require ([regex]::Match($old, $proofPattern).Success) 'Backup proof section missing.'
    Require ([regex]::Match($old, $proofPattern).Value -ceq [regex]::Match($current, $proofPattern).Value) 'Main proof changed during editorial cleanup.'
    Require ([regex]::Match($old, $theoremPattern).Value -ceq [regex]::Match($current, $theoremPattern).Value) 'Main theorem changed during editorial cleanup.'
    $proofPreserved = $true
}

$latest = ReadJson 'validation/machin_benchmark_2026-09-07.json'
$external = ReadJson 'validation/external_comparison_2026-09-07.json'
$weighted = ReadJson 'validation/weighted_cpp_performance_2026-09-07.json'
$large = ReadJson 'validation/cpp_euler_final_benchmark_2026-09-06.json'
Require ($latest.complete -and $latest.runs.Count -eq 27) 'Incomplete latest benchmark.'
Require ($external.verified -and $weighted.verified -and $large.verified) 'Earlier report not verified.'
$checkedHashes = @{}
$omittedExternal = @{}
foreach ($report in @($latest, $external, $weighted, $large)) {
    foreach ($entry in $report.source_sha256.GetEnumerator()) {
        $path = Join-Path $projectRoot $entry.Key
        if (-not (Test-Path -LiteralPath $path) -and $entry.Key -match '/bellard_pi1\.c$') {
            $omittedExternal[$entry.Key] = $entry.Value
            continue
        }
        Require ((HashFile $path) -eq $entry.Value) "Measured source changed: $($entry.Key)"
        $checkedHashes[$entry.Key] = $entry.Value
    }
}
$latestNames = @{machin='Machin, certified'; 'euler-certified'='Euler, certified'; 'euler-archived'='Euler, fixed guard'}
foreach ($summary in $latest.summary) {
    $samples = @($latest.runs | Where-Object {$_.N -eq $summary.N -and $_.mode -eq $summary.mode})
    $times = @($samples.wall_seconds | Sort-Object)
    Require ($samples.Count -eq 3 -and $summary.count -eq 3) 'Latest sample count mismatch.'
    Require ($times[1] -eq $summary.wall_median -and $times[0] -eq $summary.wall_min -and $times[2] -eq $summary.wall_max) 'Latest summary mismatch.'
    Require (($samples.peak_working_set_bytes | Measure-Object -Maximum).Maximum -eq $summary.peak_working_set_bytes) 'Latest memory summary mismatch.'
    foreach ($sample in $samples) {
        Require ($sample.digits -ceq $latest.references[[string]$sample.N].digits) 'Latest stored digit mismatch.'
        if ($sample.mode -ne 'euler-archived') { Require $sample.certified 'Missing internal certificate.' }
    }
    $row = '& {0} & {1} & {2}--{3} & {4}' -f $latestNames[$summary.mode],(F $summary.wall_median 3),(F $summary.wall_min 3),(F $summary.wall_max 3),(F ($summary.peak_working_set_bytes/1MB) 1)
    Require ($source.Contains($row)) "Missing or incorrect latest table row: $row"
}
$externalNames = @{bellard='Bellard';gourdon='Gourdon (reimplementation)';v5='2-adic Euler summation';mpfr_prefix='MPFR'}
foreach ($summary in $external.summary) {
    $samples = @($external.runs | Where-Object {$_.N -eq $summary.N -and $_.engine -eq $summary.engine})
    Require ($samples.Count -eq 3 -and @($samples.wall_seconds | Sort-Object)[1] -eq $summary.wall_median) 'External median mismatch.'
    foreach ($sample in $samples) {
        Require ($sample.digits -ceq $external.references[[string]$sample.N].Substring(0,$sample.computed_width)) 'External stored digit mismatch.'
    }
    if (-not $externalNames.ContainsKey($summary.engine)) {continue}
    $width = if ($summary.engine -eq 'v5') {12} else {9}
    $row = '& {0} & {1} & {2} & {3}--{4} & {5}' -f $externalNames[$summary.engine],$width,(F $summary.wall_median 3),(F $summary.wall_min 3),(F $summary.wall_max 3),(F ($summary.peak_bytes_max/1MB) 1)
    Require ($source.Contains($row)) "External table mismatch: $row"
}
foreach ($summary in $weighted.summary) {
    $samples = @($weighted.runs | Where-Object {$_.N -eq $summary.N})
    Require ($samples.Count -eq 3 -and @($samples.wall_seconds | Sort-Object)[1] -eq $summary.wall_median) 'Weighted median mismatch.'
    foreach ($sample in $samples) {Require ($sample.digits -ceq $weighted.references[[string]$sample.N]) 'Weighted stored digit mismatch.'}
    if ($summary.N -lt 1000000) {
        $row = '& Weighted reconstruction & 12 & {0} & {1}--{2} & {3}' -f (F $summary.wall_median 3),(F $summary.wall_min 3),(F $summary.wall_max 3),(F ($summary.peak_bytes_max/1MB) 1)
    } else {
        $row = '& Weighted reconstruction & 3 & {0} ({1}--{2}) & {3}' -f (F $summary.wall_median 2),(F $summary.wall_min 2),(F $summary.wall_max 2),(F ($summary.peak_bytes_max/1MB) 1)
    }
    Require ($source.Contains($row)) "Weighted table mismatch: $row"
}
foreach ($sample in $large.runs) {
    $name = @{v3='Direct binary splitting';v5='2-adic Euler summation'}[$sample.engine]
    $row = '& {0} & 1 & {1} & {2}' -f $name,(F $sample.wall_seconds 2),(F ($sample.peak_working_set_bytes/1MB) 1)
    Require ($source.Contains($row)) "Larger-position table mismatch: $row"
    Require ($sample.digits -ceq $large.references[[string]$sample.N].digits) 'Large stored digit mismatch.'
}

$mapping = ReadJson 'validation/manuscript_sources.json'
Require ($mapping.Count -eq @($mapping.archive | Sort-Object -Unique).Count) 'Duplicate archive paths.'
$stage = Join-Path $projectRoot ('tmp/pdfs/consolidated_' + [guid]::NewGuid().ToString('N'))
[IO.Directory]::CreateDirectory($stage) | Out-Null
$stagedZip = Join-Path $stage 'decimal_extraction_source.zip'
Add-Type -AssemblyName System.IO.Compression.FileSystem
$zip = [IO.Compression.ZipFile]::Open($stagedZip, [IO.Compression.ZipArchiveMode]::Create)
$manifest = @()
try {
    foreach ($item in ($mapping | Sort-Object archive)) {
        Require ($item.archive -notmatch '(^/|\\|(^|/)\.\.(/|$)|\.(pdf|exe|aux|log|out)$)') "Invalid archive path: $($item.archive)"
        Require ($item.archive -eq $main -or $item.archive.StartsWith('anc/')) 'Unexpected top-level file.'
        $file = Join-Path $projectRoot $item.source
        [IO.Compression.ZipFileExtensions]::CreateEntryFromFile($zip,$file,$item.archive,[IO.Compression.CompressionLevel]::Optimal) | Out-Null
        $manifest += [ordered]@{path=$item.archive;bytes=(Get-Item -LiteralPath $file).Length;sha256=(HashFile $file)}
    }
} finally {$zip.Dispose()}
$buildDir = Join-Path $stage 'isolated'
[IO.Compression.ZipFile]::ExtractToDirectory($stagedZip,$buildDir)
foreach ($entry in $manifest) {
    Require ((HashFile (Join-Path $buildDir $entry.path)) -eq $entry.sha256) "Archive byte mismatch: $($entry.path)"
}
Push-Location $buildDir
try {
    for ($pass=1; $pass -le 3; $pass++) {
        $buildText = & $PdfLaTex -no-shell-escape -interaction=nonstopmode -halt-on-error $main 2>&1
        WriteUtf8 (Join-Path $stage "build-pass-$pass.txt") ($buildText -join [char]10)
        Require ($LASTEXITCODE -eq 0) "PDFLaTeX failed; see $stage/build-pass-$pass.txt"
    }
} finally {Pop-Location}
$log = [IO.File]::ReadAllText((Join-Path $buildDir 'decimal_extraction_paper.log'))
$issues = @([regex]::Matches($log,'(?im)^(?!Package:).*?(?:warning|Overfull|Underfull|Missing character).*$') | ForEach-Object {$_.Value})
Require ($issues.Count -eq 0) ("TeX build issues: " + ($issues -join '; '))
$pdf = Join-Path $buildDir $pdfName
$info = (& $PdfInfo $pdf) -join [char]10
Require ($LASTEXITCODE -eq 0) 'pdfinfo failed.'
$pages = [int][regex]::Match($info,'(?m)^Pages:\s+(\d+)').Groups[1].Value
Require ($pages -gt 0) 'No PDF pages.'
$textPath = Join-Path $stage 'manuscript.txt'
& $PdfToText -layout $pdf $textPath
Require ($LASTEXITCODE -eq 0) 'pdftotext failed.'
$pageTexts = [IO.File]::ReadAllText($textPath).Split([char]12)
for ($index=0; $index -lt $pages; $index++) {
    Require ($pageTexts[$index].Trim().Length -gt 0 -and $pageTexts[$index] -notmatch '\?\?') "Empty page or unresolved reference: $($index+1)"
    Require ($pageTexts[$index].TrimEnd() -match ('(?:^|\s)' + ($index+1) + '$')) "Incorrect page number: $($index+1)"
}
$fonts = @(& $PdfFonts $pdf)
Require ($LASTEXITCODE -eq 0) 'pdffonts failed.'
$fontRows = @($fonts | Select-Object -Skip 2 | Where-Object {$_.Trim()})
Require ($fontRows.Count -gt 0 -and ($fonts -join ' ') -notmatch 'Type 3') 'Unexpected font inventory.'
foreach ($font in $fontRows) { Require (($font -split '\s+')[-5] -eq 'yes') 'Unembedded PDF font.' }
WriteUtf8 (Join-Path $stage 'fonts.txt') ($fonts -join [char]10)

$out = Join-Path $projectRoot 'output/arxiv'
$pdfOut = Join-Path $projectRoot 'output/pdf'
[IO.Directory]::CreateDirectory($out) | Out-Null
[IO.Directory]::CreateDirectory($pdfOut) | Out-Null
Copy-Item -LiteralPath $stagedZip -Destination (Join-Path $out 'decimal_extraction_source.zip') -Force
Copy-Item -LiteralPath $pdf -Destination (Join-Path $projectRoot $pdfName) -Force
Copy-Item -LiteralPath $pdf -Destination (Join-Path $pdfOut $pdfName) -Force
$archiveHash = HashFile $stagedZip
WriteUtf8 (Join-Path $out 'source_manifest.json') ([ordered]@{archive='decimal_extraction_source.zip';archive_sha256=$archiveHash;files=$manifest} | ConvertTo-Json -Depth 8)
$metadataPath = Join-Path $out 'metadata.txt'
$metadata = [IO.File]::ReadAllText($metadataPath)
$abstract = [regex]::Match($source,'(?s)\\begin\{abstract\}(.*?)\\end\{abstract\}').Groups[1].Value.Trim() -replace '\s+',' '
$metadata = $metadata -replace '(?m)^Comments:.*$',"Comments: $pages pages; 3 performance tables and 1 complexity table; supplementary code and benchmark data"
$metadata = ($metadata -split 'Abstract:',2)[0] + 'Abstract:' + [char]10 + $abstract
WriteUtf8 $metadataPath $metadata.TrimEnd()
$verification = [ordered]@{
    verified_at_utc=[DateTime]::UtcNow.ToString('o')
    pages=$pages;performance_tables=3;complexity_tables=1
    source_sha256=(HashFile (Join-Path $projectRoot $main));pdf_sha256=(HashFile $pdf);archive_sha256=$archiveHash
    source_archive_files=$manifest.Count
    isolated_build=@{passes=3;shell_escape=$false;only_archive_contents_required=$true;warnings=$issues;directory=$buildDir}
    source_checks=@{unique_labels=$labels.Count;all_cross_references_resolve=$true;cited_bibliography_entries=$bib.Count;all_performance_rows_match_reports=$true;main_theorem_and_proof_preserved=$proofPreserved;development_paths_absent=$true}
    pdf_checks=@{all_pages_have_text=$true;page_numbers_verified=$true;all_fonts_embedded=$true;type3_fonts=$false;font_count=$fontRows.Count}
    benchmark_checks=@{new_runs=27;historical_datasets_reused_without_new_timing=$true;source_hashes=$checkedHashes;external_source_not_locally_available=$omittedExternal}
    visual_review=@{completed=$false;note='Render and inspect every page, then record the reviewed PDF hash.'}
    scope='Editorial and reproduction verification. No new independent proof audit, benchmark run, external speed record, or arXiv submission.'
}
WriteUtf8 (Join-Path $out 'verification.json') ($verification | ConvertTo-Json -Depth 10)
Write-Host "PASS: $pages pages; $($manifest.Count) source files; all performance rows, source hashes, references, fonts, and isolated build checked."
Write-Host "Visual review pending. Build and extracted text: $stage"
