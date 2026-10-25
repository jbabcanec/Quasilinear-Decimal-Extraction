param(
    [string]$ReportPath = 'build-verification.json'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
Push-Location $PSScriptRoot
try {
    $manifest = Get-Content -Raw -LiteralPath 'lake-manifest.json' | ConvertFrom-Json
    $mathlib = $manifest.packages | Where-Object { $_.name -eq 'mathlib' }
    $actualRevision = (& git -C '.lake/packages/mathlib' rev-parse HEAD | Out-String).Trim()
    if ($LASTEXITCODE -ne 0 -or $actualRevision -ne $mathlib.rev) {
        throw 'The installed mathlib revision does not match lake-manifest.json. Run lake update and lake exe cache get.'
    }
    $dependencyChanges = @(& git -C '.lake/packages/mathlib' status --porcelain --untracked-files=no)
    if ($LASTEXITCODE -ne 0 -or $dependencyChanges.Count -gt 0) {
        throw 'The mathlib checkout has tracked modifications; this does not match the pinned dependency snapshot.'
    }
    $dependencyRevisions = [ordered]@{}
    foreach ($package in $manifest.packages) {
        if ($package.type -ne 'git') {
            throw "Unpinned dependency kind: $($package.name)"
        }
        $packagePath = Join-Path '.lake/packages' $package.name
        $revision = (& git -C $packagePath rev-parse HEAD | Out-String).Trim()
        if ($LASTEXITCODE -ne 0 -or $revision -ne $package.rev) {
            throw "Dependency revision mismatch: $($package.name)"
        }
        $changes = @(& git -C $packagePath status --porcelain --untracked-files=no)
        if ($LASTEXITCODE -ne 0 -or $changes.Count -gt 0) {
            throw "Dependency has tracked source changes: $($package.name)"
        }
        $dependencyRevisions[$package.name] = $revision
    }
    $sourceFiles = @('lean-toolchain', 'lakefile.toml', 'lake-manifest.json', 'BBP.lean', 'Audit.lean', 'verify.ps1')
    $sourceFiles += @(Get-ChildItem -LiteralPath 'BBP' -File -Filter '*.lean' |
        Sort-Object Name | ForEach-Object { 'BBP/' + $_.Name })
    $rootImports = @(Get-Content -LiteralPath 'BBP.lean' | Where-Object { $_ -match '^import BBP\.' })
    foreach ($moduleFile in Get-ChildItem -LiteralPath 'BBP' -File -Filter '*.lean') {
        $requiredImport = 'import BBP.' + $moduleFile.BaseName
        if ($rootImports -notcontains $requiredImport) {
            throw "Module lacks an explicit root import and would escape the audit: $($moduleFile.Name)"
        }
    }
    $hashesBefore = [ordered]@{}
    foreach ($file in $sourceFiles) {
        $hashesBefore[$file] = (Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash.ToLowerInvariant()
    }

    $buildOutput = @(& lake build 2>&1 | ForEach-Object { $_.ToString() })
    $buildExit = $LASTEXITCODE
    $buildOutput | Set-Content -LiteralPath '.lake/verification-build.log' -Encoding utf8
    if ($buildExit -ne 0) {
        $buildOutput | Write-Output
        throw 'Lean project build failed.'
    }
    $auditOutput = @(& lake env lean Audit.lean 2>&1 | ForEach-Object { $_.ToString() })
    $auditExit = $LASTEXITCODE
    $auditOutput | Set-Content -LiteralPath '.lake/verification-audit.log' -Encoding utf8
    $auditOutput | Write-Output
    if ($auditExit -ne 0 -or -not ($auditOutput -match 'AXIOM AUDIT PASS:')) {
        throw 'The axiom audit failed or did not produce its success marker.'
    }
    foreach ($file in $sourceFiles) {
        $hashAfter = (Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash.ToLowerInvariant()
        if ($hashAfter -ne $hashesBefore[$file]) {
            throw "Source changed during verification: $file. Rerun when edits are complete."
        }
    }
    $leanVersion = (& lake env lean --version | Out-String).Trim()
    if ($LASTEXITCODE -ne 0) { throw 'Could not obtain Lean version.' }
    $report = [ordered]@{
        verified_utc = [DateTime]::UtcNow.ToString('o')
        lean_version = $leanVersion
        mathlib_revision = $actualRevision
        dependency_revisions = $dependencyRevisions
        build_exit_code = $buildExit
        axiom_audit_exit_code = $auditExit
        axiom_audit = $auditOutput
        allowed_axioms = @('propext', 'Classical.choice', 'Quot.sound')
        scope = 'Formal components only. The complete extractor and its machine bit-complexity theorem remain unproved; see README.md.'
        source_sha256 = $hashesBefore
    }
    $report | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $ReportPath -Encoding utf8
    Write-Output "Verified formal components; report: $ReportPath"
}
finally {
    Pop-Location
}
