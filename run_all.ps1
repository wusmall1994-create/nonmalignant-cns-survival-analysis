$ErrorActionPreference = 'Stop'

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
if (-not $env:SEER_EXPORT_DIR) {
    throw 'Set SEER_EXPORT_DIR to the directory containing the SEER*Stat exports.'
}
if (-not $env:ANALYSIS_OUTPUT_DIR) {
    $env:ANALYSIS_OUTPUT_DIR = Join-Path $scriptDir 'outputs'
}

$pythonExe = if ($env:PYTHON_EXE) {
    $env:PYTHON_EXE
} else {
    Get-Command python.exe -ErrorAction Stop | Select-Object -ExpandProperty Source -First 1
}
$rCandidates = @(
    $env:RSCRIPT_EXE,
    (Get-Command Rscript.exe -ErrorAction SilentlyContinue | Select-Object -ExpandProperty Source -First 1)
) | Where-Object { $_ -and (Test-Path -LiteralPath $_) }

if (-not $rCandidates) {
    throw 'Rscript.exe was not found.'
}

& $pythonExe (Join-Path $scriptDir 'analyze.py')
if ($LASTEXITCODE -ne 0) { throw 'Python analysis failed.' }

$rscriptExe = $rCandidates | Select-Object -First 1
& $rscriptExe (Join-Path $scriptDir 'individual_analysis.R')
if ($LASTEXITCODE -ne 0) { throw 'R individual-level analysis failed.' }

& $rscriptExe (Join-Path $scriptDir 'audit_sessions.R')
if ($LASTEXITCODE -ne 0) { throw 'SEER session/export audit failed.' }

& $rscriptExe (Join-Path $scriptDir 'pohar_perme_analysis.R')
if ($LASTEXITCODE -ne 0) { throw 'Pohar-Perme sensitivity analysis failed.' }

& $rscriptExe (Join-Path $scriptDir 'advanced_robustness_analysis.R')
if ($LASTEXITCODE -ne 0) { throw 'Advanced robustness analysis failed.' }

& $rscriptExe (Join-Path $scriptDir 'make_figures.R')
if ($LASTEXITCODE -ne 0) { throw 'R figure generation failed.' }

Write-Host 'Analysis, sensitivity analyses, and figures completed successfully.'
