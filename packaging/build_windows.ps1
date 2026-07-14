param(
    [string]$Python = ""
)

$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$ReleaseDir = Join-Path $ProjectRoot "release"
$WorkDir = Join-Path $ProjectRoot "build\pyinstaller"
$SpecPath = Join-Path $PSScriptRoot "genesnap_workbench.spec"
$IconScript = Join-Path $PSScriptRoot "generate_windows_icon.py"

if (-not $Python) {
    $LocalPython = Join-Path $ProjectRoot ".venv\Scripts\python.exe"
    if (Test-Path -LiteralPath $LocalPython) {
        $Python = (Resolve-Path $LocalPython).Path
    }
    else {
        $Python = (Get-Command python -ErrorAction Stop).Source
    }
}

Push-Location $ProjectRoot
try {
    $env:PYTHONPATH = "src"
    $env:QT_QPA_PLATFORM = "offscreen"
    & $Python $IconScript
    if ($LASTEXITCODE -ne 0) {
        throw "Windows icon generation failed."
    }
    & $Python -m unittest discover -s tests -v
    if ($LASTEXITCODE -ne 0) {
        throw "Full test suite failed. Packaging stopped."
    }

    New-Item -ItemType Directory -Force -Path $ReleaseDir | Out-Null
    New-Item -ItemType Directory -Force -Path $WorkDir | Out-Null
    & $Python -m PyInstaller `
        --noconfirm `
        --clean `
        --distpath $ReleaseDir `
        --workpath $WorkDir `
        $SpecPath
    if ($LASTEXITCODE -ne 0) {
        throw "PyInstaller build failed."
    }

    $ExePath = Join-Path $ReleaseDir "GeneSnapWorkbench.exe"
    $SmokeDir = Join-Path `
        $ProjectRoot `
        ("build\packaged-smoke\" + [Guid]::NewGuid().ToString("N"))
    $SmokeDataDir = Join-Path $SmokeDir "data"
    $SmokeReportPath = Join-Path $SmokeDir "report.json"
    New-Item -ItemType Directory -Force -Path $SmokeDir | Out-Null
    $SmokeArguments = @(
        "--smoke-test",
        "--data-dir",
        ('"{0}"' -f $SmokeDataDir),
        "--smoke-report",
        ('"{0}"' -f $SmokeReportPath)
    )
    $SmokeProcess = Start-Process `
        -FilePath $ExePath `
        -ArgumentList $SmokeArguments `
        -WindowStyle Hidden `
        -Wait `
        -PassThru
    if ($SmokeProcess.ExitCode -ne 0) {
        throw "Packaged smoke test failed with exit code $($SmokeProcess.ExitCode)."
    }
    if (-not (Test-Path -LiteralPath $SmokeReportPath)) {
        throw "Packaged smoke test did not write its report."
    }
    $SmokePayload = Get-Content `
        -LiteralPath $SmokeReportPath `
        -Raw `
        -Encoding UTF8 | ConvertFrom-Json
    if (-not $SmokePayload.ok) {
        throw "Packaged smoke report indicates failure."
    }

    $ReleaseDocuments = @(
        "LICENSE",
        "NOTICE",
        "THIRD_PARTY_NOTICES.md"
    )
    foreach ($DocumentName in $ReleaseDocuments) {
        Copy-Item `
            -LiteralPath (Join-Path $ProjectRoot $DocumentName) `
            -Destination (Join-Path $ReleaseDir $DocumentName) `
            -Force
    }

    $Hash = Get-FileHash -LiteralPath $ExePath -Algorithm SHA256
    $HashLine = "$($Hash.Hash.ToLowerInvariant())  GeneSnapWorkbench.exe"
    Set-Content `
        -LiteralPath (Join-Path $ReleaseDir "GeneSnapWorkbench.sha256") `
        -Value $HashLine `
        -Encoding ASCII
    Write-Host "Windows artifact: $ExePath"
    Write-Host "Packaged smoke: OK ($SmokeReportPath)"
    Write-Host "SHA-256: $($Hash.Hash)"
}
finally {
    Pop-Location
}
