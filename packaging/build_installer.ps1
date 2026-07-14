param(
    [string]$ISCC = ""
)

$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$ExePath = Join-Path $ProjectRoot "release\GeneSnapWorkbench.exe"
$ScriptPath = Join-Path $PSScriptRoot "GeneSnapWorkbench.iss"

if (-not (Test-Path -LiteralPath $ExePath)) {
    throw "请先运行 packaging/build_windows.ps1 生成 GeneSnapWorkbench.exe"
}
if (-not $ISCC) {
    $candidates = @(
        "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe",
        "$env:ProgramFiles(x86)\Inno Setup 6\ISCC.exe",
        "$env:ProgramFiles\Inno Setup 6\ISCC.exe"
    )
    $ISCC = $candidates | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
}
if (-not $ISCC) {
    throw "没有找到 Inno Setup 6 的 ISCC.exe"
}

& $ISCC $ScriptPath
if ($LASTEXITCODE -ne 0) {
    throw "Inno Setup 构建失败"
}

$Installer = Join-Path $ProjectRoot "release\installer\GeneSnapWorkbench_Setup_v0.3.4.exe"
$Hash = Get-FileHash -LiteralPath $Installer -Algorithm SHA256
Set-Content `
    -LiteralPath (Join-Path $ProjectRoot "release\installer\GeneSnapWorkbench_Setup_v0.3.4.sha256") `
    -Value "$($Hash.Hash.ToLowerInvariant())  GeneSnapWorkbench_Setup_v0.3.4.exe" `
    -Encoding ASCII
Write-Host "Windows installer: $Installer"
Write-Host "SHA-256: $($Hash.Hash)"
