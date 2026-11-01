# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
#
# Runs kronspice's container checks on Windows (Docker Desktop), from the
# repository root. VS Code: Terminal > Run Task > "kronspice: ...".
#
#   powershell -ExecutionPolicy Bypass -File internal/golden/docker/run.ps1 -Step all
#
# Steps:
#   build   build the tools image (Go 1.24, gcc for -race, ngspice)
#   golden  regenerate testdata/golden from internal/golden/ref/ngspice/*.cir
#   twin    write testdata/golden/ngspice/*.csv: every test netlist run
#           through ngspice itself (ngspice_all.sh), for TestNgspiceTwin
#   test    go vet + go test -race ./... (CLAUDE.md §7, §8 rule 10)
#   all     build, golden, twin, test
#
# Every step writes its full output to internal/golden/docker/logs/<step>.log
# and a one-line summary to logs/summary.txt, so the results can be read
# back from the folder afterwards.
param(
    [ValidateSet('build', 'golden', 'twin', 'test', 'all')]
    [string]$Step = 'all'
)
$ErrorActionPreference = 'Stop'
# Docker prints UTF-8 (Ω, Z₀, −); decode it as such instead of the console's
# legacy code page.
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$root = Resolve-Path (Join-Path $PSScriptRoot '..\..\..')
Set-Location $root
$logs = Join-Path $PSScriptRoot 'logs'
New-Item -ItemType Directory -Force -Path $logs | Out-Null
$summary = Join-Path $logs 'summary.txt'
$image = 'kronspice-tools'
# Named volumes keep the Go module and build caches between runs.
$mounts = @('-v', "${root}:/src", '-v', 'kronspice-gomod:/go/pkg/mod', '-v', 'kronspice-gocache:/root/.cache/go-build')

function Invoke-Step([string]$name, [string[]]$dockerArgs) {
    $log = Join-Path $logs "$name.log"
    $start = Get-Date
    Write-Host "== ${name}: docker $($dockerArgs -join ' ')"
    "# $name  $(Get-Date -Format o)" | Set-Content -Encoding utf8 $log
    # Merge stderr into stdout so the log is complete; keep the exit code.
    $ErrorActionPreference = 'Continue'
    # Write the log as UTF-8 line by line (Tee-Object in Windows PowerShell
    # 5.1 writes UTF-16).
    & docker @dockerArgs 2>&1 | ForEach-Object {
        $l = "$_"
        Write-Host $l
        Add-Content -Encoding utf8 -Path $log -Value $l
    }
    $code = $LASTEXITCODE
    $ErrorActionPreference = 'Stop'
    $secs = [int]((Get-Date) - $start).TotalSeconds
    $line = "{0}  {1,-6} exit {2}  {3} s  {4}" -f (Get-Date -Format s), $name, $code, $secs, $log
    Add-Content -Encoding utf8 $summary $line
    Write-Host $line
    if ($code -ne 0) { throw "step $name failed (exit $code); see $log" }
}

# Windows PowerShell 5.1 turns a native command's stderr into an error when
# it is redirected, so relax the preference around docker calls.
$ErrorActionPreference = 'Continue'
docker version --format '{{.Server.Version}}' *> $null
$ok = $LASTEXITCODE -eq 0
$ErrorActionPreference = 'Stop'
if (-not $ok) { throw 'Docker is not running: start Docker Desktop first.' }

if ($Step -in 'build', 'all') {
    Invoke-Step 'build' @('build', '-t', $image, 'internal/golden/docker')
}
if ($Step -in 'golden', 'all') {
    Invoke-Step 'golden' (@('run', '--rm') + $mounts + @($image, 'sh', 'internal/golden/ref/ngspice_ref.sh'))
}
if ($Step -in 'twin', 'all') {
    Invoke-Step 'twin' (@('run', '--rm') + $mounts + @($image, 'sh', 'internal/golden/ref/ngspice_all.sh'))
}
if ($Step -in 'test', 'all') {
    Invoke-Step 'test' (@('run', '--rm') + $mounts + @($image, 'sh', '-c', 'go version && go vet ./... && go test -race -count=1 -timeout 60m ./...'))
}
Write-Host 'done'
