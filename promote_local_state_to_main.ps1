# Ejecutar desde la raíz del repositorio local.
# Este script convierte el working tree local completo en la rama main remota
# y elimina la única rama remota adicional.
#
# Requiere: Git, Git LFS y acceso autenticado al remoto origin.

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Invoke-Git {
    param(
        [Parameter(ValueFromRemainingArguments = $true)]
        [string[]] $GitArgs
    )

    & git @GitArgs
    if ($LASTEXITCODE -ne 0) {
        throw "Falló: git $($GitArgs -join ' ')"
    }
}

function Get-GitOutput {
    param(
        [Parameter(ValueFromRemainingArguments = $true)]
        [string[]] $GitArgs
    )

    $output = & git @GitArgs
    if ($LASTEXITCODE -ne 0) {
        throw "Falló: git $($GitArgs -join ' ')"
    }
    return @($output)
}

Write-Host ""
Write-Host "=== 1. Verificación del repositorio ===" -ForegroundColor Cyan

$insideRepo = (Get-GitOutput rev-parse --is-inside-work-tree | Select-Object -First 1).Trim()
if ($insideRepo -ne "true") {
    throw "La carpeta actual no es un repositorio Git."
}

$repoRoot = (Get-GitOutput rev-parse --show-toplevel | Select-Object -First 1).Trim()
Set-Location $repoRoot
Write-Host "Repositorio local: $repoRoot"

$originUrl = (Get-GitOutput remote get-url origin | Select-Object -First 1).Trim()
Write-Host "Origin: $originUrl"

if ($originUrl -notmatch "sermunagu/PNNN_PNGMP_DOMP_fair_comparison(\.git)?$") {
    throw "El remoto origin no parece ser sermunagu/PNNN_PNGMP_DOMP_fair_comparison."
}

& git lfs version
if ($LASTEXITCODE -ne 0) {
    throw "Git LFS no está instalado. Instálalo antes de continuar."
}
Invoke-Git lfs install --local

# Un submódulo con cambios sin commit no puede conservarse exactamente mediante
# el commit del repositorio padre. Se aborta para no perder ese estado.
$dirtySubmodules = @(
    & git submodule foreach --recursive --quiet `
        'test -z "$(git status --porcelain)" || echo "$displaypath"'
)
if ($LASTEXITCODE -ne 0) {
    throw "No se pudo comprobar el estado de los submódulos."
}
if ($dirtySubmodules.Count -gt 0) {
    Write-Host "Submódulos con cambios locales:" -ForegroundColor Yellow
    $dirtySubmodules | ForEach-Object { Write-Host "  $_" }
    throw "Hay cambios sin commit dentro de un submódulo. Deben guardarse primero dentro del propio submódulo."
}

Write-Host ""
Write-Host "Estado que se publicará íntegramente:" -ForegroundColor Cyan
Invoke-Git status --short --branch
Invoke-Git lfs status

Write-Host ""
Write-Host "=== 2. Commit de TODO el estado local ===" -ForegroundColor Cyan

$sourceBranch = (Get-GitOutput branch --show-current | Select-Object -First 1).Trim()
if ([string]::IsNullOrWhiteSpace($sourceBranch)) {
    throw "HEAD está detached. Cambia primero a la rama que contiene tu estado local."
}

Write-Host "Rama local de origen: $sourceBranch"

# El usuario ha confirmado que TODO el working tree local pertenece al estado válido.
Invoke-Git add -A

& git diff --cached --quiet
$hasStagedChanges = ($LASTEXITCODE -ne 0)

if ($hasStagedChanges) {
    Invoke-Git commit -m "Promote validated local state to main"
} else {
    Write-Host "No hay cambios nuevos que commitear; se usará el HEAD local actual."
}

$sourceCommit = (Get-GitOutput rev-parse HEAD | Select-Object -First 1).Trim()
Write-Host "Commit local que será main: $sourceCommit"

# Tras git add/commit, el estado debe ser reproducible y limpio.
$remainingStatus = @(Get-GitOutput status --porcelain)
if ($remainingStatus.Count -gt 0) {
    Write-Host "Cambios que siguen sin quedar incluidos:" -ForegroundColor Yellow
    $remainingStatus | ForEach-Object { Write-Host $_ }
    throw "El working tree no ha quedado limpio. No se publicará un estado incompleto."
}

Write-Host ""
Write-Host "=== 3. Sincronización y sustitución de main ===" -ForegroundColor Cyan

# Actualiza origin/main justo antes del force-with-lease.
Invoke-Git fetch origin --prune

if ($sourceBranch -ne "main") {
    # Crea o mueve la rama local main exactamente al commit validado.
    Invoke-Git switch -C main $sourceCommit
} else {
    Write-Host "El estado validado ya está en la rama local main."
}

# --force-with-lease permite que el local sea la fuente de verdad sin pisar
# silenciosamente un cambio remoto aparecido después del fetch.
Invoke-Git push --force-with-lease --set-upstream origin main

# Fuerza la subida de todos los objetos LFS alcanzables desde main.
Invoke-Git lfs push --all origin main

Write-Host ""
Write-Host "=== 4. Eliminación de la rama remota adicional ===" -ForegroundColor Cyan

Invoke-Git fetch origin --prune

$remoteBranches = @(
    Get-GitOutput for-each-ref "--format=%(refname:strip=3)" refs/remotes/origin |
        ForEach-Object { $_.Trim() } |
        Where-Object {
            -not [string]::IsNullOrWhiteSpace($_) -and
            $_ -ne "HEAD" -and
            $_ -ne "main"
        } |
        Sort-Object -Unique
)

if ($remoteBranches.Count -gt 1) {
    Write-Host "Se han encontrado varias ramas remotas adicionales:" -ForegroundColor Yellow
    $remoteBranches | ForEach-Object { Write-Host "  $_" }
    throw "No se borrarán varias ramas automáticamente. Elimina solo las que correspondan."
}

if ($remoteBranches.Count -eq 1) {
    $branchToDelete = $remoteBranches[0]
    Write-Host "Eliminando rama remota adicional: $branchToDelete"
    Invoke-Git push origin --delete $branchToDelete
} else {
    Write-Host "No queda ninguna rama remota adicional."
}

# Elimina también la antigua rama local desde la que se promovió el estado.
if ($sourceBranch -ne "main") {
    $localBranches = @(Get-GitOutput branch "--format=%(refname:short)")
    if ($localBranches -contains $sourceBranch) {
        Invoke-Git branch -D $sourceBranch
    }
}

Invoke-Git fetch origin --prune

Write-Host ""
Write-Host "=== 5. Verificación final ===" -ForegroundColor Cyan

$localMain = (Get-GitOutput rev-parse main | Select-Object -First 1).Trim()
$remoteMainLine = (Get-GitOutput ls-remote origin refs/heads/main | Select-Object -First 1)
if ([string]::IsNullOrWhiteSpace($remoteMainLine)) {
    throw "No se ha encontrado refs/heads/main en el remoto."
}
$remoteMain = ($remoteMainLine -split "\s+")[0].Trim()

Write-Host "main local : $localMain"
Write-Host "main remoto: $remoteMain"

if ($localMain -ne $remoteMain) {
    throw "La rama main remota no coincide con el estado local."
}

$remainingRemoteBranches = @(
    Get-GitOutput for-each-ref "--format=%(refname:strip=3)" refs/remotes/origin |
        ForEach-Object { $_.Trim() } |
        Where-Object {
            -not [string]::IsNullOrWhiteSpace($_) -and $_ -ne "HEAD"
        } |
        Sort-Object -Unique
)

Write-Host "Ramas remotas restantes:"
$remainingRemoteBranches | ForEach-Object { Write-Host "  $_" }

if (-not ($remainingRemoteBranches.Count -eq 1 -and $remainingRemoteBranches[0] -eq "main")) {
    throw "La verificación detectó ramas remotas distintas de main."
}

Invoke-Git status --short --branch
Invoke-Git log --oneline --decorate -5

Write-Host ""
Write-Host "COMPLETADO: el estado local validado es ahora main y no queda otra rama remota." -ForegroundColor Green
