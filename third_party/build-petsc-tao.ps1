[CmdletBinding()]
param(
   [switch]$Clean,
   [switch]$ConfigureOnly,
   [ValidateRange(1, 64)]
   [int]$Jobs = 4
)

$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$manifest = Import-PowerShellDataFile (Join-Path $PSScriptRoot 'manifest.psd1')
$sourceRoot = Join-Path $root 'third_party\src'
$petscSource = Join-Path $sourceRoot 'petsc-msvc'
$win32feSource = Join-Path $sourceRoot 'win32fe'
$install = Join-Path $root 'third_party\install\petsc-oneapi'
$architecture = 'arch-udar-tao-opt'

$oneApiSetVars = @(
   'C:\OPT\Intel\oneAPI\setvars.bat',
   'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
) | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
$cygwin = @('C:\cygwin64\bin\bash.exe') |
           Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
if (-not $oneApiSetVars -or -not $cygwin)
{
   throw 'PETSc/TAO requires Intel oneAPI and Cygwin64 with python3 and make.'
}

# setvars can emit both PATH and Path in a case-sensitive process environment.
# Import the uppercase value and discard the stale mixed-case duplicate.
$environmentLines = cmd.exe /d /s /c ('call "{0}" intel64 >nul && set' -f $oneApiSetVars)
foreach ($line in $environmentLines)
{
   if ($line -cmatch '^([^=]+)=(.*)$' -and $matches[1] -cne 'Path')
   {
      [Environment]::SetEnvironmentVariable($matches[1], $matches[2], 'Process')
   }
}
$cl = Join-Path $env:VCToolsInstallDir 'bin\Hostx64\x64\cl.exe'
if (-not (Test-Path -LiteralPath $cl))
{
   throw "Visual C++ compiler not found at $cl."
}
$lib = Join-Path (Split-Path $cl) 'lib.exe'
$cygpath = Join-Path (Split-Path $cygwin) 'cygpath.exe'
$env:UDAR_CL = $cl
$env:UDAR_LIB = $lib
$env:UDAR_CYGPATH = $cygpath
$env:TEMP = 'C:\cygwin64\tmp'
$env:TMP = $env:TEMP

function Convert-ToCygwinPath([string]$Path)
{
   (& $cygpath -u ([IO.Path]::GetFullPath($Path))).Trim()
}

function Convert-ToShortCygwinPath([string]$Path)
{
   $short = (& $cygpath -d ([IO.Path]::GetFullPath($Path))).Trim()
   (& $cygpath -u $short).Trim()
}

function Invoke-Cygwin([string]$Command, [string]$Name)
{
   $statusName = '{0}-{1}-{2}.status' -f $Name, $PID, [Guid]::NewGuid().ToString('N')
   $status = Join-Path $root "third_party\build\$statusName"
   $statusUnix = Convert-ToCygwinPath $status
   $runtimePath = '/usr/bin:/cygdrive/c/OPT/Intel/oneAPI/mpi/latest/bin:' +
                  '/cygdrive/c/OPT/Intel/oneAPI/compiler/latest/bin:' +
                  '/cygdrive/c/OPT/Intel/oneAPI/mkl/latest/bin'
   $wrapped = "umask 000; export PATH=$runtimePath; $Command; rc=`$?; printf '%s\n' `$rc > '$statusUnix'; exit `$rc"
   & $cygwin --noprofile --norc -c $wrapped
   $deadline = (Get-Date).AddHours(4)
   while (-not (Test-Path -LiteralPath $status) -and (Get-Date) -lt $deadline)
   {
      Start-Sleep -Seconds 2
   }
   if (-not (Test-Path -LiteralPath $status))
   {
      throw "$Name did not complete within four hours."
   }
   $code = [int](Get-Content -LiteralPath $status -Raw).Trim()
   Remove-Item -LiteralPath $status -Force
   if ($code -ne 0) { throw "$Name failed with exit code $code." }
}

function Sync-PinnedRepository([string]$Repository, [string]$Revision,
                              [string]$Directory)
{
   if (-not (Test-Path -LiteralPath (Join-Path $Directory '.git')))
   {
      git -c core.autocrlf=false clone $Repository $Directory
      if ($LASTEXITCODE -ne 0) { throw "Failed to clone $Repository." }
      git -C $Directory checkout --detach $Revision
      if ($LASTEXITCODE -ne 0) { throw "Failed to check out $Revision." }
   }
   $actual = (git -C $Directory rev-parse HEAD).Trim()
   if ($actual -ne $Revision)
   {
      throw "$Directory is at $actual; expected pinned revision $Revision."
   }
}

New-Item -ItemType Directory -Force -Path $sourceRoot,
   (Join-Path $root 'third_party\build') | Out-Null
Sync-PinnedRepository $manifest.PetscRepository $manifest.PetscCommit $petscSource
Sync-PinnedRepository $manifest.Win32feRepository $manifest.Win32feCommit $win32feSource

$win32fePatch = Join-Path $PSScriptRoot 'patches\win32fe-current-windows.patch'
git -C $win32feSource apply --check $win32fePatch 2>$null
if ($LASTEXITCODE -eq 0)
{
   git -C $win32feSource apply $win32fePatch
}
else
{
   git -C $win32feSource apply --reverse --check $win32fePatch 2>$null
   if ($LASTEXITCODE -ne 0) { throw 'The win32fe compatibility patch does not apply.' }
}

if ($Clean)
{
   $archPath = Join-Path $petscSource $architecture
   foreach ($directory in $archPath, $install)
   {
      if (Test-Path -LiteralPath $directory)
      {
         $resolved = (Resolve-Path -LiteralPath $directory).Path
         if (-not $resolved.StartsWith($root + [IO.Path]::DirectorySeparatorChar,
                                      [StringComparison]::OrdinalIgnoreCase))
         {
            throw "Refusing to remove path outside the workspace: $resolved"
         }
         Remove-Item -LiteralPath $resolved -Recurse -Force
      }
   }
}

$clUnix = Convert-ToShortCygwinPath $cl
$win32feUnix = Convert-ToCygwinPath $win32feSource
$win32feProgram = Join-Path $win32feSource 'bin\win32fe.exe'
$win32feRuntime = Join-Path $win32feSource 'bin\win32feutils.dll'
if (-not (Test-Path -LiteralPath $win32feProgram) -or
    -not (Test-Path -LiteralPath $win32feRuntime))
{
   Invoke-Cygwin "cd '$win32feUnix'; make CC='$clUnix' CXX='$clUnix'" 'win32fe-build'
}

$petscUnix = Convert-ToCygwinPath $petscSource
$installUnix = Convert-ToCygwinPath $install
$compiler = "$win32feUnix/bin/win32fe.exe cl"
$archiver = "$win32feUnix/bin/win32fe.exe lib"
Invoke-Cygwin "cd '$petscUnix'; rm -f RDict.log RDict.log.bkp configure.log configure.log.bkp" 'petsc-log-clean'
$configure = @(
   "cd '$petscUnix';",
   "export PETSC_ARCH=$architecture;",
   "./configure --prefix='$installUnix'",
   ('--with-cc="' + $compiler + '" --with-cxx="' + $compiler + '" --with-fc=0'),
   ('--with-ar="' + $archiver + '" --with-ranlib=/usr/bin/true'),
   '--with-debugging=0 --with-shared-libraries=0 --with-single-library=0',
   'CFLAGS+="-MD -Zc:preprocessor"',
   'CXXFLAGS+="-MD -Zc:preprocessor"',
   '--with-scalar-type=real --with-precision=double --with-64-bit-indices=0',
   '--with-mpi=1 --with-mpi-include=/cygdrive/c/OPT/Intel/oneAPI/mpi/latest/include',
   '--with-mpi-lib=/cygdrive/c/OPT/Intel/oneAPI/mpi/latest/lib/impi.lib',
   '--with-mpiexec=/cygdrive/c/OPT/Intel/oneAPI/mpi/latest/bin/mpiexec.exe',
   '--with-blaslapack-lib="-L/cygdrive/c/OPT/Intel/oneAPI/mkl/latest/lib mkl_intel_lp64_dll.lib mkl_sequential_dll.lib mkl_core_dll.lib"',
   '--known-64-bit-blas-indices=0',
   '--with-x=0 --with-windows-graphics=0'
) -join ' '
Invoke-Cygwin $configure 'petsc-configure'

$petscHeader = Join-Path $petscSource "$architecture\include\petscconf.h"
$petscConfiguration = Get-Content -LiteralPath $petscHeader -Raw
if ($petscConfiguration -notmatch 'PETSC_USE_REAL_DOUBLE' -or
    $petscConfiguration -notmatch 'PETSC_HAVE_CXX' -or
    $petscConfiguration -match 'PETSC_USE_64BIT_INDICES' -or
    $petscConfiguration -match 'PETSC_HAVE_64BIT_BLAS_INDICES' -or
    $petscConfiguration -match 'PETSC_USE_DEBUG')
{
   throw 'PETSc configuration does not satisfy the UDAR C++/ABI contract.'
}

if (-not $ConfigureOnly)
{
   Invoke-Cygwin "cd '$petscUnix'; make -j$Jobs PETSC_DIR='$petscUnix' PETSC_ARCH=$architecture all" 'petsc-build'
   Invoke-Cygwin "cd '$petscUnix'; make PETSC_DIR='$petscUnix' PETSC_ARCH=$architecture install" 'petsc-install'
}

Write-Host "PETSc/TAO installed at $install"
Write-Host 'ABI: MSVC /MD, Intel MPI, oneMKL LP64, real double, 32-bit PetscInt, Release.'
