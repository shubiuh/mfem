[CmdletBinding()]
param(
   [switch]$Clean,
   [switch]$ConfigureOnly
)

$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$oneApiCandidates = @(
   'C:\OPT\Intel\oneAPI\setvars.bat',
   'C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
)
$oneApiSetVars = $oneApiCandidates | Where-Object {
   Test-Path -LiteralPath $_
} | Select-Object -First 1
if (-not $oneApiSetVars)
{
   throw 'Intel oneAPI HPC Toolkit is not installed in C:\\OPT\\Intel\\oneAPI or C:\\Program Files (x86)\\Intel\\oneAPI. Install MPI, ifx, and oneMKL and rerun this script.'
}

# Import the environment emitted by setvars.bat into this PowerShell process.
$envLines = cmd.exe /d /s /c ('call "{0}" intel64 >nul && set' -f $oneApiSetVars)
foreach ($line in $envLines)
{
   if ($line -match '^([^=]+)=(.*)$')
   {
      [Environment]::SetEnvironmentVariable($matches[1], $matches[2], 'Process')
   }
}

# On this Windows oneAPI installation, icx is the MSVC-compatible driver for
# both C and C++; it keeps C, C++, and ifx Fortran objects ABI-compatible while
# retaining C99 support needed by hypre.
foreach ($command in 'icx', 'ifx', 'mpiexec', 'cmake', 'ninja', 'git')
{
   if (-not (Get-Command $command -ErrorAction SilentlyContinue))
   {
      throw "oneAPI environment is incomplete: '$command' is not available. Ensure the HPC Toolkit components are installed."
   }
}

$manifest = Import-PowerShellDataFile (Join-Path $PSScriptRoot 'manifest.psd1')

$src = Join-Path $root 'third_party\src'
$build = Join-Path $root 'third_party\build\oneapi'
$install = Join-Path $root 'third_party\install\oneapi'
if ($Clean)
{
   foreach ($directory in $build, $install)
   {
      if (Test-Path -LiteralPath $directory) { Remove-Item -LiteralPath $directory -Recurse -Force }
   }
}
New-Item -ItemType Directory -Force -Path $src, $build, $install | Out-Null

function Sync-PinnedRepository([string]$Repository, [string]$Tag, [string]$Directory)
{
   if (-not (Test-Path -LiteralPath (Join-Path $Directory '.git')))
   {
      git clone --branch $Tag --depth 1 $Repository $Directory
   }
   else
   {
      git -C $Directory fetch --depth 1 origin $Tag
      git -C $Directory checkout --detach FETCH_HEAD
   }
   if ($LASTEXITCODE -ne 0)
   {
      throw "Failed to fetch $Repository at $Tag."
   }
}

Sync-PinnedRepository $manifest.HypreRepository $manifest.HypreTag (Join-Path $src 'hypre')
Sync-PinnedRepository $manifest.MumpsRepository $manifest.MumpsTag (Join-Path $src 'mumps-superbuild')

$metisSource = Join-Path $src 'metis-src'
if (-not (Test-Path -LiteralPath (Join-Path $metisSource 'CMakeLists.txt')))
{
   New-Item -ItemType Directory -Force -Path $metisSource | Out-Null
   $metisArchive = Join-Path $src 'metis-5.1.0.tar.gz'
   Invoke-WebRequest -Uri $manifest.MetisArchive -OutFile $metisArchive
   $actualHash = (Get-FileHash -LiteralPath $metisArchive -Algorithm SHA256).Hash.ToLowerInvariant()
   if ($actualHash -ne $manifest.MetisSha256)
   {
      throw "METIS archive checksum mismatch: expected $($manifest.MetisSha256), got $actualHash."
   }
   Push-Location $src
   try
   {
      tar.exe -xzf $metisArchive -C $metisSource --strip-components=1
      if ($LASTEXITCODE -ne 0) { throw 'Failed to extract the METIS archive.' }
   }
   finally
   {
      Pop-Location
   }
}

# METIS 5.1.0 predates modern MSVC's C99 math declarations. Its legacy rint
# macro corrupts the declaration in current Windows SDK headers.
$metisPatch = Join-Path $PSScriptRoot 'patches\metis-5.1.0-modern-msvc.patch'
Push-Location $root
try
{
   git apply --check --directory='third_party/src/metis-src' $metisPatch 2>$null
   if ($LASTEXITCODE -eq 0)
   {
      git apply --directory='third_party/src/metis-src' $metisPatch
      if ($LASTEXITCODE -ne 0) { throw 'Failed to patch METIS for modern MSVC.' }
   }
   else
   {
      git apply --reverse --check --directory='third_party/src/metis-src' $metisPatch 2>$null
      if ($LASTEXITCODE -ne 0) { throw 'The METIS MSVC compatibility patch does not apply.' }
   }
}
finally
{
   Pop-Location
}

$common = @(
   '-G', 'Ninja',
   "-DCMAKE_INSTALL_PREFIX=$install",
   '-DCMAKE_BUILD_TYPE=Release',
   '-DCMAKE_C_COMPILER=icx',
   '-DCMAKE_CXX_COMPILER=icx',
   '-DCMAKE_Fortran_COMPILER=ifx',
   '-DBUILD_SHARED_LIBS=OFF'
)

# METIS is a C library. Build it with MSVC so it remains compatible with the
# cl-built MFEM library even when the installed icx is older than the current
# Visual Studio headers.
$metisBuild = Join-Path $build 'metis'
cmake -S $metisSource -B $metisBuild -G Ninja `
   "-DCMAKE_INSTALL_PREFIX=$install" `
   '-DCMAKE_BUILD_TYPE=Release' `
   '-DCMAKE_C_COMPILER=cl' `
   '-DCMAKE_POLICY_VERSION_MINIMUM=3.5' `
   '-DSHARED=OFF'
if ($LASTEXITCODE -ne 0) { throw 'METIS configuration failed.' }
if (-not $ConfigureOnly)
{
   cmake --build $metisBuild --parallel
   if ($LASTEXITCODE -ne 0) { throw 'METIS build failed.' }
   # METIS 5.1.0 disables its install rules under MSVC, so install the two
   # artifacts MFEM consumes explicitly.
   cmake -E make_directory (Join-Path $install 'include') (Join-Path $install 'lib')
   cmake -E copy_if_different (Join-Path $metisSource 'include\metis.h') (Join-Path $install 'include\metis.h')
   cmake -E copy_if_different (Join-Path $metisBuild 'libmetis\metis.lib') (Join-Path $install 'lib\metis.lib')
   if ($LASTEXITCODE -ne 0) { throw 'METIS installation failed.' }
}

$hypreSource = Join-Path $src 'hypre\src'
$hypreBuild = Join-Path $build 'hypre'
cmake -S $hypreSource -B $hypreBuild @common '-DHYPRE_ENABLE_MPI=ON' '-DHYPRE_ENABLE_OPENMP=OFF' '-DHYPRE_BUILD_EXAMPLES=OFF' '-DHYPRE_BUILD_TESTS=OFF'
if ($LASTEXITCODE -ne 0) { throw 'Hypre configuration failed.' }
if (-not $ConfigureOnly)
{
   cmake --build $hypreBuild --parallel
   if ($LASTEXITCODE -ne 0) { throw 'Hypre build failed.' }
   cmake --install $hypreBuild
   if ($LASTEXITCODE -ne 0) { throw 'Hypre installation failed.' }
}

$mumpsSource = Join-Path $src 'mumps-superbuild'
$mumpsBuild = Join-Path $build 'mumps'
cmake -S $mumpsSource -B $mumpsBuild @common '-DMUMPS_parallel=ON' '-DBUILD_SINGLE=OFF' '-DBUILD_DOUBLE=ON' '-DBUILD_COMPLEX=OFF' '-DBUILD_COMPLEX16=ON' '-DMUMPS_scalapack=OFF'
if ($LASTEXITCODE -ne 0) { throw 'MUMPS configuration failed.' }
if (-not $ConfigureOnly)
{
   cmake --build $mumpsBuild --parallel
   if ($LASTEXITCODE -ne 0) { throw 'MUMPS build failed.' }
   cmake --install $mumpsBuild
   if ($LASTEXITCODE -ne 0) { throw 'MUMPS installation failed.' }
}

Write-Host "Dependencies installed at $install"
Write-Host 'Next: cmake --preset oneapi-mpi-mumps'
