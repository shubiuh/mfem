# Windows oneAPI third-party dependencies

This directory provides one ABI-compatible dependency stack for the parallel
UDAR example and MFEM:

```
Intel MPI + oneAPI compilers + oneMKL
              |              |
              |              +-- MUMPS (Fortran/C/C++)
              +----------------- hypre (C/C++)
              +----------------- METIS (C)
                                      |
                                      +-- MFEM (C++)
```

All dependency components and MFEM must be configured in the same initialized
oneAPI shell.
Do not mix MS-MPI, MinGW, or a different MSVC runtime with this stack. MUMPS
uses `ifx` for Fortran and hypre uses the MSVC-compatible `icx` driver for
C/C++. Every component links Intel MPI and oneMKL from the same oneAPI release.

The optional `build-petsc-tao.ps1` driver adds PETSc/TAO to the same ABI family.
It builds PETSc with MSVC `/MD`, Intel MPI, oneMKL LP64, real double scalars,
and 32-bit indices. The static PETSc libraries are linked into `udar_tao.dll`,
so Python still crosses only the small UDAR C ABI.

The current Visual Studio 2026 STL requires Clang 20 or newer, while the
installed oneAPI compiler is 2024.1. MFEM is therefore compiled with `cl`; it
uses the same MSVC ABI as the `icx`-built hypre and MUMPS libraries.

## Prerequisite

Install the oneAPI HPC Toolkit, including Intel MPI, Intel Fortran, and
oneMKL, plus the Visual Studio C++ build tools. The script discovers oneAPI
from `C:\\OPT\\Intel\\oneAPI` first, then the default Intel installation
location.

## Build

From the **Intel oneAPI command prompt for Intel 64 for Visual Studio** (the
script additionally initializes the oneAPI environment):

```powershell
.\third_party\build-oneapi.ps1
```

The script clones pinned source revisions into `third_party/src`, builds and
installs METIS, hypre, and the scivision MUMPS superbuild under
`third_party/install/oneapi`, and then configures MFEM with the resulting
prefix. Use `-Clean` only to remove the generated build and install trees.
METIS 5.1.0 is downloaded from MFEM's TPL mirror and SHA-256 verified. The
script supplies the modern CMake policy floor, applies the included current-
MSVC math-header compatibility patch, and copies `metis.h` and `metis.lib`
explicitly because that legacy release disables CMake install rules on MSVC.

The MUMPS build requests both `dmumps` and `zmumps`: MFEM's double-precision
UDAR implementation uses `zmumps` when native complex MUMPS is available.

## Configure MFEM after the dependencies

Open the Intel oneAPI command prompt for Intel 64 for Visual Studio, or
initialize the environment from PowerShell:

```powershell
cmd /d /s /c 'call "C:\OPT\Intel\oneAPI\setvars.bat" intel64 && powershell'
```

Then configure and build the library and all standard examples:

```powershell
cmake --preset oneapi-mpi-mumps
cmake --build --preset oneapi-mpi-mumps --target mfem
cmake --build --preset oneapi-mpi-mumps --target examples --parallel
```

The generated MFEM build tree is `build/oneapi-mpi-mumps`; the static archive
is `build/oneapi-mpi-mumps/mfem.lib`. METIS is enabled because `ParMesh`
partitions serial input meshes even for a one-rank MPI run; this is required by
`ex1p`, `udar3p`, and similar parallel examples. ParMETIS is not required for
this workflow. MUMPS itself uses MPI and oneMKL LAPACK/BLAS in this stack.

## Smoke tests

Run these from the initialized oneAPI environment after building. Change to
the relevant executable directory before launching parallel programs because
their default mesh paths are relative to the process working directory. The
UDAR miniapp uses `../../data/inline-hex.mesh` from its build directory:

```powershell
.\build\oneapi-mpi-mumps\examples\ex1.exe -no-vis
```

```powershell
Push-Location .\build\oneapi-mpi-mumps\examples
mpiexec -n 2 .\ex1p.exe -no-vis
Pop-Location
Push-Location .\build\oneapi-mpi-mumps\miniapps\udar
mpiexec -n 2 .\udar3p.exe `
  -mumps -freq 2000 -ns 1 -sp 0.08 -out udar3p_smoke.csv
Pop-Location
mpiexec -n 2 .\udar3p.exe `
  -gmres -freq 2000 -ns 1 -sp 0.08 -out udar3p_gmres_smoke.csv
Pop-Location
```

The parallel UDAR example can be built with
`cmake --build --preset oneapi-mpi-mumps --target udar3p` and run through
`mpiexec` from the same shell. The commands above write their CSV files to
`build/oneapi-mpi-mumps/examples`. To launch `udar3p` from the repository root
instead, supply the mesh explicitly with
`-m .\build\oneapi-mpi-mumps\data\inline-hex.mesh`.
