# UDAR inversion: Python interface with C++ ROL and PETSc/TAO cores

This is the first executable inversion workflow for `examples/udar3p.cpp`.
Python handles configuration, CSV observations and the forward adapter. The
optimization algorithm, objective, regularization, finite-difference derivatives
and bound handling run in C++ using either Trilinos/ROL or PETSc/TAO. No SciPy
optimizer is used.

```text
Python API / JSON configuration
          |
          v
udar_rol or udar_tao shared library: bounded limited-memory quasi-Newton
          |
          v
synchronous prediction callback -> persistent udar3p server -> MFEM + MUMPS/hypre
```

The current workflow inverts one or both isotropic conductivities in S/m:
`sigma_background` and `sigma_target`. The target parameter is available for
slabs and boxes. Geometry, anisotropy, source/receiver layout and mesh stay fixed.
Each frequency/station/polarization is still solved by the existing forward code.

By default, the wrapper launches one long-lived `udar3p` process for an inversion.
MPI, the mesh, and the Nedelec finite-element space are initialized once; each
request updates the scalar isotropic background/target conductivity and writes a
new tensor CSV. Material-dependent forms and numeric factorizations are rebuilt
for every request. Set `forward.persistent` to `false` to retain the one-process-
per-evaluation compatibility path. The server log and each request's output make
the optimizer-to-forward interaction auditable.

## Build the C++ cores

Build `udar3p` first using the existing oneAPI/MUMPS instructions in
`examples/udar/README.md`. The inversion DLL is a separate CMake project and does
not rebuild or link MFEM. Source dependencies are pinned in
`third_party/manifest.psd1`; download once from the repository root:

```sh
git clone --depth 1 --filter=blob:none --sparse --branch trilinos-release-17-2-1 https://github.com/trilinos/Trilinos.git third_party/src/trilinos-rol
git -C third_party/src/trilinos-rol sparse-checkout set packages/rol
git clone --depth 1 --branch v1.14 https://github.com/zeux/pugixml.git third_party/src/pugixml-rol
```

Tested revisions: Trilinos `a0de2fe38d6b68f88fc4557abb003ff379561933` and
pugixml `db78afc2b7d8f043b4bc6b185635d949ea2ed2a8`.

In an initialized Windows oneAPI/Visual C++ command prompt:

```sh
cmake -S inversion/cpp -B build/udar-inversion -G Ninja -DCMAKE_BUILD_TYPE=Release -DUDAR_BLAS_LIBRARIES=C:/OPT/Intel/oneAPI/mkl/latest/lib/mkl_rt.lib
cmake --build build/udar-inversion --parallel 4
```

If Ninja is not on PATH, set `CMAKE_MAKE_PROGRAM` to its executable, or use an
installed CMake generator. With a multi-configuration generator pass
`--config Release` and specify the resulting DLL using `--library`.

On Linux with a BLAS/LAPACK development installation:

```sh
cmake -S inversion/cpp -B build/udar-inversion -DCMAKE_BUILD_TYPE=Release
cmake --build build/udar-inversion -j4
```

The Windows oneAPI configuration has been executed and tested locally; the Linux
build recipe is provided but has not been executed here. Edit the example JSON
executable paths for Linux. CMake accepts `ROL_SOURCE_DIR`, `PUGIXML_SOURCE_DIR`
and `UDAR_BLAS_LIBRARIES` for alternative installations. The adapter uses ROL's
standalone compatibility components, including its bundled Teuchos-lite support.
It deliberately avoids a complete Trilinos build. Other Trilinos releases are
not yet certified against this interface.

### PETSc/TAO on Windows

The TAO backend uses PETSc 3.25.5 `TAOBLMVM`, built as static Release component
libraries and linked into `udar_tao.dll`. The component layout avoids Cygwin's
command-line limit when MSVC archives the complete PETSc object set. Run the
dependency driver from PowerShell:

```powershell
.\third_party\build-petsc-tao.ps1 -Jobs 4
```

The driver requires Cygwin64 with `python3` and `make`, Visual C++ build tools,
and the same Intel oneAPI installation used for MFEM. It downloads the pinned
PETSc and official `win32fe` sources, applies
`third_party/patches/win32fe-current-windows.patch`, and installs PETSc under
`third_party/install/petsc-oneapi`. The patch fixes current Windows ACL handling
and supplies explicit Cygwin, compiler, and librarian paths. It changes the
build helper only; it does not patch PETSc numerical code.

The ABI is fixed to Visual C++ `/MD`, Release, real double scalars, 32-bit
`PetscInt`, Intel MPI, and oneMKL LP64. PETSc is static, which avoids another
C++ runtime boundary and leaves only the stable UDAR C ABI visible to Python.
Do not substitute a MinGW PETSc library in this build because it is not
ABI-compatible with the MSVC MFEM stack.

After PETSc is installed, build both inversion libraries in the initialized
oneAPI/Visual C++ environment:

```powershell
cmake -S inversion/cpp -B build/udar-inversion -G Ninja `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_CXX_COMPILER=cl `
  -DUDAR_BUILD_ROL=ON -DUDAR_BUILD_TAO=ON `
  -DUDAR_BLAS_LIBRARIES=C:/OPT/Intel/oneAPI/mkl/latest/lib/mkl_rt.lib `
  -DPETSC_DIR=third_party/install/petsc-oneapi
cmake --build build/udar-inversion --parallel 4
```

On another Windows installation, set `UDAR_PETSC_EXTERNAL_LIBRARIES` if Intel
MPI or oneMKL are outside `C:/OPT/Intel/oneAPI`. The list must match PETSc's
installed `petscvariables` file.

## Python setup and runs

Use Python 3.10 or newer. The wrapper uses only the Python standard library.
Either install with `python -m pip install -e inversion` or, from PowerShell:

```powershell
$env:PYTHONPATH = Join-Path $PWD 'inversion/python'
python -m udar_inversion synthetic inversion/configs/uniform.json --output build/uniform-inversion
python -m udar_inversion synthetic inversion/configs/box.json --output build/box-inversion
python -m udar_inversion synthetic inversion/configs/uniform.json --backend tao --output build/uniform-tao-inversion
```

Use a new output directory for each run; existing runs are never overwritten.
In Windows, initialize oneAPI before launching Python so the forward child and
oneMKL can locate their runtime libraries. The wrapper registers `MKLROOT/bin`
for DLL loading. Additional DLL directories can be supplied through the
semicolon-separated `UDAR_DLL_DIRS` variable. Specify a non-default core with
`--library path/to/udar_rol.dll` or `UDAR_ROL_LIBRARY`.
For TAO, use `--backend tao`; its corresponding override is
`UDAR_TAO_LIBRARY`.

For observed data:

```powershell
python -m udar_inversion invert inversion/configs/box.json --observed measurements.csv --output build/measured-inversion
```

The public low-level API is `udar_inversion.minimize(predict, initial, lower,
upper, observed, stddev, backend="rol" | "tao", ...)`; parameter vectors are **natural-log
conductivities**. The JSON workflow takes physical conductivities and converts
them automatically. `udar_inversion.workflow.run_inversion` and
`synthetic_recovery` are also callable from Python.

## Objective and data contract

For model `m = log(sigma)` the C++ objective is

```text
Phi(m) = 0.5 sum_i ((d_pred_i(m) - d_obs_i) / stddev_i)^2
         + 0.5 beta sum_j (m_j - m_prior_j)^2
```

Bounds are supplied in S/m and transformed to log space. `prior` defaults to the
initial conductivity; `optimizer.regularization` is beta (default zero).
The C++ gradient differentiates residuals using central finite differences with
absolute log-space step `difference_step`; near bounds it uses second-order
inward differences. Callback failures and non-finite predictions terminate the
run, rather than producing an apparently valid objective. Exceptions are caught
on each side of the C ABI and propagated back to Python.

Input observations use the complete `udar3p` tensor CSV format: station,
measured depth, frequency, spacing, Tx/Rx coordinates and all 18 real/imaginary
H-tensor columns. Select the fitted columns with `channels`; the default is all
18. Data is sorted by `(station, frequency_hz, spacing_m)` and then by the
configured channel order. Row order is irrelevant. Duplicate keys, missing
channels, non-finite values and differing acquisition keys/coordinates are
rejected. Frequency/spacing keys must match numerically exactly. Tool frame,
moment normalization, SI units and phasor convention must also match: the CSV
alone cannot verify these, so retain the matching forward configuration.

Weights are fixed from the observations, with a positive absolute noise floor:

```text
stddev(Hab_real) = stddev(Hab_imag)
                 = sqrt((relative * abs(Hab_observed))^2 + absolute^2)
```

The initial implementation supports these diagonal weights, not correlated
measurement errors. Example noise levels are numerical test weights, not an
instrument noise specification. The synthetic examples are noiseless and fit
imaginary diagonal components that are conductivity-sensitive at the selected
frequencies. They deliberately use the tiny existing mesh and PEC outer boundary
to test the inversion plumbing; they do not establish field accuracy or a
realistic open-domain UDAR response.

## Results and failure reporting

Each inversion directory contains `configuration.json`, `observed.csv`,
`predicted.csv`, `result.json`, `residuals.json`, `rol.log` or `tao.log`, and individual forward
CSVs plus a server log and evaluation trace. Every forward evaluation also writes
ParaView VTK data under `forward/vtk/evaluation_00001/`: `model/model.pvd`
contains `Rh_ohm_m` and `Rv_ohm_m`, while each `fields_fNNN_stationNNN/*.pvd`
contains real and imaginary `E` and `H` for the three transmitter axes. The
`iteration_models.csv` file maps every evaluated inversion model to its VTK
directory; finite-difference evaluations are included so the result is fully
auditable. Synthetic runs write the reference model under `true_model/vtk/`.
Open a `.pvd` file in ParaView; it loads the associated parallel `.pvtu` and
per-rank `.vtu` files automatically. Results include conductivity,
resistivity, weighted RMS, objective reduction, projected-gradient norm,
iterations, forward counts, and persistent-server starts. `converged` is true
only when the selected C++ backend reports convergence. An iteration limit or
tiny step is recorded separately and the command exits with code 2.
Exceptions write `failure.json` when raised during inversion. A synthetic run
also writes `recovery.json`; its check requires convergence and <1% relative
parameter error. It is a same-discretization software recovery check.

## Checks

With the compiled core discoverable and the runtime initialized:

```powershell
python -m unittest discover -s inversion/tests -v
```

These tests exercise each available ROL and TAO library: coupled least squares,
an active bound, regularization, callback failure, invalid data sizes, non-finite data,
non-convergence reporting, and acquisition alignment. The two synthetic commands
above additionally exercise the real MFEM executable.

Verified locally on 2026-09-10 with Windows, Visual C++, oneMKL, and the existing
single-rank MUMPS forward executable:

| Recovery check | Initial conductivity (S/m) | Recovered conductivity (S/m) | ROL iterations | Forward evaluations | Server starts |
|---|---|---|---|---|---|
| Uniform | background 0.1 | background 0.3000000000000001 (truth 0.3) | 9 | 32 | 1 |
| 3D box | background 0.1, target 0.1 | background 0.1999999999976, target 0.05000000000549 (truth 0.2, 0.05) | 13 | 73 | 1 |

All 11 automated tests passed. These results establish end-to-end software
consistency, not identifiability with noisy field data or physical accuracy.

## Next stages

1. Refactor the persistent executable loop into a reusable C++ forward context,
   cache invariant source/acquisition work, and measure whether symbolic solver
   analysis can be reused safely across material updates.
2. Add discrete tangent/adjoint operations and test Jacobian-adjoint dot products
   and Taylor convergence. The current finite differences cost approximately
   `2 * number_of_parameters` extra forward evaluations per gradient and are
   intended for small parameter counts. They are not suitable for cellwise 3D
   inversion.
3. Add a field parameterization, spatial regularization and matrix-free
   Gauss-Newton/ROL SimOpt interfaces. Preserve the existing prediction/data
   contract while extending the C ABI with derivative callbacks.
4. Add geometry derivatives only with a validated interface treatment. A positive
   `interface_width` smooths the forward model but does not itself validate its
   derivatives or physical bias.

Set `forward.launcher` to an argument list such as `["mpiexec", "-n", "2"]` to
request MPI forward solves. Python and ROL should run once, outside `mpiexec`.
Optimization is serial in this milestone; distributed control vectors and
parallel source/frequency groups are future work. MPI forward launches have not
been validated by the single-rank recovery checks. A launcher timeout may need
MPI-runtime-specific cleanup on a cluster.

Upstream references: [ROL](https://trilinos.github.io/rol.html),
[ROL source](https://github.com/trilinos/Trilinos/tree/trilinos-release-17-2-1/packages/rol),
[pugixml](https://github.com/zeux/pugixml/tree/v1.14), and
[PETSc/TAO](https://petsc.org/release/manual/tao/).
