# MFEM 3D UDAR forward solver

The Python workflow with a C++ Trilinos/ROL inversion core is documented in
[inversion/README.md](../../inversion/README.md). It supports bounded recovery
of background/target isotropic conductivities and synthetic recovery checks.
The staged performance roadmap based on a review of Palace through v0.18.0 is
documented in [PALACE_UPGRADE_PLAN.md](PALACE_UPGRADE_PLAN.md).

`udar3p` is a first parallel forward-modeling application for ultra-deep
azimuthal resistivity (UDAR) studies. It uses MFEM Nedelec elements for the
frequency-domain electric field and recovers the full complex 3x3 magnetic
transfer tensor in the tool coordinate system. Receiver locations are placed
behind each transmitter along the negative tool axis by the requested spacing.

The miniapp is self-contained under `miniapps/udar`: `udar3p.cpp` owns command
line handling and the survey workflow, while the adaptation, PML, distributed
point-location, and linear-solver implementations live in separate
`udar3p_*.hpp/.cpp` modules. It has dedicated CMake and make build files and
keeps the executable name `udar3p` for compatibility with the inversion
wrapper.

## Capabilities

- Multiple frequencies and transmitter-receiver spacings in one run.
- Multiple stations along a straight well trajectory.
- Three orthogonal magnetic-dipole transmitters and receivers (`Hxx` through
  `Hzz`), using an ideal point dipole or a resolved finite-radius loop.
- Full-tensor, diagonal-anisotropic, or isotropic conductivity values.
- Uniform, dipping-slab, and bounded 3D-geobody formation models.
- Optional smooth formation interfaces for differentiable geometry parameters.
- Native complex MUMPS when `zmumps` is available, real-block MUMPS otherwise,
  or iterative flexible GMRES with shifted Hypre AMS and distributed ILU.
- Tidy CSV results plus 3x3 amplitude and phase curve panels.
- ParaView VTK output for the Rh/Rv model and complex electric and magnetic
  fields.

The CSV stores `H(receiver axis, transmitter axis)` normalized by transmitter
moment. With SI mesh coordinates its units are `m^-3`. Real and imaginary
parts use the `exp(+i omega t)` phasor convention. The H tensor axes use the
tool frame; conductivity tensors supplied on the command line use global mesh
axes.

## Build

Configure MFEM with MPI and Hypre. Enable MUMPS to use the direct solver:

```sh
cmake -S . -B build \
  -DMFEM_USE_MPI=ON \
  -DMFEM_USE_METIS=ON \
  -DMFEM_USE_MUMPS=ON \
  -DMETIS_DIR=/path/to/metis \
  -DMUMPS_DIR=/path/to/mumps
cmake --build build --target udar3p -j
```

The MUMPS installation should include `zmumps` for the native complex path.
If it does not, `udar3p` automatically uses MFEM's real 2x2 representation
with the real MUMPS solver. METIS is needed because the application constructs
a parallel mesh by partitioning a serial mesh file.

## Run

The default mesh-scale survey is small enough for a workflow check:

```sh
mpirun -np 4 ./build/miniapps/udar/udar3p \
  -m data/inline-hex.mesh -rs 1 -rp 1 \
  -freq 2000,6000,12000,48000,96000 \
  -sp 0.08,0.16,0.24 -gmres \
  -out udar_tensor.csv
```

Use `-mumps` for the direct solver. Run `udar3p -h` to see all geometry,
material, trajectory, source, boundary, and iterative-solver controls.
MUMPS symbolic ordering is reused by default across persistent-server material
updates on the same mesh. Use `-no-mumps-reuse` only when debugging a suspected
sparsity-pattern change.
The `-gmres` path uses FGMRES with distributed ILU by default. Add
`-pc ams` for the positive shifted H(curl) operator with Hypre AMS. The AMS
defaults follow Palace's driven-solver choices: cycle type 14, two AMS cycles,
and one smoothing sweep. `-ams-mode complex` applies FGMRES to MFEM's complex
operator without assembling the full 2-by-2 matrix; `-ams-mode real` uses the
explicit real block for diagnosis and comparison. Use `-ams-cycle`, `-ams-it`,
`-ams-smooth`, and `-ams-print` to tune Hypre. A failed AMS solve retries with
distributed ILU by default; disable this with `-no-ams-fallback` when measuring
AMS robustness. The default relative tolerance is `1e-4`. MUMPS remains the
reference solver.

### Tool-region octree mesh

Add `-octree` to retain the coarse hexahedral mesh away from the tool while
applying graded nonconforming octree refinement around every transmitter and
receiver:

```sh
mpirun -np 4 ./build/miniapps/udar/udar3p \
  -m data/inline-hex.mesh -freq 2000,48000,96000 -sp 0.08,0.16 \
  -octree -octree-l 3 -octree-r 0.08 -octree-pad 0.02 \
  -mumps -out udar_octree.csv
```

`-octree-l` is the number of local levels; `-octree-r` is the radius of the
finest transmitter-receiver tube; and `-octree-pad` extends that tube radially
and beyond the first/last tool position. The first selected level is wider
than the next, creating coarse transition rings instead of refining the entire
formation. The base mesh must contain hexahedra. Start with a conservative
radius that resolves the loop or point-source neighborhood, then verify both
the tensor response and the outer-boundary error as the radius, level count,
and domain size change.

Use `-si 12,13,14` (or `--station-indices`) to calculate selected original
survey indices only. The CSV retains the original station and measured-depth
values, so it can be merged with earlier station data without renumbering.
This is useful for incremental inversions; it does not change the physics at a
selected station.

### ParaView output

Add `-vtk DIRECTORY` to write VTK files. For each CSV output named
`NAME.csv`, the writer creates `DIRECTORY/NAME/model/model.pvd` with
per-rank legacy `.vtk` pieces containing `Rh_ohm_m` and `Rv_ohm_m` as
**cell data**, and one `fields_fNNN_stationNNN/*.pvd` collection per frequency
and transmitter station. The field collections contain `E_real_x`, `E_imag_x`,
`E_magnitude_x`, `H_real_x`, `H_imag_x`, `H_magnitude_x`, and corresponding
y/z transmitter fields. Open a `.pvd` file in ParaView; the model collection
locates its per-rank `.vtk` pieces and field collections locate their `.pvtu`
and `.vtu` pieces. Rh is the average xx/yy component of the inverse
conductivity tensor and Rv is its zz component, so these labels use the global
mesh axes.

### Tetrahedra and PML

`udar3p` accepts both hexahedral and tetrahedral 3D meshes. A Cartesian PML
does not require hexahedra: its stretch tensor is evaluated in physical
coordinates at quadrature points. Axis-aligned structured hex layers are
usually more efficient and give cleaner resolution. A tetrahedral PML is valid
when its elements are well shaped, size transitions are gradual, and the layer
is adequately resolved. The `-octree`
option performs MFEM nonconforming nested refinement around the tool path for
either element type; its historical option name is retained for compatibility.
For example, `-m ../../data/inline-tet.mesh -octree -octree-l 3` starts from the
moderate tetrahedral mesh and refines only the transmitter-receiver region.

Use `-pml` to add a Cartesian, polynomially graded coordinate-stretch PML at
all six outer bounding-box faces. The implementation transforms both the
curl-curl tensor and the complete complex mass tensor, including anisotropic
formation conductivity. `-pml-t`, `-pml-s`, and `-pml-o` set its thickness,
maximum artificial stretch conductivity in S/m, and polynomial order. The
default thickness is 15% of the smallest domain span and the default maximum
stretch conductivity is ten times the background xx conductivity. Keep the
PML outside the source, receivers, and region of interest, and use a PEC outer
boundary (`-bc 0`) for the usual terminating configuration.

The model VTK files contain cell arrays `PML_region` (0 or 1), `PML_profile`
(normalized grading), `PML_sigma_S_m` (artificial stretch conductivity), and
`cell_size`. In ParaView, color by `PML_region`, apply a Threshold filter with
value 1, and enable surface edges to inspect the PML elements. Rh and Rv remain
the physical formation properties and do not include the artificial PML.
Use at least four elements, preferably six to ten, through the PML thickness;
the `cell_size` array makes an under-resolved absorbing layer easy to identify.

For repeated conductivity evaluations, `-server DIRECTORY` starts a persistent
stdin/stdout protocol used by the Python inversion wrapper. It initializes MPI,
the mesh, and the finite-element space once. Rank 0 prints
`UDAR_SERVER_READY 1`; send
`E REQUEST_ID SIGMA_BACKGROUND SIGMA_TARGET` for each scalar isotropic model.
The response CSV is `DIRECTORY/evaluation_REQUEST_ID.csv`, with a five-digit
zero-padded request ID, and `UDAR_SERVER_DONE REQUEST_ID` signals that the file
is complete. Send `Q` for an orderly shutdown. All ranks must participate in
the same server launch, and the output directory must already exist. The matrix
and its numeric factorization are rebuilt after each conductivity update.

### Windows oneAPI MPI smoke test

Do not start `udar3p.exe` directly from Explorer or an ordinary terminal. It
requires `impi.dll` and the oneAPI compiler/MKL runtimes. The repository
launcher initializes oneAPI and uses Intel MPI's local-root mode so the child
process inherits the runtime path:

```powershell
$env:UDAR_MPI_NP = 2
.\miniapps\udar\run_udar3p_oneapi.cmd `
  -m .\data\inline-tet.mesh -pml -bc 0 -pc ams `
  -freq 48000 -si 0 -ns 2 -out udar_ams.csv
```

The launcher uses `C:\OPT\Intel\oneAPI\setvars.bat` by default. Set
`ONEAPI_ROOT` first if oneAPI is installed elsewhere. Pass `-m` explicitly;
paths are resolved from the directory where the launcher is invoked.

For an interactive oneAPI shell, initialize it from a regular PowerShell
session:

```powershell
cmd /d /s /c 'call "C:\OPT\Intel\oneAPI\setvars.bat" --force && powershell'
```

In the initialized PowerShell, start at the MFEM repository root and run:

```powershell
Push-Location .\build\oneapi-mpi-mumps\miniapps\udar
mpiexec -localroot -genv PATH $env:PATH -n 2 .\udar3p.exe `
  -mumps -freq 2000 -ns 1 -sp 0.08 -out udar3p_smoke.csv
Pop-Location
```

This is the two-rank command used for the Windows smoke test. Entering the
built examples directory is important: the default mesh
`../../data/inline-hex.mesh` is resolved relative to the process working
directory. The output is
`build/oneapi-mpi-mumps/miniapps/udar/udar3p_smoke.csv`. Alternatively, launch from
the repository root and pass the built mesh explicitly:

```powershell
mpiexec -localroot -genv PATH $env:PATH -n 2 .\build\oneapi-mpi-mumps\miniapps\udar\udar3p.exe `
  -m .\build\oneapi-mpi-mumps\data\inline-hex.mesh `
  -mumps -freq 2000 -ns 1 -sp 0.08 -out udar3p_smoke.csv
```

For a physical model, coordinates and spacings must use meters and conductivity
must use S/m. The default `-src 0` is a COMSOL-style ideal magnetic point dipole
implemented as a weak load, so the transmitter must lie inside a mesh element.
Use `-src 1 -loop-r R -wire-r A` for a smooth circular loop with radius `R` and
Gaussian wire radius `A`; several elements must resolve `A`. Field accuracy
should be checked with local refinement around every transmitter and receiver.
The outer boundary should be several skin depths from the tool and target. The
first-order absorbing condition (`-bc 1`) is useful for prototyping; production
models should be checked against domain enlargement and eventually a PML.

The slab and box models use sharp interfaces by default. Set `-iw WIDTH` to use
a smooth hyperbolic-tangent transition when developing inversion of interface
position, slab thickness, or dip. Choose `WIDTH` relative to the local element
size, and verify that both the mesh and the transition width are converged; a
positive width changes the modeled formation rather than resolving a sharp
interface exactly.

## Plot curves

```sh
python miniapps/udar/plot_udar_tensor.py udar_tensor.csv \
  --prefix udar_run
```

This creates `udar_run_amplitude.png` and `udar_run_phase.png`, each with a 3x3
panel of tensor curves versus measured depth. Curves are grouped by frequency
and spacing in a style suitable for comparing shallow and deep UDAR responses.

These plots show absolute transfer-tensor amplitude and phase. Proprietary
GeoSphere-style calibrated attenuation, phase-shift, and directional products
require exact antenna geometry, receiver calibration, bucking/ratio rules, and
measurement definitions; those should be added as a separate instrument model
after the forward fields are benchmarked.
