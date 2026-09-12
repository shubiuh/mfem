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
point-location, frequency-system assembly, electromagnetic postprocessing, and
linear-solver implementations live in separate `udar3p_*.hpp/.cpp` modules.
`UDARFrequencySystem` owns the Maxwell form, coefficient lifetimes, assembled
operator, real-block conversion, and solver for one frequency. These objects
persist across inversion requests, allowing MUMPS symbolic state to be reused
while material-dependent operators are rebuilt safely. The miniapp has
dedicated CMake and make build files and keeps the executable name `udar3p` for
compatibility with the inversion wrapper.

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
y/z transmitter fields. They also contain `electric_energy_J_m3_x`,
`magnetic_energy_J_m3_x`, `ohmic_loss_W_m3_x`, and the active
`poynting_W_m2_x` vector for each transmitter axis. Energy uses the configured
isotropic permittivity and permeability. Ohmic loss uses only the physical
formation conductivity and excludes the artificial PML stretch. Open a `.pvd`
file in ParaView; the model collection
locates its per-rank `.vtk` pieces and field collections locate their `.pvtu`
and `.vtu` pieces. Rh is the average xx/yy component of the inverse
conductivity tensor and Rv is its zz component, so these labels use the global
mesh axes.

Every response CSV has a sibling `NAME_timing.csv`. Each frequency row records
the global element and complex-unknown counts, solver and preconditioner modes,
maximum MPI-rank wall times for operator assembly, solver setup, right-hand
side construction, solves, field output, and receiver sampling, plus aggregate
Krylov iterations and ILU recoveries. These measurements make frequency-mesh
and solver choices directly comparable during inversion development.

### Measuring derived-field output before device implementation

Add `-vtk-profile` to a representative forward run with `-vtk DIRECTORY` to
measure host derived-coefficient evaluation inside ParaView `Save()`. The
timing CSV now includes `vtk_save_s`, `derived_eval_s`, and `vtk_profiled`.
`vtk_save_s` includes mesh/field serialization and file writing;
`derived_eval_s` includes host FE interpolation and coefficient algebra for
magnitudes, energy, loss and Poynting fields. It is a subset of Save time,
which is itself within `field_output_s`; do not add these nested measurements
when computing total runtime. Times are maxima across MPI ranks, so different
columns may refer to different ranks. A zero `derived_eval_s` with
`vtk_profiled=0` means unmeasured, not free. Profiling adds two clock reads per
coefficient call; compare with an unprofiled run to assess overhead.

Build the explicit CMake target `udar_output_benchmark` for an output-only
scaling study using the production writer. Example arguments:

```text
-n 32 -o 2 -r 3 -out build/udar-output-n32
```

It constructs manufactured ND fields and curl-derived RT fields for all three
transmitter orientations, alternates raw E/H and complete output, then runs a
separate instrumented pass. Each repeated mode overwrites its own collection;
all timings are retained in `timing.csv`. The current production writer uses
one visualization subdivision (eight vertex samples per hex), independently
of FE order. This benchmark estimates output cost at large FE sizes but does
not establish its share of a production Maxwell solve. See
[OUTPUT_PERFORMANCE.md](OUTPUT_PERFORMANCE.md) for measured results and the
remaining device-validation gate. The production output backend remains host.

### Partial assembly

`-assembly auto` (the default) selects partial assembly only for constant-order
hex meshes with Nedelec order at least 2, complex FGMRES with AMS
(`-gmres -pc ams -ams-mode complex`) with AMS cycle 13 (`-ams-cycle 13`),
disabled ILU fallback
(`-no-ams-fallback`), PEC boundaries (`-bc 0`), and no PML. Locally refined
nonconforming hex meshes are supported. `-assembly full` always uses the
assembled path. `-assembly partial` requests PA but obeys the same guards;
an incompatible request uses full assembly and prints the reason. Solver,
boundary and fallback choices are never changed implicitly.

Only the Maxwell operator uses partial assembly. The positive shifted AMS
auxiliary operator remains fully assembled, so this mode does not remove all
high-order matrix storage. MUMPS, ILU, explicit real-block mode, order-1,
tetrahedral/mixed meshes, variable order, PML and absorbing-boundary runs retain
full assembly. PA reduces Maxwell matrix storage; it does not guarantee faster
solves or improve AMS convergence. With ILU fallback disabled, a failed AMS
solve reports nonconvergence as usual.

This initial PA gate admits only the tested AMS cycle 13. Cycle 14 reached
the iteration limit in a higher-order forward check using the assembled
operator as well; its high-order convergence needs separate investigation.
An example selection, appended to the normal survey arguments, is:

```text
-o 2 -gmres -pc ams -ams-mode complex -ams-cycle 13 -no-ams-fallback -bc 0 -no-pml -assembly auto
```

The startup log reports the actual selection for every frequency, and each
timing CSV includes an `assembly_mode` column. Material updates rebuild PA
quadrature data and the AMS auxiliary matrix. The explicit CMake target
`udar_assembly_test` compares full and partial operator actions on locally
refined order-2 and order-3 hex meshes, including anisotropic conductivity,
PEC rows and repeated material updates. Run MPI workers with the initialized
oneAPI environment; on Windows pass `-genv PATH "%PATH%"` to `mpiexec`.

### Mesh hierarchy and PML

The forward solver retains a `ParFiniteElementSpaceHierarchy` throughout the
survey and persistent inversion requests. Serial uniform refinements (`-rs`)
define its coarse mesh. The mesh is then partitioned once; each tool-region
refinement pass (`-octree-l`) and each parallel uniform refinement (`-rp`)
adds a separate mesh and Nedelec space. Refinement copies the parent mesh, so
coarse levels remain intact. No repartitioning is performed between levels.
Both hex and tet tool refinement use MFEM nonconforming meshes.

Each adjacent pair retains an MFEM `TrueTransferOperator` with hanging-node
constraints and edge orientations handled by the finite element spaces.
`Mult` prolongs coarse true DOFs and `MultTranspose` restricts residuals.
`GetFESpaceAtLevel` exposes each level for auxiliary Maxwell assembly and
level-specific essential-boundary DOFs. The existing forward solver uses
`GetFinestFESpace`; retaining the hierarchy does not yet select a geometric
multigrid cycle. Level operators, smoothers and the coarse solve are the next
implementation stage. Keeping these meshes, spaces and transfers increases
memory usage, and the inherited coarse partition may need load-balancing work
before large distributed local-refinement runs.

Build the explicit CMake target `udar_hierarchy_test` to check hex and tet
hierarchies with two local refinement passes followed by a uniform pass.
Run it in the initialized oneAPI environment, directly or through MPI. It
checks coarse-mesh preservation, interpolation of a representable field with
nonzero curl, and the prolongation/restriction transpose identity.

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

When PML is enabled, the miniapp analyzes all six Cartesian sides before the
first solve. For each side it reports the number of intersected elements, the
largest element span normal to the boundary, a conservative layer count
`floor(PML thickness / maximum normal span)`, and the largest normalized
grading change contained in one element. A layer count below four or a grading
jump above 0.35 produces a warning. The conservative definition avoids
overstating resolution on tetrahedral, nonuniform, and locally refined meshes.

The complete report is written to `pml_resolution.csv` in the ParaView or
persistent-server directory. Without either directory it is written beside the
response as `NAME_pml_resolution.csv`. Its `max_conductivity_jump_S_m` column
converts the normalized grading jump to the corresponding artificial stretch-
conductivity change. Compare this report between frequency meshes together
with receiver-response convergence; the diagnostic measures mesh resolution
of the chosen PML profile and does not by itself establish boundary accuracy.

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
