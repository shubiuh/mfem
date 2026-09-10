# MFEM 3D UDAR forward solver

`udar3p` is a first parallel forward-modeling application for ultra-deep
azimuthal resistivity (UDAR) studies. It uses MFEM Nedelec elements for the
frequency-domain electric field and recovers the full complex 3x3 magnetic
transfer tensor in the tool coordinate system. Receiver locations are placed
behind each transmitter along the negative tool axis by the requested spacing.

## Capabilities

- Multiple frequencies and transmitter-receiver spacings in one run.
- Multiple stations along a straight well trajectory.
- Three orthogonal magnetic-dipole transmitters and receivers (`Hxx` through
  `Hzz`), using an ideal point dipole or a resolved finite-radius loop.
- Full-tensor, diagonal-anisotropic, or isotropic conductivity values.
- Uniform, dipping-slab, and bounded 3D-geobody formation models.
- Optional smooth formation interfaces for differentiable geometry parameters.
- Native complex MUMPS when `zmumps` is available, real-block MUMPS otherwise,
  or iterative flexible GMRES with a distributed ILU preconditioner.
- Tidy CSV results plus 3x3 amplitude and phase curve panels.

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
mpirun -np 4 ./build/examples/udar3p \
  -m data/inline-hex.mesh -rs 1 -rp 1 \
  -freq 2000,6000,12000,48000,96000 \
  -sp 0.08,0.16,0.24 -gmres \
  -out udar_tensor.csv
```

Use `-mumps` for the direct solver. Run `udar3p -h` to see all geometry,
material, trajectory, source, boundary, and iterative-solver controls.
The `-gmres` path uses FGMRES with distributed ILU and defaults to a relative
residual tolerance of `1e-4`; tighten `-rtol` only after checking convergence
for the chosen mesh and material contrast. MUMPS is preferred for reference
solutions.

### Windows oneAPI MPI smoke test

First initialize oneAPI from a regular PowerShell session:

```powershell
cmd /d /s /c 'call "C:\OPT\Intel\oneAPI\setvars.bat" intel64 && powershell'
```

In the initialized PowerShell, start at the MFEM repository root and run:

```powershell
Push-Location .\build\oneapi-mpi-mumps\examples
mpiexec -n 2 .\udar3p.exe `
  -mumps -freq 2000 -ns 1 -sp 0.08 -out udar3p_smoke.csv
Pop-Location
```

This is the two-rank command used for the Windows smoke test. Entering the
built examples directory is important: the default mesh
`../data/inline-hex.mesh` is resolved relative to the process working
directory. The output is
`build/oneapi-mpi-mumps/examples/udar3p_smoke.csv`. Alternatively, launch from
the repository root and pass the built mesh explicitly:

```powershell
mpiexec -n 2 .\build\oneapi-mpi-mumps\examples\udar3p.exe `
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
python examples/udar/plot_udar_tensor.py udar_tensor.csv \
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
