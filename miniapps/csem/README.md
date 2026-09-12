# CSEM spectral-element forward solver and inversion

`csem` is a distributed-memory MFEM miniapp implementing the numerical chain
described in the three Weiss papers under `doc/paper`:

- the frequency-domain total-electric-field equation with the `exp(+i omega t)`
  convention, displacement current, selectable PEC or first-order absorbing
  boundaries, and a finite-length grounded-wire source with configurable
  position, azimuth, dip, length, radius, and current;
- high-order curl-conforming Nedelec elements on hexahedra, using
  Gauss-Lobatto transverse interpolation and the paper's Gauss-point
  tangential basis family;
- the real two-by-two block formulation and the PRESB inner-outer algorithm;
- multifrequency survey assembly, arbitrarily oriented electric-field
  receiver responses on configurable survey lines, and magnetic-field
  recovery from Faraday's law;
- matrix-free directional tangent and adjoint sensitivities, a streamed
  centered finite-difference audit, and a bounded, regularized Gauss-Newton
  inversion in log horizontal/vertical resistivity and dip.

The implementation normalizes mesh coordinates by the longest physical domain
dimension `L`. The assembled coefficients contain the corresponding `L^2`
factors, so the scaled problem is algebraically equivalent to the physical
equation while avoiding the most severe unit-driven conditioning.

## Mathematical formulation

After multiplication by free-space permeability and coordinate scaling, the
discrete complex system is

```text
(B + i M_sigma) E = b,
B = K - M_epsilon,
b = -i omega mu0 J_s.
```

Writing `E = E_R + i E_I` and substituting `-E_I` gives the Weiss system

```text
[ M_sigma  -B ] [ E_R  ] = [ Im(b) ]
[ B         M ] [-E_I  ]   [ Re(b) ].
```

For `H = M_sigma + B`, every PRESB application follows Algorithm 1 in the
solver paper:

```text
H g = f1 + f2
q   = f1 - M_sigma g
H h = q
w1  = g + h,  w2 = -h.
```

Both `H` systems are solved by FGMRES preconditioned with hypre AMS. For
high-order spaces, AMS acts on MFEM's spectrally equivalent low-order-refined
operator. The outer method is FGMRES because the inexact inner solves make the
preconditioner variable. `--solver presb-mumps` retains the outer PRESB method
but solves its two inner systems directly, matching the second iterative
configuration evaluated in the papers. `--solver mumps` selects the fully
direct comparison path when MFEM was built with MUMPS.

The inversion parameters are

```text
m = [log(Rh_0), log(Rv_0), dip_0, ..., log(Rh_N), log(Rv_N), dip_N].
```

Each model unit is transversely isotropic. With `dip = 0`, its symmetry axis
is vertical. Positive dip tilts that axis toward +x in the x-z plane:

```text
n = [sin(dip), 0, cos(dip)]
Sigma = sigma_h I + (sigma_v - sigma_h) n n^T,
sigma_h = 1/Rh,  sigma_v = 1/Rv.
```

`-mp two` retains the geometric background/target model. `-mp region` assigns
one triplet to each positive mesh element attribute. `-mp cell` assigns one
triplet to every global volume-element number. True and starting properties
are read with `-mf` and `-imf` using `id,Rh,Rv,dip_degrees` rows. Region IDs
are one-based MFEM attributes; cell IDs and the two-model background/target
IDs are zero-based.

For a model-space direction `v`, the code contracts the conductivity tensor
derivatives locally and performs one directional tangent solve per survey
task:

```text
A delta_E = -(dA[m;v]) E,
J v = R delta_E.
```

No matrix with `number_of_data * number_of_parameters` entries is formed.
The transpose action first maps arbitrary real/imaginary receiver weights
through the electric or magnetic receiver transpose, solves the transposed
real block system, and accumulates all model derivatives in one element and
quadrature traversal:

```text
A_real_block^T z = R^T w,
J^T w = contraction(z, dA/dm, E).
```

For the complex-symmetric CSEM matrix, the transposed real-block solve is
performed with the existing forward solver by conjugating its right-hand side
and solution. Magnetic receiver transpose actions include the curl transpose,
Faraday scaling, and the required real/imaginary sign change.

The Gauss-Newton operator therefore applies

```text
J^T W^2 J + lambda L^T L + damping I
```

inside conjugate gradients without storing either `J` or the normal matrix.
Its tolerance is controlled by `-gnct` and defaults to `1e-6`, consistent with
the accuracy floor of the inversion tangent/adjoint PDE solves. `-gnd` controls
the separate diagonal damping. Bound projection and backtracking make the
compact example robust enough for synthetic recovery and observed-data
studies. Rh/Rv are bounded in log space using `-rmin/-rmax`; dip is bounded
internally in radians using the degree-valued `-dmin/-dmax` options.

### Spatial regularization

The former diagonal reference-model penalty has been replaced by a sparse
spatial graph operator. Model-unit centroids are computed collectively from
the distributed mesh. Each unit is connected to its `-rkn` nearest units, and
duplicate undirected edges are removed. For an edge `(i,j)`, `L` contains
three first-difference rows:

```text
s_h [log(Rh_j)-log(Rh_i)]
s_v [log(Rv_j)-log(Rv_i)]
s_d [dip_j-dip_i]
```

The scales are configured with `-rrhs`, `-rrvs`, and `-rds`. The objective,
gradient, and Gauss-Newton approximation now use

```text
Phi_reg = 0.5 lambda ||L (m-m_ref)||^2
g_reg   = lambda L^T L (m-m_ref)
H_reg x = lambda L^T L x.
```

Levenberg damping remains a separate `damping I` term and is not reported as
regularization. `-rkn 0` disables graph edges. A one-unit model naturally has
zero spatial edges. The inversion history reports data and spatial objectives
separately.

### Tagged production meshes and topography

`-mat materials.csv` switches material evaluation from the demonstration
`z < 0` air split to volume-attribute tags on an imported hexahedral mesh. The
mesh coordinates are used exactly as imported before the common normalization
by `L`; no flat-surface test is applied. Air, seawater, earth, infrastructure,
and terrain-following interfaces can therefore be represented by conforming
mesh regions with arbitrary topography.

The material table schema is

```text
attribute,Rh_ohm_m,Rv_ohm_m,dip_degrees,relative_permittivity,invertible
1,100000000,100000000,0,1,0
2,1000,1500,15,10,1
```

Every volume attribute present in the mesh must have exactly one row.
`invertible=0` keeps the tagged material fixed in forward, tangent, adjoint,
and inversion calculations. `invertible=1` makes the material active. With
`-mp region`, each active attribute produces one Rh/Rv/dip triplet. With
`-mp cell`, only cells belonging to active attributes produce triplets; fixed
air or seawater cells are omitted from the model vector. Cell model IDs remain
the original global mesh element numbers even though the internal vector is
compact. Tagged materials require `-mp region` or `-mp cell`.

The material table initializes both the true and starting active model. Use
`-mf` and `-imf` to override the true and starting Rh/Rv/dip values. The
permittivity column always follows the material tag and is not inverted.

The included `sample_topography.mesh` is a two-layer terrain-following mesh:
attribute 1 is fixed air, attribute 2 is invertible earth, and boundary tags
1–6 identify the exterior sides. Its companion files are
`sample_topography_materials.csv` and `sample_topography_initial.csv`.

### Optional first-order absorbing boundary

The default `-bc pec` retains homogeneous tangential electric conditions on
all exterior boundaries. `-bc absorbing` applies a first-order impedance
condition to the attributes selected by `-ab`. An empty `-ab` selection means
all boundary attributes. Boundary attributes not selected by `-ab` remain
PEC, allowing, for example, absorbing side/top boundaries with a PEC bottom.

For the `exp(+i omega t)` convention, the weak boundary coefficient is

```text
epsilon_c = epsilon_r - i sigma/(omega epsilon_0)
k_ext     = k_0 sqrt(mu_r epsilon_c)
gamma     = i L k_ext / mu_r.
```

The complex tangential boundary mass `gamma <E_t,F_t>` is included in the
forward, tangent, and adjoint systems. The PRESB shifted auxiliary operator
also receives the corresponding positive boundary shift. Exterior properties
are set with `-abcs`, `-abce`, and `-abcm`.

Example with absorbing side and air-top tags 1–5 while retaining a PEC bottom
tag 6:

```sh
mpiexec -n 2 build/oneapi-mpi-mumps/miniapps/csem/csem \
  -m miniapps/csem/sample_topography.mesh \
  -mat miniapps/csem/sample_topography_materials.csv \
  -mp region -imf miniapps/csem/sample_topography_initial.csv \
  -bc absorbing -ab "1,2,3,4,5" -abcs 0 -abce 1 -abcm 1 \
  -vis -p csem_topography
```

This boundary is a local first-order truncation, not a perfectly matched
layer. Boundaries should still be placed sufficiently far from sources,
receivers, and targets, especially for diffusive low-frequency fields.

For the `exp(+i omega t)` convention, the magnetic field written to VTK is

```text
H_R = -curl(E_I)/(omega mu),
H_I =  curl(E_R)/(omega mu).
```

The coordinate normalization is included in this recovery, so magnetic values
are in A/m and electric values are in V/m.

## Build and run

Configure MFEM with MPI and hypre, then build the target:

```sh
cmake --build build/oneapi-mpi-mumps --target csem
```

Forward response with the paper-style PRESB tolerances:

```sh
mpiexec -n 2 build/oneapi-mpi-mumps/miniapps/csem/csem -no-vis
```

Multifrequency response and one ParaView collection per frequency:

```sh
mpiexec -n 2 build/oneapi-mpi-mumps/miniapps/csem/csem \
  -fs "1,10,100" -vis -p csem_multifrequency
```

Rotate the grounded wire by 30 degrees and measure the horizontal electric
field along a perpendicular receiver line:

```sh
mpiexec -n 2 build/oneapi-mpi-mumps/miniapps/csem/csem \
  -saz 30 -sdip 0 -raz 120 -rdip 0 -rlaz 120 -nr 11 -no-vis
```

For irregular surveys, `-rf receivers.csv` replaces the generated line. Each
non-header row contains either a location that uses the global receiver axis,
or a location followed by its own electric measurement-axis angles:

```text
x_m,y_m,z_m
-1200,650,0

x_m,y_m,z_m,azimuth_deg,dip_deg
300,900,25,90,0
```

`-rc E` (the default) projects the electric field onto each receiver axis.
`-rc H` applies the same receiver geometry to the magnetic field recovered by
Faraday's law. Both modes are used consistently for forward responses,
analytic tangent sensitivities, Jacobian checks, and inversion.

The default `-rm gaussian` evaluates a normalized volume receiver functional,
which has a well-defined transpose and is therefore used by inversion.
`-rm point` evaluates the finite-element field directly at each receiver and
is provided for forward comparisons with the point samples in the paper.
Point receivers are intentionally rejected for inversion and adjoint checks.

`-sf transmitter.csv` replaces the generated straight wire with a connected
polyline. Each non-header row is one physical-coordinate vertex and current
flows in row order; the global `--source-current` and `--source-radius` apply
to every segment:

```text
x_m,y_m,z_m
-500,0,0
0,250,0
500,0,0
```

Ready-to-run receiver and transmitter-catalog examples are included as
`sample_receivers.csv` and `sample_transmitters.csv`.

For multiple independent straight transmitter stations, use
`-sc transmitters.csv`. Each station is solved separately at every frequency
and stacked into the forward vector and Jacobian in frequency/source/receiver
order:

```text
x_m,y_m,z_m,length_m,azimuth_deg,dip_deg,current_A,radius_m
-600,0,0,800,0,0,1,350
600,0,0,800,180,0,1,350
```

The final radius column is optional and defaults to `--source-radius`. The
included `sample_transmitters.csv` is ready to run. `--source-file` and
`--source-catalog` are mutually exclusive because a source file represents one
connected polyline transmitter while a catalog represents separate stations.

Two source discretizations are available. `-sm gaussian` uses the regularized
volume current and is robust on deliberately coarse demonstration meshes.
`-sm line` assembles the wire distribution directly as weighted point deltas
from Gauss quadrature on every segment. Increase `-sqo` when a receiver is
close to the transmitter or the finite-element order is high.

Synthetic inversion:

```sh
mpiexec -n 2 build/oneapi-mpi-mumps/miniapps/csem/csem -inv -ni 5
```

Audit every matrix-free basis-direction tangent action using centered finite
differences, followed by an independent `<Jv,w> = <v,J^T w>` dot test:

```sh
mpiexec -n 2 build/oneapi-mpi-mumps/miniapps/csem/csem \
  -fs "3,10,30" -vj -jfd 1e-4 -p csem_jacobian_check
```

Direct-solver comparison, if MUMPS is enabled:

```sh
mpiexec -n 2 build/oneapi-mpi-mumps/miniapps/csem/csem -s mumps -inv
```

PRESB with direct MUMPS inner solves:

```sh
mpiexec -n 2 build/oneapi-mpi-mumps/miniapps/csem/csem -s presb-mumps
```

The forward run writes `csem_response.csv`. Synthetic inversion writes the
generated observations, the per-iteration history, recovered responses, and
one ParaView field collection per frequency. Every collection contains
`E_real`, `E_imaginary`, `H_real`, `H_imaginary`, `Rh_ohm_m`, `Rv_ohm_m`,
`sigma_h_S_per_m`, `sigma_v_S_per_m`, `dip_degrees`, and
`relative_permittivity`. `material_attribute` records the original volume tag
and `invertible_material` is 1 on active cells and 0 on fixed cells. The mesh
coordinates in these files are normalized by the printed length scale `L`;
multiply displayed coordinates by `L` to recover metres. Field values retain
physical units.

## Layered-earth reference and convergence tests

`paper_scenario2_reference.csv` contains independent 1-D responses for the
paper's Scenario 2 at 1000 Hz. The model is a 10,000 ohm-m half-space with a
100 ohm-m layer from 500 to 1000 m depth and 100,000,000 ohm-m air. A 200 m
x-directed grounded cable carries 0.5 A, giving the stated 100 A-m source
moment. Ex, Hy, and Hz are sampled from 1 to 5 km in 250 m cross-line steps.
The complex responses were generated with `empymod` 2.6.0 using finite-bipole
Gauss integration and its QWE Hankel-transform method. Routine regressions use
the committed CSV and have no Python-package dependency beyond the standard
library.

`csem_convergence.py` constructs the 40 x 40 x 40 km tagged hexahedral model,
places the two 100 m source segments on mesh edges, runs p- and h-refinement
cases, and writes every response, solver log, and a combined
`scenario2_convergence.csv`. The summary contains mean, RMS, and maximum
relative errors for amplitude and for the real and imaginary parts, following
equation (25) of the paper. With `--suite paper`, it also checks the expected
14 x 14 x 14 topology and the Table 9 H(curl) degree-of-freedom counts. The
paper's printed p=2 count (70,466) differs from the tensor-product conformity
count (70,644); the harness documents and accepts both, while MFEM produces
70,644.

The practical CTest regression covers Ex, Hy, and Hz. It uses p=1,2 for the p
trend and one global p=1 mesh refinement for the h trend in one invocation,
reusing their common baseline:

```sh
ctest --test-dir build/oneapi-mpi-mumps -R csem_layered --output-on-failure
```

They are registered when testing, MPI, MUMPS, and a Python interpreter are
available. MUMPS is used to keep solver error out of the discretization test.
The full paper-order study is explicit because Table 9 reports 1,058,610
unknowns and roughly 45 GB of factorization memory at p=5:

```sh
python miniapps/csem/csem_convergence.py \
  --executable build/oneapi-mpi-mumps/miniapps/csem/csem \
  --suite paper --orders 1,2,3,4,5 --h-levels 0,1 --h-order 2 \
  --components Ex,Hy,Hz --solver mumps \
  --output-dir csem-paper-convergence
```

Additional global refinements can be requested with, for example,
`--h-levels 0,1,2`; their memory grows rapidly. `--suite quick` substitutes a
smaller mesh for development runs while retaining the same physical layered
model and reference response. Use `--report-only` to collect curves without
applying regression gates.

Maintainers can reproduce the committed table in an environment containing
`empymod==2.6.0` and NumPy:

```sh
python miniapps/csem/generate_layered_reference.py
```

The observed-data CSV schema is:

```text
frequency_hz,source,receiver,x_m,y_m,z_m,direction_x,direction_y,direction_z,real,imag,stddev_real,stddev_imag
```

Receiver coordinates, line azimuth, and measurement-axis azimuth/dip are
controlled by command-line geometry options or `--receiver-file`. The three
direction columns are the unit vector onto which the electric field is
projected, or the magnetic field when `--receiver-field H` is selected. An
observed CSV must have the same field choice and frequency/source/receiver row
count and ordering. For compatibility, the earlier twelve-column one-source,
nine-column multifrequency, and eight-column single-frequency CSV formats are
also accepted.

With multiple sources, visualization directories include `_s000`, `_s001`,
and so on before the frequency tag. This prevents one transmitter's field
collection from overwriting another's.

Independent frequency/source tasks can be distributed across MPI subgroups
with `-sg N`. Ranks are divided as evenly as possible into `N` communicators;
task number `frequency_index * source_count + source_index` is assigned
round-robin to a group. Each group performs its own distributed H(curl) solve,
then one group root contributes its receiver block to a world-communicator
sum. This preserves exactly the serial frequency/source/receiver ordering for
forward vectors and Jacobians. `N` must not exceed the number of MPI ranks.

The default forward tolerances are `1e-8` outer and `1e-2` inner. The default
inversion sensitivities use `1e-4` outer and `1e-1` inner, following the
performance/accuracy choice reported by Weiss, Rochlitz, and Guenther. Use
`-h` to list all options. At startup the executable also prints the exact
lossy-medium plane-wave skin depth and the paper's recommended five-skin-depth
boundary distance. The lightweight default mesh is deliberately smaller; use
`--mesh` with a graded physical-coordinate hexahedral mesh for publication
accuracy and place every truncation boundary roughly five skin depths from the
source, receivers, and target at the lowest simulated frequency.

## Scope

This is a focused, executable miniapp rather than a survey-scale replacement
for custEM/pyGIMLi. It supports smooth Gaussian and line-delta approximations
to straight or polyline grounded wires. Gaussian receiver functionals keep
inversion and adjoint evaluation stable on coarse demonstration meshes, while
point sampling supports reference forward studies. The generated source and electric
receiver axes can have arbitrary azimuth and dip. Receivers may be generated on
an arbitrarily oriented horizontal line or read at independent 3-D locations
and orientations from a CSV file. The forward operator, PRESB
factorization, multifrequency assembly, matrix-free tangent/adjoint
sensitivities, and Gauss-Newton objective support multiple transmitters and
geometric, mesh-region, or cellwise anisotropic resistivity parameterizations.

The current matrix-free implementation prioritizes parameter-space memory:
storage is `O(number_of_data + number_of_parameters)` rather than
`O(number_of_data * number_of_parameters)`. Each operator action currently
reassembles and resolves its assigned frequency/source systems. Reusing cached
linearizations or checkpointed forward fields is a future performance
optimization and does not change the tangent/adjoint algebra.
