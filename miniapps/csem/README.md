# CSEM spectral-element forward solver and inversion

`csem` is a distributed-memory MFEM miniapp implementing the numerical chain
described in the three Weiss papers under `doc/paper`:

- the frequency-domain total-electric-field equation with the `exp(+i omega t)`
  convention, displacement current, homogeneous tangential electric boundary
  conditions, and a finite-length grounded-wire source with configurable
  position, azimuth, dip, length, radius, and current;
- high-order curl-conforming Nedelec elements on hexahedra, using
  Gauss-Lobatto interpolation and integrated-GLL basis families;
- the real two-by-two block formulation and the PRESB inner-outer algorithm;
- multifrequency survey assembly, arbitrarily oriented electric-field
  receiver responses on configurable survey lines, and magnetic-field
  recovery from Faraday's law;
- an explicit sensitivity matrix from discrete tangent solves, a centered
  finite-difference Jacobian audit, and a bounded, regularized Gauss-Newton
  inversion in natural-log conductivity.

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
m = [log(sigma_background), log(sigma_target)].
```

For parameter `m_j`, the code differentiates the assembled system and solves

```text
A dE/dm_j = -(dA/dm_j) E.
```

Receiver functionals applied to these tangent fields form the explicit
Jacobian. A matrix-free normal operator applies

```text
J^T W^2 J + lambda I + damping I
```

inside conjugate gradients. Bound projection and backtracking make the compact
example robust enough for synthetic recovery and small observed-data studies.

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

Ready-to-run versions of both formats are included as
`sample_receivers.csv` and `sample_transmitter.csv`.

Synthetic two-property inversion:

```sh
mpiexec -n 2 build/oneapi-mpi-mumps/miniapps/csem/csem -inv -ni 5
```

Audit the two analytic tangent-sensitivity columns using centered finite
differences:

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
`E_real`, `E_imaginary`, `H_real`, `H_imaginary`, conductivity, log10
conductivity, resistivity, and relative permittivity. The mesh coordinates in
these files are normalized by the printed length scale `L`; multiply displayed
coordinates by `L` to recover metres. Field values retain physical units.

The observed-data CSV schema is:

```text
frequency_hz,receiver,x_m,y_m,z_m,direction_x,direction_y,direction_z,real,imag,stddev_real,stddev_imag
```

Receiver coordinates, line azimuth, and measurement-axis azimuth/dip are
controlled by command-line geometry options or `--receiver-file`. The three
direction columns are the unit vector onto which the electric field is
projected, or the magnetic field when `--receiver-field H` is selected. An
observed CSV must have the same field choice, frequency/receiver row count, and
ordering. For compatibility, the earlier nine-column multifrequency and
eight-column single-frequency CSV formats are also accepted.

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
for custEM/pyGIMLi. It uses a smooth Gaussian approximation to a straight or
polyline grounded wire and Gaussian receiver functionals to keep source and
receiver evaluation stable on coarse demonstration meshes. The generated
source and electric
receiver axes can have arbitrary azimuth and dip. Receivers may be generated on
an arbitrarily oriented horizontal line or read at independent 3-D locations
and orientations from a CSV file. The forward operator, PRESB
factorization, multifrequency assembly, tangent sensitivities, and Gauss-Newton
objective are the same discrete methods that extend to many transmitter
stations and a cellwise conductivity parameterization.
