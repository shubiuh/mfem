# Paper-method implementation map

This file tracks how the focused MFEM `csem` miniapp maps to the three Weiss
papers in `doc/paper`. It separates implemented numerical methods from the
survey-scale features that still belong in later increments.

## Stage 1: high-order forward formulation — implemented

- Total electric-field frequency-domain Maxwell equation with the
  `exp(+i omega t)` convention.
- Conductivity and displacement-current terms with spatially variable
  conductivity/permittivity.
- Homogeneous tangential electric boundaries or a selectable first-order
  complex impedance boundary on any subset of mesh boundary attributes.
- Passive lossy-exterior impedance from conductivity, relative permittivity,
  and relative permeability, with the matching boundary term in the PRESB
  shifted auxiliary operator.
- Curl-conforming Nedelec elements on hexahedra with the paper's
  Gauss-Lobatto closed and Gauss-Legendre open basis-point families.
- Configurable polynomial order and imported physical-coordinate hexahedral
  meshes; MFEM handles local/global edge orientation.
- Attribute-driven production material tables with fixed and invertible
  regions. Imported terrain-following geometry supplies topography directly,
  so no global flat air/earth test is used for tagged meshes.
- Compact cellwise parameter vectors exclude fixed tagged cells while retaining
  original global cell IDs in model input and output.
- Finite-length grounded-wire source with configurable center, azimuth, dip,
  length, current, and normalized Gaussian regularization that is stable on
  coarse demonstration meshes.
- Connected polyline transmitters from physical-coordinate vertex CSV files;
  the current follows the listed segment order.
- Selectable Gaussian-volume or quadrature-assembled line-delta source
  discretization on every straight/polyline segment.
- Multiple independent straight transmitter stations from a catalog, stacked
  consistently across forward responses, tangent Jacobians, inversion data,
  and source-tagged VTK collections.
- Independent frequency/source tasks distributed round-robin over balanced
  MPI communicator groups, followed by deterministic world-communicator
  assembly of forward vectors and Jacobians.
- Arbitrarily oriented electric receiver axes on configurable horizontal
  survey lines or at independent 3-D locations read from a geometry CSV;
  receiver orientation is recorded in the response CSV.
- Multifrequency survey assembly and frequency-tagged response CSV output.
- Magnetic-field recovery from Faraday's law and ParaView output of real and
  imaginary E/H vectors plus the constitutive-property fields.
- Electric or magnetic receiver projections use the same arbitrary 3-D
  positions and axes in forward modeling, tangent sensitivities, Jacobian
  verification, and inversion.
- Exact lossy-medium plane-wave skin-depth reporting with the paper's
  five-skin-depth boundary guideline.

## Stage 4: layered-earth verification and convergence — implemented

- Independent complex Scenario 2 layered-earth responses at 1000 Hz for Ex,
  Hy, and Hz, sampled at the paper's 1--5 km cross-line receiver offsets.
- Reproducible reference generator pinned to `empymod` 2.6.0 with a finite
  200 m, 0.5 A bipole and QWE Hankel transforms; routine tests consume the
  committed CSV without requiring `empymod`.
- A 40 x 40 x 40 km, 14 x 14 x 14 tagged mesh generator with two equal cells
  above the conductive layer, three inside it, and two graded basement cells,
  as specified for the paper's coarse Scenario 2 model.
- The 200 m grounded cable is split into the paper's two 100 m, edge-aligned
  electric-dipole segments. Receivers use direct forward point interpolation.
- Systematic configurable p=1,...,5 and global h-refinement runs, with mean,
  RMS, and maximum amplitude/real/imaginary errors and per-case solver logs.
- Regression checks for error reduction, final error bounds, and the paper's
  H(curl) degree-of-freedom sequence, including a documented allowance for
  the apparent p=2 digit transposition in Table 9.
- A practical combined CTest p- and h-regression for all three reported
  components when
  MPI, MUMPS, and Python are available; million-DOF full-paper cases remain an
  explicit opt-in run because of their published memory cost.

## Stage 2: Weiss solver framework — implemented

- Real two-by-two square-block form.
- Literal PRESB Algorithm 1: two solves with the shifted H(curl) operator, one
  conductivity-mass product, and the stated vector combinations.
- Outer FGMRES with inner FGMRES preconditioned by low-order-refined hypre AMS.
- Alternative PRESB configuration with direct MUMPS inner solves
  (`--solver presb-mumps`).
- Fully direct complex/system MUMPS reference (`--solver mumps`).
- Paper-style forward tolerances (`1e-8` outer, `1e-2` inner) and inversion
  sensitivity tolerances (`1e-4` outer, `1e-1` inner).
- True residual reporting and reusable operator/factorization within each
  forward-plus-tangent solve group.

## Stage 3: anisotropic region/cell inversion — implemented

- Geometric two-unit, mesh-attribute regionwise, or global-element cellwise
  parameterization.
- Transversely isotropic conductivity tensor per unit from Rh, Rv, and dip,
  with analytic derivatives with respect to log(Rh), log(Rv), and dip.
- Positive log-resistivity bounds and independent dip-angle bounds.
- Analytic directional tangent solves implementing matrix-free `J v` actions.
- Receiver-transpose adjoint solves and elementwise contraction of all local
  conductivity derivatives implementing matrix-free `J^T w` actions.
- Electric and magnetic receiver adjoints, including the magnetic curl
  transpose, Faraday scaling, and real/imaginary sign convention.
- Multifrequency and multisource tangent/adjoint assembly across survey MPI
  groups without storing a dense survey Jacobian.
- Error-weighted, spatially regularized Gauss-Newton objective relative to a
  reference model.
- Sparse centroid-neighbor graph with separate first-difference rows for
  log(Rh), log(Rv), and dip; the resulting L-transpose-L term is used in the
  objective, gradient, and normal operator.
- Matrix-free `J^T W^2 J + lambda L^T L + damping I` products inside conjugate
  gradients for an arbitrary number of model parameters.
- Damping, bound projection, backtracking, and convergence checks.
- Synthetic-data creation, optional noise, observed-data CSV input, iteration
  history, recovered responses, and recovered-model VTK output.
- Streamed centered finite-difference audit of every basis-direction tangent
  action plus an independent tangent/adjoint dot-product test.

## Later survey-scale increments

The executable is intentionally small enough to remain an MFEM miniapp. These
extensions are not claimed by the current code and are the natural next steps:

1. Reuse or checkpoint frequency/source forward states and assembled solver
   objects across repeated Krylov actions. The current implementation removes
   dense-Jacobian storage but recomputes each assigned forward linearization.
2. Automated graded production-mesh generation satisfying the five-skin-depth
   boundary criterion, plus higher-order absorbing conditions or PML where the
   first-order impedance truncation is insufficient.

The default remains a compact geometric background/target model, but the same
forward operator, PRESB factorization, tangent equations, and Gauss-Newton
algebra now accept three anisotropic parameters for every mesh attribute or
global cell. Both the survey sensitivity and spatial regularization paths are
now matrix-free or sparse in parameter space; production-scale performance can
be improved further through the forward-state reuse listed above.
