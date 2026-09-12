# Paper-method implementation map

This file tracks how the focused MFEM `csem` miniapp maps to the three Weiss
papers in `doc/paper`. It separates implemented numerical methods from the
survey-scale features that still belong in later increments.

## Stage 1: high-order forward formulation — implemented

- Total electric-field frequency-domain Maxwell equation with the
  `exp(+i omega t)` convention.
- Conductivity and displacement-current terms, spatially variable
  conductivity/permittivity, and homogeneous tangential electric boundary
  conditions.
- Curl-conforming Nedelec elements on hexahedra with Gauss-Lobatto and
  integrated-GLL basis families.
- Configurable polynomial order and imported physical-coordinate hexahedral
  meshes; MFEM handles local/global edge orientation.
- Finite-length grounded-wire source with configurable center, azimuth, dip,
  length, current, and normalized Gaussian regularization that is stable on
  coarse demonstration meshes.
- Connected polyline transmitters from physical-coordinate vertex CSV files;
  the current follows the listed segment order.
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

## Stage 3: inversion — implemented reduced problem

- Natural-log conductivity parameterization with positive bounds.
- Analytic discrete tangent solves and an explicit real-valued Jacobian.
- Multifrequency Jacobian stacking.
- Error-weighted, Tikhonov-regularized Gauss-Newton objective.
- Matrix-free normal-operator products inside conjugate gradients.
- Damping, bound projection, backtracking, and convergence checks.
- Synthetic-data creation, optional noise, observed-data CSV input, iteration
  history, recovered responses, and recovered-model VTK output.
- Centered finite-difference audit of every analytic Jacobian entry.

## Later survey-scale increments

The executable is intentionally small enough to remain an MFEM miniapp. These
extensions are not claimed by the current code and are the natural next steps:

1. Exact edge-segment delta sources (the implemented straight/polyline source
   uses Gaussian cross-section and end regularization).
2. Many transmitter stations and distribution of independent frequencies and
   sources across MPI communicator groups.
3. Cellwise or regionwise high-dimensional conductivity parameterization,
   sparse regularization operators, and matrix-free J/J-transpose actions that
   avoid storing a survey-scale dense Jacobian.
4. Topography/material models read from tagged production meshes, absorbing
   layers where desired, and automated graded-mesh generation satisfying the
   five-skin-depth boundary criterion.
5. Semi-analytical layered-earth regression data and systematic h/p convergence
   studies matching the numerical examples in the spectral-element paper.

The current reduced inverse has two parameters (background and target), one
straight transmitter, and one configurable electric projection per Gaussian
receiver. Those reductions do not change the implemented forward operator,
PRESB factorization, tangent equations, or Gauss-Newton algebra; they keep
verification runs practical.
