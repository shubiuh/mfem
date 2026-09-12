# Palace-informed UDAR upgrade plan

This review uses the local Palace checkout at commit `a0655fbf` (after the
0.18.0 release, 2026-09-09). Palace is Apache-2.0 while MFEM is BSD-3-Clause.
The UDAR changes therefore adopt algorithms and public MFEM/HYPRE interfaces;
they do not copy Palace implementation files into MFEM.

## Relevant Palace evolution

| Palace release/change | Useful lesson for UDAR | UDAR action |
|---|---|---|
| 0.6.0: complex Hypre operators and h/p geometric multigrid with AMS | Keep complex algebra separate from the real auxiliary preconditioner and preserve mesh hierarchies | Add a shifted H(curl) preconditioner, then an optional geometric hierarchy |
| 0.9.0: SuperLU_DIST/STRUMPACK, solver controls, CGS2 | Direct solvers need selectable ordering and symbolic reuse; iterative solvers need robust orthogonalization | Expose MUMPS ordering/reuse first, then add solver backends only when their ABI is validated |
| 0.10.0: real/imaginary phasor and energy output | Complex output needs explicit phase and derived quantities | Retain E/H real, imaginary, and magnitude fields; add energy and Poynting fields later |
| 0.12.0: complex operators independent of PETSc, libCEED partial assembly, flux error estimates, AMR | PETSc should remain optional; high-order/GPU execution and goal-oriented refinement require operator-level separation | Refactor assembly behind a forward-operator class before adding partial assembly and error-driven AMR |
| 0.13.0: partial assembly by default at high order and device-aware execution | Full sparse assembly suits the present order-1 direct solve; high-order forward meshes need a matrix-free path | Select full or partial assembly from order and solver compatibility |
| 0.14.0: complex coarse solve and consolidated postprocessing | A real approximation is cheaper, but complex coarse solves can be more robust in conductive media | Support both preconditioner modes and record the selected mode |
| 0.16.0: tunable AMS iterations | One fixed AMS setup is not robust across order, conductivity, and frequency | Default AMS work to finite-element order and expose cycle/smoothing controls |
| 0.17.0: cached GSLIB point locator and pre-partition region refinement | Survey points and prescribed refinements should be prepared once per mesh | Cache receiver locations across frequencies and inversion evaluations |
| 0.18.0: batched reductions, device ParaView fields, operator/preconditioner reuse | Reduce collectives, avoid unused output work, and reuse setup while sparsity is unchanged | Batch receiver sampling, skip unused H projection, then reuse solver analysis |
| commits `75b4f6d8`, `a0655fbf`: orthogonalization and reduced weighted-CGS/operator costs | Global reductions dominate at scale | Add reduction and iteration counters before selecting a low-reduction Krylov method |

## Implementation order

1. **Cached, batched survey sampling (implemented).** Locate all receivers once per station,
   evaluate all three transmitter fields in one pass, and use one MPI reduction.
   Do not project H into an RT space unless field visualization is requested.
2. **Solver reuse and controls (symbolic reuse and timings implemented).**
   Reuse MUMPS symbolic analysis while the matrix pattern is unchanged and
   record assembly, setup, solve, field-output, sampling, Krylov-iteration,
   and fallback metrics. Next expose ordering, BLR, and out-of-core controls.
3. **H(curl) preconditioning (implemented).** Assemble a positive shifted
   auxiliary operator and configure Hypre AMS. Complex-operator and explicit
   real-block modes are selectable, with lazy distributed-ILU recovery. Final
   multi-rank performance validation remains part of the Windows MPI checks.
4. **Mesh hierarchy and AMR (PML diagnostics and hierarchy implemented).** Report a
   conservative element-layer count, maximum normal element span, and maximum
   polynomial-grading jump independently on all six Cartesian PML sides.
   Retain parallel local and uniform refinement levels, Nedelec spaces and
   true-DOF prolongation/transpose restriction for geometric multigrid. Next
   add level-wise auxiliary Maxwell operators, smoothers and a coarse solver,
   then use curl/flux recovery plus Dörfler marking for response-driven refinement. Frequency
   meshes remain fixed during each inversion objective/gradient evaluation.
5. **Forward-operator refactor (frequency systems implemented).**
   `UDARFrequencySystem` now owns frequency-dependent coefficients, Maxwell
   assembly, real-block conversion, and direct or iterative solver state.
   Mesh/material configuration, sources, receiver evaluation, and output have
   explicit remaining boundaries. Guarded native partial assembly is now
   available for order >= 2 hex meshes with complex AMS cycle 13, PEC, no PML and no
   ILU fallback. The AMS auxiliary matrix remains assembled. Validate GPU
   execution and expand boundary/element support before extending this gate.
6. **Postprocessing (field diagnostics implemented).** Output electric and
   magnetic energy density, physical ohmic loss density, active Poynting
   vector, mesh size, and solver timing metadata. Next add integrated power
   balances. Output profiling and a production-writer scaling benchmark now
   isolate host derived-field evaluation. Introduce device output after a
   representative survey confirms a material end-to-end share, and validate
   transfer-inclusive speedup and numerical agreement on an accelerator build;
   the current oneAPI build is CPU-only. See OUTPUT_PERFORMANCE.md.
7. **Multi-frequency reduction.** Introduce a projection-based adaptive
   frequency sweep for line search and early inversion iterations only after
   full-order solves and derivatives are stable, with periodic validation.

Every performance change should preserve the complex 3x3 receiver tensor to a
documented tolerance. Boundary/PML convergence and mesh convergence remain the
accuracy gates; a faster unconverged model is unsuitable for inversion.
