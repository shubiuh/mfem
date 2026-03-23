// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

//                    MFEM Example: Octree EM (ex_octree_em)
//
// Compile with: make ex_octree_em
//
// Sample runs:
//   ex_octree_em
//   ex_octree_em -o 1 -r 2 -amr 4
//   ex_octree_em -o 2 -r 1 -amr 6 -e 1e-4
//   ex_octree_em -m ../data/inline-hex.mesh -o 1 -r 2 -amr 4
//   ex_octree_em -no-vis
//
// Description:
//   This example solves the definite Maxwell (curl-curl) equation
//
//       curl ( mu^{-1} curl E ) + sigma E = f
//
//   on an all-hexahedral domain using Nedelec (H(curl)) finite elements
//   of order >= 1.  The mesh is adaptively refined using an octree-style
//   nonconforming (NC) refinement strategy: hanging nodes are strictly
//   limited to at most ONE level (the 2:1 balance rule of standard octree
//   meshes).
//
//   The NC DOF constraints arising at the hanging faces/edges are enforced
//   automatically by MFEM through the conforming prolongation matrix P
//   obtained from FiniteElementSpace::GetConformingProlongation().
//   FormLinearSystem() applies these constraints transparently.
//
//   An element-wise L2-norm error estimator drives the refinement; the
//   ThresholdRefiner enforces nc_limit = 1.  After each AMR cycle the
//   GridFunction is updated (prolongated) onto the new mesh and the
//   problem is re-solved.
//
//   Exact solution (for verification):
//       E(x,y,z) = ( sin(kappa*y), sin(kappa*z), sin(kappa*x) )
//   with kappa = freq * pi.
//
//   The corresponding right-hand side is:
//       f = (1 + kappa^2) * E   (using mu^{-1}=1, sigma=1).
//
//   We recommend viewing Examples 3 and 15 before this example.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// ── Global parameters ────────────────────────────────────────────────────────
real_t freq  = 1.0;   // oscillation frequency for the exact solution
real_t kappa;          // = freq * pi

// ── Forward declarations ─────────────────────────────────────────────────────
void E_exact(const Vector &x, Vector &E);
void f_exact(const Vector &x, Vector &f);

// ── Helper: rebuild FES, reassemble, solve, return L2 error ─────────────────
// (used inside the AMR loop; returns error norm for reporting)
static real_t SolveOnCurrentMesh(Mesh                 &mesh,
                                  FiniteElementSpace   &fes,
                                  GridFunction         &x,
                                  bool                  visualization,
                                  socketstream         *sout,
                                  int                   amr_it);

// =============================================================================
int main(int argc, char *argv[])
{
   // ── 1. Parse command-line options ─────────────────────────────────────────
   const char *mesh_file     = "";   // empty → build Cartesian hex mesh
   int         order         = 1;
   int         ref_levels    = 1;   // initial uniform refinements
   int         amr_iter      = 4;   // number of AMR cycles
   real_t      max_err       = 5e-4; // local error threshold for refinement
   real_t      hysteresis    = 0.25; // derefinement safety factor
   bool        visualization = true;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file,  "-m",    "--mesh",
                  "Hex mesh file (empty = built-in 4x4x4 Cartesian hex).");
   args.AddOption(&order,      "-o",    "--order",
                  "Nedelec FE order (>= 1).");
   args.AddOption(&ref_levels, "-r",    "--ref-levels",
                  "Number of initial uniform refinements.");
   args.AddOption(&amr_iter,   "-amr",  "--amr-iter",
                  "Number of AMR refinement cycles.");
   args.AddOption(&max_err,    "-e",    "--max-err",
                  "Maximum element error for AMR threshold.");
   args.AddOption(&hysteresis, "-y",    "--hysteresis",
                  "Derefinement safety coefficient (fraction of max_err).");
   args.AddOption(&freq,       "-f",    "--frequency",
                  "Frequency parameter for the exact solution.");
   args.AddOption(&visualization, "-vis", "--visualization",
                  "-no-vis", "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   kappa = freq * M_PI;

   // ── 2. Build or load mesh ──────────────────────────────────────────────────
   // We work exclusively with hexahedral meshes.  If no file is given, build
   // a simple 4x4x4 unit-cube hex mesh.
   Mesh mesh;
   if (string(mesh_file).empty())
   {
      cout << "Building 4x4x4 Cartesian hexahedral mesh.\n";
      mesh = Mesh::MakeCartesian3D(4, 4, 4, Element::HEXAHEDRON,
                                   1.0, 1.0, 1.0);
   }
   else
   {
      mesh = Mesh(mesh_file, 1, 1);
   }

   int dim  = mesh.Dimension();
   MFEM_VERIFY(dim == 3, "This example requires a 3-D hexahedral mesh.");

   // Activate NCMesh tracking (mandatory for hanging-node hex refinement).
   // false = don't force simplices into NC mode (irrelevant for pure hex).
   mesh.EnsureNCMesh(false);

   // ── 3. Uniform initial refinement ─────────────────────────────────────────
   for (int l = 0; l < ref_levels; l++)
   {
      // UniformRefinement on an NCMesh-tracked hex mesh still preserves
      // the nc_limit = 0 (conforming) state because all elements split.
      mesh.UniformRefinement();
   }
   cout << "Initial mesh: " << mesh.GetNE() << " elements.\n";

   // ── 4. FE space (Nedelec H(curl)) ─────────────────────────────────────────
   ND_FECollection       fec(order, dim);
   FiniteElementSpace    fes(&mesh, &fec);

   // Solution GridFunction initialised to zero.
   GridFunction x(&fes);
   x = 0.0;

   // ── 5. GLVis socket (optional) ────────────────────────────────────────────
   socketstream sout;
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      sout.open(vishost, visport);
      if (!sout) { visualization = false; }
      else       { sout.precision(8); }
   }

   // ── 6. AMR loop ───────────────────────────────────────────────────────────
   //
   // Error indicator: per-element L2 norm of (E_exact - E_h).
   // Because we know the exact solution we can use it directly as a reliable
   // AMR indicator.  The nc_limit = 1 argument to RefineByError /
   // DerefineByError enforces the octree 2:1 balance rule: no face may border
   // elements that differ by more than one refinement level.
   //
   // In a production solver without a known exact solution, replace this with
   // a residual-based or jump-penalty estimator computed from curl
   // discontinuities across interior faces.

   for (int it = 0; it <= amr_iter; it++)
   {
      cout << "\n=== AMR iteration " << it
           << "  |  elements: " << mesh.GetNE()
           << "  |  true DOFs: " << fes.GetTrueVSize() << " ===\n";

      // Solve and report global L2 error.
      real_t global_err = SolveOnCurrentMesh(mesh, fes, x,
                                              visualization,
                                              visualization ? &sout : nullptr,
                                              it);
      cout << "|| E_h - E ||_{L2} = " << global_err << "\n";

      if (it == amr_iter) { break; }

      // ── 6a. Per-element error indicator ───────────────────────────────────
      // ComputeElementL2Errors fills a Vector with the L2 error norm on
      // each element; this serves as our refinement indicator.
      VectorFunctionCoefficient E_coef(mesh.SpaceDimension(), E_exact);
      Vector elem_errors(mesh.GetNE());
      x.ComputeElementL2Errors(E_coef, elem_errors);

      // ── 6b. Refine elements above the threshold ────────────────────────────
      //   nc_limit = 1 means MFEM calls NCMesh::LimitNCLevel(1) after
      //   refinement, adding extra splits to neighbours as needed to maintain
      //   the 2:1 octree balance.
      bool refined = mesh.RefineByError(elem_errors, max_err,
                                         /*nonconforming=*/-1,
                                         /*nc_limit=*/1);
      if (!refined)
      {
         cout << "All elements satisfy the error threshold – stopping.\n";
         break;
      }

      // ── 6c. Update FES and prolongate solution ─────────────────────────────
      //   fes.Update()  rebuilds the DOF table and the conforming P for the
      //   new NC mesh (master/slave edge constraints at hanging faces).
      //   x.Update()    uses the RefinementOperator (P_ref) to prolongate
      //   the current solution onto the refined mesh.
      fes.Update();
      x.Update();

      // ── 6d. Optional derefinement ──────────────────────────────────────────
      // Elements with error well below the threshold are coarsened back.
      // nc_limit = 1 prevents coarsening that would violate the 2:1 balance.
      bool derefined = mesh.DerefineByError(elem_errors,
                                             hysteresis * max_err,
                                             /*nc_limit=*/1);
      if (derefined)
      {
         cout << "   Derefinement applied.\n";
         fes.Update();
         x.Update();
      }
   }

   // ── 10. Save final mesh and solution ──────────────────────────────────────
   {
      ofstream mesh_ofs("octree_em.mesh");
      mesh_ofs.precision(8);
      mesh.Print(mesh_ofs);

      ofstream sol_ofs("octree_em_sol.gf");
      sol_ofs.precision(8);
      x.Save(sol_ofs);

      cout << "\nFinal mesh saved to octree_em.mesh\n"
           << "Solution saved to octree_em_sol.gf\n"
           << "(View with: glvis -m octree_em.mesh -g octree_em_sol.gf)\n";
   }

   return 0;
}

// =============================================================================
// SolveOnCurrentMesh
// Assembles and solves  (curl mu^{-1} curl + sigma) E = f  on the current mesh.
// Returns the L2 error against the exact solution.
// =============================================================================
static real_t SolveOnCurrentMesh(Mesh               &mesh,
                                   FiniteElementSpace &fes,
                                   GridFunction       &x,
                                   bool                visualization,
                                   socketstream       *sout,
                                   int                 amr_it)
{
   int dim  = mesh.Dimension();
   int sdim = mesh.SpaceDimension();

   // ── Essential BCs: PEC (E x n = 0) on all boundary faces ─────────────────
   // Translates to zero tangential E ↔ ess_tdof on all boundary edges.
   Array<int> ess_bdr(mesh.bdr_attributes.Max());
   ess_bdr = 1;
   Array<int> ess_tdof_list;
   fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

   // Coefficients: mu^{-1} = 1, sigma = 1.
   ConstantCoefficient muinv(1.0);
   ConstantCoefficient sigma(1.0);

   // ── Bilinear form:  a(u,v) = (mu^{-1} curl u, curl v) + (sigma u, v) ─────
   BilinearForm a(&fes);
   a.AddDomainIntegrator(new CurlCurlIntegrator(muinv));
   a.AddDomainIntegrator(new VectorFEMassIntegrator(sigma));
   a.Assemble();

   // ── Linear form:  b(v) = (f, v) ──────────────────────────────────────────
   VectorFunctionCoefficient f_coef(sdim, f_exact);
   LinearForm b(&fes);
   b.AddDomainIntegrator(new VectorFEDomainLFIntegrator(f_coef));
   b.Assemble();

   // ── Project exact solution to set non-homogeneous BC on boundary edges ────
   VectorFunctionCoefficient E_coef(sdim, E_exact);
   x.ProjectCoefficient(E_coef);

   // ── Form reduced system:  P^T A P X = P^T b ──────────────────────────────
   // FormLinearSystem automatically applies the conforming prolongation P
   // (built from the NCMesh master/slave lists) and eliminates essential DOFs.
   OperatorPtr A;
   Vector      B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

   cout << "   System size (true DOFs after NC constraints): "
        << A->Height() << "\n";

   // ── Solve ─────────────────────────────────────────────────────────────────
   // GSSmoother + PCG is adequate for small meshes.
   // For large 3-D problems replace with HypreAMS (see ex3p.cpp).
#ifndef MFEM_USE_SUITESPARSE
   GSSmoother M((SparseMatrix &)(*A));
   PCG(*A, M, B, X, 1, 1000, 1e-12, 0.0);
#else
   UMFPackSolver umf_solver;
   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
   umf_solver.SetOperator(*A);
   umf_solver.Mult(B, X);
#endif

   // ── Recover FEM solution from the reduced system ──────────────────────────
   a.RecoverFEMSolution(X, b, x);

   // ── Visualise ─────────────────────────────────────────────────────────────
   if (visualization && sout && sout->is_open())
   {
      *sout << "solution\n" << mesh << x
            << "window_title 'AMR it " << amr_it << "'\n" << flush;
   }

   // ── Return L2 error ───────────────────────────────────────────────────────
   return x.ComputeL2Error(E_coef);
}

// =============================================================================
// Exact solution  E = (sin(kappa y), sin(kappa z), sin(kappa x))
// =============================================================================
void E_exact(const Vector &x, Vector &E)
{
   E(0) = sin(kappa * x(1));
   E(1) = sin(kappa * x(2));
   E(2) = sin(kappa * x(0));
}

// =============================================================================
// RHS   f = (curl curl + sigma) E  with mu^{-1}=1, sigma=1
//   curl curl E = kappa^2 * E  for this choice, so f = (1 + kappa^2) * E
// =============================================================================
void f_exact(const Vector &x, Vector &f)
{
   f(0) = (1.0 + kappa * kappa) * sin(kappa * x(1));
   f(1) = (1.0 + kappa * kappa) * sin(kappa * x(2));
   f(2) = (1.0 + kappa * kappa) * sin(kappa * x(0));
}
