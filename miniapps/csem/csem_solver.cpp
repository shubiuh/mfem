// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#include "csem_solver.hpp"

namespace mfem
{
namespace csem
{

WeissBlockOperator::WeissBlockOperator(const Operator &B_, const Operator &M_)
   : Operator(2 * B_.Height()), B(B_), M(M_), block_size(B_.Height()),
     work(block_size)
{
}

void WeissBlockOperator::Mult(const Vector &x, Vector &y) const
{
   x1.MakeRef(const_cast<Vector &>(x), 0, block_size);
   x2.MakeRef(const_cast<Vector &>(x), block_size, block_size);
   y1.MakeRef(y, 0, block_size);
   y2.MakeRef(y, block_size, block_size);

   M.Mult(x1, y1);
   B.Mult(x2, work);
   y1 -= work;
   B.Mult(x1, y2);
   M.Mult(x2, work);
   y2 += work;
}

PRESBPreconditioner::PRESBPreconditioner(const Operator &mass_,
                                         Solver &inner_solver_)
   : Solver(2 * mass_.Height()), mass(mass_), inner_solver(inner_solver_),
     iterative_inner_solver(dynamic_cast<IterativeSolver *>(&inner_solver_)),
     block_size(mass_.Height()), g_rhs(block_size), g(block_size),
     mass_g(block_size), h_rhs(block_size), h(block_size)
{
}

void PRESBPreconditioner::SetOperator(const Operator &)
{
   MFEM_ABORT("PRESBPreconditioner is bound to one CSEM operator.");
}

void PRESBPreconditioner::Mult(const Vector &f, Vector &w) const
{
   f1.MakeRef(const_cast<Vector &>(f), 0, block_size);
   f2.MakeRef(const_cast<Vector &>(f), block_size, block_size);
   w1.MakeRef(w, 0, block_size);
   w2.MakeRef(w, block_size, block_size);

   g_rhs = f1;
   g_rhs += f2;
   g = 0.0;
   inner_solver.Mult(g_rhs, g);
   if (iterative_inner_solver)
   {
      MFEM_VERIFY(iterative_inner_solver->GetConverged(),
                  "The first inner PRESB solve did not converge.");
      inner_iterations += iterative_inner_solver->GetNumIterations();
   }

   mass.Mult(g, mass_g);
   h_rhs = f1;
   h_rhs -= mass_g;
   h = 0.0;
   inner_solver.Mult(h_rhs, h);
   if (iterative_inner_solver)
   {
      MFEM_VERIFY(iterative_inner_solver->GetConverged(),
                  "The second inner PRESB solve did not converge.");
      inner_iterations += iterative_inner_solver->GetNumIterations();
   }

   w1 = g;
   w1 += h;
   w2 = h;
   w2 *= -1.0;
   applications++;
}

} // namespace csem
} // namespace mfem
