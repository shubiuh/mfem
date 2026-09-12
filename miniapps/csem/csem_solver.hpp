// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#ifndef MFEM_MINIAPPS_CSEM_SOLVER_HPP
#define MFEM_MINIAPPS_CSEM_SOLVER_HPP

#include "mfem.hpp"

namespace mfem
{
namespace csem
{

/// The Weiss real system [M, -B; B, M] for B = K - M_epsilon.
class WeissBlockOperator : public Operator
{
private:
   const Operator &B;
   const Operator &M;
   int block_size;
   mutable Vector x1, x2, y1, y2, work;

public:
   WeissBlockOperator(const Operator &B_, const Operator &M_);
   void Mult(const Vector &x, Vector &y) const override;
};

/// Algorithm 1 from Weiss et al.: two H solves and one mass product.
class PRESBPreconditioner : public Solver
{
private:
   const Operator &mass;
   Solver &inner_solver;
   IterativeSolver *iterative_inner_solver;
   int block_size;
   mutable Vector f1, f2, w1, w2;
   mutable Vector g_rhs, g, mass_g, h_rhs, h;
   mutable long long applications = 0;
   mutable long long inner_iterations = 0;

public:
   PRESBPreconditioner(const Operator &mass_, Solver &inner_solver_);
   void SetOperator(const Operator &) override;
   void Mult(const Vector &f, Vector &w) const override;

   long long Applications() const { return applications; }
   long long InnerIterations() const { return inner_iterations; }
};

} // namespace csem
} // namespace mfem

#endif
