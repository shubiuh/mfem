// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#ifndef MFEM_UDAR3P_SOLVER_HPP
#define MFEM_UDAR3P_SOLVER_HPP

#include "mfem.hpp"
#include "udar3p_pml.hpp"

#include <memory>
#include <string>

namespace mfem
{

/** @brief Controls the iterative Maxwell solver and its preconditioner. */
struct UDARIterativeSolverOptions
{
   /** @brief Use @c "ilu" or the shifted @c "ams" preconditioner. */
   std::string preconditioner = "ilu";
   /** @brief Apply AMS to the complex operator or its assembled real block. */
   std::string system_mode = "complex";
   int maximum_iterations = 2000;
   int krylov_dimension = 200;
   real_t relative_tolerance = 1.0e-4;
   real_t absolute_tolerance = 0.0;
   int ams_cycle_type = 14;
   int ams_iterations = 2;
   int ams_smoothing_iterations = 1;
   int ams_print_level = 0;
   bool ilu_fallback = true;
};

/** @brief Iterative solver for the packed real representation of UDAR fields.

    The AMS path assembles the positive auxiliary Maxwell operator
    @f$\nabla\times\mu_r^{-1}\nabla\times +
    (k_0^2\epsilon_r I+\omega\mu_0\sigma)@f$. Two copies of its inverse
    precondition the real and imaginary field blocks. With the Hermitian MFEM
    convention, the imaginary copy receives the required minus sign.

    The @c complex system mode applies FGMRES directly to
    ComplexOperator and avoids materializing its full real matrix. A partially
    assembled operator requires AMS, complex mode and disabled ILU fallback;
    the positive AMS auxiliary matrix is still fully assembled. The
    @c real mode uses the assembled distributed real 2-by-2 matrix. If AMS
    fails to converge and fallback is enabled, the solve is retried with the
    existing distributed Schur-corrected ILU preconditioner. */
class UDARIterativeSolver : public Solver
{
 private:
   const ComplexOperator &complex_operator;
   UDARIterativeSolverOptions options;
   mutable std::unique_ptr<HypreParMatrix> real_block_operator;
   std::unique_ptr<HypreParMatrix> shifted_operator;
   std::unique_ptr<HypreAMS> ams;
   std::unique_ptr<ScaledOperator> scaled_ams;
   std::unique_ptr<BlockDiagonalPreconditioner> block_preconditioner;
   mutable std::unique_ptr<HypreILU> ilu;
   std::unique_ptr<FGMRESSolver> primary_solver;
   mutable std::unique_ptr<FGMRESSolver> fallback_solver;
   mutable bool used_fallback = false;
   mutable int last_iterations = 0;
   mutable real_t last_final_norm = 0.0;

   void ConfigureKrylov(FGMRESSolver &solver) const;
   void BuildILUFallback() const;

 public:
   UDARIterativeSolver(const ComplexOperator &operator_,
                       ParFiniteElementSpace &edge_space,
                       MatrixCoefficient &conductivity,
                       const UDARCartesianPML *pml, real_t omega,
                       real_t epsilon0, real_t mu0, real_t eps_r, real_t mu_r,
                       const Array<int> &essential_tdofs,
                       Array<int> *absorbing_boundary, real_t absorbing_mass,
                       const UDARIterativeSolverOptions &options_);

   void SetOperator(const Operator &op) override;
   void Mult(const Vector &right_hand_side, Vector &solution) const override;

   /** @brief Return whether the most recent solve required ILU recovery. */
   bool UsedFallback() const { return used_fallback; }

   /** @brief Return Krylov iterations used by the most recent solve.

       When AMS recovery is required, this includes both the failed primary
       solve and the successful ILU solve. */
   int GetLastIterations() const { return last_iterations; }

   /** @brief Return the final residual norm from the successful solver. */
   real_t GetLastFinalNorm() const { return last_final_norm; }

   /** @brief Return a concise description for logs. */
   std::string Description() const;
};

} // namespace mfem

#endif
