// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#ifndef MFEM_UDAR3P_FREQUENCY_HPP
#define MFEM_UDAR3P_FREQUENCY_HPP

#include "mfem.hpp"
#include "udar3p_pml.hpp"
#include "udar3p_solver.hpp"

#include <memory>
#include <string>

namespace mfem
{

/** @brief Fixed configuration for one frequency-domain Maxwell system. */
struct UDARFrequencySystemOptions
{
   real_t frequency = 0.0;
   real_t epsilon0 = 0.0;
   real_t mu0 = 0.0;
   real_t relative_permittivity = 1.0;
   real_t relative_permeability = 1.0;
   int boundary_condition = 0;
   bool use_mumps = false;
   bool reuse_mumps_ordering = true;
   /** @brief Select @c auto, @c full, or guarded @c partial assembly. */
   std::string assembly = "auto";
   UDARIterativeSolverOptions iterative;
};

/** @brief Statistics accumulated while solving one set of source vectors. */
struct UDARFrequencySolveStats
{
   int krylov_iterations = 0;
   int ilu_fallbacks = 0;
};

/** @brief Own one frequency's Maxwell operator and linear solver.

    The object persists across inversion evaluations. AssembleOperator()
    rebuilds material-dependent coefficients, forms, and matrices, while a
    MUMPS instance remains alive so it can reuse symbolic ordering for an
    unchanged mesh and sparsity pattern. Iterative solver and preconditioner
    state are rebuilt with the operator because their shifted coefficients
    depend on conductivity. */
class UDARFrequencySystem
{
 private:
   ParFiniteElementSpace &space;
   MatrixCoefficient &conductivity;
   const UDARCartesianPML *pml;
   Array<int> essential_tdofs;
   UDARFrequencySystemOptions options;
   real_t omega;
   real_t wave_number;
   real_t absorbing_mass = 0.0;
   bool partial_assembly = false;
   std::string assembly_reason;

   std::unique_ptr<Coefficient> inverse_permeability;
   std::unique_ptr<Coefficient> real_mass;
   std::unique_ptr<MatrixCoefficient> imaginary_mass;
   std::unique_ptr<UDARPMLCurlCoefficient> pml_curl_real;
   std::unique_ptr<UDARPMLCurlCoefficient> pml_curl_imaginary;
   std::unique_ptr<UDARPMLMassCoefficient> pml_mass_real;
   std::unique_ptr<UDARPMLMassCoefficient> pml_mass_imaginary;
   std::unique_ptr<Coefficient> absorbing_real;
   std::unique_ptr<Coefficient> absorbing_imaginary;
   Array<int> absorbing_boundary;
   std::unique_ptr<ParSesquilinearForm> form;
   OperatorPtr system_operator;
   std::unique_ptr<HypreParMatrix> real_block_operator;
   std::unique_ptr<Solver> direct_solver;
   std::unique_ptr<UDARIterativeSolver> iterative_solver;

   void ResetOperator();
   void CreateDirectSolver();

 public:
   UDARFrequencySystem(ParFiniteElementSpace &space_,
                       MatrixCoefficient &conductivity_,
                       const UDARCartesianPML *pml_,
                       const Array<int> &essential_tdofs_,
                       const UDARFrequencySystemOptions &options_);

   /** @brief Assemble the Maxwell operator for the current conductivity.
       @param background_conductivity Conductivity used by the absorbing
       boundary approximation, in S/m. */
   void AssembleOperator(real_t background_conductivity);

   /** @brief Configure or update the selected solver for the assembled
       matrix. */
   void ConfigureSolver();

   /** @brief Solve all supplied source vectors.
       @return Aggregate iterative-solver statistics; zeros for MUMPS. */
   UDARFrequencySolveStats Solve(const Array<const Vector *> &right_hand_sides,
                                 Array<Vector *> &solutions) const;

   /** @brief Return a concise solver description for logs and metadata. */
   std::string SolverDescription() const;

   real_t Frequency() const { return options.frequency; }
   real_t Omega() const { return omega; }
   bool UsesMUMPS() const { return options.use_mumps; }
   /** @brief Actual Maxwell assembly mode; AMS remains fully assembled. */
   const char *AssemblyMode() const
   { return partial_assembly ? "partial" : "full"; }
   /** @brief Explain why full assembly was selected. */
   const std::string &AssemblyReason() const { return assembly_reason; }
   /** @brief Read-only assembled or partially assembled operator view.
       Invalidated by the next AssembleOperator() call. */
   const ComplexOperator &GetOperator() const
   {
      MFEM_VERIFY(system_operator.Ptr(), "The operator has not been assembled.");
      return *system_operator.As<ComplexOperator>();
   }
   const std::string &Preconditioner() const
   { return options.iterative.preconditioner; }
   const std::string &SystemMode() const
   { return options.iterative.system_mode; }
};

} // namespace mfem

#endif
