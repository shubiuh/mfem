// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#include "udar3p_solver.hpp"

#include <cmath>

namespace mfem
{

namespace
{

class ShiftedMassCoefficient : public MatrixCoefficient
{
 private:
   MatrixCoefficient &conductivity;
   real_t dielectric_mass;
   real_t conductivity_scale;

 public:
   ShiftedMassCoefficient(MatrixCoefficient &conductivity_,
                          real_t dielectric_mass_, real_t conductivity_scale_)
       : MatrixCoefficient(3), conductivity(conductivity_),
         dielectric_mass(dielectric_mass_),
         conductivity_scale(conductivity_scale_)
   {
   }

   void Eval(DenseMatrix &value, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      conductivity.Eval(value, T, ip);
      value *= conductivity_scale;
      for (int d = 0; d < 3; d++) { value(d, d) += dielectric_mass; }
   }
};

} // namespace

UDARIterativeSolver::UDARIterativeSolver(
    const ComplexOperator &operator_, ParFiniteElementSpace &edge_space,
    MatrixCoefficient &conductivity, const UDARCartesianPML *pml, real_t omega,
    real_t epsilon0, real_t mu0, real_t eps_r, real_t mu_r,
    const Array<int> &essential_tdofs, Array<int> *absorbing_boundary,
    real_t absorbing_mass, const UDARIterativeSolverOptions &options_)
    : Solver(operator_.Height()), complex_operator(operator_), options(options_)
{
   MFEM_VERIFY(dynamic_cast<const ComplexHypreParMatrix *>(&operator_) ||
                   (options.preconditioner == "ams" &&
                    options.system_mode == "complex" && !options.ilu_fallback),
               "A matrix-free operator requires complex AMS without ILU fallback.");
   MFEM_VERIFY(options.preconditioner == "ilu" ||
                   options.preconditioner == "ams",
               "Iterative preconditioner must be 'ilu' or 'ams'.");
   MFEM_VERIFY(options.system_mode == "complex" ||
                   options.system_mode == "real",
               "Iterative system mode must be 'complex' or 'real'.");

   if (options.preconditioner == "ilu")
   {
      BuildILUFallback();
      primary_solver = std::move(fallback_solver);
      return;
   }

   int rank = 0;
   MPI_Comm_rank(edge_space.GetComm(), &rank);
   if (rank == 0 && options.ams_print_level > 0)
   {
      out << "AMS: assembling positive shifted H(curl) operator ...\n"
          << std::flush;
   }

   const real_t wave_number = omega * std::sqrt(mu0 * epsilon0);
   ConstantCoefficient inverse_mu_r(1.0 / mu_r);
   ShiftedMassCoefficient shifted_mass(
       conductivity, wave_number * wave_number * eps_r, omega * mu0);
   std::unique_ptr<UDARPMLShiftedCurlCoefficient> pml_curl;
   std::unique_ptr<UDARPMLShiftedMassCoefficient> pml_mass;
   ParBilinearForm shifted_form(&edge_space);
   if (pml)
   {
      pml_curl.reset(
          new UDARPMLShiftedCurlCoefficient(*pml, omega, epsilon0, 1.0 / mu_r));
      pml_mass.reset(new UDARPMLShiftedMassCoefficient(
          conductivity, *pml, omega, epsilon0,
          wave_number * wave_number * eps_r, omega * mu0));
      shifted_form.AddDomainIntegrator(new CurlCurlIntegrator(*pml_curl));
      shifted_form.AddDomainIntegrator(new VectorFEMassIntegrator(*pml_mass));
   }
   else
   {
      shifted_form.AddDomainIntegrator(new CurlCurlIntegrator(inverse_mu_r));
      shifted_form.AddDomainIntegrator(
          new VectorFEMassIntegrator(shifted_mass));
   }
   ConstantCoefficient boundary_mass(absorbing_mass);
   if (absorbing_boundary && absorbing_mass > 0.0)
   {
      shifted_form.AddBoundaryIntegrator(
          new VectorFEMassIntegrator(boundary_mass), *absorbing_boundary);
   }
   shifted_form.Assemble();
   OperatorPtr shifted_operator_pointer;
   shifted_form.FormSystemMatrix(essential_tdofs, shifted_operator_pointer);
   // FormSystemMatrix returns a non-owning view of the matrix retained by the
   // bilinear form. Keep an independent matrix because shifted_form is local.
   shifted_operator.reset(
       new HypreParMatrix(*shifted_operator_pointer.As<HypreParMatrix>()));

   if (rank == 0 && options.ams_print_level > 0)
   {
      out << "AMS: constructing auxiliary-space interpolation ...\n"
          << std::flush;
   }

   ams.reset(new HypreAMS(*shifted_operator, &edge_space));
   HYPRE_Solver hypre_ams = static_cast<HYPRE_Solver>(*ams);
   HYPRE_AMSSetCycleType(hypre_ams, options.ams_cycle_type);
   HYPRE_AMSSetMaxIter(hypre_ams, options.ams_iterations);
   HYPRE_AMSSetTol(hypre_ams, 0.0);
   HYPRE_AMSSetSmoothingOptions(hypre_ams, 2, options.ams_smoothing_iterations,
                                1.0, 1.0);
   ams->SetPrintLevel(options.ams_print_level);

   if (rank == 0 && options.ams_print_level > 0)
   {
      out << "AMS: configuring block preconditioner ...\n" << std::flush;
   }

   Array<int> offsets(3);
   offsets[0] = 0;
   offsets[1] = edge_space.GetTrueVSize();
   offsets[2] = edge_space.GetTrueVSize();
   offsets.PartialSum();
   scaled_ams.reset(new ScaledOperator(ams.get(), -1.0));
   block_preconditioner.reset(new BlockDiagonalPreconditioner(offsets));
   block_preconditioner->SetDiagonalBlock(0, ams.get());
   block_preconditioner->SetDiagonalBlock(1, scaled_ams.get());

   primary_solver.reset(new FGMRESSolver(edge_space.GetComm()));
   ConfigureKrylov(*primary_solver);
   if (options.system_mode == "real")
   {
      real_block_operator.reset(
          dynamic_cast<const ComplexHypreParMatrix &>(complex_operator)
          .GetSystemMatrix());
      primary_solver->SetOperator(*real_block_operator);
   }
   else
   {
      primary_solver->SetOperator(complex_operator);
   }
   primary_solver->SetPreconditioner(*block_preconditioner);
   if (rank == 0 && options.ams_print_level > 0)
   {
      out << "AMS: setup complete.\n" << std::flush;
   }
}

void UDARIterativeSolver::ConfigureKrylov(FGMRESSolver &solver) const
{
   solver.SetPrintLevel(0);
   solver.SetKDim(options.krylov_dimension);
   solver.SetMaxIter(options.maximum_iterations);
   solver.SetRelTol(options.relative_tolerance);
   solver.SetAbsTol(options.absolute_tolerance);
}

void UDARIterativeSolver::BuildILUFallback() const
{
   if (!real_block_operator)
   {
      const auto *matrix =
          dynamic_cast<const ComplexHypreParMatrix *>(&complex_operator);
      MFEM_VERIFY(matrix, "ILU requires an assembled complex matrix.");
      real_block_operator.reset(matrix->GetSystemMatrix());
   }
   ilu.reset(new HypreILU);
   ilu->SetType(10);
   ilu->SetLevelOfFill(2);
   ilu->SetOperator(*real_block_operator);
   fallback_solver.reset(new FGMRESSolver(real_block_operator->GetComm()));
   ConfigureKrylov(*fallback_solver);
   fallback_solver->SetOperator(*real_block_operator);
   fallback_solver->SetPreconditioner(*ilu);
}

void UDARIterativeSolver::SetOperator(const Operator &)
{
   MFEM_ABORT("UDARIterativeSolver is constructed for one frequency operator.");
}

void UDARIterativeSolver::Mult(const Vector &right_hand_side,
                               Vector &solution) const
{
   used_fallback = false;
   primary_solver->Mult(right_hand_side, solution);
   last_iterations = primary_solver->GetNumIterations();
   last_final_norm = primary_solver->GetFinalNorm();
   if (primary_solver->GetConverged()) { return; }

   MFEM_VERIFY(options.preconditioner == "ams" && options.ilu_fallback,
               "FGMRES failed to converge with " << Description() << '.');
   if (!fallback_solver) { BuildILUFallback(); }
   solution = 0.0;
   fallback_solver->Mult(right_hand_side, solution);
   used_fallback = true;
   last_iterations += fallback_solver->GetNumIterations();
   last_final_norm = fallback_solver->GetFinalNorm();
   MFEM_VERIFY(fallback_solver->GetConverged(),
               "FGMRES failed with both shifted AMS and distributed ILU.");
}

std::string UDARIterativeSolver::Description() const
{
   if (options.preconditioner == "ilu") { return "FGMRES+ILU"; }
   return "FGMRES+shifted-AMS(" + options.system_mode + ")";
}

} // namespace mfem
