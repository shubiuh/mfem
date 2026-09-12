// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#include "udar3p_frequency.hpp"

#include <cmath>
#include <complex>

namespace mfem
{

namespace
{

class FrequencyScaledMatrixCoefficient : public MatrixCoefficient
{
 private:
   MatrixCoefficient &coefficient;
   real_t scale;

 public:
   FrequencyScaledMatrixCoefficient(MatrixCoefficient &coefficient_,
                                    real_t scale_)
      : MatrixCoefficient(3), coefficient(coefficient_), scale(scale_)
   {
   }

   void Eval(DenseMatrix &value, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      coefficient.Eval(value, T, ip);
      value *= scale;
   }
};

} // namespace

UDARFrequencySystem::UDARFrequencySystem(
    ParFiniteElementSpace &space_, MatrixCoefficient &conductivity_,
    const UDARCartesianPML *pml_, const Array<int> &essential_tdofs_,
    const UDARFrequencySystemOptions &options_)
   : space(space_), conductivity(conductivity_), pml(pml_),
     essential_tdofs(essential_tdofs_), options(options_),
     omega(2.0 * M_PI * options.frequency),
     wave_number(omega * std::sqrt(options.mu0 * options.epsilon0))
{
   MFEM_VERIFY(options.frequency > 0.0, "Frequency must be positive.");
   MFEM_VERIFY(options.epsilon0 > 0.0 && options.mu0 > 0.0,
               "Vacuum material constants must be positive.");
   MFEM_VERIFY(options.relative_permittivity > 0.0 &&
                   options.relative_permeability > 0.0,
               "Relative material constants must be positive.");
   MFEM_VERIFY(options.boundary_condition == 0 ||
                   options.boundary_condition == 1,
               "Boundary condition must be PEC or first-order absorbing.");
   MFEM_VERIFY(options.assembly == "auto" || options.assembly == "full" ||
                   options.assembly == "partial",
               "Assembly must be 'auto', 'full', or 'partial'.");
   // Check every rank: one unsupported element rules out the PA kernels.
   int local_supported = !space.IsVariableOrder() &&
                         space.FEColl()->GetOrder() >= 2;
   for (int element = 0; element < space.GetNE(); element++)
   {
      if (space.GetMesh()->GetElementBaseGeometry(element) != Geometry::CUBE)
      {
         local_supported = 0;
      }
   }
   int supported = 0;
   MPI_Allreduce(&local_supported, &supported, 1, MPI_INT, MPI_MIN,
                 space.GetComm());
   if (options.assembly == "full")
   { assembly_reason = "full assembly requested"; }
   else if (options.use_mumps)
   { assembly_reason = "MUMPS requires an assembled matrix"; }
   else if (options.iterative.preconditioner != "ams")
   { assembly_reason = "partial assembly requires AMS"; }
   else if (options.iterative.system_mode != "complex")
   { assembly_reason = "the real-block mode requires an assembled matrix"; }
   else if (options.iterative.ilu_fallback)
   { assembly_reason = "ILU fallback requires an assembled matrix"; }
   else if (options.iterative.ams_cycle_type != 13)
   { assembly_reason = "partial assembly currently requires validated AMS cycle 13"; }
   else if (!supported)
   { assembly_reason = "partial assembly requires uniform order >= 2 hex elements"; }
   else if (pml)
   { assembly_reason = "PML partial assembly is not enabled"; }
   else if (options.boundary_condition != 0)
   { assembly_reason = "partial assembly currently requires PEC boundaries"; }
   else
   { partial_assembly = true; }
   if (options.use_mumps) { CreateDirectSolver(); }
}

void UDARFrequencySystem::ResetOperator()
{
   iterative_solver.reset();
   real_block_operator.reset();
   system_operator.Clear();
   form.reset();
   absorbing_imaginary.reset();
   absorbing_real.reset();
   pml_mass_imaginary.reset();
   pml_mass_real.reset();
   pml_curl_imaginary.reset();
   pml_curl_real.reset();
   imaginary_mass.reset();
   real_mass.reset();
   inverse_permeability.reset();
   absorbing_boundary.SetSize(0);
   absorbing_mass = 0.0;
}

void UDARFrequencySystem::CreateDirectSolver()
{
#ifdef MFEM_USE_COMPLEX_MUMPS
   auto *mumps = new ComplexMUMPSSolver(space.GetComm());
   mumps->SetPrintLevel(0);
   mumps->SetMatrixSymType(ComplexMUMPSSolver::UNSYMMETRIC);
   mumps->SetReorderingReuse(options.reuse_mumps_ordering);
   direct_solver.reset(mumps);
#elif defined(MFEM_USE_MUMPS)
   auto *mumps = new MUMPSSolver(space.GetComm());
   mumps->SetPrintLevel(0);
   mumps->SetMatrixSymType(MUMPSSolver::UNSYMMETRIC);
   mumps->SetReorderingReuse(options.reuse_mumps_ordering);
   direct_solver.reset(mumps);
#else
   MFEM_ABORT("MUMPS was requested, but MFEM_USE_MUMPS is disabled.");
#endif
}

void UDARFrequencySystem::AssembleOperator(real_t background_conductivity)
{
   ResetOperator();
   const ComplexOperator::Convention convention = ComplexOperator::HERMITIAN;
   inverse_permeability.reset(
       new ConstantCoefficient(1.0 / options.relative_permeability));
   real_mass.reset(new ConstantCoefficient(
       -wave_number * wave_number * options.relative_permittivity));
   imaginary_mass.reset(new FrequencyScaledMatrixCoefficient(
       conductivity, omega * options.mu0));
   form.reset(new ParSesquilinearForm(&space, convention));
   if (partial_assembly) { form->SetAssemblyLevel(AssemblyLevel::PARTIAL); }

   if (pml)
   {
      pml_curl_real.reset(new UDARPMLCurlCoefficient(
          *pml, omega, options.epsilon0,
          1.0 / options.relative_permeability, false));
      pml_curl_imaginary.reset(new UDARPMLCurlCoefficient(
          *pml, omega, options.epsilon0,
          1.0 / options.relative_permeability, true));
      pml_mass_real.reset(new UDARPMLMassCoefficient(
          conductivity, *pml, omega, options.epsilon0,
          -wave_number * wave_number * options.relative_permittivity,
          omega * options.mu0, false));
      pml_mass_imaginary.reset(new UDARPMLMassCoefficient(
          conductivity, *pml, omega, options.epsilon0,
          -wave_number * wave_number * options.relative_permittivity,
          omega * options.mu0, true));
      form->AddDomainIntegrator(
          new CurlCurlIntegrator(*pml_curl_real),
          new CurlCurlIntegrator(*pml_curl_imaginary));
      form->AddDomainIntegrator(
          new VectorFEMassIntegrator(*pml_mass_real),
          new VectorFEMassIntegrator(*pml_mass_imaginary));
   }
   else
   {
      form->AddDomainIntegrator(
          new CurlCurlIntegrator(*inverse_permeability), nullptr);
      form->AddDomainIntegrator(new VectorFEMassIntegrator(*real_mass),
                                new VectorFEMassIntegrator(*imaginary_mass));
   }

   if (options.boundary_condition == 1 &&
       space.GetParMesh()->bdr_attributes.Size())
   {
      absorbing_boundary.SetSize(
          space.GetParMesh()->bdr_attributes.Max());
      absorbing_boundary = 1;
      const std::complex<real_t> imaginary_unit(0.0, 1.0);
      const std::complex<real_t> background_epsilon =
          options.relative_permittivity -
          imaginary_unit * background_conductivity /
              (omega * options.epsilon0);
      const std::complex<real_t> boundary_wave_number =
          wave_number *
          std::sqrt(options.relative_permeability * background_epsilon);
      // exp(+i omega t), outgoing exp(-i k r): the weak boundary term is
      // +i k/mu_r times the tangential-field mass matrix.
      const std::complex<real_t> gamma =
          imaginary_unit * boundary_wave_number /
          options.relative_permeability;
      absorbing_mass = std::abs(gamma);
      absorbing_real.reset(new ConstantCoefficient(gamma.real()));
      absorbing_imaginary.reset(new ConstantCoefficient(gamma.imag()));
      form->AddBoundaryIntegrator(
          new VectorFEMassIntegrator(*absorbing_real),
          new VectorFEMassIntegrator(*absorbing_imaginary),
          absorbing_boundary);
   }

   form->Assemble();
   form->FormSystemMatrix(essential_tdofs, system_operator);
}

void UDARFrequencySystem::ConfigureSolver()
{
   MFEM_VERIFY(system_operator.Ptr(),
               "AssembleOperator must be called before ConfigureSolver.");
   if (options.use_mumps)
   {
#ifdef MFEM_USE_COMPLEX_MUMPS
      direct_solver->SetOperator(
          *system_operator.As<ComplexHypreParMatrix>());
#elif defined(MFEM_USE_MUMPS)
      real_block_operator.reset(
          system_operator.As<ComplexHypreParMatrix>()->GetSystemMatrix());
      direct_solver->SetOperator(*real_block_operator);
#endif
      return;
   }

   iterative_solver.reset(new UDARIterativeSolver(
       *system_operator.As<ComplexOperator>(), space, conductivity, pml,
       omega, options.epsilon0, options.mu0, options.relative_permittivity,
       options.relative_permeability, essential_tdofs,
       absorbing_boundary.Size() ? &absorbing_boundary : nullptr,
       absorbing_mass, options.iterative));
}

UDARFrequencySolveStats UDARFrequencySystem::Solve(
    const Array<const Vector *> &right_hand_sides,
    Array<Vector *> &solutions) const
{
   MFEM_VERIFY(right_hand_sides.Size() == solutions.Size(),
               "Right-hand-side and solution counts must match.");
   UDARFrequencySolveStats stats;
   if (options.use_mumps)
   {
      MFEM_VERIFY(direct_solver, "The direct solver is not configured.");
      direct_solver->ArrayMult(right_hand_sides, solutions);
      return stats;
   }

   MFEM_VERIFY(iterative_solver, "The iterative solver is not configured.");
   for (int source = 0; source < right_hand_sides.Size(); source++)
   {
      iterative_solver->Mult(*right_hand_sides[source], *solutions[source]);
      stats.krylov_iterations += iterative_solver->GetLastIterations();
      stats.ilu_fallbacks += iterative_solver->UsedFallback() ? 1 : 0;
   }
   return stats;
}

std::string UDARFrequencySystem::SolverDescription() const
{
   return options.use_mumps ? "MUMPS" : iterative_solver->Description();
}

} // namespace mfem
