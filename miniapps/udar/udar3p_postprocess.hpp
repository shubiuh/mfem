// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#ifndef MFEM_UDAR3P_POSTPROCESS_HPP
#define MFEM_UDAR3P_POSTPROCESS_HPP

#include "mfem.hpp"

#include <memory>
#include <string>
#include <vector>

namespace mfem
{

/** @brief Local-rank output timings; evaluation is measured only on request.
    Derived evaluation includes host FE interpolation and coefficient algebra,
    but excludes serialization. Profiling adds two clock reads per evaluation;
    compare with an unprofiled run before drawing performance conclusions. */
struct UDARVTKOutputStats
{
   double save_seconds = 0.0;
   double derived_evaluation_seconds = 0.0;
   long long derived_evaluations = 0;
};

/** @brief Write the production ParaView E/H and derived-field collection.
    The optional profile measures host derived-coefficient evaluation inside
    Save(). The derived_fields switch is for controlled output benchmarks.
    Returned timings are local; callers must reduce them across MPI ranks. */
UDARVTKOutputStats WriteUDARVTKFields(
    const std::string &directory, const std::string &evaluation_name,
    int frequency_index, int station, real_t frequency, ParMesh &mesh,
    MatrixCoefficient &conductivity, real_t permittivity, real_t permeability,
    const std::vector<std::unique_ptr<ParComplexGridFunction>> &electric_fields,
    const std::vector<std::unique_ptr<ParGridFunction>> &magnetic_real,
    const std::vector<std::unique_ptr<ParGridFunction>> &magnetic_imaginary,
    bool profile = false, bool derived_fields = true);

/** @brief Time-average electric energy density for a complex phasor field.

    This coefficient evaluates @f$\epsilon |E|^2/4@f$ in J/m^3 for a
    nondispersive isotropic permittivity. */
class UDARElectricEnergyDensityCoefficient : public Coefficient
{
 private:
   VectorCoefficient &real_part;
   VectorCoefficient &imaginary_part;
   real_t permittivity;

 public:
   UDARElectricEnergyDensityCoefficient(VectorCoefficient &real_part_,
                                        VectorCoefficient &imaginary_part_,
                                        real_t permittivity_);

   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override;
};

/** @brief Time-average magnetic energy density for a complex phasor field.

    This coefficient evaluates @f$\mu |H|^2/4@f$ in J/m^3 for an isotropic
    permeability. */
class UDARMagneticEnergyDensityCoefficient : public Coefficient
{
 private:
   VectorCoefficient &real_part;
   VectorCoefficient &imaginary_part;
   real_t permeability;

 public:
   UDARMagneticEnergyDensityCoefficient(VectorCoefficient &real_part_,
                                        VectorCoefficient &imaginary_part_,
                                        real_t permeability_);

   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override;
};

/** @brief Time-average conductive loss density for a complex electric field.

    For a real conductivity tensor this evaluates
    @f$\left(E_r^T\sigma E_r + E_i^T\sigma E_i\right)/2@f$ in W/m^3. */
class UDAROhmicLossDensityCoefficient : public Coefficient
{
 private:
   VectorCoefficient &real_part;
   VectorCoefficient &imaginary_part;
   MatrixCoefficient &conductivity;

 public:
   UDAROhmicLossDensityCoefficient(VectorCoefficient &real_part_,
                                   VectorCoefficient &imaginary_part_,
                                   MatrixCoefficient &conductivity_);

   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override;
};

/** @brief Time-average active Poynting vector for complex E and H fields.

    The returned vector is
    @f$\mathrm{Re}(E\times H^*)/2 =
    (E_r\times H_r + E_i\times H_i)/2@f$ in W/m^2. */
class UDARTimeAveragePoyntingCoefficient : public VectorCoefficient
{
 private:
   VectorCoefficient &electric_real;
   VectorCoefficient &electric_imaginary;
   VectorCoefficient &magnetic_real;
   VectorCoefficient &magnetic_imaginary;

 public:
   UDARTimeAveragePoyntingCoefficient(
       VectorCoefficient &electric_real_,
       VectorCoefficient &electric_imaginary_,
       VectorCoefficient &magnetic_real_,
       VectorCoefficient &magnetic_imaginary_);

   void Eval(Vector &value, ElementTransformation &T,
             const IntegrationPoint &ip) override;
};

} // namespace mfem

#endif
