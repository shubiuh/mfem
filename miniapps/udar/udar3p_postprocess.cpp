// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#include "udar3p_postprocess.hpp"

#include <cmath>
#include <chrono>
#include <iomanip>
#include <sstream>

using namespace std;

namespace mfem
{

namespace
{

/** @brief Instrument host evaluation without changing coefficient values. */
class TimedScalarCoefficient : public Coefficient
{
   Coefficient &coefficient;
   UDARVTKOutputStats &stats;
 public:
   TimedScalarCoefficient(Coefficient &coefficient_, UDARVTKOutputStats &stats_)
      : coefficient(coefficient_), stats(stats_) { }
   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
   {
      const auto start = std::chrono::steady_clock::now();
      const real_t value = coefficient.Eval(T, ip);
      stats.derived_evaluation_seconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - start).count();
      stats.derived_evaluations++;
      return value;
   }
};

class TimedVectorCoefficient : public VectorCoefficient
{
   VectorCoefficient &coefficient;
   UDARVTKOutputStats &stats;
 public:
   TimedVectorCoefficient(VectorCoefficient &coefficient_,
                          UDARVTKOutputStats &stats_)
      : VectorCoefficient(coefficient_.GetVDim()),
        coefficient(coefficient_), stats(stats_) { }
   void Eval(Vector &value, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      const auto start = std::chrono::steady_clock::now();
      coefficient.Eval(value, T, ip);
      stats.derived_evaluation_seconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - start).count();
      stats.derived_evaluations++;
   }
};

/** @brief Magnitude of a complex vector field represented by real and imaginary
   parts.

    The coefficient evaluates @f$\sqrt{\mathbf{u}_r\cdot\mathbf{u}_r +
    \mathbf{u}_i\cdot\mathbf{u}_i}@f$, which is the magnitude of a complex
    phasor vector at the requested physical point. */
class ComplexVectorMagnitudeCoefficient : public Coefficient
{
 private:
   VectorCoefficient &real_part;
   VectorCoefficient &imaginary_part;

 public:
   ComplexVectorMagnitudeCoefficient(VectorCoefficient &real_part_,
                                     VectorCoefficient &imaginary_part_)
       : real_part(real_part_), imaginary_part(imaginary_part_)
   {
      MFEM_VERIFY(real_part.GetVDim() == imaginary_part.GetVDim(),
                  "Complex vector parts must have the same dimension.");
   }

   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
   {
      Vector real_value(real_part.GetVDim());
      Vector imaginary_value(imaginary_part.GetVDim());
      real_part.Eval(real_value, T, ip);
      imaginary_part.Eval(imaginary_value, T, ip);
      return sqrt(InnerProduct(real_value, real_value) +
                  InnerProduct(imaginary_value, imaginary_value));
   }
};

void CrossProduct3(const Vector &a, const Vector &b, Vector &value)
{
   value.SetSize(3);
   value[0] = a[1] * b[2] - a[2] * b[1];
   value[1] = a[2] * b[0] - a[0] * b[2];
   value[2] = a[0] * b[1] - a[1] * b[0];
}

void VerifyComplexParts(VectorCoefficient &real_part,
                        VectorCoefficient &imaginary_part)
{
   MFEM_VERIFY(real_part.GetVDim() == 3 && imaginary_part.GetVDim() == 3,
               "UDAR postprocessing requires three-component fields.");
}

} // namespace

UDARElectricEnergyDensityCoefficient::
UDARElectricEnergyDensityCoefficient(VectorCoefficient &real_part_,
                                     VectorCoefficient &imaginary_part_,
                                     real_t permittivity_)
   : real_part(real_part_), imaginary_part(imaginary_part_),
     permittivity(permittivity_)
{
   VerifyComplexParts(real_part, imaginary_part);
   MFEM_VERIFY(permittivity > 0.0, "Permittivity must be positive.");
}

real_t UDARElectricEnergyDensityCoefficient::Eval(
    ElementTransformation &T, const IntegrationPoint &ip)
{
   Vector real_value(3), imaginary_value(3);
   real_part.Eval(real_value, T, ip);
   imaginary_part.Eval(imaginary_value, T, ip);
   return 0.25 * permittivity *
          (InnerProduct(real_value, real_value) +
           InnerProduct(imaginary_value, imaginary_value));
}

UDARMagneticEnergyDensityCoefficient::
UDARMagneticEnergyDensityCoefficient(VectorCoefficient &real_part_,
                                     VectorCoefficient &imaginary_part_,
                                     real_t permeability_)
   : real_part(real_part_), imaginary_part(imaginary_part_),
     permeability(permeability_)
{
   VerifyComplexParts(real_part, imaginary_part);
   MFEM_VERIFY(permeability > 0.0, "Permeability must be positive.");
}

real_t UDARMagneticEnergyDensityCoefficient::Eval(
    ElementTransformation &T, const IntegrationPoint &ip)
{
   Vector real_value(3), imaginary_value(3);
   real_part.Eval(real_value, T, ip);
   imaginary_part.Eval(imaginary_value, T, ip);
   return 0.25 * permeability *
          (InnerProduct(real_value, real_value) +
           InnerProduct(imaginary_value, imaginary_value));
}

UDAROhmicLossDensityCoefficient::UDAROhmicLossDensityCoefficient(
    VectorCoefficient &real_part_, VectorCoefficient &imaginary_part_,
    MatrixCoefficient &conductivity_)
   : real_part(real_part_), imaginary_part(imaginary_part_),
     conductivity(conductivity_)
{
   VerifyComplexParts(real_part, imaginary_part);
   MFEM_VERIFY(conductivity.GetHeight() == 3 && conductivity.GetWidth() == 3,
               "UDAR loss density requires a 3-by-3 conductivity tensor.");
}

real_t UDAROhmicLossDensityCoefficient::Eval(
    ElementTransformation &T, const IntegrationPoint &ip)
{
   Vector real_value(3), imaginary_value(3), current(3);
   DenseMatrix sigma(3);
   real_part.Eval(real_value, T, ip);
   imaginary_part.Eval(imaginary_value, T, ip);
   conductivity.Eval(sigma, T, ip);
   sigma.Mult(real_value, current);
   real_t loss = InnerProduct(real_value, current);
   sigma.Mult(imaginary_value, current);
   loss += InnerProduct(imaginary_value, current);
   return 0.5 * loss;
}

UDARTimeAveragePoyntingCoefficient::UDARTimeAveragePoyntingCoefficient(
    VectorCoefficient &electric_real_,
    VectorCoefficient &electric_imaginary_,
    VectorCoefficient &magnetic_real_,
    VectorCoefficient &magnetic_imaginary_)
   : VectorCoefficient(3), electric_real(electric_real_),
     electric_imaginary(electric_imaginary_), magnetic_real(magnetic_real_),
     magnetic_imaginary(magnetic_imaginary_)
{
   VerifyComplexParts(electric_real, electric_imaginary);
   VerifyComplexParts(magnetic_real, magnetic_imaginary);
}

void UDARTimeAveragePoyntingCoefficient::Eval(
    Vector &value, ElementTransformation &T, const IntegrationPoint &ip)
{
   Vector electric_r(3), electric_i(3), magnetic_r(3), magnetic_i(3);
   Vector real_cross(3), imaginary_cross(3);
   electric_real.Eval(electric_r, T, ip);
   electric_imaginary.Eval(electric_i, T, ip);
   magnetic_real.Eval(magnetic_r, T, ip);
   magnetic_imaginary.Eval(magnetic_i, T, ip);
   CrossProduct3(electric_r, magnetic_r, real_cross);
   CrossProduct3(electric_i, magnetic_i, imaginary_cross);
   value.SetSize(3);
   value = real_cross;
   value += imaginary_cross;
   value *= 0.5;
}

/** @brief Write complex E and H fields for one frequency and tool station.
    @param directory Root directory selected with @c -vtk.
    @param evaluation_name Stem of the CSV response file.
    @param frequency_index Zero-based index in the requested frequency list.
    @param station Zero-based transmitter station index.
    @param frequency Frequency in Hz, stored as the ParaView time value.
    @param mesh Parallel mesh.
    @param conductivity Physical formation conductivity tensor.
    @param permittivity Physical formation permittivity in F/m.
    @param permeability Physical formation permeability in H/m.
    @param electric_fields Three complex E fields for x, y, and z transmitters.
    @param magnetic_real Real H fields in the RT space.
    @param magnetic_imaginary Imaginary H fields in the RT space.

    Field names use transmitter-axis suffixes: for example @c E_real_x and
    @c H_imag_z. The scalar fields @c E_magnitude_x and @c H_magnitude_x are
    the magnitudes of their corresponding complex vector phasors. The output
    also contains electric and magnetic energy density, physical ohmic loss
    density, and the time-average active Poynting vector for each transmitter.
    Values use the exp(+i omega t) phasor convention. */
UDARVTKOutputStats WriteUDARVTKFields(
    const string &directory, const string &evaluation_name, int frequency_index,
    int station, real_t frequency, ParMesh &mesh,
    MatrixCoefficient &conductivity, real_t permittivity,
    real_t permeability,
    const vector<unique_ptr<ParComplexGridFunction>> &electric_fields,
    const vector<unique_ptr<ParGridFunction>> &magnetic_real,
    const vector<unique_ptr<ParGridFunction>> &magnetic_imaginary,
    bool profile, bool derived_fields)
{
   UDARVTKOutputStats stats;
   if (directory.empty()) { return stats; }

   ostringstream name;
   name << "fields_f" << setw(3) << setfill('0') << frequency_index
        << "_station" << setw(3) << station;
   ParaViewDataCollection collection(name.str(), &mesh);
   collection.SetPrefixPath(directory + "/" + evaluation_name + "/");
   collection.SetDataFormat(VTKFormat::BINARY);
   collection.SetHighOrderOutput(true);
   collection.SetCycle(station);
   collection.SetTime(frequency);
   std::vector<std::unique_ptr<Coefficient>> profiled_scalars;
   std::vector<std::unique_ptr<VectorCoefficient>> profiled_vectors;
   auto register_scalar = [&](const std::string &name, Coefficient *coefficient)
   {
      if (!derived_fields) { return; }
      if (profile)
      {
         profiled_scalars.emplace_back(
             new TimedScalarCoefficient(*coefficient, stats));
         coefficient = profiled_scalars.back().get();
      }
      collection.RegisterCoeffField(name, coefficient);
   };
   auto register_vector = [&](const std::string &name,
                              VectorCoefficient *coefficient)
   {
      if (!derived_fields) { return; }
      if (profile)
      {
         profiled_vectors.emplace_back(
             new TimedVectorCoefficient(*coefficient, stats));
         coefficient = profiled_vectors.back().get();
      }
      collection.RegisterVCoeffField(name, coefficient);
   };
   const char axis[3] = {'x', 'y', 'z'};
   vector<unique_ptr<VectorGridFunctionCoefficient>> field_coefficients;
   vector<unique_ptr<ComplexVectorMagnitudeCoefficient>> magnitude_coefficients;
   vector<unique_ptr<UDARElectricEnergyDensityCoefficient>> electric_energy;
   vector<unique_ptr<UDARMagneticEnergyDensityCoefficient>> magnetic_energy;
   vector<unique_ptr<UDAROhmicLossDensityCoefficient>> ohmic_loss;
   vector<unique_ptr<UDARTimeAveragePoyntingCoefficient>> poynting;
   field_coefficients.reserve(12);
   magnitude_coefficients.reserve(6);
   electric_energy.reserve(3);
   magnetic_energy.reserve(3);
   ohmic_loss.reserve(3);
   poynting.reserve(3);
   for (int transmitter_axis = 0; transmitter_axis < 3; transmitter_axis++)
   {
      const string suffix(1, axis[transmitter_axis]);
      collection.RegisterField("E_real_" + suffix,
                               &electric_fields[transmitter_axis]->real());
      collection.RegisterField("E_imag_" + suffix,
                               &electric_fields[transmitter_axis]->imag());
      collection.RegisterField("H_real_" + suffix,
                               magnetic_real[transmitter_axis].get());
      collection.RegisterField("H_imag_" + suffix,
                               magnetic_imaginary[transmitter_axis].get());
      field_coefficients.emplace_back(
          make_unique<VectorGridFunctionCoefficient>(
              &electric_fields[transmitter_axis]->real()));
      field_coefficients.emplace_back(
          make_unique<VectorGridFunctionCoefficient>(
              &electric_fields[transmitter_axis]->imag()));
      field_coefficients.emplace_back(
          make_unique<VectorGridFunctionCoefficient>(
              magnetic_real[transmitter_axis].get()));
      field_coefficients.emplace_back(
          make_unique<VectorGridFunctionCoefficient>(
              magnetic_imaginary[transmitter_axis].get()));
      magnitude_coefficients.emplace_back(
          make_unique<ComplexVectorMagnitudeCoefficient>(
              *field_coefficients[field_coefficients.size() - 4],
              *field_coefficients[field_coefficients.size() - 3]));
      magnitude_coefficients.emplace_back(
          make_unique<ComplexVectorMagnitudeCoefficient>(
              *field_coefficients[field_coefficients.size() - 2],
              *field_coefficients[field_coefficients.size() - 1]));
      register_scalar(
          "E_magnitude_" + suffix,
          magnitude_coefficients[magnitude_coefficients.size() - 2].get());
      register_scalar("H_magnitude_" + suffix,
                      magnitude_coefficients.back().get());
      VectorCoefficient &electric_real =
          *field_coefficients[field_coefficients.size() - 4];
      VectorCoefficient &electric_imaginary =
          *field_coefficients[field_coefficients.size() - 3];
      VectorCoefficient &magnetic_real_part =
          *field_coefficients[field_coefficients.size() - 2];
      VectorCoefficient &magnetic_imaginary_part =
          *field_coefficients[field_coefficients.size() - 1];
      electric_energy.emplace_back(
          make_unique<UDARElectricEnergyDensityCoefficient>(
              electric_real, electric_imaginary, permittivity));
      magnetic_energy.emplace_back(
          make_unique<UDARMagneticEnergyDensityCoefficient>(
              magnetic_real_part, magnetic_imaginary_part, permeability));
      ohmic_loss.emplace_back(make_unique<UDAROhmicLossDensityCoefficient>(
          electric_real, electric_imaginary, conductivity));
      poynting.emplace_back(make_unique<UDARTimeAveragePoyntingCoefficient>(
          electric_real, electric_imaginary, magnetic_real_part,
          magnetic_imaginary_part));
      register_scalar("electric_energy_J_m3_" + suffix,
                      electric_energy.back().get());
      register_scalar("magnetic_energy_J_m3_" + suffix,
                      magnetic_energy.back().get());
      register_scalar("ohmic_loss_W_m3_" + suffix,
                      ohmic_loss.back().get());
      register_vector("poynting_W_m2_" + suffix,
                      poynting.back().get());
   }
   StopWatch timer;
   timer.Start();
   collection.Save();
   timer.Stop();
   stats.save_seconds = timer.RealTime();
   return stats;
}

} // namespace mfem
