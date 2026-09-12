// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#include "udar3p_pml.hpp"

#include <algorithm>
#include <cmath>

namespace mfem
{
namespace
{

void GetPhysicalPoint(ElementTransformation &T, const IntegrationPoint &ip,
                      Vector &point)
{
   point.SetSize(3);
   T.Transform(ip, point);
}

} // namespace

UDARCartesianPML::UDARCartesianPML(const Vector &minimum_,
                                   const Vector &maximum_, real_t thickness_,
                                   real_t maximum_conductivity_,
                                   int polynomial_order_)
    : minimum(minimum_), maximum(maximum_), thickness(thickness_),
      maximum_conductivity(maximum_conductivity_),
      polynomial_order(polynomial_order_)
{
   MFEM_VERIFY(minimum.Size() == 3 && maximum.Size() == 3,
               "The PML domain must be three-dimensional.");
   MFEM_VERIFY(thickness > 0.0 && maximum_conductivity > 0.0,
               "PML thickness and maximum conductivity must be positive.");
   MFEM_VERIFY(polynomial_order >= 1,
               "PML polynomial order must be at least one.");
   for (int d = 0; d < 3; d++)
   {
      MFEM_VERIFY(maximum[d] - minimum[d] > 2.0 * thickness,
                  "PML thickness leaves no physical interior along an axis.");
   }
}

real_t UDARCartesianPML::AxisProfile(const Vector &point, int axis) const
{
   MFEM_ASSERT(axis >= 0 && axis < 3, "Invalid PML coordinate axis.");
   const real_t distance =
       std::min(point[axis] - minimum[axis], maximum[axis] - point[axis]);
   if (distance >= thickness) { return 0.0; }
   const real_t coordinate =
       std::max(real_t(0.0), (thickness - distance) / thickness);
   return pow(coordinate, polynomial_order);
}

real_t UDARCartesianPML::Profile(const Vector &point) const
{
   real_t profile = 0.0;
   for (int d = 0; d < 3; d++)
   {
      profile = std::max(profile, AxisProfile(point, d));
   }
   return profile;
}

void UDARCartesianPML::Stretch(const Vector &point, real_t omega,
                               real_t epsilon0,
                               std::complex<real_t> stretch[3]) const
{
   MFEM_ASSERT(omega > 0.0 && epsilon0 > 0.0,
               "PML stretch requires positive frequency and permittivity.");
   const std::complex<real_t> imaginary_unit(0.0, 1.0);
   for (int d = 0; d < 3; d++)
   {
      const real_t sigma = maximum_conductivity * AxisProfile(point, d);
      stretch[d] = 1.0 - imaginary_unit * sigma / (omega * epsilon0);
   }
}

UDARPMLCurlCoefficient::UDARPMLCurlCoefficient(const UDARCartesianPML &pml_,
                                               real_t omega_, real_t epsilon0_,
                                               real_t inverse_mu_r_,
                                               bool imaginary_part_)
    : MatrixCoefficient(3), pml(pml_), omega(omega_), epsilon0(epsilon0_),
      inverse_mu_r(inverse_mu_r_), imaginary_part(imaginary_part_)
{
}

void UDARPMLCurlCoefficient::Eval(DenseMatrix &value, ElementTransformation &T,
                                  const IntegrationPoint &ip)
{
   Vector point;
   GetPhysicalPoint(T, ip, point);
   std::complex<real_t> stretch[3];
   pml.Stretch(point, omega, epsilon0, stretch);
   const std::complex<real_t> determinant =
       stretch[0] * stretch[1] * stretch[2];

   value.SetSize(3);
   value = 0.0;
   for (int d = 0; d < 3; d++)
   {
      const std::complex<real_t> entry =
          inverse_mu_r * stretch[d] * stretch[d] / determinant;
      value(d, d) = imaginary_part ? entry.imag() : entry.real();
   }
}

UDARPMLMassCoefficient::UDARPMLMassCoefficient(MatrixCoefficient &formation_,
                                               const UDARCartesianPML &pml_,
                                               real_t omega_, real_t epsilon0_,
                                               real_t dielectric_mass_,
                                               real_t conductivity_scale_,
                                               bool imaginary_part_)
    : MatrixCoefficient(3), formation(formation_), pml(pml_), omega(omega_),
      epsilon0(epsilon0_), dielectric_mass(dielectric_mass_),
      conductivity_scale(conductivity_scale_), imaginary_part(imaginary_part_)
{
}

void UDARPMLMassCoefficient::Eval(DenseMatrix &value, ElementTransformation &T,
                                  const IntegrationPoint &ip)
{
   DenseMatrix sigma;
   formation.Eval(sigma, T, ip);
   Vector point;
   GetPhysicalPoint(T, ip, point);
   std::complex<real_t> stretch[3];
   pml.Stretch(point, omega, epsilon0, stretch);
   const std::complex<real_t> determinant =
       stretch[0] * stretch[1] * stretch[2];

   value.SetSize(3);
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         std::complex<real_t> physical(0.0, conductivity_scale * sigma(i, j));
         if (i == j) { physical += dielectric_mass; }
         const std::complex<real_t> transformed =
             determinant * physical / (stretch[i] * stretch[j]);
         value(i, j) = imaginary_part ? transformed.imag() : transformed.real();
      }
   }
}

UDARPMLShiftedCurlCoefficient::UDARPMLShiftedCurlCoefficient(
    const UDARCartesianPML &pml_, real_t omega_, real_t epsilon0_,
    real_t inverse_mu_r_)
    : MatrixCoefficient(3), pml(pml_), omega(omega_), epsilon0(epsilon0_),
      inverse_mu_r(inverse_mu_r_)
{
}

void UDARPMLShiftedCurlCoefficient::Eval(DenseMatrix &value,
                                         ElementTransformation &T,
                                         const IntegrationPoint &ip)
{
   Vector point;
   GetPhysicalPoint(T, ip, point);
   std::complex<real_t> stretch[3];
   pml.Stretch(point, omega, epsilon0, stretch);
   const real_t determinant = std::abs(stretch[0] * stretch[1] * stretch[2]);

   value.SetSize(3);
   value = 0.0;
   for (int d = 0; d < 3; d++)
   {
      value(d, d) = inverse_mu_r * std::norm(stretch[d]) / determinant;
   }
}

UDARPMLShiftedMassCoefficient::UDARPMLShiftedMassCoefficient(
    MatrixCoefficient &formation_, const UDARCartesianPML &pml_, real_t omega_,
    real_t epsilon0_, real_t dielectric_mass_, real_t conductivity_scale_)
    : MatrixCoefficient(3), formation(formation_), pml(pml_), omega(omega_),
      epsilon0(epsilon0_), dielectric_mass(dielectric_mass_),
      conductivity_scale(conductivity_scale_)
{
}

void UDARPMLShiftedMassCoefficient::Eval(DenseMatrix &value,
                                         ElementTransformation &T,
                                         const IntegrationPoint &ip)
{
   formation.Eval(value, T, ip);
   value *= conductivity_scale;
   for (int d = 0; d < 3; d++) { value(d, d) += dielectric_mass; }

   Vector point;
   GetPhysicalPoint(T, ip, point);
   std::complex<real_t> stretch[3];
   pml.Stretch(point, omega, epsilon0, stretch);
   const real_t determinant = std::abs(stretch[0] * stretch[1] * stretch[2]);
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         value(i, j) *=
             determinant / (std::abs(stretch[i]) * std::abs(stretch[j]));
      }
   }
}

} // namespace mfem
