// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#include "udar3p_pml.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>

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
   return std::max(SideProfile(point, axis, false),
                   SideProfile(point, axis, true));
}

real_t UDARCartesianPML::SideProfile(const Vector &point, int axis,
                                     bool upper) const
{
   MFEM_ASSERT(axis >= 0 && axis < 3, "Invalid PML coordinate axis.");
   const real_t distance = upper ? maximum[axis] - point[axis]
                                 : point[axis] - minimum[axis];
   if (distance >= thickness) { return 0.0; }
   const real_t coordinate = std::min(
       real_t(1.0),
       std::max(real_t(0.0), (thickness - distance) / thickness));
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

std::string UDARPMLSideDiagnostics::SideName() const
{
   static const char axis_name[3] = {'x', 'y', 'z'};
   return std::string(1, axis_name[axis]) + (upper ? "_max" : "_min");
}

std::vector<UDARPMLSideDiagnostics> AnalyzeUDARPMLResolution(
    ParMesh &mesh, const UDARCartesianPML &pml)
{
   std::vector<UDARPMLSideDiagnostics> diagnostics(6);
   for (int axis = 0; axis < 3; axis++)
   {
      for (int side = 0; side < 2; side++)
      {
         UDARPMLSideDiagnostics &entry = diagnostics[2 * axis + side];
         entry.axis = axis;
         entry.upper = side == 1;
         entry.boundary_coordinate =
             entry.upper ? pml.Maximum()[axis] : pml.Minimum()[axis];
         entry.interface_coordinate =
             entry.boundary_coordinate +
             (entry.upper ? -pml.Thickness() : pml.Thickness());

         long long local_elements = 0;
         real_t local_maximum_span = 0.0;
         real_t local_maximum_jump = 0.0;
         for (int element = 0; element < mesh.GetNE(); element++)
         {
            Array<int> vertices;
            mesh.GetElementVertices(element, vertices);
            real_t element_minimum =
                std::numeric_limits<real_t>::infinity();
            real_t element_maximum =
                -std::numeric_limits<real_t>::infinity();
            real_t profile_minimum =
                std::numeric_limits<real_t>::infinity();
            real_t profile_maximum = 0.0;
            for (int vertex : vertices)
            {
               Vector point(3);
               for (int d = 0; d < 3; d++)
               {
                  point[d] = mesh.GetVertex(vertex)[d];
               }
               element_minimum = std::min(element_minimum, point[axis]);
               element_maximum = std::max(element_maximum, point[axis]);
               const real_t profile =
                   pml.SideProfile(point, axis, entry.upper);
               profile_minimum = std::min(profile_minimum, profile);
               profile_maximum = std::max(profile_maximum, profile);
            }
            const bool intersects =
                entry.upper
                    ? element_maximum > entry.interface_coordinate
                    : element_minimum < entry.interface_coordinate;
            if (!intersects) { continue; }
            local_elements++;
            local_maximum_span =
                std::max(local_maximum_span,
                         element_maximum - element_minimum);
            local_maximum_jump =
                std::max(local_maximum_jump,
                         profile_maximum - profile_minimum);
         }

         MPI_Allreduce(&local_elements, &entry.intersected_elements, 1,
                       MPI_LONG_LONG, MPI_SUM, mesh.GetComm());
         MPI_Allreduce(&local_maximum_span, &entry.maximum_normal_span, 1,
                       MFEM_MPI_REAL_T, MPI_MAX, mesh.GetComm());
         MPI_Allreduce(&local_maximum_jump, &entry.maximum_profile_jump, 1,
                       MFEM_MPI_REAL_T, MPI_MAX, mesh.GetComm());
         if (entry.maximum_normal_span > 0.0)
         {
            const real_t ratio =
                pml.Thickness() / entry.maximum_normal_span;
            entry.conservative_layers =
                static_cast<int>(std::floor(ratio + 1.0e-10));
         }
         entry.maximum_conductivity_jump =
             pml.MaximumConductivity() * entry.maximum_profile_jump;
      }
   }
   return diagnostics;
}

void WriteUDARPMLResolutionCSV(
    const std::string &path,
    const std::vector<UDARPMLSideDiagnostics> &diagnostics,
    const UDARCartesianPML &pml)
{
   std::ofstream output(path);
   MFEM_VERIFY(output.good(), "Unable to open PML diagnostics " << path);
   output << std::setprecision(9) << std::scientific
          << "side,axis,boundary_coordinate,interface_coordinate,"
          << "thickness,polynomial_order,intersected_elements,"
          << "conservative_layers,max_normal_span,max_profile_jump,"
          << "max_conductivity_jump_S_m\n";
   for (const UDARPMLSideDiagnostics &entry : diagnostics)
   {
      output << entry.SideName() << ',' << entry.axis << ','
             << entry.boundary_coordinate << ',' << entry.interface_coordinate
             << ',' << pml.Thickness() << ',' << pml.PolynomialOrder() << ','
             << entry.intersected_elements << ',' << entry.conservative_layers
             << ',' << entry.maximum_normal_span << ','
             << entry.maximum_profile_jump << ','
             << entry.maximum_conductivity_jump << '\n';
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
