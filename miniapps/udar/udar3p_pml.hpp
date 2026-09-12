// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#ifndef MFEM_UDAR3P_PML_HPP
#define MFEM_UDAR3P_PML_HPP

#include "mfem.hpp"

#include <complex>
#include <string>
#include <vector>

namespace mfem
{

/** @brief Cartesian, polynomially graded perfectly matched layer.

    The coordinate stretch uses @f$s_d=1-i\sigma_d/(\omega\epsilon_0)@f$,
    which is the attenuating sign for the @f$\exp(+i\omega t)@f$ convention.
    The supplied maximum conductivity controls the strength of the artificial
    stretch and is not added to the physical formation conductivity. */
class UDARCartesianPML
{
 private:
   Vector minimum;
   Vector maximum;
   real_t thickness;
   real_t maximum_conductivity;
   int polynomial_order;

 public:
   UDARCartesianPML(const Vector &minimum_, const Vector &maximum_,
                    real_t thickness_, real_t maximum_conductivity_,
                    int polynomial_order_);

   /** @brief Return the normalized grading on one coordinate axis. */
   real_t AxisProfile(const Vector &point, int axis) const;

   /** @brief Return the normalized grading for one side of an axis.
       @param point Physical coordinate.
       @param axis Cartesian axis in the range [0, 2].
       @param upper Select the maximum-coordinate side when true. */
   real_t SideProfile(const Vector &point, int axis, bool upper) const;

   /** @brief Return the largest directional grading at a point. */
   real_t Profile(const Vector &point) const;

   /** @brief Return whether a point is in any PML slab. */
   bool Contains(const Vector &point) const { return Profile(point) > 0.0; }

   /** @brief Return the artificial conductivity represented by the profile. */
   real_t ArtificialConductivity(const Vector &point) const
   { return maximum_conductivity * Profile(point); }

   /** @brief Evaluate the three diagonal coordinate-stretch factors. */
   void Stretch(const Vector &point, real_t omega, real_t epsilon0,
                std::complex<real_t> stretch[3]) const;

   const Vector &Minimum() const { return minimum; }
   const Vector &Maximum() const { return maximum; }
   real_t Thickness() const { return thickness; }
   real_t MaximumConductivity() const { return maximum_conductivity; }
   int PolynomialOrder() const { return polynomial_order; }
};

/** @brief Mesh-resolution diagnostics for one Cartesian PML side. */
struct UDARPMLSideDiagnostics
{
   int axis = 0;
   bool upper = false;
   real_t boundary_coordinate = 0.0;
   real_t interface_coordinate = 0.0;
   long long intersected_elements = 0;
   int conservative_layers = 0;
   real_t maximum_normal_span = 0.0;
   real_t maximum_profile_jump = 0.0;
   real_t maximum_conductivity_jump = 0.0;

   /** @brief Return a stable label such as @c x_min or @c z_max. */
   std::string SideName() const;
};

/** @brief Analyze element resolution independently on all six PML sides.

    The conservative layer count is floor(thickness / maximum normal element
    span). The grading jump is the largest variation of the side-specific
    normalized PML profile across the vertices of any intersecting element.
    Returned values are global across the parallel mesh. */
std::vector<UDARPMLSideDiagnostics> AnalyzeUDARPMLResolution(
    ParMesh &mesh, const UDARCartesianPML &pml);

/** @brief Write one CSV row per Cartesian PML side. */
void WriteUDARPMLResolutionCSV(
    const std::string &path,
    const std::vector<UDARPMLSideDiagnostics> &diagnostics,
    const UDARCartesianPML &pml);

/** @brief Real or imaginary PML-transformed inverse-permeability tensor. */
class UDARPMLCurlCoefficient : public MatrixCoefficient
{
 private:
   const UDARCartesianPML &pml;
   real_t omega;
   real_t epsilon0;
   real_t inverse_mu_r;
   bool imaginary_part;

 public:
   UDARPMLCurlCoefficient(const UDARCartesianPML &pml_, real_t omega_,
                          real_t epsilon0_, real_t inverse_mu_r_,
                          bool imaginary_part_);

   void Eval(DenseMatrix &value, ElementTransformation &T,
             const IntegrationPoint &ip) override;
};

/** @brief Real or imaginary PML-transformed Maxwell mass tensor.

    This coefficient transforms the complete physical tensor
    @f$-k_0^2\epsilon_r I+i\omega\mu_0\sigma@f$, including anisotropic
    formation conductivity, by @f$\det(J)J^{-1}(.)J^{-T}@f$. */
class UDARPMLMassCoefficient : public MatrixCoefficient
{
 private:
   MatrixCoefficient &formation;
   const UDARCartesianPML &pml;
   real_t omega;
   real_t epsilon0;
   real_t dielectric_mass;
   real_t conductivity_scale;
   bool imaginary_part;

 public:
   UDARPMLMassCoefficient(MatrixCoefficient &formation_,
                          const UDARCartesianPML &pml_, real_t omega_,
                          real_t epsilon0_, real_t dielectric_mass_,
                          real_t conductivity_scale_, bool imaginary_part_);

   void Eval(DenseMatrix &value, ElementTransformation &T,
             const IntegrationPoint &ip) override;
};

/** @brief Positive PML curl tensor used by the shifted AMS operator.

    The diagonal entries are the magnitudes of the complex coordinate-
    transformed inverse-permeability tensor. */
class UDARPMLShiftedCurlCoefficient : public MatrixCoefficient
{
 private:
   const UDARCartesianPML &pml;
   real_t omega;
   real_t epsilon0;
   real_t inverse_mu_r;

 public:
   UDARPMLShiftedCurlCoefficient(const UDARCartesianPML &pml_, real_t omega_,
                                 real_t epsilon0_, real_t inverse_mu_r_);

   void Eval(DenseMatrix &value, ElementTransformation &T,
             const IntegrationPoint &ip) override;
};

/** @brief Positive PML mass tensor used by the shifted AMS operator.

    If @f$D=k_0^2\epsilon_r I+\omega\mu_0\sigma@f$ is the positive physical
    mass tensor, this coefficient returns
    @f$\widehat D_{ij}=|\det J|D_{ij}/(|s_i||s_j|)@f$. This is a diagonal
    congruence of @f$D@f$ and therefore remains symmetric positive definite. */
class UDARPMLShiftedMassCoefficient : public MatrixCoefficient
{
 private:
   MatrixCoefficient &formation;
   const UDARCartesianPML &pml;
   real_t omega;
   real_t epsilon0;
   real_t dielectric_mass;
   real_t conductivity_scale;

 public:
   UDARPMLShiftedMassCoefficient(MatrixCoefficient &formation_,
                                 const UDARCartesianPML &pml_, real_t omega_,
                                 real_t epsilon0_, real_t dielectric_mass_,
                                 real_t conductivity_scale_);

   void Eval(DenseMatrix &value, ElementTransformation &T,
             const IntegrationPoint &ip) override;
};

} // namespace mfem

#endif
