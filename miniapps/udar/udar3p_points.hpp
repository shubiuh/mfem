// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#ifndef MFEM_UDAR3P_POINTS_HPP
#define MFEM_UDAR3P_POINTS_HPP

#include "mfem.hpp"

#include <vector>

namespace mfem
{

/** @brief Cache element locations for a fixed set of receiver points.

    Point locations are independent of frequency and material values, so the
    cache can be reused by every forward evaluation on the same mesh. Points
    on a shared face retain all matching elements and are sampled by the same
    arithmetic average used by the original UDAR point evaluator. */
class UDARBatchedPointEvaluator
{
private:
   ParMesh &mesh;
   std::vector<std::vector<int>> local_elements;
   std::vector<std::vector<IntegrationPoint>> local_points;
   std::vector<int> global_matches;

public:
   /** @brief Locate all receiver points and cache their reference coordinates.
       @param mesh_ Parallel mesh used for all subsequent evaluations.
       @param points Physical receiver coordinates.
       @param tolerance Physical inverse-transformation tolerance. */
   UDARBatchedPointEvaluator(ParMesh &mesh_, const std::vector<Vector> &points,
                             real_t tolerance);

   /** @brief Return the number of cached receiver points. */
   int Size() const { return static_cast<int>(global_matches.size()); }

   /** @brief Return whether receiver @a point was found in the global mesh. */
   bool Found(int point) const { return global_matches[point] > 0; }

   /** @brief Evaluate H for all receivers and transmitter solutions at once.
       @param electric_fields Complex electric fields, one per transmitter.
       @param omega Angular frequency in radians per second.
       @param mu0 Vacuum permeability.
       @param mu_r Relative permeability.
       @param field_real Per-point matrices whose columns are transmitter
       fields and rows are global Cartesian components.
       @param field_imaginary Imaginary counterparts of @a field_real.

       The method performs a single MPI reduction for every receiver and
       transmitter combination. */
   void EvaluateMagneticFields(
       const std::vector<const ParComplexGridFunction *> &electric_fields,
       real_t omega, real_t mu0, real_t mu_r,
       std::vector<DenseMatrix> &field_real,
       std::vector<DenseMatrix> &field_imaginary) const;
};

} // namespace mfem

#endif
