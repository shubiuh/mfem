// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#include "udar3p_adapt.hpp"

#include <algorithm>
#include <cmath>

namespace mfem
{
namespace
{

real_t Dot3(const Vector &a, const Vector &b)
{ return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

real_t PointSegmentDistance(const Vector &point, const Vector &start,
                            const Vector &end)
{
   Vector direction(end);
   direction -= start;
   const real_t length_squared = Dot3(direction, direction);
   real_t fraction = 0.0;
   if (length_squared > 0.0)
   {
      Vector offset(point);
      offset -= start;
      fraction = std::max(real_t(0.0), std::min(real_t(1.0),
                                                 Dot3(offset, direction) /
                                                     length_squared));
   }
   Vector closest(start);
   closest.Add(fraction, direction);
   closest -= point;
   return closest.Norml2();
}

} // namespace

int RefineUDARToolRegion(Mesh &mesh, const Vector &survey_start,
                         const Vector &survey_end, const Vector &tool_axis,
                         real_t maximum_spacing, real_t core_radius,
                         real_t padding, int levels)
{
   MFEM_VERIFY(levels >= 1, "Tool-region refinement needs at least one level.");
   MFEM_VERIFY(core_radius > 0.0 && padding >= 0.0,
               "Tool-region radius must be positive and padding nonnegative.");
   for (int element = 0; element < mesh.GetNE(); element++)
   {
      const Geometry::Type geometry = mesh.GetElementBaseGeometry(element);
      MFEM_VERIFY(geometry == Geometry::CUBE || geometry == Geometry::TETRAHEDRON,
                  "Tool-region refinement requires hexahedral or tetrahedral "
                  "elements.");
   }

   Vector refinement_start(survey_start);
   refinement_start.Add(-(maximum_spacing + padding), tool_axis);
   Vector refinement_end(survey_end);
   refinement_end.Add(padding, tool_axis);
   mesh.EnsureNCMesh();

   int passes = 0;
   for (int level = 0; level < levels; level++)
   {
      const real_t level_radius =
          core_radius * pow(2.0, static_cast<real_t>(levels - level - 1)) +
          padding;
      Array<int> elements;
      for (int element = 0; element < mesh.GetNE(); element++)
      {
         Vector center;
         mesh.GetElementCenter(element, center);
         // The enclosing-sphere bound is deliberately conservative for both
         // element types, so every cell touching the refinement tube is kept.
         const real_t element_radius =
             0.5 * sqrt(3.0) * mesh.GetElementSize(element, 2);
         if (PointSegmentDistance(center, refinement_start, refinement_end) <=
             level_radius + element_radius)
         {
            elements.Append(element);
         }
      }
      if (elements.Size() == 0) { break; }
      mesh.GeneralRefinement(elements, 1);
      passes++;
   }
   return passes;
}

} // namespace mfem
