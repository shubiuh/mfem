// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#ifndef MFEM_UDAR3P_ADAPT_HPP
#define MFEM_UDAR3P_ADAPT_HPP

#include "mfem.hpp"

namespace mfem
{

/** @brief Refine hexahedral or tetrahedral cells near a UDAR tool trajectory.

    The serial mesh is converted to an MFEM nonconforming mesh and refined in
    nested distance bands. This works for both cube and tetrahedron elements;
    MFEM supplies the appropriate child-element topology. */
int RefineUDARToolRegion(Mesh &mesh, const Vector &survey_start,
                         const Vector &survey_end, const Vector &tool_axis,
                         real_t maximum_spacing, real_t core_radius,
                         real_t padding, int levels);

} // namespace mfem

#endif
