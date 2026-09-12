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

/** @brief Retain parallel local-refinement levels and H(curl) transfers.

    The hierarchy must start with a nonconforming parallel mesh and a Nedelec
    space. Each pass copies the finest mesh before refining it, preserving its
    parent level and partition. All ranks participate, including ranks with no
    marked elements. The hierarchy owns the added meshes, spaces and true-DOF
    prolongations; its finite element collection must outlive the hierarchy.
    No repartitioning is performed between levels. */
int RefineUDARToolRegion(ParFiniteElementSpaceHierarchy &hierarchy,
                         const Vector &survey_start, const Vector &survey_end,
                         const Vector &tool_axis, real_t maximum_spacing,
                         real_t core_radius, real_t padding, int levels);

} // namespace mfem

#endif
