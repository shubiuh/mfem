// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#include "../udar3p_adapt.hpp"

#include <algorithm>
#include <cmath>

using namespace mfem;

/** @brief Check nested ND interpolation and transpose restriction on hex/tet.

    Run with one or more MPI ranks. The localized refinement also exercises
    collective participation by ranks without locally marked cells. */
int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   for (Element::Type type : {Element::HEXAHEDRON, Element::TETRAHEDRON})
   {
      Mesh serial = Mesh::MakeCartesian3D(4, 4, 4, type);
      serial.EnsureNCMesh(true);
      ND_FECollection collection(1, 3);
      auto *mesh = new ParMesh(MPI_COMM_WORLD, serial);
      auto *space = new ParFiniteElementSpace(mesh, &collection);
      ParFiniteElementSpaceHierarchy hierarchy(mesh, space, true, true);
      const auto coarse_elements = mesh->GetGlobalNE();
      Vector start({0.03, 0.04, 0.05}), end({0.08, 0.04, 0.05});
      Vector axis({1.0, 0.0, 0.0});
      const int passes = RefineUDARToolRegion(
          hierarchy, start, end, axis, 0.01, 0.015, 0.0, 2);
      MFEM_VERIFY(passes == 2, "Expected two retained local-refinement levels.");
      MFEM_VERIFY(mesh->GetGlobalNE() == coarse_elements,
                  "Local refinement modified the retained coarse mesh.");
      hierarchy.AddUniformlyRefinedLevel();
      MFEM_VERIFY(hierarchy.GetNumLevels() == 4, "Missing hierarchy levels.");
      VectorFunctionCoefficient field(3, [](const Vector &x, Vector &value)
      {
         value.SetSize(3);
         value[0] = 1.0 - x[1];
         value[1] = 2.0 + x[0];
         value[2] = 3.0;
      });
      for (int level = 0; level + 1 < hierarchy.GetNumLevels(); level++)
      {
         auto &coarse = hierarchy.GetFESpaceAtLevel(level);
         auto &fine = hierarchy.GetFESpaceAtLevel(level + 1);
         Operator &transfer = *hierarchy.GetProlongationAtLevel(level);
         ParGridFunction coarse_field(&coarse), fine_field(&fine);
         coarse_field.ProjectCoefficient(field);
         Vector x, y(fine.GetTrueVSize());
         coarse_field.GetTrueDofs(x);
         transfer.Mult(x, y);
         fine_field.SetFromTrueDofs(y);
         const real_t error = fine_field.ComputeL2Error(field);
         MFEM_VERIFY(error < 1.0e-10, "ND transfer lost a representable field.");

         Vector probe(y.Size()), restricted(x.Size());
         for (int i = 0; i < probe.Size(); i++)
         {
            probe[i] = std::sin(real_t(i + 1 + Mpi::WorldRank()));
         }
         transfer.MultTranspose(probe, restricted);
         const real_t lhs = InnerProduct(MPI_COMM_WORLD, y, probe);
         const real_t rhs = InnerProduct(MPI_COMM_WORLD, x, restricted);
         MFEM_VERIFY(std::abs(lhs - rhs) <
                     1.0e-10 * std::max(real_t(1.0), std::abs(lhs)),
                     "Hierarchy restriction is not the transpose of prolongation.");
         if (Mpi::Root())
         {
            out << (type == Element::HEXAHEDRON ? "hex" : "tet")
                << " level " << level << " -> " << level + 1
                << ": field error " << error
                << ", transpose error " << std::abs(lhs - rhs) << '\n';
         }
      }
   }
   return 0;
}
