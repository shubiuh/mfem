// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#include "../udar3p_frequency.hpp"

#include <cmath>
#include <string>

using namespace mfem;

/** @brief Compare full and partial complex operator actions, including PEC
    rows, anisotropic conductivity, material updates and hanging edges. */
int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   for (int order : {2, 3})
   {
      Mesh serial = Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON);
      serial.EnsureNCMesh(true);
      Array<int> marked(1);
      marked[0] = 0;
      serial.GeneralRefinement(marked, 1);
      ParMesh mesh(MPI_COMM_WORLD, serial);
      ND_FECollection collection(order, 3);
      ParFiniteElementSpace space(&mesh, &collection);
      Array<int> boundary(mesh.bdr_attributes.Max()), essential;
      boundary = 1;
      space.GetEssentialTrueDofs(boundary, essential);
      DenseMatrix sigma(3);
      sigma = 0.0;
      sigma(0, 0) = 0.1;
      sigma(1, 1) = 0.2;
      sigma(2, 2) = 0.3;
      sigma(0, 1) = sigma(1, 0) = 0.02;
      MatrixFunctionCoefficient conductivity(3,
          [&sigma](const Vector &, DenseMatrix &value) { value = sigma; });
      UDARFrequencySystemOptions options;
      options.frequency = 48000.0;
      options.epsilon0 = 8.854187817e-12;
      options.mu0 = 4.0e-7 * M_PI;
      options.iterative.preconditioner = "ams";
      options.iterative.ams_cycle_type = 13;
      options.iterative.ilu_fallback = false;
      options.assembly = "partial";
      UDARFrequencySystem partial(space, conductivity, nullptr, essential, options);
      MFEM_VERIFY(std::string(partial.AssemblyMode()) == "partial",
                  "Compatible high-order configuration must use PA.");
      options.assembly = "full";
      UDARFrequencySystem full(space, conductivity, nullptr, essential, options);
      for (int update = 0; update < 2; update++)
      {
         sigma *= 2.0;
         partial.AssembleOperator(0.1);
         full.AssembleOperator(0.1);
         MFEM_VERIFY(!dynamic_cast<const ComplexHypreParMatrix *>(
                         &partial.GetOperator()),
                     "PA unexpectedly materialized a complex sparse matrix.");
         Vector x(2 * space.GetTrueVSize()), y(x.Size()), reference(x.Size());
         for (int i = 0; i < x.Size(); i++)
         { x[i] = std::sin(real_t(i + 1 + Mpi::WorldRank())); }
         partial.GetOperator().Mult(x, y);
         full.GetOperator().Mult(x, reference);
         y -= reference;
         const real_t relative_error = std::sqrt(
             InnerProduct(MPI_COMM_WORLD, y, y) /
             InnerProduct(MPI_COMM_WORLD, reference, reference));
         MFEM_VERIFY(relative_error < 1.0e-11,
                     "Full and partial Maxwell actions differ.");
         if (Mpi::Root())
         { out << "Order " << order << ", update " << update
               << ": relative operator error " << relative_error << '\n'; }
      }
      for (int incompatible = 0; incompatible < 5; incompatible++)
      {
         auto guarded = options;
         guarded.assembly = "partial";
         if (incompatible == 0) { guarded.iterative.preconditioner = "ilu"; }
         if (incompatible == 1) { guarded.iterative.system_mode = "real"; }
         if (incompatible == 2) { guarded.iterative.ilu_fallback = true; }
         if (incompatible == 3) { guarded.boundary_condition = 1; }
         if (incompatible == 4) { guarded.iterative.ams_cycle_type = 14; }
         UDARFrequencySystem system(space, conductivity, nullptr, essential,
                                    guarded);
         MFEM_VERIFY(std::string(system.AssemblyMode()) == "full" &&
                     !system.AssemblyReason().empty(),
                     "Incompatible configuration bypassed the assembly guard.");
      }
      Vector minimum({0.0, 0.0, 0.0}), maximum({1.0, 1.0, 1.0});
      UDARCartesianPML pml(minimum, maximum, 0.15, 1.0, 2);
      options.assembly = "partial";
      UDARFrequencySystem pml_system(space, conductivity, &pml, essential, options);
      MFEM_VERIFY(std::string(pml_system.AssemblyMode()) == "full",
                  "PML must retain full assembly.");
      ND_FECollection low_collection(1, 3);
      ParFiniteElementSpace low_space(&mesh, &low_collection);
      UDARFrequencySystem low_system(low_space, conductivity, nullptr,
                                     Array<int>(), options);
      MFEM_VERIFY(std::string(low_system.AssemblyMode()) == "full",
                  "Order one must retain full assembly.");
      Mesh tet_serial = Mesh::MakeCartesian3D(2, 2, 2, Element::TETRAHEDRON);
      tet_serial.EnsureNCMesh(true);
      ParMesh tet_mesh(MPI_COMM_WORLD, tet_serial);
      ParFiniteElementSpace tet_space(&tet_mesh, &collection);
      UDARFrequencySystem tet_system(tet_space, conductivity, nullptr,
                                     Array<int>(), options);
      MFEM_VERIFY(std::string(tet_system.AssemblyMode()) == "full",
                  "Tetrahedra must retain full assembly.");
   }
   return 0;
}
