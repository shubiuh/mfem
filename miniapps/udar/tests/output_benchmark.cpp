// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#include "../udar3p_postprocess.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>

using namespace mfem;

/** @brief Measure the production writer using manufactured ND/RT fields.
    This isolates output cost, not the fraction of an actual Maxwell solve.
    Raw/full runs are alternated to reduce ordering bias. Profiling is measured
    separately so clock overhead does not affect the unprofiled comparison. */
int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   int cells = 16, order = 2, repeats = 3;
   const char *directory = "udar-output-benchmark";
   OptionsParser args(argc, argv);
   args.AddOption(&cells, "-n", "--cells-per-axis", "Hex cells per axis.");
   args.AddOption(&order, "-o", "--order", "Nedelec order.");
   args.AddOption(&repeats, "-r", "--repeats", "Unprofiled repetitions.");
   args.AddOption(&directory, "-out", "--output", "Benchmark output directory.");
   args.Parse();
   if (!args.Good()) { if (Mpi::Root()) { args.PrintUsage(out); } return 1; }
   MFEM_VERIFY(cells > 0 && order > 0 && repeats > 0, "Invalid benchmark size.");
   Mesh serial = Mesh::MakeCartesian3D(cells, cells, cells, Element::HEXAHEDRON);
   serial.EnsureNCMesh(true);
   ParMesh mesh(MPI_COMM_WORLD, serial);
   serial.Clear();
   ND_FECollection nd(order, 3);
   RT_FECollection rt(order - 1, 3);
   ParFiniteElementSpace edge(&mesh, &nd), face(&mesh, &rt);
   ParDiscreteLinearOperator curl(&edge, &face);
   curl.AddDomainInterpolator(new CurlInterpolator);
   curl.Assemble();
   curl.Finalize();
   std::vector<std::unique_ptr<ParComplexGridFunction>> electric;
   std::vector<std::unique_ptr<ParGridFunction>> magnetic_r, magnetic_i;
   for (int axis = 0; axis < 3; axis++)
   {
      electric.emplace_back(new ParComplexGridFunction(&edge));
      VectorFunctionCoefficient real_field(3,
          [axis](const Vector &x, Vector &value)
      {
         value.SetSize(3);
         for (int d = 0; d < 3; d++)
         { value[d] = std::sin((axis + d + 1) * x[d] + x[(d + 1) % 3]); }
      });
      VectorFunctionCoefficient imaginary_field(3,
          [axis](const Vector &x, Vector &value)
      {
         value.SetSize(3);
         for (int d = 0; d < 3; d++)
         { value[d] = std::cos((axis + d + 1) * x[d] + x[(d + 1) % 3]); }
      });
      electric.back()->real().ProjectCoefficient(real_field);
      electric.back()->imag().ProjectCoefficient(imaginary_field);
      magnetic_r.emplace_back(new ParGridFunction(&face));
      magnetic_i.emplace_back(new ParGridFunction(&face));
      curl.Mult(electric.back()->imag(), *magnetic_r.back());
      curl.Mult(electric.back()->real(), *magnetic_i.back());
      *magnetic_r.back() *= -1.0 / (2.0 * M_PI * 48000.0 * 4.0e-7 * M_PI);
      *magnetic_i.back() *= 1.0 / (2.0 * M_PI * 48000.0 * 4.0e-7 * M_PI);
   }
   MatrixFunctionCoefficient sigma(3, [](const Vector &x, DenseMatrix &value)
   {
      value.SetSize(3);
      value = 0.0;
      value(0, 0) = 0.1 + 0.05 * x[0];
      value(1, 1) = 0.2;
      value(2, 2) = 0.3;
      value(0, 1) = value(1, 0) = 0.02;
   });
   const auto elements = mesh.GetGlobalNE();
   const auto unknowns = edge.GlobalTrueVSize();
   std::ofstream report;
   if (Mpi::Root())
   {
      std::filesystem::create_directories(directory);
      report.open(std::string(directory) + "/timing.csv");
      MFEM_VERIFY(report.good(), "Cannot open benchmark report.");
      report << "mode,repeat,ranks,elements,order,complex_unknowns,save_s,"
             << "derived_eval_s,derived_calls\n" << std::setprecision(12);
   }
   MPI_Barrier(mesh.GetComm());
   auto run = [&](const std::string &name, int repetition, bool profile,
                  bool derived)
   {
      MPI_Barrier(mesh.GetComm());
      const auto stats = WriteUDARVTKFields(
          directory, name, 0, 0, 48000.0, mesh, sigma, 8.854187817e-12,
          4.0e-7 * M_PI, electric, magnetic_r, magnetic_i, profile, derived);
      const double local[2] = {stats.save_seconds, stats.derived_evaluation_seconds};
      double maximum[2];
      long long evaluations;
      MPI_Reduce(local, maximum, 2, MPI_DOUBLE, MPI_MAX, 0, mesh.GetComm());
      MPI_Reduce(&stats.derived_evaluations, &evaluations, 1, MPI_LONG_LONG,
                 MPI_SUM, 0, mesh.GetComm());
      if (Mpi::Root())
      {
         report << name << ',' << repetition << ',' << Mpi::WorldSize() << ','
                << elements << ',' << order << ',' << unknowns << ','
                << maximum[0] << ',' << maximum[1] << ',' << evaluations << '\n';
         report.flush();
         out << name << " repeat " << repetition << ": Save " << maximum[0]
             << " s, derived Eval " << maximum[1] << " s\n" << std::flush;
      }
   };
   for (int repeat = 0; repeat < repeats; repeat++)
   {
      if (repeat % 2 == 0) { run("raw", repeat, false, false); }
      run("full", repeat, false, true);
      if (repeat % 2 != 0) { run("raw", repeat, false, false); }
   }
   run("profiled", 0, true, true);
   return 0;
}
