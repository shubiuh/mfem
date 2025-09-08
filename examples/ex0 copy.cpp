#include <iostream>

#include <mfem.hpp>
#include "mkl.h"
#include <cstdlib> // For exit()

// Custom assert-like macro
#define MY_REQUIRE(cond)                                                                           \
  if (!(cond)) {                                                                                   \
    std::cerr << "Check failed: " << #cond << " at " << __FILE__ << ":" << __LINE__ << std::endl;  \
    std::exit(EXIT_FAILURE);                                                                       \
  } else {                                                                                         \
    std::cout << "Check passed: " << #cond << " at " << __FILE__ << ":" << __LINE__ << std::endl;  \
  }

using namespace mfem;

double uexact(const Vector &x) {
  double u;
  switch (x.Size()) {
  case 1:
    u = 3.0 + 2.0 * x(0) - 0.5 * x(0) * x(0);
    break;
  case 2:
    u = 1.0 + 0.2 * x(0) - 0.9 * x(0) * x(1) + x(1) * x(1) * x(0);
    break;
  default:
    u = x(2) * x(2) * x(2) - 5.0 * x(0) * x(0) * x(1) * x(2);
    break;
  }
  return u;
}

void gradexact(const Vector &x, Vector &grad) {
  grad.SetSize(x.Size());
  switch (x.Size()) {
  case 1:
    grad[0] = 2.0 - x(0);
    break;
  case 2:
    grad[0] = 0.2 - 0.9 * x(1) + x(1) * x(1);
    grad[1] = -0.9 * x(0) + 2.0 * x(0) * x(1);
    break;
  default:
    grad[0] = -10.0 * x(0) * x(1) * x(2);
    grad[1] = -5.0 * x(0) * x(0) * x(2);
    grad[2] = 3.0 * x(2) * x(2) - 5.0 * x(0) * x(0) * x(1);
    break;
  }
}

double d2uexact(const Vector &x) // returns \Delta u
{
  double d2u;
  switch (x.Size()) {
  case 1:
    d2u = -1.0;
    break;
  case 2:
    d2u = 2.0 * x(0);
    break;
  default:
    d2u = -10.0 * x(1) * x(2) + 6.0 * x(2);
    break;
  }
  return d2u;
}

double fexact(const Vector &x) // returns -\Delta u
{
  double d2u = d2uexact(x);
  return -d2u;
}

int main(int argc, char *argv[]) {
  std::cout << mkl_get_max_threads() << std::endl; // in VS2022, check properties->intel library for
                                                   // OneAPI->Use oneMKL (Parallel)
  mkl_set_dynamic(0);
  mkl_set_num_threads(12);
  std::cout << mkl_get_max_threads() << std::endl; // in VS2022, check properties->intel library for
                                                   // OneAPI->Use oneMKL (Parallel)
  const int ne = 2;
  StopWatch sw; // Start the timer
  for (int dim = 1; dim <= 3; ++dim) {
    sw.Start();
    Mesh mesh;
    if (dim == 1) {
      mesh = Mesh::MakeCartesian1D(ne, 1.0);
    } else if (dim == 2) {
      mesh = Mesh::MakeCartesian2D(ne, ne, Element::QUADRILATERAL, 1, 1.0, 1.0);
    } else {
      mesh = Mesh::MakeCartesian3D(ne, ne, ne, Element::HEXAHEDRON, 1.0, 1.0, 1.0);
    }
    for (unsigned int i = 0; i < 3; ++i) {
      mesh.UniformRefinement();
      out << "Refined mesh has " << mesh.GetNE() << " elements." << std::endl;
    }
    int order = 2;
    H1_FECollection fec(order, dim);
    FiniteElementSpace fespace(&mesh, &fec);
    Array<int> ess_tdof_list, ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 1;
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    out << "Number of finite element unknowns: " << fespace.GetTrueVSize() << std::endl;

    FunctionCoefficient f(fexact);
    LinearForm b(&fespace);
    b.AddDomainIntegrator(new DomainLFIntegrator(f));
    b.Assemble();

    BilinearForm a(&fespace);
    ConstantCoefficient one(1.0);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));
    a.Assemble();

    GridFunction x(&fespace);
    FunctionCoefficient uex(uexact);
    x = 0.0;
    x.ProjectBdrCoefficient(uex, ess_bdr);

    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

// #ifdef MFEM_USE_SUITESPARSE
//       {
//          UMFPackSolver umf_solver;
//          umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
//          umf_solver.SetOperator(*A);
//          umf_solver.Mult(B, X);
//
//          Vector Y(X.Size());
//          A->Mult(X, Y);
//          Y -= B;
//          assert(Y.Norml2() < 1.e-12);
//
//          a.RecoverFEMSolution(X, b, x);
//          VectorFunctionCoefficient grad(dim, gradexact);
//          double error = x.ComputeH1Error(&uex, &grad);
//          MY_REQUIRE(error < 1.e-12);
//       }
// #endif
#ifdef MFEM_USE_MKL_PARDISO
    {
      PardisoSolver pardiso_solver;
      pardiso_solver.SetPrintLevel(1);
      pardiso_solver.SetOperator(*A);
      pardiso_solver.Mult(B, X);

      Vector Y(X.Size());
      A->Mult(X, Y);
      Y -= B;
      MY_REQUIRE(Y.Norml2() < 1.e-12);

      a.RecoverFEMSolution(X, b, x);
      VectorFunctionCoefficient grad(dim, gradexact);
      double error = x.ComputeH1Error(&uex, &grad);
      MY_REQUIRE(error < 1.e-12);
    }
#endif
    sw.Stop();
    out << "Elapsed time for dimension " << dim << ": " << sw.RealTime() << " seconds."
        << std::endl;
  }
}