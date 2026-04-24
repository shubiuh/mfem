#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

template <typename T> T pow2(const T &x) { return x*x; }

enum bc_type
{
   pec_bc  = 0,
   abc1_bc = 1
};

// ---------------------------------------------------------------------------
// Problem data
// ---------------------------------------------------------------------------
static int dim = 0;
static real_t mu = 1.0;
static real_t epsilon = 1.0;
static real_t sigma = 0.0;
static real_t omega = 0.0;

// phasor_sign = +1  -> Engineer convention: + i*omega*sigma
// phasor_sign = -1  -> e^{-i omega t} convention (physics): - i*omega*sigma
static int phasor_sign = -1;

// Background properties used by the 1st-order absorbing boundary condition.
static real_t mu_bdr = 1.0;
static real_t epsilon_bdr = 1.0;

void source(const Vector &x, Vector &f)
{
   // Simple localized electric-field RHS for testing / prototyping.
   // If you want a physical impressed electric current J under the e^{-iwt}
   // convention, replace this by f = -i*omega*J (or +i*omega*J for e^{+iwt}).
   Vector center(dim);
   center = 0.0;

   real_t r2 = 0.0;
   for (int i = 0; i < dim; i++)
   {
      r2 += pow2(x(i) - center(i));
   }

   const real_t a = 16.0; // Gaussian width parameter
   const real_t amp = exp(-a * r2);

   f.SetSize(dim);
   f = 0.0;
   f(0) = amp;
}

int main(int argc, char *argv[])
{
   const char *mesh_file = "../data/inline-hex.mesh";
   int order = 1;
   int ref_levels = 2;
   real_t freq = 1.0;
   int bc = (int)abc1_bc;
   bool herm_conv = true;
   bool umf_solver = false;
   bool visualization = true;
   bool pa = false;
   const char *device_config = "cpu";

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&ref_levels, "-ref", "--refinements",
                  "Number of uniform refinements.");
   args.AddOption(&mu, "-mu", "--permeability", "Magnetic permeability.");
   args.AddOption(&epsilon, "-eps", "--permittivity", "Electric permittivity.");
   args.AddOption(&sigma, "-sigma", "--conductivity", "Electric conductivity.");
   args.AddOption(&mu_bdr, "-mu-bdr", "--boundary-permeability",
                  "Boundary/background permeability used by the 1st-order ABC.");
   args.AddOption(&epsilon_bdr, "-eps-bdr", "--boundary-permittivity",
                  "Boundary/background permittivity used by the 1st-order ABC.");
   args.AddOption(&freq, "-f", "--frequency", "Frequency in Hz.");
   args.AddOption(&bc, "-bc", "--boundary-condition",
                  "Boundary condition: 0 = PEC, 1 = 1st-order absorbing (Silver-Muller).");
   args.AddOption(&phasor_sign, "-phasor", "--phasor-sign",
                  "Imaginary-term sign: +1 for +i*omega*(...), -1 for -i*omega*(...).");
   args.AddOption(&herm_conv, "-herm", "--hermitian", "-no-herm", "--no-hermitian",
                  "Use convention for Hermitian operators.");
#ifdef MFEM_USE_SUITESPARSE
   args.AddOption(&umf_solver, "-umf", "--umfpack", "-no-umf", "--no-umfpack",
                  "Use the UMFPack solver.");
#endif
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis", "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa", "--no-partial-assembly",
                  "Enable Partial Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   MFEM_VERIFY(phasor_sign == 1 || phasor_sign == -1,
               "phasor_sign must be +1 or -1.");
   MFEM_VERIFY(bc == (int)pec_bc || bc == (int)abc1_bc,
               "bc must be 0 (PEC) or 1 (ABC1).");

   Device device(device_config);
   device.Print();

   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   dim = mesh->Dimension();
   omega = 2.0 * M_PI * freq;

   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }

   FiniteElementCollection *fec = new ND_FECollection(order, dim);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << endl;

   Array<int> ess_bdr;
   Array<int> ess_tdof_list;
   if (mesh->bdr_attributes.Size())
   {
      ess_bdr.SetSize(mesh->bdr_attributes.Max());
      ess_bdr = 0;
      if (bc == (int)pec_bc)
      {
         ess_bdr = 1;
      }
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   ComplexOperator::Convention conv =
      herm_conv ? ComplexOperator::HERMITIAN : ComplexOperator::BLOCK_SYMMETRIC;

   VectorFunctionCoefficient f_rhs(dim, source);
   ComplexLinearForm b(fespace, conv);
   b.AddDomainIntegrator(nullptr, new VectorFEDomainLFIntegrator(f_rhs));
   b.Assemble();

   ComplexGridFunction x(fespace);
   x = 0.0;

   // ------------------------------------------------------------------------
   // curl(mu^{-1} curl E) - omega^2 * epsilon * E
   //     + i * phasor_sign * omega * sigma * E = f
   // with either:
   //   PEC:  n x E = 0
   // or
   //   ABC1: n x H = Y * n x (n x E),  Y = sqrt(epsilon_bdr/mu_bdr)
   //         which yields an additional boundary loss term
   //         + i * phasor_sign * omega * Y * <E_t, F_t>_Gamma
   // ------------------------------------------------------------------------
   const real_t imag_loss = phasor_sign * omega * sigma;
   const real_t eta_inv = sqrt(epsilon_bdr / mu_bdr);
   const real_t imag_abc = phasor_sign * omega * eta_inv;

   ConstantCoefficient muinv(1.0 / mu);
   ConstantCoefficient mass_re(-pow2(omega) * epsilon);
   ConstantCoefficient mass_im(imag_loss);
   ConstantCoefficient abc_im(imag_abc);

   SesquilinearForm a(fespace, conv);
   a.AddDomainIntegrator(new CurlCurlIntegrator(muinv), nullptr);
   a.AddDomainIntegrator(new VectorFEMassIntegrator(mass_re),
                         new VectorFEMassIntegrator(mass_im));

   Array<int> abc_bdr;
   if (bc == (int)abc1_bc && mesh->bdr_attributes.Size())
   {
      abc_bdr.SetSize(mesh->bdr_attributes.Max());
      abc_bdr = 1;
      a.AddBoundaryIntegrator(nullptr,
                              new VectorFEMassIntegrator(abc_im),
                              abc_bdr);
   }

   if (pa) { a.SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   a.Assemble(0);

   OperatorPtr A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

#ifdef MFEM_USE_SUITESPARSE
   if (!pa && umf_solver)
   {
      ComplexUMFPackSolver csolver(*A.As<ComplexSparseMatrix>());
      csolver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
      csolver.SetPrintLevel(1);
      csolver.Mult(B, X);
   }
   else
#endif
   {
      ConstantCoefficient abs_mass(pow2(omega) * epsilon + std::abs(imag_loss));
      ConstantCoefficient abs_abc(std::abs(imag_abc));

      BilinearForm prec(fespace);
      prec.AddDomainIntegrator(new CurlCurlIntegrator(muinv));
      prec.AddDomainIntegrator(new VectorFEMassIntegrator(abs_mass));
      if (bc == (int)abc1_bc && mesh->bdr_attributes.Size())
      {
         prec.AddBoundaryIntegrator(new VectorFEMassIntegrator(abs_abc), abc_bdr);
      }
      if (pa) { prec.SetAssemblyLevel(AssemblyLevel::PARTIAL); }
      prec.Assemble();

      Array<int> offsets(3);
      offsets[0] = 0;
      offsets[1] = fespace->GetTrueVSize();
      offsets[2] = fespace->GetTrueVSize();
      offsets.PartialSum();

      std::unique_ptr<Operator> pc_r;
      std::unique_ptr<Operator> pc_i;
      const real_t s = (conv == ComplexOperator::HERMITIAN) ? -1.0 : 1.0;

      if (pa)
      {
         pc_r.reset(new OperatorJacobiSmoother(prec, ess_tdof_list));
         pc_i.reset(new ScaledOperator(pc_r.get(), s));
      }
      else
      {
         OperatorPtr P;
         prec.SetDiagonalPolicy(Operator::DIAG_ONE);
         prec.FormSystemMatrix(ess_tdof_list, P);
         pc_r.reset(new GSSmoother(*P.As<SparseMatrix>()));
         pc_i.reset(new ScaledOperator(pc_r.get(), s));
      }

      BlockDiagonalPreconditioner BDP(offsets);
      BDP.SetDiagonalBlock(0, pc_r.get());
      BDP.SetDiagonalBlock(1, pc_i.get());

      GMRESSolver gmres;
      gmres.SetPrintLevel(1);
      gmres.SetKDim(200);
      gmres.SetMaxIter(pa ? 5000 : 2000);
      gmres.SetRelTol(1e-6);
      gmres.SetAbsTol(0.0);
      gmres.SetOperator(*A);
      gmres.SetPreconditioner(BDP);
      gmres.Mult(B, X);
   }

   a.RecoverFEMSolution(X, b, x);

   {
      ofstream mesh_ofs("ex25_subsurface.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      ofstream sol_r_ofs("ex25_subsurface-sol_r.gf");
      ofstream sol_i_ofs("ex25_subsurface-sol_i.gf");
      sol_r_ofs.precision(8);
      sol_i_ofs.precision(8);
      x.real().Save(sol_r_ofs);
      x.imag().Save(sol_i_ofs);
   }

   if (visualization)
   {
      string keys = (dim == 3) ? "keys macF\n" : "keys amrRljcUUuu\n";
      char vishost[] = "localhost";
      int visport = 19916;

      socketstream sol_sock_re(vishost, visport);
      sol_sock_re.precision(8);
      sol_sock_re << "solution\n"
                  << *mesh << x.real() << keys
                  << "window_title 'E real part'" << flush;

      socketstream sol_sock_im(vishost, visport);
      sol_sock_im.precision(8);
      sol_sock_im << "solution\n"
                  << *mesh << x.imag() << keys
                  << "window_title 'E imag part'" << flush;
   }

   delete fespace;
   delete fec;
   delete mesh;
   return 0;
}
