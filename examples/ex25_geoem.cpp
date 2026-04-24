// MFEM Example 25-GeoEM: conductive curl-curl E formulation.
// Time convention: exp(+j omega t).
// Equation:
//   curl mu_r^{-1} curl E - k0^2 (eps_r - j sigma/(omega eps0)) E
//     = -j omega mu0 Js.

#include "mfem.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace std;

static const real_t MU0 = 4.0e-7 * M_PI;
static const real_t EPS0 = 8.854187817e-12;

static real_t sqr(real_t x) { return x*x; }

static void GetBox(const Mesh &mesh, Vector &mn, Vector &mx)
{
   const int dim = mesh.Dimension();
   mn.SetSize(dim); mx.SetSize(dim);
   mn = numeric_limits<real_t>::infinity();
   mx = -numeric_limits<real_t>::infinity();
   for (int i = 0; i < mesh.GetNV(); i++)
   {
      const real_t *v = mesh.GetVertex(i);
      for (int d = 0; d < dim; d++)
      {
         mn[d] = std::min(mn[d], v[d]);
         mx[d] = std::max(mx[d], v[d]);
      }
   }
}

static void ParseTensor(const char *txt, DenseMatrix &M, int dim, real_t diag)
{
   M.SetSize(dim); M = 0.0;
   for (int i = 0; i < dim; i++) { M(i,i) = diag; }
   if (!txt || string(txt).empty()) { return; }
   string s(txt); replace(s.begin(), s.end(), ',', ' ');
   stringstream ss(s); vector<real_t> a; real_t x;
   while (ss >> x) { a.push_back(x); }
   if ((int)a.size() == dim*dim)
   {
      int k = 0;
      for (int i = 0; i < dim; i++)
         for (int j = 0; j < dim; j++) { M(i,j) = a[k++]; }
   }
   else if (dim == 2 && a.size() == 9)
   {
      M(0,0) = a[0]; M(0,1) = a[1];
      M(1,0) = a[3]; M(1,1) = a[4];
   }
   else
   {
      out << "Warning: sigma tensor parse failed; using isotropic sigma.\n";
   }
}

class SigmaCoef : public MatrixCoefficient
{
   int dim, model, ti_axis;
   real_t sig, sigh, sigv, sigx, sigy, sigz;
   DenseMatrix tensor;
public:
   SigmaCoef(int dim_, int model_, int ti_axis_, real_t sig_, real_t sigh_,
             real_t sigv_, real_t sigx_, real_t sigy_, real_t sigz_,
             const DenseMatrix &tensor_)
      : MatrixCoefficient(dim_), dim(dim_), model(model_), ti_axis(ti_axis_),
        sig(sig_), sigh(sigh_), sigv(sigv_), sigx(sigx_), sigy(sigy_),
        sigz(sigz_)
   {
      if (ti_axis < 0 || ti_axis >= dim) { ti_axis = dim - 1; }
      tensor.SetSize(dim); tensor = 0.0;
      for (int i = 0; i < dim; i++)
         for (int j = 0; j < dim; j++) { tensor(i,j) = tensor_(i,j); }
   }

   void Eval(DenseMatrix &M, ElementTransformation &, const IntegrationPoint &) override
   {
      M.SetSize(dim); M = 0.0;
      if (model == 0)
      {
         for (int i = 0; i < dim; i++) { M(i,i) = sig; }
      }
      else if (model == 1)
      {
         for (int i = 0; i < dim; i++) { M(i,i) = sigh; }
         M(ti_axis, ti_axis) = sigv;
      }
      else if (model == 2)
      {
         M(0,0) = sigx;
         if (dim > 1) { M(1,1) = sigy; }
         if (dim > 2) { M(2,2) = sigz; }
      }
      else
      {
         for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++) { M(i,j) = tensor(i,j); }
      }
   }
};

class ScaleMatCoef : public MatrixCoefficient
{
   int dim;
   MatrixCoefficient &A;
   real_t c;
public:
   ScaleMatCoef(int dim_, MatrixCoefficient &A_, real_t c_)
      : MatrixCoefficient(dim_), dim(dim_), A(A_), c(c_) { }
   void Eval(DenseMatrix &M, ElementTransformation &T, const IntegrationPoint &ip) override
   {
      A.Eval(M, T, ip);
      M *= c;
   }
};

class PosMassCoef : public MatrixCoefficient
{
   int dim;
   MatrixCoefficient &sigma;
   real_t k0, omega, epsr;
public:
   PosMassCoef(int dim_, MatrixCoefficient &sigma_, real_t k0_, real_t omega_, real_t epsr_)
      : MatrixCoefficient(dim_), dim(dim_), sigma(sigma_), k0(k0_), omega(omega_), epsr(epsr_) { }
   void Eval(DenseMatrix &M, ElementTransformation &T, const IntegrationPoint &ip) override
   {
      DenseMatrix S; sigma.Eval(S, T, ip);
      M.SetSize(dim); M = 0.0;
      for (int i = 0; i < dim; i++)
      {
         real_t rowsum = 0.0;
         for (int j = 0; j < dim; j++) { rowsum += fabs(S(i,j)); }
         M(i,i) = k0*k0*fabs(epsr) + omega*MU0*rowsum;
      }
   }
};

class GaussianJ : public VectorCoefficient
{
   int dim;
   Vector xc, p;
   real_t width, amp;
public:
   GaussianJ(int dim_, const Vector &xc_, const Vector &p_, real_t width_, real_t amp_)
      : VectorCoefficient(dim_), dim(dim_), xc(xc_), p(p_), width(width_), amp(amp_) { }
   void Eval(Vector &J, ElementTransformation &T, const IntegrationPoint &ip) override
   {
      Vector x(dim); T.Transform(ip, x);
      real_t r2 = 0.0;
      for (int i = 0; i < dim; i++) { r2 += sqr(x[i] - xc[i]); }
      J.SetSize(dim);
      const real_t g = amp*exp(-r2/(width*width));
      for (int i = 0; i < dim; i++) { J[i] = g*p[i]; }
   }
};

int main(int argc, char *argv[])
{
   const char *mesh_file = "../data/inline-hex.mesh";
   const char *sig_tensor = "";
   int order = 1, ref = 0, smodel = 0, ti_axis = -1, bc = 0;
   real_t freq = 1000.0, mur = 1.0, epsr = 10.0;
   real_t sig = 1e-2, sigh = 1e-2, sigv = 1e-3;
   real_t sigx = 1e-2, sigy = 1e-2, sigz = 1e-3;
   real_t abc_sig = -1.0, abc_mur = -1.0, abc_epsr = -1.0;
   real_t sx = numeric_limits<real_t>::quiet_NaN();
   real_t sy = numeric_limits<real_t>::quiet_NaN();
   real_t sz = numeric_limits<real_t>::quiet_NaN();
   real_t jx = 1.0, jy = 0.0, jz = 0.0, sw = -1.0, swr = 0.05, amp = 1.0;
   bool vis = true, herm = true, umf = false;
   const char *device_config = "cpu";

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file.");
   args.AddOption(&order, "-o", "--order", "Nedelec order.");
   args.AddOption(&ref, "-ref", "--refinements", "Uniform refinements.");
   args.AddOption(&freq, "-f", "--frequency", "Frequency [Hz].");
   args.AddOption(&mur, "-mur", "--relative-permeability", "mu_r.");
   args.AddOption(&epsr, "-epsr", "--relative-permittivity", "epsilon_r.");
   args.AddOption(&smodel, "-sm", "--sigma-model", "0 iso, 1 TI, 2 biaxial, 3 tensor.");
   args.AddOption(&sig, "-sig", "--sigma", "Isotropic sigma [S/m].");
   args.AddOption(&sigh, "-sigh", "--sigma-horizontal", "TI horizontal sigma [S/m].");
   args.AddOption(&sigv, "-sigv", "--sigma-vertical", "TI axis sigma [S/m].");
   args.AddOption(&ti_axis, "-tiaxis", "--ti-axis", "TI axis index; default last coordinate.");
   args.AddOption(&sigx, "-sigx", "--sigma-x", "Biaxial sigma_x [S/m].");
   args.AddOption(&sigy, "-sigy", "--sigma-y", "Biaxial sigma_y [S/m].");
   args.AddOption(&sigz, "-sigz", "--sigma-z", "Biaxial sigma_z [S/m].");
   args.AddOption(&sig_tensor, "-sigt", "--sigma-tensor", "Full sigma tensor, row-major.");
   args.AddOption(&bc, "-bc", "--boundary-condition", "0 PEC, 1 first-order ABC.");
   args.AddOption(&abc_sig, "-abc-sig", "--abc-sigma", "ABC exterior sigma. Default sigma.");
   args.AddOption(&abc_mur, "-abc-mur", "--abc-mu-r", "ABC exterior mu_r. Default mur.");
   args.AddOption(&abc_epsr, "-abc-epsr", "--abc-eps-r", "ABC exterior eps_r. Default epsr.");
   args.AddOption(&sx, "-srcx", "--source-x", "Source x, default center.");
   args.AddOption(&sy, "-srcy", "--source-y", "Source y, default center.");
   args.AddOption(&sz, "-srcz", "--source-z", "Source z, default center.");
   args.AddOption(&jx, "-jx", "--current-x", "Current moment x.");
   args.AddOption(&jy, "-jy", "--current-y", "Current moment y.");
   args.AddOption(&jz, "-jz", "--current-z", "Current moment z.");
   args.AddOption(&sw, "-sw", "--source-width", "Gaussian width, default relative.");
   args.AddOption(&swr, "-swr", "--source-width-relative", "Width/domain diagonal.");
   args.AddOption(&amp, "-sa", "--source-amplitude", "Current amplitude.");
   args.AddOption(&herm, "-herm", "--hermitian", "-no-herm", "--no-hermitian", "Hermitian complex convention.");
#ifdef MFEM_USE_SUITESPARSE
   args.AddOption(&umf, "-umf", "--umfpack", "-no-umf", "--no-umfpack", "Use ComplexUMFPack.");
#endif
   args.AddOption(&vis, "-vis", "--visualization", "-no-vis", "--no-visualization", "GLVis output.");
   args.AddOption(&device_config, "-d", "--device", "Device string.");
   args.Parse();
   if (!args.Good()) { args.PrintUsage(cout); return 1; }
   args.PrintOptions(cout);
   if (freq <= 0.0) { mfem::err << "Frequency must be positive.\n"; return 2; }

   smodel = std::max(0, std::min(3, smodel));
   bc = std::max(0, std::min(1, bc));

   Device device(device_config); device.Print();
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   const int dim = mesh->Dimension();
   if (dim < 2 || dim > 3) { mfem::err << "2D or 3D mesh required.\n"; return 3; }
   for (int l = 0; l < ref; l++) { mesh->UniformRefinement(); }

   const real_t omega = 2.0*M_PI*freq;
   const real_t k0 = omega*sqrt(MU0*EPS0);
   if (ti_axis < 0 || ti_axis >= dim) { ti_axis = dim - 1; }

   DenseMatrix tensor; ParseTensor(sig_tensor, tensor, dim, sig);
   Vector mn, mx; GetBox(*mesh, mn, mx);
   real_t diag2 = 0.0;
   for (int i = 0; i < dim; i++) { diag2 += sqr(mx[i] - mn[i]); }
   const real_t dbox = sqrt(diag2);

   Vector xc(dim), p(dim); p = 0.0;
   for (int i = 0; i < dim; i++) { xc[i] = 0.5*(mn[i] + mx[i]); }
   if (isfinite(sx)) { xc[0] = sx; }
   if (dim > 1 && isfinite(sy)) { xc[1] = sy; }
   if (dim > 2 && isfinite(sz)) { xc[2] = sz; }
   p[0] = jx; if (dim > 1) { p[1] = jy; } if (dim > 2) { p[2] = jz; }
   if (sw <= 0.0) { sw = std::max(swr*dbox, 1e-12); }

   if (abc_sig < 0.0) { abc_sig = sig; }
   if (abc_mur < 0.0) { abc_mur = mur; }
   if (abc_epsr < 0.0) { abc_epsr = epsr; }

   ND_FECollection fec(order, dim);
   FiniteElementSpace fes(mesh, &fec);
   out << "Unknowns: " << fes.GetTrueVSize() << "\n";
   out << "Convention exp(+j omega t), f=" << freq << " Hz, k0=" << k0 << "\n";
   out << "Solving curl mur^-1 curl E - k0^2(epsr - j*sigma/(omega*eps0))E = -j*omega*mu0*Js\n";

   Array<int> ess_tdof, ess_bdr;
   if (mesh->bdr_attributes.Size())
   {
      ess_bdr.SetSize(mesh->bdr_attributes.Max());
      ess_bdr = (bc == 0) ? 1 : 0;
      fes.GetEssentialTrueDofs(ess_bdr, ess_tdof);
   }

   ComplexOperator::Convention conv = herm ? ComplexOperator::HERMITIAN
                                           : ComplexOperator::BLOCK_SYMMETRIC;
   SigmaCoef sigma(dim, smodel, ti_axis, sig, sigh, sigv, sigx, sigy, sigz, tensor);
   ScaleMatCoef im_mass(dim, sigma, omega*MU0);
   PosMassCoef prec_mass(dim, sigma, k0, omega, epsr);
   ConstantCoefficient mur_inv(1.0/mur), re_mass(-k0*k0*epsr);

   GaussianJ Js(dim, xc, p, sw, amp);
   ConstantCoefficient rhs_scale(-omega*MU0);
   ScalarVectorProductCoefficient rhs_im(rhs_scale, Js);

   ComplexLinearForm b(&fes, conv);
   b.AddDomainIntegrator(NULL, new VectorFEDomainLFIntegrator(rhs_im));
   b.Assemble();

   ComplexGridFunction x(&fes); x = complex<real_t>(0.0, 0.0);
   SesquilinearForm a(&fes, conv);
   a.AddDomainIntegrator(new CurlCurlIntegrator(mur_inv), NULL);
   a.AddDomainIntegrator(new VectorFEMassIntegrator(re_mass), new VectorFEMassIntegrator(im_mass));

   ConstantCoefficient *abc_re = NULL, *abc_im = NULL, *abc_abs = NULL;
   Array<int> abc_bdr;
   if (bc == 1 && mesh->bdr_attributes.Size())
   {
      abc_bdr.SetSize(mesh->bdr_attributes.Max()); abc_bdr = 1;
      const complex<real_t> j(0.0, 1.0);
      const complex<real_t> epsc = abc_epsr - j*abc_sig/(omega*EPS0);
      const complex<real_t> kb = k0*sqrt(abc_mur*epsc);
      const complex<real_t> gamma = -j*kb/abc_mur;
      out << "ABC gamma = " << gamma.real() << " + j " << gamma.imag() << "\n";
      abc_re = new ConstantCoefficient(gamma.real());
      abc_im = new ConstantCoefficient(gamma.imag());
      abc_abs = new ConstantCoefficient(abs(gamma));
      a.AddBoundaryIntegrator(new VectorFEMassIntegrator(*abc_re), new VectorFEMassIntegrator(*abc_im), abc_bdr);
   }

   a.Assemble(0);
   OperatorPtr A; Vector B, X;
   a.FormLinearSystem(ess_tdof, x, b, A, X, B);

   bool solved = false;
#ifdef MFEM_USE_SUITESPARSE
   if (umf)
   {
      ComplexUMFPackSolver solver(*A.As<ComplexSparseMatrix>());
      solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
      solver.SetPrintLevel(1);
      solver.Mult(B, X);
      solved = true;
   }
#endif
   if (!solved)
   {
      BilinearForm pc(&fes);
      pc.AddDomainIntegrator(new CurlCurlIntegrator(mur_inv));
      pc.AddDomainIntegrator(new VectorFEMassIntegrator(prec_mass));
      if (bc == 1 && abc_abs) { pc.AddBoundaryIntegrator(new VectorFEMassIntegrator(*abc_abs), abc_bdr); }
      pc.Assemble();
      OperatorPtr PCOp; pc.SetDiagonalPolicy(Operator::DIAG_ONE); pc.FormSystemMatrix(ess_tdof, PCOp);

      Array<int> offsets(3);
      offsets[0] = 0; offsets[1] = fes.GetTrueVSize(); offsets[2] = fes.GetTrueVSize();
      offsets.PartialSum();
      const real_t s = (conv == ComplexOperator::HERMITIAN) ? -1.0 : 1.0;
      GSSmoother pc_r(*PCOp.As<SparseMatrix>());
      ScaledOperator pc_i(&pc_r, s);
      BlockDiagonalPreconditioner block_pc(offsets);
      block_pc.SetDiagonalBlock(0, &pc_r);
      block_pc.SetDiagonalBlock(1, &pc_i);

      GMRESSolver gmres;
      gmres.SetPrintLevel(1); gmres.SetKDim(200); gmres.SetMaxIter(2000);
      gmres.SetRelTol(1e-8); gmres.SetAbsTol(0.0);
      gmres.SetOperator(*A); gmres.SetPreconditioner(block_pc); gmres.Mult(B, X);
   }

   a.RecoverFEMSolution(X, b, x);
   ofstream mesh_ofs("ex25_geoem.mesh"); mesh_ofs.precision(8); mesh->Print(mesh_ofs);
   ofstream solr("ex25_geoem-sol_r.gf"), soli("ex25_geoem-sol_i.gf");
   solr.precision(8); soli.precision(8); x.real().Save(solr); x.imag().Save(soli);
   out << "Saved ex25_geoem.mesh, ex25_geoem-sol_r.gf, ex25_geoem-sol_i.gf\n";

   if (vis)
   {
      char host[] = "localhost"; int port = 19916;
      string keys = (dim == 3) ? "keys macF\n" : "keys amrRljcUUuu\n";
      socketstream sr(host, port); sr.precision(8);
      sr << "solution\n" << *mesh << x.real() << keys << "window_title 'GeoEM E real'" << flush;
      socketstream si(host, port); si.precision(8);
      si << "solution\n" << *mesh << x.imag() << keys << "window_title 'GeoEM E imag'" << flush;
   }

   delete abc_re; delete abc_im; delete abc_abs; delete mesh;
   return 0;
}
