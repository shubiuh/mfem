//                       MFEM UDAR 3D Forward Solver
//
// Compile with: make udar3p
//
// Sample runs:
//   mpirun -np 4 udar3p -m ../data/inline-hex.mesh -rs 1 -rp 1
//   mpirun -np 4 udar3p -mumps -freq 2000,6000,12000,48000,96000
//
// Description:
//   This application computes frequency-domain ultra-deep azimuthal
//   resistivity (UDAR) responses for a straight triaxial logging tool.  At
//   every survey station it solves for three orthogonal magnetic-dipole
//   transmitters and samples three receiver components at every requested
//   transmitter-receiver spacing.  The result is the full complex 3x3
//   magnetic-field transfer tensor.
//
//   The exp(+i omega t) convention is used.  The electric-field equation is
//
//     curl(mu_r^-1 curl(E)) - k0^2 eps_r E + i omega mu0 sigma E
//       = -i omega mu0 J_s,
//
//   The source can be an ideal point magnetic dipole, introduced weakly
//   through m dot curl(v), or a resolved finite-radius current loop.  The
//   magnetic field is recovered from
//   H = i curl(E)/(omega mu0 mu_r).

#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace std;

namespace
{

const real_t mu0 = 4.0e-7 * M_PI;
const real_t eps0 = 8.854187817e-12;

real_t Dot3(const Vector &a, const Vector &b)
{
   return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void Cross3(const Vector &a, const Vector &b, Vector &c)
{
   c.SetSize(3);
   c[0] = a[1] * b[2] - a[2] * b[1];
   c[1] = a[2] * b[0] - a[0] * b[2];
   c[2] = a[0] * b[1] - a[1] * b[0];
}

void Normalize3(Vector &v, const char *name)
{
   const real_t n = v.Norml2();
   MFEM_VERIFY(n > 1.0e-14, name << " must have nonzero length.");
   v /= n;
}

vector<real_t> ParseList(const char *text)
{
   string input(text ? text : "");
   replace(input.begin(), input.end(), ',', ' ');
   stringstream stream(input);
   vector<real_t> values;
   real_t value;
   while (stream >> value) { values.push_back(value); }
   return values;
}

bool ParseVector3(const char *text, Vector &v)
{
   const vector<real_t> values = ParseList(text);
   if (values.empty()) { return false; }
   MFEM_VERIFY(values.size() == 3, "Expected three comma-separated values.");
   v.SetSize(3);
   for (int i = 0; i < 3; i++) { v[i] = values[i]; }
   return true;
}

void ParseTensor3(const char *text, DenseMatrix &tensor)
{
   const vector<real_t> values = ParseList(text);
   tensor.SetSize(3);
   tensor = 0.0;

   if (values.size() == 1)
   {
      for (int i = 0; i < 3; i++) { tensor(i, i) = values[0]; }
   }
   else if (values.size() == 3)
   {
      for (int i = 0; i < 3; i++) { tensor(i, i) = values[i]; }
   }
   else
   {
      MFEM_VERIFY(values.size() == 9,
                  "Conductivity must contain 1, 3, or 9 values.");
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++) { tensor(i, j) = values[3 * i + j]; }
      }
   }

   real_t tensor_scale = 0.0;
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         MFEM_VERIFY(isfinite(tensor(i, j)),
                     "Conductivity entries must be finite.");
         tensor_scale = max(tensor_scale, fabs(tensor(i, j)));
      }
   }
   MFEM_VERIFY(tensor_scale > 0.0,
               "Conductivity tensor must be positive definite.");

   const real_t symmetry_tolerance =
      100.0 * numeric_limits<real_t>::epsilon() * tensor_scale;
   for (int i = 0; i < 3; i++)
   {
      for (int j = i + 1; j < 3; j++)
      {
         MFEM_VERIFY(fabs(tensor(i, j) - tensor(j, i)) <=
                     symmetry_tolerance,
                     "Conductivity tensor must be symmetric.");
         const real_t symmetric_value =
            0.5 * (tensor(i, j) + tensor(j, i));
         tensor(i, j) = tensor(j, i) = symmetric_value;
      }
   }

   DenseMatrix normalized(tensor);
   normalized *= 1.0 / tensor_scale;
   const real_t leading_minor_1 = normalized(0, 0);
   const real_t leading_minor_2 =
      normalized(0, 0) * normalized(1, 1) -
      normalized(0, 1) * normalized(1, 0);
   const real_t determinant =
      normalized(0, 0) * (normalized(1, 1) * normalized(2, 2) -
                          normalized(1, 2) * normalized(2, 1)) -
      normalized(0, 1) * (normalized(1, 0) * normalized(2, 2) -
                          normalized(1, 2) * normalized(2, 0)) +
      normalized(0, 2) * (normalized(1, 0) * normalized(2, 1) -
                          normalized(1, 1) * normalized(2, 0));
   MFEM_VERIFY(leading_minor_1 > 0.0 && leading_minor_2 > 0.0 &&
               determinant > 0.0,
               "Conductivity tensor must be positive definite.");
}

void BuildToolFrame(const Vector &start, const Vector &end,
                    const Vector &up, DenseMatrix &frame)
{
   Vector z(end);
   z -= start;
   Normalize3(z, "Survey trajectory");

   Vector x(up);
   x.Add(-Dot3(up, z), z);
   if (x.Norml2() < 1.0e-10)
   {
      x = 0.0;
      x[(fabs(z[0]) < 0.8) ? 0 : 1] = 1.0;
      x.Add(-Dot3(x, z), z);
   }
   Normalize3(x, "Tool high-side reference");

   Vector y(3);
   Cross3(z, x, y);
   Normalize3(y, "Tool transverse axis");

   frame.SetSize(3);
   for (int d = 0; d < 3; d++)
   {
      frame(d, 0) = x[d];
      frame(d, 1) = y[d];
      frame(d, 2) = z[d];
   }
}

void GetFrameAxis(const DenseMatrix &frame, int axis, Vector &v)
{
   v.SetSize(3);
   for (int d = 0; d < 3; d++) { v[d] = frame(d, axis); }
}

class FormationSigmaCoefficient : public MatrixCoefficient
{
private:
   int model;
   DenseMatrix background;
   DenseMatrix target;
   Vector slab_normal;
   Vector slab_center;
   real_t slab_thickness;
   real_t interface_width;
   Vector box_minimum;
   Vector box_maximum;

public:
   FormationSigmaCoefficient(int model_, const DenseMatrix &background_,
                             const DenseMatrix &target_,
                             const Vector &slab_normal_,
                             const Vector &slab_center_,
                             real_t slab_thickness_, real_t interface_width_,
                             const Vector &box_minimum_,
                             const Vector &box_maximum_)
      : MatrixCoefficient(3), model(model_), background(background_),
        target(target_), slab_normal(slab_normal_), slab_center(slab_center_),
        slab_thickness(slab_thickness_), interface_width(interface_width_),
        box_minimum(box_minimum_), box_maximum(box_maximum_)
   {
      Normalize3(slab_normal, "Slab normal");
   }

   void Eval(DenseMatrix &sigma, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      Vector x(3);
      T.Transform(ip, x);
      real_t target_weight = 0.0;

      if (model == 1)
      {
         Vector delta(x);
         delta -= slab_center;
         const real_t distance = Dot3(delta, slab_normal);
         if (interface_width > 0.0)
         {
            const real_t half_thickness = 0.5 * slab_thickness;
            target_weight = 0.5 *
                            (tanh((half_thickness + distance) /
                                  interface_width) +
                             tanh((half_thickness - distance) /
                                  interface_width));
         }
         else
         {
            target_weight =
               fabs(distance) <= 0.5 * slab_thickness ? 1.0 : 0.0;
         }
      }
      else if (model == 2)
      {
         target_weight = 1.0;
         if (interface_width > 0.0)
         {
            for (int d = 0; d < 3; d++)
            {
               target_weight *= 0.5 *
                  (tanh((x[d] - box_minimum[d]) / interface_width) -
                   tanh((x[d] - box_maximum[d]) / interface_width));
            }
         }
         else
         {
            for (int d = 0; d < 3; d++)
            {
               target_weight *=
                  x[d] >= box_minimum[d] && x[d] <= box_maximum[d] ?
                  1.0 : 0.0;
            }
         }
      }

      sigma.SetSize(3);
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            sigma(i, j) = (1.0 - target_weight) * background(i, j) +
                          target_weight * target(i, j);
         }
      }
   }
};

class ScaledMatrixCoefficient : public MatrixCoefficient
{
private:
   MatrixCoefficient &coefficient;
   real_t scale;

public:
   ScaledMatrixCoefficient(MatrixCoefficient &coefficient_, real_t scale_)
      : MatrixCoefficient(3), coefficient(coefficient_), scale(scale_) { }

   void Eval(DenseMatrix &value, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      coefficient.Eval(value, T, ip);
      value *= scale;
   }
};

class FiniteLoopCurrentCoefficient : public VectorCoefficient
{
private:
   Vector center;
   Vector normal;
   real_t loop_radius;
   real_t wire_radius;
   real_t current_scale;

public:
   FiniteLoopCurrentCoefficient(const Vector &center_, const Vector &moment,
                                real_t loop_radius_, real_t wire_radius_)
      : VectorCoefficient(3), center(center_), normal(moment),
        loop_radius(loop_radius_), wire_radius(wire_radius_)
   {
      MFEM_VERIFY(loop_radius > 0.0 && wire_radius > 0.0,
                  "Loop and wire radii must be positive.");
      const real_t moment_magnitude = normal.Norml2();
      MFEM_VERIFY(moment_magnitude > 0.0,
                  "Finite loop requires a nonzero magnetic moment.");
      normal /= moment_magnitude;
      const real_t current = moment_magnitude / (M_PI * loop_radius * loop_radius);
      current_scale = current / (M_PI * wire_radius * wire_radius);
   }

   void Eval(Vector &current_density, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      Vector relative(3), radial(3), tangent(3);
      T.Transform(ip, relative);
      relative -= center;
      const real_t axial = Dot3(relative, normal);
      radial = relative;
      radial.Add(-axial, normal);
      const real_t radius = radial.Norml2();

      current_density.SetSize(3);
      current_density = 0.0;
      if (radius < 1.0e-14) { return; }

      Cross3(normal, radial, tangent);
      tangent /= radius;
      const real_t tube_distance2 =
         (radius - loop_radius) * (radius - loop_radius) + axial * axial;
      const real_t amplitude =
         current_scale * exp(-tube_distance2 / (wire_radius * wire_radius));
      current_density.Add(amplitude, tangent);
   }
};

bool AddPointMagneticDipoleRHS(ParMesh &mesh, ParFiniteElementSpace &space,
                               const Vector &point, const Vector &moment,
                               real_t omega, real_t tolerance,
                               ParComplexLinearForm &linear_form)
{
   int local_element = -1;
   IntegrationPoint local_ip;
   for (int element = 0; element < mesh.GetNE(); element++)
   {
      ElementTransformation *T = mesh.GetElementTransformation(element);
      IntegrationPoint ip;
      if (T->TransformBack(point, ip, tolerance) ==
          InverseElementTransformation::Inside)
      {
         local_element = element;
         local_ip = ip;
         break;
      }
   }

   const int no_owner = Mpi::WorldSize();
   const int candidate_owner = (local_element >= 0) ? Mpi::WorldRank() : no_owner;
   int owner = no_owner;
   MPI_Allreduce(&candidate_owner, &owner, 1, MPI_INT, MPI_MIN,
                 mesh.GetComm());
   if (owner == no_owner) { return false; }

   if (Mpi::WorldRank() == owner)
   {
      ElementTransformation *T = mesh.GetElementTransformation(local_element);
      T->SetIntPoint(&local_ip);
      const FiniteElement *finite_element = space.GetFE(local_element);
      DenseMatrix curl_shape(finite_element->GetDof(), 3);
      finite_element->CalcPhysCurlShape(*T, curl_shape);
      Vector element_vector(finite_element->GetDof());
      for (int i = 0; i < finite_element->GetDof(); i++)
      {
         real_t value = 0.0;
         for (int d = 0; d < 3; d++)
         {
            value += moment[d] * curl_shape(i, d);
         }
         element_vector[i] = -omega * mu0 * value;
      }
      Array<int> vdofs;
      DofTransformation dof_transformation;
      space.GetElementVDofs(local_element, vdofs, dof_transformation);
      dof_transformation.TransformDual(element_vector);
      linear_form.imag().AddElementVector(vdofs, element_vector);
   }
   return true;
}

bool PointIsInside(ParMesh &mesh, const Vector &point, real_t tolerance)
{
   int local_found = 0;
   for (int element = 0; element < mesh.GetNE() && !local_found; element++)
   {
      ElementTransformation *T = mesh.GetElementTransformation(element);
      IntegrationPoint ip;
      local_found =
         (T->TransformBack(point, ip, tolerance) ==
          InverseElementTransformation::Inside);
   }
   int global_found = 0;
   MPI_Allreduce(&local_found, &global_found, 1, MPI_INT, MPI_MAX,
                 mesh.GetComm());
   return global_found != 0;
}

bool SampleMagneticField(ParMesh &mesh,
                         const ParComplexGridFunction &electric_field,
                         const Vector &point, real_t omega, real_t mu_r,
                         real_t tolerance, Vector &field_real,
                         Vector &field_imag)
{
   Vector local_curl_real(3), local_curl_imag(3);
   local_curl_real = 0.0;
   local_curl_imag = 0.0;
   int local_found = 0;

   for (int element = 0; element < mesh.GetNE(); element++)
   {
      ElementTransformation *T = mesh.GetElementTransformation(element);
      IntegrationPoint ip;
      if (T->TransformBack(point, ip, tolerance) !=
          InverseElementTransformation::Inside)
      {
         continue;
      }
      T->SetIntPoint(&ip);
      Vector element_curl_real(3), element_curl_imag(3);
      electric_field.real().GetCurl(*T, element_curl_real);
      electric_field.imag().GetCurl(*T, element_curl_imag);
      local_curl_real += element_curl_real;
      local_curl_imag += element_curl_imag;
      local_found++;
   }

   Vector curl_real(3), curl_imag(3);
   int found = 0;
   MPI_Allreduce(local_curl_real.GetData(), curl_real.GetData(), 3,
                 MFEM_MPI_REAL_T, MPI_SUM, mesh.GetComm());
   MPI_Allreduce(local_curl_imag.GetData(), curl_imag.GetData(), 3,
                 MFEM_MPI_REAL_T, MPI_SUM, mesh.GetComm());
   MPI_Allreduce(&local_found, &found, 1, MPI_INT, MPI_SUM, mesh.GetComm());

   if (!found) { return false; }
   curl_real /= found;
   curl_imag /= found;

   // exp(+i omega t): H = i curl(E)/(omega mu0 mu_r).
   const real_t scale = 1.0 / (omega * mu0 * mu_r);
   field_real.SetSize(3);
   field_imag.SetSize(3);
   for (int d = 0; d < 3; d++)
   {
      field_real[d] = -scale * curl_imag[d];
      field_imag[d] =  scale * curl_real[d];
   }
   return true;
}

void SetHomogeneousEssentialRHS(const Array<int> &essential_tdofs, Vector &rhs,
                                int true_size)
{
   for (int i = 0; i < essential_tdofs.Size(); i++)
   {
      const int tdof = essential_tdofs[i];
      rhs[tdof] = 0.0;
      rhs[tdof + true_size] = 0.0;
   }
}

void WriteCSVHeader(ostream &output)
{
   output << "station,md_m,frequency_hz,spacing_m,"
          << "tx_x,tx_y,tx_z,rx_x,rx_y,rx_z";
   const char axis[3] = {'x', 'y', 'z'};
   for (int receiver = 0; receiver < 3; receiver++)
   {
      for (int transmitter = 0; transmitter < 3; transmitter++)
      {
         output << ",H" << axis[receiver] << axis[transmitter] << "_real"
                << ",H" << axis[receiver] << axis[transmitter] << "_imag";
      }
   }
   output << '\n';
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   const int myid = Mpi::WorldRank();

   const char *mesh_file = "../data/inline-hex.mesh";
   const char *frequency_text = "2000,12000,48000";
   const char *spacing_text = "";
   const char *survey_start_text = "";
   const char *survey_end_text = "";
   const char *up_text = "0,0,1";
   const char *sigma_background_text = "0.1";
   const char *sigma_target_text = "0.01";
   const char *slab_center_text = "";
   const char *slab_normal_text = "0,0,1";
   const char *box_minimum_text = "";
   const char *box_maximum_text = "";
   const char *output_file = "udar_tensor.csv";
   const char *device_config = "cpu";
   int order = 1;
   int serial_refinements = 0;
   int parallel_refinements = 0;
   int stations = 7;
   int formation_model = 2;
   int boundary_condition = 1;
   int source_type = 0;
   int max_iterations = 2000;
   int krylov_dimension = 200;
   real_t eps_r = 1.0;
   real_t mu_r = 1.0;
   real_t slab_thickness = -1.0;
   real_t interface_width = 0.0;
   real_t transmitter_moment = 1.0;
   real_t loop_radius = -1.0;
   real_t wire_radius = -1.0;
   real_t measured_depth_origin = 0.0;
   real_t relative_tolerance = 1.0e-4;
   real_t absolute_tolerance = 0.0;
   real_t point_tolerance = 1.0e-10;
   bool use_mumps = false;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "3D mesh file.");
   args.AddOption(&order, "-o", "--order", "Nedelec finite element order.");
   args.AddOption(&serial_refinements, "-rs", "--refinements-serial",
                  "Number of serial uniform refinements.");
   args.AddOption(&parallel_refinements, "-rp", "--refinements-parallel",
                  "Number of parallel uniform refinements.");
   args.AddOption(&frequency_text, "-freq", "--frequencies",
                  "Comma-separated frequencies in Hz.");
   args.AddOption(&spacing_text, "-sp", "--spacings",
                  "Comma-separated transmitter-receiver spacings in mesh units.");
   args.AddOption(&stations, "-ns", "--num-stations",
                  "Number of transmitter stations on the straight survey.");
   args.AddOption(&survey_start_text, "-p0", "--survey-start",
                  "Survey start x,y,z; default is based on the mesh bounds.");
   args.AddOption(&survey_end_text, "-p1", "--survey-end",
                  "Survey end x,y,z; default is based on the mesh bounds.");
   args.AddOption(&measured_depth_origin, "-md0", "--measured-depth-origin",
                  "Measured depth assigned to the first station.");
   args.AddOption(&up_text, "-up", "--tool-up",
                  "High-side reference vector used to construct the tool frame.");
   args.AddOption(&transmitter_moment, "-mom", "--transmitter-moment",
                  "Magnetic dipole moment in A m^2.");
   args.AddOption(&source_type, "-src", "--source-type",
                  "0 ideal point magnetic dipole, 1 finite-radius loop.");
   args.AddOption(&loop_radius, "-loop-r", "--loop-radius",
                  "Finite-loop radius in mesh units.");
   args.AddOption(&wire_radius, "-wire-r", "--wire-radius",
                  "Gaussian wire cross-section radius in mesh units.");
   args.AddOption(&formation_model, "-model", "--formation-model",
                  "0 uniform, 1 dipping slab, 2 bounded 3D geobody.");
   args.AddOption(&sigma_background_text, "-sig-bg", "--sigma-background",
                  "Background conductivity: 1, 3 diagonal, or 9 tensor values [S/m].");
   args.AddOption(&sigma_target_text, "-sig-target", "--sigma-target",
                  "Target conductivity: 1, 3 diagonal, or 9 tensor values [S/m].");
   args.AddOption(&slab_center_text, "-slab-c", "--slab-center",
                  "Point on the center plane of the target slab.");
   args.AddOption(&slab_normal_text, "-slab-n", "--slab-normal",
                  "Target slab normal vector.");
   args.AddOption(&slab_thickness, "-slab-t", "--slab-thickness",
                  "Target slab thickness.");
   args.AddOption(&interface_width, "-iw", "--interface-width",
                  "Smooth interface width; zero gives a sharp interface.");
   args.AddOption(&box_minimum_text, "-box-min", "--box-minimum",
                  "Minimum x,y,z corner of the target geobody.");
   args.AddOption(&box_maximum_text, "-box-max", "--box-maximum",
                  "Maximum x,y,z corner of the target geobody.");
   args.AddOption(&eps_r, "-epsr", "--relative-permittivity",
                  "Relative permittivity.");
   args.AddOption(&mu_r, "-mur", "--relative-permeability",
                  "Relative permeability.");
   args.AddOption(&boundary_condition, "-bc", "--boundary-condition",
                  "0 PEC far boundary, 1 first-order absorbing boundary.");
   args.AddOption(&use_mumps, "-mumps", "--mumps-solver", "-gmres",
                  "--iterative-solver", "Use MUMPS or iterative GMRES+AMS.");
   args.AddOption(&max_iterations, "-maxit", "--maximum-iterations",
                  "Maximum GMRES iterations.");
   args.AddOption(&krylov_dimension, "-kdim", "--krylov-dimension",
                  "Restart dimension for GMRES.");
   args.AddOption(&relative_tolerance, "-rtol", "--relative-tolerance",
                  "GMRES relative tolerance.");
   args.AddOption(&absolute_tolerance, "-atol", "--absolute-tolerance",
                  "GMRES absolute tolerance.");
   args.AddOption(&point_tolerance, "-ptol", "--point-tolerance",
                  "Physical point-location tolerance.");
   args.AddOption(&output_file, "-out", "--output", "Tensor CSV output file.");
   args.AddOption(&device_config, "-d", "--device", "Device configuration.");
   args.Parse();
   if (!args.Good())
   {
      if (myid == 0) { args.PrintUsage(cout); }
      return 1;
   }
   if (myid == 0) { args.PrintOptions(cout); }

   MFEM_VERIFY(order >= 1, "Finite element order must be at least one.");
   MFEM_VERIFY(stations >= 1, "At least one survey station is required.");
   MFEM_VERIFY(formation_model >= 0 && formation_model <= 2,
               "Formation model must be 0, 1, or 2.");
   MFEM_VERIFY(boundary_condition == 0 || boundary_condition == 1,
               "Boundary condition must be 0 or 1.");
   MFEM_VERIFY(source_type == 0 || source_type == 1,
               "Source type must be 0 or 1.");
   MFEM_VERIFY(eps_r > 0.0 && mu_r > 0.0,
               "Relative permittivity and permeability must be positive.");
   MFEM_VERIFY(interface_width >= 0.0,
               "Formation-interface width must be nonnegative.");
   MFEM_VERIFY(transmitter_moment != 0.0,
               "Transmitter moment must be nonzero.");

   vector<real_t> frequencies = ParseList(frequency_text);
   MFEM_VERIFY(!frequencies.empty(), "At least one frequency is required.");
   for (real_t frequency : frequencies)
   {
      MFEM_VERIFY(frequency > 0.0, "All frequencies must be positive.");
   }

   Device device(device_config);
   if (myid == 0) { device.Print(); }

   Mesh serial_mesh(mesh_file, 1, 1);
   MFEM_VERIFY(serial_mesh.Dimension() == 3, "UDAR requires a 3D mesh.");
   for (int level = 0; level < serial_refinements; level++)
   {
      serial_mesh.UniformRefinement();
   }

   ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
   serial_mesh.Clear();
   for (int level = 0; level < parallel_refinements; level++)
   {
      mesh.UniformRefinement();
   }

   Vector bounds_minimum(3), bounds_maximum(3);
   mesh.GetBoundingBox(bounds_minimum, bounds_maximum);
   Vector span(bounds_maximum);
   span -= bounds_minimum;
   const real_t min_span = min(span[0], min(span[1], span[2]));
   MFEM_VERIFY(min_span > 0.0, "Mesh bounding box must have positive volume.");

   vector<real_t> spacings = ParseList(spacing_text);
   if (spacings.empty())
   {
      spacings.push_back(0.08 * min_span);
      spacings.push_back(0.16 * min_span);
   }
   for (real_t spacing : spacings)
   {
      MFEM_VERIFY(spacing > 0.0, "All spacings must be positive.");
   }

   Vector survey_start(3), survey_end(3);
   const bool have_start = ParseVector3(survey_start_text, survey_start);
   const bool have_end = ParseVector3(survey_end_text, survey_end);
   if (!have_start || !have_end)
   {
      for (int d = 0; d < 3; d++)
      {
         survey_start[d] = 0.5 * (bounds_minimum[d] + bounds_maximum[d]);
         survey_end[d] = survey_start[d];
      }
      survey_start[0] = bounds_minimum[0] + 0.35 * span[0];
      survey_end[0] = bounds_minimum[0] + 0.73 * span[0];
      survey_start[1] = survey_end[1] =
         bounds_minimum[1] + 0.47 * span[1];
      survey_start[2] = survey_end[2] =
         bounds_minimum[2] + 0.53 * span[2];
   }

   Vector up(3);
   MFEM_VERIFY(ParseVector3(up_text, up), "Tool-up vector is required.");
   DenseMatrix tool_frame(3);
   BuildToolFrame(survey_start, survey_end, up, tool_frame);
   Vector tool_axis(3);
   GetFrameAxis(tool_frame, 2, tool_axis);

   DenseMatrix sigma_background(3), sigma_target(3);
   ParseTensor3(sigma_background_text, sigma_background);
   ParseTensor3(sigma_target_text, sigma_target);

   Vector slab_center(3), slab_normal(3), box_minimum(3), box_maximum(3);
   if (!ParseVector3(slab_center_text, slab_center))
   {
      slab_center = bounds_minimum;
      slab_center += bounds_maximum;
      slab_center *= 0.5;
   }
   MFEM_VERIFY(ParseVector3(slab_normal_text, slab_normal),
               "Slab normal is required.");
   if (slab_thickness <= 0.0) { slab_thickness = 0.20 * min_span; }

   if (!ParseVector3(box_minimum_text, box_minimum) ||
       !ParseVector3(box_maximum_text, box_maximum))
   {
      for (int d = 0; d < 3; d++)
      {
         box_minimum[d] = bounds_minimum[d] + 0.30 * span[d];
         box_maximum[d] = bounds_minimum[d] + 0.70 * span[d];
      }
      box_minimum[0] = bounds_minimum[0] + 0.48 * span[0];
      box_maximum[0] = bounds_minimum[0] + 0.68 * span[0];
   }
   for (int d = 0; d < 3; d++)
   {
      MFEM_VERIFY(box_minimum[d] < box_maximum[d],
                  "Each box minimum must be less than its maximum.");
   }
   if (loop_radius <= 0.0) { loop_radius = 0.03 * min_span; }
   if (wire_radius <= 0.0) { wire_radius = 0.20 * loop_radius; }
   if (myid == 0)
   {
      cout << "UDAR unknowns: ";
   }
   ND_FECollection finite_elements(order, 3);
   ParFiniteElementSpace space(&mesh, &finite_elements);
   if (myid == 0)
   {
      cout << space.GlobalTrueVSize() << " complex ("
           << 2 * space.GlobalTrueVSize() << " real block unknowns)\n";
      cout << "Frequencies: ";
      for (real_t frequency : frequencies) { cout << frequency << ' '; }
      cout << "Hz\nSpacings: ";
      for (real_t spacing : spacings) { cout << spacing << ' '; }
      cout << "mesh units\n";
      cout << "Source: "
           << (source_type == 0 ? "point magnetic dipole" : "finite current loop")
           << '\n';
   }

   FormationSigmaCoefficient sigma(formation_model, sigma_background,
                                   sigma_target, slab_normal, slab_center,
                                   slab_thickness, interface_width,
                                   box_minimum, box_maximum);

   Array<int> essential_tdofs;
   if (boundary_condition == 0 && mesh.bdr_attributes.Size())
   {
      Array<int> essential_boundary(mesh.bdr_attributes.Max());
      essential_boundary = 1;
      space.GetEssentialTrueDofs(essential_boundary, essential_tdofs);
   }

   ofstream output;
   if (myid == 0)
   {
      output.open(output_file);
      MFEM_VERIFY(output.good(), "Unable to open output file " << output_file);
      output << setprecision(16) << scientific;
      WriteCSVHeader(output);
   }

   Vector trajectory(survey_end);
   trajectory -= survey_start;
   const real_t survey_length = trajectory.Norml2();
   const int true_size = space.GetTrueVSize();
   const ComplexOperator::Convention convention = ComplexOperator::HERMITIAN;

   for (real_t frequency : frequencies)
   {
      const real_t omega = 2.0 * M_PI * frequency;
      const real_t k0 = omega * sqrt(mu0 * eps0);
      ConstantCoefficient inverse_mu_r(1.0 / mu_r);
      ConstantCoefficient real_mass(-k0 * k0 * eps_r);
      ScaledMatrixCoefficient imaginary_mass(sigma, omega * mu0);

      unique_ptr<ConstantCoefficient> abc_real;
      unique_ptr<ConstantCoefficient> abc_imaginary;
      Array<int> absorbing_boundary;

      ParSesquilinearForm form(&space, convention);
      form.AddDomainIntegrator(new CurlCurlIntegrator(inverse_mu_r), nullptr);
      form.AddDomainIntegrator(new VectorFEMassIntegrator(real_mass),
                               new VectorFEMassIntegrator(imaginary_mass));

      if (boundary_condition == 1 && mesh.bdr_attributes.Size())
      {
         absorbing_boundary.SetSize(mesh.bdr_attributes.Max());
         absorbing_boundary = 1;
         const complex<real_t> imaginary_unit(0.0, 1.0);
         const complex<real_t> background_epsilon =
            eps_r - imaginary_unit * sigma_background(0, 0) / (omega * eps0);
         const complex<real_t> boundary_wave_number =
            k0 * sqrt(mu_r * background_epsilon);
         // exp(+i omega t), outgoing exp(-i k r): the weak boundary term is
         // +i k/mu_r times the tangential-field mass matrix.
         const complex<real_t> gamma =
            imaginary_unit * boundary_wave_number / mu_r;
         abc_real.reset(new ConstantCoefficient(gamma.real()));
         abc_imaginary.reset(new ConstantCoefficient(gamma.imag()));
         form.AddBoundaryIntegrator(new VectorFEMassIntegrator(*abc_real),
                                    new VectorFEMassIntegrator(*abc_imaginary),
                                    absorbing_boundary);
      }

      form.Assemble();
      OperatorPtr system_operator;
      form.FormSystemMatrix(essential_tdofs, system_operator);

      unique_ptr<HypreParMatrix> real_block_operator;
      unique_ptr<HypreILU> ilu;
      unique_ptr<Solver> solver;

      if (use_mumps)
      {
#ifdef MFEM_USE_COMPLEX_MUMPS
         ComplexMUMPSSolver *mumps = new ComplexMUMPSSolver(MPI_COMM_WORLD);
         mumps->SetPrintLevel(0);
         mumps->SetMatrixSymType(ComplexMUMPSSolver::UNSYMMETRIC);
         mumps->SetOperator(*system_operator.As<ComplexHypreParMatrix>());
         solver.reset(mumps);
#elif defined(MFEM_USE_MUMPS)
         real_block_operator.reset(
            system_operator.As<ComplexHypreParMatrix>()->GetSystemMatrix());
         MUMPSSolver *mumps = new MUMPSSolver(MPI_COMM_WORLD);
         mumps->SetPrintLevel(0);
         mumps->SetMatrixSymType(MUMPSSolver::UNSYMMETRIC);
         mumps->SetOperator(*real_block_operator);
         solver.reset(mumps);
#else
         MFEM_ABORT("MUMPS was requested, but MFEM_USE_MUMPS is disabled.");
#endif
      }
      else
      {
         // Apply ILU to the complete real 2x2 representation. A diagonal AMS
         // approximation leaves the conductivity-dominated gradient subspace
         // poorly conditioned and can stagnate even on the smoke-test mesh.
         real_block_operator.reset(
            system_operator.As<ComplexHypreParMatrix>()->GetSystemMatrix());
         ilu.reset(new HypreILU);
         // Include a global Schur-complement correction; type 0 is only
         // block-Jacobi ILU and stagnates once the mesh spans MPI ranks.
         ilu->SetType(10);
         ilu->SetLevelOfFill(2);
         ilu->SetOperator(*real_block_operator);

         FGMRESSolver *fgmres = new FGMRESSolver(MPI_COMM_WORLD);
         fgmres->SetPrintLevel(0);
         fgmres->SetKDim(krylov_dimension);
         fgmres->SetMaxIter(max_iterations);
         fgmres->SetRelTol(relative_tolerance);
         fgmres->SetAbsTol(absolute_tolerance);
         fgmres->SetOperator(*real_block_operator);
         fgmres->SetPreconditioner(*ilu);
         solver.reset(fgmres);
      }

      if (myid == 0)
      {
         cout << "Solving " << frequency << " Hz with "
              << (use_mumps ? "MUMPS" : "FGMRES+ILU") << "\n";
      }

      for (int station = 0; station < stations; station++)
      {
         const real_t fraction =
            (stations == 1) ? 0.0 : real_t(station) / real_t(stations - 1);
         Vector transmitter(survey_start);
         transmitter.Add(fraction, trajectory);
         MFEM_VERIFY(PointIsInside(mesh, transmitter, point_tolerance),
                     "Transmitter station " << station << " lies outside the mesh.");
         vector<Vector> rhs(3), solution(3);
         for (int transmitter_axis = 0; transmitter_axis < 3;
              transmitter_axis++)
         {
            Vector moment(3);
            GetFrameAxis(tool_frame, transmitter_axis, moment);
            moment *= transmitter_moment;
            unique_ptr<FiniteLoopCurrentCoefficient> loop_current;
            unique_ptr<ConstantCoefficient> loop_rhs_scale;
            unique_ptr<ScalarVectorProductCoefficient> loop_rhs;
            ParComplexLinearForm linear_form(&space, convention);
            if (source_type == 0)
            {
               linear_form.Assemble();
               MFEM_VERIFY(AddPointMagneticDipoleRHS(
                              mesh, space, transmitter, moment, omega,
                              point_tolerance, linear_form),
                           "Transmitter station " << station
                           << " lies outside the mesh.");
            }
            else
            {
               loop_current.reset(new FiniteLoopCurrentCoefficient(
                                     transmitter, moment, loop_radius,
                                     wire_radius));
               loop_rhs_scale.reset(new ConstantCoefficient(-omega * mu0));
               loop_rhs.reset(new ScalarVectorProductCoefficient(
                                 *loop_rhs_scale, *loop_current));
               linear_form.AddDomainIntegrator(
                  nullptr, new VectorFEDomainLFIntegrator(*loop_rhs));
               linear_form.Assemble();
            }

            rhs[transmitter_axis].SetSize(2 * true_size);
            linear_form.ParallelAssemble(rhs[transmitter_axis]);
            SetHomogeneousEssentialRHS(essential_tdofs,
                                       rhs[transmitter_axis], true_size);
            solution[transmitter_axis].SetSize(2 * true_size);
            solution[transmitter_axis] = 0.0;
         }

         Array<const Vector *> rhs_array(3);
         Array<Vector *> solution_array(3);
         for (int axis = 0; axis < 3; axis++)
         {
            rhs_array[axis] = &rhs[axis];
            solution_array[axis] = &solution[axis];
         }
         if (use_mumps)
         {
            // ComplexMUMPSSolver and MUMPSSolver both use one batched solve for
            // the three source polarizations and retain the frequency's factors.
            solver->ArrayMult(rhs_array, solution_array);
         }
         else
         {
            IterativeSolver *iterative =
               dynamic_cast<IterativeSolver *>(solver.get());
            MFEM_VERIFY(iterative, "Internal iterative solver setup failed.");
            for (int axis = 0; axis < 3; axis++)
            {
               iterative->Mult(rhs[axis], solution[axis]);
               MFEM_VERIFY(iterative->GetConverged(),
                           "GMRES failed at frequency " << frequency
                           << " Hz, station " << station
                           << ", transmitter axis " << axis << '.');
            }
         }

         vector<unique_ptr<ParComplexGridFunction>> electric_fields(3);
         for (int axis = 0; axis < 3; axis++)
         {
            electric_fields[axis].reset(new ParComplexGridFunction(&space));
            Vector solution_real, solution_imaginary;
            solution_real.MakeRef(solution[axis], 0, true_size);
            solution_imaginary.MakeRef(solution[axis], true_size, true_size);
            electric_fields[axis]->real().SetFromTrueDofs(solution_real);
            electric_fields[axis]->imag().SetFromTrueDofs(solution_imaginary);
         }

         for (real_t spacing : spacings)
         {
            Vector receiver(transmitter);
            receiver.Add(-spacing, tool_axis);
            complex<real_t> response[3][3];
            bool receiver_found = true;

            for (int transmitter_axis = 0; transmitter_axis < 3;
                 transmitter_axis++)
            {
               Vector field_real(3), field_imaginary(3);
               const bool found = SampleMagneticField(
                  mesh, *electric_fields[transmitter_axis], receiver, omega,
                  mu_r, point_tolerance, field_real, field_imaginary);
               receiver_found = receiver_found && found;

               for (int receiver_axis = 0; receiver_axis < 3; receiver_axis++)
               {
                  Vector receiver_direction(3);
                  GetFrameAxis(tool_frame, receiver_axis, receiver_direction);
                  if (found)
                  {
                     response[receiver_axis][transmitter_axis] =
                        complex<real_t>(Dot3(receiver_direction, field_real),
                                        Dot3(receiver_direction,
                                             field_imaginary)) /
                        transmitter_moment;
                  }
                  else
                  {
                     const real_t nan =
                        numeric_limits<real_t>::quiet_NaN();
                     response[receiver_axis][transmitter_axis] =
                        complex<real_t>(nan, nan);
                  }
               }
            }

            if (!receiver_found && myid == 0)
            {
               cerr << "Warning: receiver for station " << station
                    << ", spacing " << spacing << " lies outside the mesh.\n";
            }

            if (myid == 0)
            {
               output << station << ','
                      << measured_depth_origin + fraction * survey_length << ','
                      << frequency << ',' << spacing;
               for (int d = 0; d < 3; d++) { output << ',' << transmitter[d]; }
               for (int d = 0; d < 3; d++) { output << ',' << receiver[d]; }
               for (int receiver_axis = 0; receiver_axis < 3; receiver_axis++)
               {
                  for (int transmitter_axis = 0; transmitter_axis < 3;
                       transmitter_axis++)
                  {
                     output << ',' << response[receiver_axis][transmitter_axis].real()
                            << ',' << response[receiver_axis][transmitter_axis].imag();
                  }
               }
               output << '\n';
            }
         }

         if (myid == 0)
         {
            cout << "  station " << station + 1 << '/' << stations << "\n";
         }
      }
   }

   if (myid == 0)
   {
      cout << "Wrote full complex 3x3 H tensor responses to "
           << output_file << '\n';
   }
   return 0;
}
