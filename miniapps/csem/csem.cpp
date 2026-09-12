// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

//                 MFEM Miniapp: CSEM (csem)
//
// This miniapp implements the total-electric-field controlled-source
// electromagnetic formulation, the Weiss PRESB inner-outer solver, and a
// compact two-parameter Gauss-Newton inversion.  It requires MPI and hypre.

#include "mfem.hpp"
#include "csem_solver.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::csem;
using namespace std;

namespace
{

constexpr real_t mu0 = 4.0e-7 * M_PI;
constexpr real_t epsilon0 = 8.8541878128e-12;

struct CSEMOptions
{
   const char *mesh_file = "";
   const char *solver = "presb";
   const char *observed_file = "";
   const char *receiver_file = "";
   const char *receiver_field = "E";
   const char *source_file = "";
   const char *frequencies = "";
   const char *output_prefix = "csem";
   const char *device = "cpu";

   int order = 2;
   int refinements = 0;
   int nx = 8;
   int ny = 8;
   int nz = 6;
   int receivers = 7;
   int outer_max_iterations = 300;
   int outer_kdim = 50;
   int inner_max_iterations = 100;
   int inner_kdim = 30;
   int inversion_iterations = 5;
   int line_search_steps = 8;
   int random_seed = 1729;

   real_t size_x = 6000.0;
   real_t size_y = 6000.0;
   real_t air_thickness = 1000.0;
   real_t earth_depth = 4000.0;
   real_t frequency = 10.0;
   real_t relative_permeability = 1.0;
   real_t earth_relative_permittivity = 10.0;
   real_t air_conductivity = 1.0e-8;
   real_t background_conductivity = 1.0e-3;
   real_t target_conductivity = 5.0e-2;
   real_t initial_background_conductivity = 3.0e-3;
   real_t initial_target_conductivity = 1.0e-2;
   real_t minimum_conductivity = 1.0e-6;
   real_t maximum_conductivity = 10.0;

   real_t target_x = 0.0;
   real_t target_y = 800.0;
   real_t target_z = 1200.0;
   real_t target_size_x = 1400.0;
   real_t target_size_y = 1400.0;
   real_t target_size_z = 900.0;

   real_t source_x = 0.0;
   real_t source_y = 0.0;
   real_t source_z = 0.0;
   real_t source_length = 1000.0;
   real_t source_radius = 350.0;
   real_t source_current = 1.0;
   real_t source_azimuth = 0.0;
   real_t source_dip = 0.0;

   real_t receiver_x = 0.0;
   real_t receiver_y = 800.0;
   real_t receiver_z = 0.0;
   real_t receiver_span = 3600.0;
   real_t receiver_radius = 450.0;
   real_t receiver_azimuth = 0.0;
   real_t receiver_dip = 0.0;
   real_t receiver_line_azimuth = 0.0;

   real_t forward_outer_tolerance = 1.0e-8;
   real_t forward_inner_tolerance = 1.0e-2;
   real_t inversion_outer_tolerance = 1.0e-4;
   real_t inversion_inner_tolerance = 1.0e-1;
   real_t data_relative_error = 0.05;
   real_t data_error_floor = 1.0e-3;
   real_t synthetic_noise = 0.0;
   real_t regularization = 1.0e-2;
   real_t gauss_newton_damping = 1.0e-6;
   real_t inversion_gradient_tolerance = 1.0e-5;
   real_t inversion_step_tolerance = 1.0e-4;
   real_t jacobian_fd_step = 1.0e-4;

   bool inversion = false;
   bool visualization = false;
   bool verify_jacobian = false;
};

struct SolverSettings
{
   real_t outer_tolerance;
   real_t inner_tolerance;
   bool verbose;
};

Vector UnitVectorFromAngles(real_t azimuth_degrees, real_t dip_degrees)
{
   const real_t azimuth = azimuth_degrees * M_PI / 180.0;
   const real_t dip = dip_degrees * M_PI / 180.0;
   Vector direction(3);
   direction[0] = cos(dip) * cos(azimuth);
   direction[1] = cos(dip) * sin(azimuth);
   direction[2] = sin(dip);
   return direction;
}

void GetBoundingBox(const Mesh &mesh, Vector &minimum, Vector &maximum)
{
   const int dim = mesh.SpaceDimension();
   minimum.SetSize(dim);
   maximum.SetSize(dim);
   minimum = numeric_limits<real_t>::infinity();
   maximum = -numeric_limits<real_t>::infinity();
   for (int vertex = 0; vertex < mesh.GetNV(); vertex++)
   {
      const real_t *x = mesh.GetVertex(vertex);
      for (int d = 0; d < dim; d++)
      {
         minimum[d] = min(minimum[d], x[d]);
         maximum[d] = max(maximum[d], x[d]);
      }
   }
}

class CSEMConductivity : public Coefficient
{
private:
   const CSEMOptions &options;
   real_t length_scale;
   real_t background;
   real_t target;

public:
   CSEMConductivity(const CSEMOptions &options_, real_t length_scale_)
      : options(options_), length_scale(length_scale_),
        background(options.background_conductivity),
        target(options.target_conductivity)
   {
   }

   void SetLogModel(const Vector &model)
   {
      MFEM_VERIFY(model.Size() == 2, "The CSEM model has two parameters.");
      background = exp(model[0]);
      target = exp(model[1]);
   }

   bool InTarget(const Vector &scaled_x) const
   {
      const real_t x = length_scale * scaled_x[0];
      const real_t y = length_scale * scaled_x[1];
      const real_t z = length_scale * scaled_x[2];
      return fabs(x - options.target_x) <= 0.5 * options.target_size_x &&
             fabs(y - options.target_y) <= 0.5 * options.target_size_y &&
             fabs(z - options.target_z) <= 0.5 * options.target_size_z;
   }

   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override
   {
      Vector x;
      T.Transform(ip, x);
      if (length_scale * x[2] < 0.0) { return options.air_conductivity; }
      return InTarget(x) ? target : background;
   }

   real_t ParameterDerivative(const Vector &scaled_x, int parameter) const
   {
      if (length_scale * scaled_x[2] < 0.0) { return 0.0; }
      const bool target_point = InTarget(scaled_x);
      if (parameter == 0) { return target_point ? 0.0 : background; }
      return target_point ? target : 0.0;
   }
};

class RelativePermittivity : public Coefficient
{
private:
   const CSEMOptions &options;
   real_t length_scale;

public:
   RelativePermittivity(const CSEMOptions &options_, real_t length_scale_)
      : options(options_), length_scale(length_scale_)
   {
   }

   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override
   {
      Vector x;
      T.Transform(ip, x);
      return length_scale * x[2] < 0.0 ? 1.0
                                      : options.earth_relative_permittivity;
   }
};

class ScaledCoefficient : public Coefficient
{
private:
   Coefficient &coefficient;
   real_t scale;

public:
   ScaledCoefficient(Coefficient &coefficient_, real_t scale_)
      : coefficient(coefficient_), scale(scale_)
   {
   }

   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override
   {
      return scale * coefficient.Eval(T, ip);
   }
};

class ShiftedMassCoefficient : public Coefficient
{
private:
   CSEMConductivity &conductivity;
   RelativePermittivity &permittivity;
   real_t conductivity_scale;
   real_t permittivity_scale;

public:
   ShiftedMassCoefficient(CSEMConductivity &conductivity_,
                          RelativePermittivity &permittivity_,
                          real_t conductivity_scale_,
                          real_t permittivity_scale_)
      : conductivity(conductivity_), permittivity(permittivity_),
        conductivity_scale(conductivity_scale_),
        permittivity_scale(permittivity_scale_)
   {
   }

   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override
   {
      return conductivity_scale * conductivity.Eval(T, ip) -
             permittivity_scale * permittivity.Eval(T, ip);
   }
};

class ParameterDerivativeCoefficient : public Coefficient
{
private:
   CSEMConductivity &conductivity;
   int parameter;
   real_t scale;

public:
   ParameterDerivativeCoefficient(CSEMConductivity &conductivity_,
                                  int parameter_, real_t scale_)
      : conductivity(conductivity_), parameter(parameter_), scale(scale_)
   {
   }

   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override
   {
      Vector x;
      T.Transform(ip, x);
      return scale * conductivity.ParameterDerivative(x, parameter);
   }
};

class ReciprocalCoefficient : public Coefficient
{
private:
   Coefficient &coefficient;

public:
   explicit ReciprocalCoefficient(Coefficient &coefficient_)
      : coefficient(coefficient_)
   {
   }

   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override
   {
      return 1.0 / coefficient.Eval(T, ip);
   }
};

class Log10Coefficient : public Coefficient
{
private:
   Coefficient &coefficient;

public:
   explicit Log10Coefficient(Coefficient &coefficient_)
      : coefficient(coefficient_)
   {
   }

   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override
   {
      return log10(coefficient.Eval(T, ip));
   }
};

class GaussianWireCurrent : public VectorCoefficient
{
private:
   struct Segment
   {
      Vector center;
      Vector direction;
      real_t half_length;
   };

   const CSEMOptions &options;
   real_t length_scale;
   real_t right_hand_side_scale;
   vector<Segment> segments;

   void AddSegment(const Vector &start, const Vector &end)
   {
      Segment segment;
      segment.center.SetSize(3);
      segment.direction.SetSize(3);
      segment.center = start;
      segment.center += end;
      segment.center *= 0.5;
      segment.direction = end;
      segment.direction -= start;
      const real_t length = segment.direction.Norml2();
      MFEM_VERIFY(length > 0.0,
                  "Consecutive transmitter vertices must be distinct.");
      segment.direction /= length;
      segment.half_length = 0.5 * length;
      segments.push_back(segment);
   }

public:
   GaussianWireCurrent(const CSEMOptions &options_, real_t length_scale_,
                       real_t omega)
      : VectorCoefficient(3), options(options_), length_scale(length_scale_),
        right_hand_side_scale(-omega * mu0 * length_scale_ * length_scale_)
   {
      vector<Vector> vertices;
      if (!string(options.source_file).empty())
      {
         ifstream input(options.source_file);
         MFEM_VERIFY(input, "Unable to open transmitter geometry file: "
                     << options.source_file);
         string line;
         while (getline(input, line))
         {
            replace(line.begin(), line.end(), ',', ' ');
            istringstream row(line);
            Vector vertex(3);
            if (row >> vertex[0] >> vertex[1] >> vertex[2])
            {
               vertices.push_back(vertex);
            }
         }
         MFEM_VERIFY(vertices.size() >= 2,
                     "A transmitter polyline needs at least two vertices.");
      }
      else
      {
         const Vector direction =
            UnitVectorFromAngles(options.source_azimuth, options.source_dip);
         Vector center(3), start(3), end(3);
         center[0] = options.source_x;
         center[1] = options.source_y;
         center[2] = options.source_z;
         start = center;
         end = center;
         start.Add(-0.5 * options.source_length, direction);
         end.Add(0.5 * options.source_length, direction);
         vertices.push_back(start);
         vertices.push_back(end);
      }
      for (int vertex = 0;
           vertex + 1 < static_cast<int>(vertices.size()); vertex++)
      {
         AddSegment(vertices[vertex], vertices[vertex + 1]);
      }
   }

   void Eval(Vector &value, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      Vector scaled_x;
      T.Transform(ip, scaled_x);
      const real_t radius2 = options.source_radius * options.source_radius;
      value.SetSize(3);
      value = 0.0;
      for (const Segment &segment : segments)
      {
         Vector displacement(3);
         for (int d = 0; d < 3; d++)
         {
            displacement[d] = length_scale * scaled_x[d] - segment.center[d];
         }
         const real_t axial_distance = displacement * segment.direction;
         const real_t transverse2 = max(
            displacement * displacement - axial_distance * axial_distance,
            real_t(0.0));
         // Smooth each segment end over one source radius. Connected segments
         // then form a stable regularized line current even on coarse meshes.
         const real_t longitudinal =
            0.5 * (erf((axial_distance + segment.half_length) /
                       options.source_radius) -
                   erf((axial_distance - segment.half_length) /
                       options.source_radius));
         const real_t profile =
            longitudinal * exp(-transverse2 / radius2) / (M_PI * radius2);
         value.Add(right_hand_side_scale * options.source_current * profile,
                   segment.direction);
      }
   }
};

class GaussianReceiver : public VectorCoefficient
{
private:
   Vector center;
   real_t length_scale;
   real_t radius;
   real_t normalization;
   Vector direction;

public:
   GaussianReceiver(const Vector &center_, real_t length_scale_,
                    real_t radius_, const Vector &direction_)
      : VectorCoefficient(3), center(center_), length_scale(length_scale_),
        radius(radius_), direction(direction_)
   {
      normalization = pow(length_scale, 3) /
                      (pow(M_PI, 1.5) * radius * radius * radius);
   }

   void Eval(Vector &value, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      Vector scaled_x;
      T.Transform(ip, scaled_x);
      real_t distance2 = 0.0;
      for (int d = 0; d < 3; d++)
      {
         const real_t delta = length_scale * scaled_x[d] - center[d];
         distance2 += delta * delta;
      }
      value.SetSize(3);
      value = direction;
      value *= normalization * exp(-distance2 / (radius * radius));
   }
};

class GaussNewtonNormalOperator : public Operator
{
private:
   const DenseMatrix &jacobian;
   const Vector &standard_deviation;
   real_t diagonal_shift;
   mutable Vector weighted_action;

public:
   GaussNewtonNormalOperator(const DenseMatrix &jacobian_,
                             const Vector &standard_deviation_,
                             real_t diagonal_shift_)
      : Operator(jacobian_.Width()), jacobian(jacobian_),
        standard_deviation(standard_deviation_),
        diagonal_shift(diagonal_shift_),
        weighted_action(jacobian_.Height())
   {
   }

   void Mult(const Vector &x, Vector &y) const override
   {
      jacobian.Mult(x, weighted_action);
      for (int i = 0; i < weighted_action.Size(); i++)
      {
         const real_t variance =
            standard_deviation[i] * standard_deviation[i];
         weighted_action[i] /= variance;
      }
      jacobian.MultTranspose(weighted_action, y);
      y.Add(diagonal_shift, x);
   }
};

class CSEMProblem
{
private:
   MPI_Comm communicator;
   int rank;
   const CSEMOptions &options;
   real_t length_scale = 1.0;
   unique_ptr<ParMesh> mesh;
   unique_ptr<ND_FECollection> collection;
   unique_ptr<ParFiniteElementSpace> space;
   unique_ptr<ParComplexGridFunction> field;
   unique_ptr<RT_FECollection> magnetic_collection;
   unique_ptr<ParFiniteElementSpace> magnetic_space;
   unique_ptr<ParDiscreteLinearOperator> curl;
   unique_ptr<CSEMConductivity> conductivity;
   unique_ptr<RelativePermittivity> permittivity;
   Array<int> essential_tdofs;
   vector<Vector> receiver_centers;
   vector<Vector> receiver_directions;
   vector<unique_ptr<GaussianReceiver>> receiver_coefficients;
   vector<unique_ptr<ParLinearForm>> receiver_functionals;

   bool MagneticReceivers() const
   {
      return string(options.receiver_field) == "H" ||
             string(options.receiver_field) == "h";
   }

   void BuildMesh()
   {
      unique_ptr<Mesh> serial_mesh;
      if (string(options.mesh_file).empty())
      {
         serial_mesh.reset(new Mesh(Mesh::MakeCartesian3D(
            options.nx, options.ny, options.nz, Element::HEXAHEDRON,
            options.size_x, options.size_y,
            options.air_thickness + options.earth_depth)));
         const real_t shift_x = -0.5 * options.size_x;
         const real_t shift_y = -0.5 * options.size_y;
         const real_t shift_z = -options.air_thickness;
         serial_mesh->Transform(
            [=](const Vector &input, Vector &output)
            {
               output = input;
               output[0] += shift_x;
               output[1] += shift_y;
               output[2] += shift_z;
            });
      }
      else
      {
         serial_mesh.reset(new Mesh(options.mesh_file, 1, 1));
      }

      MFEM_VERIFY(serial_mesh->Dimension() == 3,
                  "The CSEM miniapp requires a three-dimensional mesh.");
      for (int element = 0; element < serial_mesh->GetNE(); element++)
      {
         MFEM_VERIFY(serial_mesh->GetElementBaseGeometry(element) ==
                        Geometry::CUBE,
                     "The spectral CSEM miniapp requires hexahedra.");
      }
      for (int level = 0; level < options.refinements; level++)
      {
         serial_mesh->UniformRefinement();
      }

      Vector minimum, maximum;
      GetBoundingBox(*serial_mesh, minimum, maximum);
      length_scale = 0.0;
      for (int d = 0; d < 3; d++)
      {
         length_scale = max(length_scale, maximum[d] - minimum[d]);
      }
      MFEM_VERIFY(length_scale > 0.0, "The mesh has zero physical extent.");

      const real_t scale = length_scale;
      serial_mesh->Transform(
         [=](const Vector &input, Vector &output)
         {
            output = input;
            output /= scale;
         });
      mesh.reset(new ParMesh(communicator, *serial_mesh));
   }

   void BuildReceivers()
   {
      auto add_receiver = [&](const Vector &center, const Vector &direction)
      {
         receiver_centers.push_back(center);
         receiver_directions.push_back(direction);
         receiver_coefficients.emplace_back(
            new GaussianReceiver(center, length_scale,
                                 options.receiver_radius, direction));
         ParFiniteElementSpace *receiver_space =
            MagneticReceivers() ? magnetic_space.get() : space.get();
         receiver_functionals.emplace_back(new ParLinearForm(receiver_space));
         receiver_functionals.back()->AddDomainIntegrator(
            new VectorFEDomainLFIntegrator(*receiver_coefficients.back()));
         receiver_functionals.back()->Assemble();
      };

      if (!string(options.receiver_file).empty())
      {
         ifstream input(options.receiver_file);
         MFEM_VERIFY(input, "Unable to open receiver geometry file: "
                     << options.receiver_file);
         string line;
         while (getline(input, line))
         {
            replace(line.begin(), line.end(), ',', ' ');
            istringstream row(line);
            vector<real_t> values;
            real_t value;
            while (row >> value) { values.push_back(value); }
            if (values.empty()) { continue; }
            MFEM_VERIFY(values.size() == 3 || values.size() == 5,
                        "Each receiver row must contain x,y,z or "
                        "x,y,z,azimuth,dip.");
            Vector center(3);
            center[0] = values[0];
            center[1] = values[1];
            center[2] = values[2];
            const Vector direction = values.size() == 5
               ? UnitVectorFromAngles(values[3], values[4])
               : UnitVectorFromAngles(options.receiver_azimuth,
                                      options.receiver_dip);
            add_receiver(center, direction);
         }
         MFEM_VERIFY(!receiver_centers.empty(),
                     "The receiver geometry file contains no receivers.");
         return;
      }

      const Vector direction = UnitVectorFromAngles(options.receiver_azimuth,
                                                    options.receiver_dip);
      const real_t line_azimuth =
         options.receiver_line_azimuth * M_PI / 180.0;
      const real_t line_x = cos(line_azimuth);
      const real_t line_y = sin(line_azimuth);
      for (int receiver = 0; receiver < options.receivers; receiver++)
      {
         Vector center(3);
         const real_t fraction = options.receivers == 1
                                    ? 0.5
                                    : real_t(receiver) /
                                         real_t(options.receivers - 1);
         const real_t offset = -0.5 * options.receiver_span +
                               fraction * options.receiver_span;
         center[0] = options.receiver_x + offset * line_x;
         center[1] = options.receiver_y + offset * line_y;
         center[2] = options.receiver_z;
         add_receiver(center, direction);
      }
   }

   void ExtractResponses(const ParComplexGridFunction &solution,
                         real_t frequency,
                         Vector &responses) const
   {
      responses.SetSize(2 * receiver_functionals.size());
      unique_ptr<ParGridFunction> magnetic_real;
      unique_ptr<ParGridFunction> magnetic_imaginary;
      if (MagneticReceivers())
      {
         magnetic_real.reset(new ParGridFunction(magnetic_space.get()));
         magnetic_imaginary.reset(new ParGridFunction(magnetic_space.get()));
         curl->Mult(solution.imag(), *magnetic_real);
         curl->Mult(solution.real(), *magnetic_imaginary);
         const real_t magnetic_scale =
            1.0 / (2.0 * M_PI * frequency * mu0 *
                   options.relative_permeability * length_scale);
         *magnetic_real *= -magnetic_scale;
         *magnetic_imaginary *= magnetic_scale;
      }
      for (int receiver = 0;
           receiver < static_cast<int>(receiver_functionals.size());
           receiver++)
      {
         if (MagneticReceivers())
         {
            responses[2 * receiver] =
               (*receiver_functionals[receiver])(*magnetic_real);
            responses[2 * receiver + 1] =
               (*receiver_functionals[receiver])(*magnetic_imaginary);
         }
         else
         {
            responses[2 * receiver] =
               (*receiver_functionals[receiver])(solution.real());
            responses[2 * receiver + 1] =
               (*receiver_functionals[receiver])(solution.imag());
         }
      }
   }

public:
   CSEMProblem(MPI_Comm communicator_, const CSEMOptions &options_)
      : communicator(communicator_), options(options_)
   {
      MPI_Comm_rank(communicator, &rank);
      BuildMesh();
      collection.reset(new ND_FECollection(
         options.order, 3, BasisType::GaussLobatto,
         BasisType::IntegratedGLL));
      space.reset(new ParFiniteElementSpace(mesh.get(), collection.get()));
      field.reset(new ParComplexGridFunction(space.get()));
      if (MagneticReceivers())
      {
         magnetic_collection.reset(new RT_FECollection(options.order - 1, 3));
         magnetic_space.reset(
            new ParFiniteElementSpace(mesh.get(), magnetic_collection.get()));
         curl.reset(new ParDiscreteLinearOperator(space.get(),
                                                  magnetic_space.get()));
         curl->AddDomainInterpolator(new CurlInterpolator);
         curl->Assemble();
         curl->Finalize();
      }
      conductivity.reset(new CSEMConductivity(options, length_scale));
      permittivity.reset(new RelativePermittivity(options, length_scale));

      Array<int> essential_boundary(mesh->bdr_attributes.Max());
      essential_boundary = 1;
      space->GetEssentialTrueDofs(essential_boundary, essential_tdofs);
      BuildReceivers();

      if (rank == 0)
      {
         cout << "CSEM mesh: " << mesh->GetGlobalNE() << " hexahedra, "
              << space->GlobalTrueVSize() << " true H(curl) DOFs\n"
              << "Length scaling L = " << length_scale << " m\n";
      }
   }

   int Rank() const { return rank; }
   int DataSize() const { return 2 * static_cast<int>(receiver_centers.size()); }
   const Vector &ReceiverCenter(int receiver) const
   {
      return receiver_centers[receiver];
   }
   const Vector &ReceiverDirection(int receiver) const
   {
      return receiver_directions[receiver];
   }

   void Solve(const Vector &log_model, real_t frequency, Vector &responses,
              DenseMatrix *jacobian, const SolverSettings &settings)
   {
      conductivity->SetLogModel(log_model);
      const real_t omega = 2.0 * M_PI * frequency;
      const real_t scaled_conductivity =
         omega * mu0 * length_scale * length_scale;
      const real_t wave_number = omega * sqrt(mu0 * epsilon0);
      const real_t scaled_permittivity =
         wave_number * wave_number * length_scale * length_scale;

      ConstantCoefficient inverse_relative_permeability(
         1.0 / options.relative_permeability);
      ScaledCoefficient conduction(*conductivity, scaled_conductivity);
      ScaledCoefficient dielectric(*permittivity, -scaled_permittivity);
      ShiftedMassCoefficient shifted_mass(
         *conductivity, *permittivity, scaled_conductivity,
         scaled_permittivity);

      GaussianWireCurrent source(options, length_scale, omega);
      ParComplexLinearForm right_hand_side(
         space.get(), ComplexOperator::HERMITIAN);
      right_hand_side.AddDomainIntegrator(
         nullptr, new VectorFEDomainLFIntegrator(source));
      right_hand_side.Assemble();

      *field = complex<real_t>(0.0, 0.0);
      ParSesquilinearForm form(space.get(), ComplexOperator::HERMITIAN);
      form.AddDomainIntegrator(
         new CurlCurlIntegrator(inverse_relative_permeability), nullptr);
      form.AddDomainIntegrator(new VectorFEMassIntegrator(dielectric),
                               new VectorFEMassIntegrator(conduction));
      form.Assemble();

      OperatorHandle system_operator;
      Vector X, B;
      form.FormLinearSystem(essential_tdofs, *field, right_hand_side,
                            system_operator, X, B);
      ComplexHypreParMatrix &complex_matrix =
         *system_operator.As<ComplexHypreParMatrix>();

      unique_ptr<Solver> direct_solver;
      unique_ptr<HypreParMatrix> real_system_matrix;
      unique_ptr<ParBilinearForm> shifted_form;
      OperatorHandle shifted_operator;
      unique_ptr<LORSolver<HypreAMS>> ams;
      unique_ptr<FGMRESSolver> inner_solver;
      unique_ptr<Solver> direct_inner_solver;
      unique_ptr<WeissBlockOperator> weiss_operator;
      unique_ptr<PRESBPreconditioner> presb;
      unique_ptr<FGMRESSolver> outer_solver;

      if (string(options.solver) == "mumps")
      {
#ifdef MFEM_USE_COMPLEX_MUMPS
         auto *mumps = new ComplexMUMPSSolver(communicator);
         mumps->SetPrintLevel(0);
         mumps->SetMatrixSymType(ComplexMUMPSSolver::UNSYMMETRIC);
         mumps->SetOperator(complex_matrix);
         direct_solver.reset(mumps);
#elif defined(MFEM_USE_MUMPS)
         real_system_matrix.reset(complex_matrix.GetSystemMatrix());
         auto *mumps = new MUMPSSolver(communicator);
         mumps->SetPrintLevel(0);
         mumps->SetMatrixSymType(MUMPSSolver::UNSYMMETRIC);
         mumps->SetOperator(*real_system_matrix);
         direct_solver.reset(mumps);
#else
         MFEM_ABORT("The requested MUMPS solver is not enabled.");
#endif
      }
      else
      {
         MFEM_VERIFY(string(options.solver) == "presb" ||
                        string(options.solver) == "presb-mumps",
                     "Solver must be 'presb', 'presb-mumps', or 'mumps'.");
         shifted_form.reset(new ParBilinearForm(space.get()));
         shifted_form->AddDomainIntegrator(
            new CurlCurlIntegrator(inverse_relative_permeability));
         shifted_form->AddDomainIntegrator(
            new VectorFEMassIntegrator(shifted_mass));
         shifted_form->Assemble();
         shifted_form->FormSystemMatrix(essential_tdofs, shifted_operator);
         HypreParMatrix &H = *shifted_operator.As<HypreParMatrix>();

         Solver *presb_inner_solver = nullptr;
         if (string(options.solver) == "presb-mumps")
         {
#ifdef MFEM_USE_MUMPS
            auto *mumps = new MUMPSSolver(communicator);
            mumps->SetPrintLevel(0);
            mumps->SetMatrixSymType(
               MUMPSSolver::SYMMETRIC_POSITIVE_DEFINITE);
            mumps->SetOperator(H);
            direct_inner_solver.reset(mumps);
            presb_inner_solver = direct_inner_solver.get();
#else
            MFEM_ABORT("The requested MUMPS inner solver is not enabled.");
#endif
         }
         else
         {
            // AMS is applied on the low-order-refined auxiliary operator.
            // This retains the paper's auxiliary-space method while making
            // the high-order GLL discretization practical.
            ams.reset(new LORSolver<HypreAMS>(*shifted_form,
                                              essential_tdofs));
            HypreAMS &hypre_ams = ams->GetSolver();
            HYPRE_Solver hypre_solver = static_cast<HYPRE_Solver>(hypre_ams);
            HYPRE_AMSSetCycleType(hypre_solver, 14);
            HYPRE_AMSSetMaxIter(hypre_solver, 2);
            HYPRE_AMSSetTol(hypre_solver, 0.0);
            HYPRE_AMSSetSmoothingOptions(hypre_solver, 2, 1, 1.0, 1.0);
            hypre_ams.SetPrintLevel(0);
            inner_solver.reset(new FGMRESSolver(communicator));
            inner_solver->SetOperator(H);
            inner_solver->SetPreconditioner(*ams);
            inner_solver->SetRelTol(settings.inner_tolerance);
            inner_solver->SetAbsTol(0.0);
            inner_solver->SetMaxIter(options.inner_max_iterations);
            inner_solver->SetKDim(options.inner_kdim);
            inner_solver->SetPrintLevel(0);
            presb_inner_solver = inner_solver.get();
         }

         weiss_operator.reset(new WeissBlockOperator(
            complex_matrix.real(), complex_matrix.imag()));
         presb.reset(new PRESBPreconditioner(complex_matrix.imag(),
                                             *presb_inner_solver));
         outer_solver.reset(new FGMRESSolver(communicator));
         outer_solver->SetOperator(*weiss_operator);
         outer_solver->SetPreconditioner(*presb);
         outer_solver->SetRelTol(settings.outer_tolerance);
         outer_solver->SetAbsTol(0.0);
         outer_solver->SetMaxIter(options.outer_max_iterations);
         outer_solver->SetKDim(options.outer_kdim);
         outer_solver->SetPrintLevel(settings.verbose && rank == 0 ? 1 : 0);
      }

      int linear_solve_number = 0;
      auto solve_linear_system = [&](const Vector &complex_rhs,
                                     Vector &complex_solution)
      {
         complex_solution.SetSize(complex_rhs.Size());
         complex_solution = 0.0;
         if (direct_solver)
         {
            direct_solver->Mult(complex_rhs, complex_solution);
         }
         else
         {
            const int block_size = complex_rhs.Size() / 2;
            Vector rhs_real, rhs_imaginary, weiss_rhs(2 * block_size);
            Vector weiss_first, weiss_second, weiss_solution(2 * block_size);
            rhs_real.MakeRef(const_cast<Vector &>(complex_rhs), 0, block_size);
            rhs_imaginary.MakeRef(const_cast<Vector &>(complex_rhs),
                                  block_size, block_size);
            weiss_first.MakeRef(weiss_rhs, 0, block_size);
            weiss_second.MakeRef(weiss_rhs, block_size, block_size);
            weiss_first = rhs_imaginary;
            weiss_second = rhs_real;
            weiss_solution = 0.0;
            outer_solver->Mult(weiss_rhs, weiss_solution);
            MFEM_VERIFY(outer_solver->GetConverged(),
                        "The outer PRESB solve did not converge.");

            Vector solution_real, solution_imaginary;
            Vector z1, z2;
            solution_real.MakeRef(complex_solution, 0, block_size);
            solution_imaginary.MakeRef(complex_solution, block_size,
                                       block_size);
            z1.MakeRef(weiss_solution, 0, block_size);
            z2.MakeRef(weiss_solution, block_size, block_size);
            solution_real = z1;
            solution_imaginary = z2;
            solution_imaginary *= -1.0;
         }
         linear_solve_number++;
      };

      solve_linear_system(B, X);
      form.RecoverFEMSolution(X, right_hand_side, *field);
      field->Sync();
      ExtractResponses(*field, frequency, responses);

      Vector AX(B.Size());
      system_operator->Mult(X, AX);
      AX -= B;
      const real_t rhs_norm =
         sqrt(max(InnerProduct(communicator, B, B), real_t(0.0)));
      const real_t residual_norm =
         sqrt(max(InnerProduct(communicator, AX, AX), real_t(0.0)));
      if (settings.verbose && rank == 0)
      {
         cout << "Forward relative residual: "
              << residual_norm / max(rhs_norm, numeric_limits<real_t>::min())
              << '\n';
      }

      if (jacobian)
      {
         jacobian->SetSize(DataSize(), 2);
         const int block_size = space->GetTrueVSize();
         Vector electric_real, electric_imaginary;
         electric_real.MakeRef(X, 0, block_size);
         electric_imaginary.MakeRef(X, block_size, block_size);

         for (int parameter = 0; parameter < 2; parameter++)
         {
            ParameterDerivativeCoefficient derivative_coefficient(
               *conductivity, parameter, scaled_conductivity);
            ParBilinearForm derivative_form(space.get());
            derivative_form.SetDiagonalPolicy(Matrix::DIAG_ZERO);
            derivative_form.AddDomainIntegrator(
               new VectorFEMassIntegrator(derivative_coefficient));
            derivative_form.Assemble();
            OperatorHandle derivative_operator;
            derivative_form.FormSystemMatrix(essential_tdofs,
                                              derivative_operator);

            Vector tangent_rhs(2 * block_size);
            Vector tangent_real, tangent_imaginary;
            tangent_real.MakeRef(tangent_rhs, 0, block_size);
            tangent_imaginary.MakeRef(tangent_rhs, block_size, block_size);
            derivative_operator->Mult(electric_imaginary, tangent_real);
            derivative_operator->Mult(electric_real, tangent_imaginary);
            tangent_imaginary *= -1.0;

            Vector tangent_solution;
            solve_linear_system(tangent_rhs, tangent_solution);
            ParComplexGridFunction tangent_field(space.get());
            tangent_field = complex<real_t>(0.0, 0.0);
            form.RecoverFEMSolution(tangent_solution, right_hand_side,
                                    tangent_field);
            tangent_field.Sync();
            Vector tangent_responses;
            ExtractResponses(tangent_field, frequency, tangent_responses);
            for (int datum = 0; datum < DataSize(); datum++)
            {
               (*jacobian)(datum, parameter) = tangent_responses[datum];
            }
         }
      }

      if (settings.verbose && rank == 0 && presb)
      {
         cout << "PRESB linear solves: " << linear_solve_number
              << ", applications: " << presb->Applications()
              << ", accumulated inner iterations: "
              << presb->InnerIterations() << '\n';
      }
   }

   void SaveField(const string &prefix, real_t frequency)
   {
      RT_FECollection curl_collection(options.order - 1, 3);
      ParFiniteElementSpace curl_space(mesh.get(), &curl_collection);
      ParDiscreteLinearOperator curl(space.get(), &curl_space);
      curl.AddDomainInterpolator(new CurlInterpolator);
      curl.Assemble();
      curl.Finalize();

      ParGridFunction magnetic_real(&curl_space);
      ParGridFunction magnetic_imaginary(&curl_space);
      curl.Mult(field->imag(), magnetic_real);
      curl.Mult(field->real(), magnetic_imaginary);
      const real_t magnetic_scale =
         1.0 / (2.0 * M_PI * frequency * mu0 *
                options.relative_permeability * length_scale);
      magnetic_real *= -magnetic_scale;
      magnetic_imaginary *= magnetic_scale;

      L2_FECollection material_collection(0, 3);
      ParFiniteElementSpace material_space(mesh.get(), &material_collection);
      ParGridFunction sigma(&material_space);
      ParGridFunction log10_sigma(&material_space);
      ParGridFunction resistivity(&material_space);
      ParGridFunction epsilon_r(&material_space);
      ReciprocalCoefficient reciprocal_sigma(*conductivity);
      Log10Coefficient logarithmic_sigma(*conductivity);
      sigma.ProjectCoefficient(*conductivity);
      log10_sigma.ProjectCoefficient(logarithmic_sigma);
      resistivity.ProjectCoefficient(reciprocal_sigma);
      epsilon_r.ProjectCoefficient(*permittivity);

      ParaViewDataCollection collection_output(prefix.c_str(), mesh.get());
      collection_output.SetDataFormat(VTKFormat::BINARY);
      collection_output.SetHighOrderOutput(true);
      collection_output.SetLevelsOfDetail(options.order);
      collection_output.SetCycle(0);
      collection_output.SetTime(frequency);
      collection_output.RegisterField("E_real", &field->real());
      collection_output.RegisterField("E_imaginary", &field->imag());
      collection_output.RegisterField("H_real", &magnetic_real);
      collection_output.RegisterField("H_imaginary", &magnetic_imaginary);
      collection_output.RegisterField("conductivity_S_per_m", &sigma);
      collection_output.RegisterField("log10_conductivity", &log10_sigma);
      collection_output.RegisterField("resistivity_ohm_m", &resistivity);
      collection_output.RegisterField("relative_permittivity", &epsilon_r);
      collection_output.Save();
   }
};

vector<real_t> ParseFrequencies(const CSEMOptions &options)
{
   vector<real_t> result;
   string text(options.frequencies);
   replace(text.begin(), text.end(), ',', ' ');
   replace(text.begin(), text.end(), ';', ' ');
   istringstream input(text);
   real_t frequency;
   while (input >> frequency)
   {
      MFEM_VERIFY(frequency > 0.0, "Frequencies must be positive.");
      result.push_back(frequency);
   }
   if (result.empty()) { result.push_back(options.frequency); }
   return result;
}

real_t PlaneWaveSkinDepth(const CSEMOptions &options, real_t frequency)
{
   const real_t omega = 2.0 * M_PI * frequency;
   const real_t permeability = mu0 * options.relative_permeability;
   const real_t permittivity =
      epsilon0 * options.earth_relative_permittivity;
   const real_t ratio = options.background_conductivity /
                        (omega * permittivity);
   const real_t attenuation = omega * sqrt(0.5 * permeability * permittivity) *
      sqrt(sqrt(1.0 + ratio * ratio) - 1.0);
   return 1.0 / attenuation;
}

string FrequencyFieldPrefix(const string &base, int frequency_index,
                            int frequency_count, real_t frequency)
{
   if (frequency_count == 1) { return base; }
   ostringstream frequency_text;
   frequency_text << setprecision(8) << frequency;
   string label = frequency_text.str();
   replace(label.begin(), label.end(), '.', 'p');
   replace(label.begin(), label.end(), '+', 'p');
   ostringstream name;
   name << base << "_f" << setw(3) << setfill('0') << frequency_index
        << '_' << label << "Hz";
   return name.str();
}

void SolveSurvey(CSEMProblem &problem, const vector<real_t> &frequencies,
                 const Vector &log_model, Vector &responses,
                 DenseMatrix *jacobian, const SolverSettings &settings,
                 const string &field_prefix = "")
{
   const int block_size = problem.DataSize();
   const int frequency_count = static_cast<int>(frequencies.size());
   responses.SetSize(block_size * frequency_count);
   if (jacobian) { jacobian->SetSize(responses.Size(), 2); }

   for (int frequency_index = 0;
         frequency_index < frequency_count;
        frequency_index++)
   {
      const real_t frequency = frequencies[frequency_index];
      if (settings.verbose && problem.Rank() == 0)
      {
         cout << "Solving frequency " << frequency << " Hz ("
              << frequency_index + 1 << '/' << frequencies.size() << ")\n";
      }
      Vector frequency_responses;
      DenseMatrix frequency_jacobian;
      problem.Solve(log_model, frequency, frequency_responses,
                    jacobian ? &frequency_jacobian : nullptr, settings);
      for (int datum = 0; datum < block_size; datum++)
      {
         const int global_datum = frequency_index * block_size + datum;
         responses[global_datum] = frequency_responses[datum];
         if (jacobian)
         {
            for (int parameter = 0; parameter < 2; parameter++)
            {
               (*jacobian)(global_datum, parameter) =
                  frequency_jacobian(datum, parameter);
            }
         }
      }
      if (!field_prefix.empty())
      {
         problem.SaveField(FrequencyFieldPrefix(
                              field_prefix, frequency_index,
                              frequency_count, frequency),
                           frequency);
      }
   }
}

real_t Objective(const Vector &prediction, const Vector &observed,
                 const Vector &standard_deviation, const Vector &model,
                 const Vector &reference_model, real_t regularization)
{
   real_t value = 0.0;
   for (int datum = 0; datum < prediction.Size(); datum++)
   {
      const real_t residual =
         (prediction[datum] - observed[datum]) / standard_deviation[datum];
      value += residual * residual;
   }
   for (int parameter = 0; parameter < model.Size(); parameter++)
   {
      const real_t difference = model[parameter] - reference_model[parameter];
      value += regularization * difference * difference;
   }
   return 0.5 * value;
}

void WriteData(const string &file_name, const CSEMProblem &problem,
               const vector<real_t> &frequencies, const Vector &data,
               const Vector &standard_deviation)
{
   if (problem.Rank() != 0) { return; }
   ofstream output(file_name);
   MFEM_VERIFY(output, "Unable to open output data file: " << file_name);
   output << "frequency_hz,receiver,x_m,y_m,z_m,direction_x,direction_y,"
             "direction_z,real,imag,"
             "stddev_real,stddev_imag\n";
   output << setprecision(16);
   const int frequency_block_size = problem.DataSize();
   MFEM_VERIFY(data.Size() == frequency_block_size * frequencies.size(),
               "Survey data size does not match the frequency list.");
   for (int frequency_index = 0;
        frequency_index < static_cast<int>(frequencies.size());
        frequency_index++)
   {
      for (int receiver = 0; receiver < frequency_block_size / 2; receiver++)
      {
         const int offset = frequency_index * frequency_block_size +
                            2 * receiver;
         const Vector &center = problem.ReceiverCenter(receiver);
         const Vector &direction = problem.ReceiverDirection(receiver);
         output << frequencies[frequency_index] << ',' << receiver << ','
                << center[0] << ',' << center[1] << ',' << center[2] << ','
                << direction[0] << ',' << direction[1] << ','
                << direction[2] << ',' << data[offset] << ','
                << data[offset + 1] << ','
                << standard_deviation[offset] << ','
                << standard_deviation[offset + 1] << '\n';
      }
   }
}

void ReadObservedData(const string &file_name, int expected_rows,
                      Vector &observed, Vector &standard_deviation)
{
   ifstream input(file_name);
   MFEM_VERIFY(input, "Unable to open observed data file: " << file_name);
   vector<real_t> data_values;
   vector<real_t> error_values;
   string line;
   while (getline(input, line))
   {
      replace(line.begin(), line.end(), ',', ' ');
      istringstream row(line);
      vector<real_t> values;
      real_t value;
      while (row >> value) { values.push_back(value); }
      if (values.size() == 12)
      {
         data_values.push_back(values[8]);
         data_values.push_back(values[9]);
         error_values.push_back(values[10]);
         error_values.push_back(values[11]);
      }
      else if (values.size() == 9)
      {
         // Backward-compatible reader for pre-orientation multifrequency CSV.
         data_values.push_back(values[5]);
         data_values.push_back(values[6]);
         error_values.push_back(values[7]);
         error_values.push_back(values[8]);
      }
      else if (values.size() == 8)
      {
         // Backward-compatible reader for the original single-frequency CSV.
         data_values.push_back(values[4]);
         data_values.push_back(values[5]);
         error_values.push_back(values[6]);
         error_values.push_back(values[7]);
      }
   }
   MFEM_VERIFY(int(data_values.size()) == 2 * expected_rows,
               "Observed CSV row count does not match frequencies times "
               "--receivers.");
   observed.SetSize(data_values.size());
   standard_deviation.SetSize(error_values.size());
   for (int i = 0; i < observed.Size(); i++)
   {
      observed[i] = data_values[i];
      standard_deviation[i] = error_values[i];
      MFEM_VERIFY(standard_deviation[i] > 0.0,
                  "Observed standard deviations must be positive.");
   }
}

void BuildSyntheticData(CSEMProblem &problem, const CSEMOptions &options,
                        const vector<real_t> &frequencies,
                        const SolverSettings &forward_settings,
                        Vector &observed, Vector &standard_deviation)
{
   Vector true_model(2);
   true_model[0] = log(options.background_conductivity);
   true_model[1] = log(options.target_conductivity);
   SolveSurvey(problem, frequencies, true_model, observed, nullptr,
               forward_settings);

   real_t maximum_amplitude = 0.0;
   for (int datum = 0; datum < observed.Size(); datum++)
   {
      maximum_amplitude = max(maximum_amplitude, fabs(observed[datum]));
   }
   standard_deviation.SetSize(observed.Size());
   mt19937 generator(options.random_seed);
   normal_distribution<real_t> normal(0.0, 1.0);
   for (int datum = 0; datum < observed.Size(); datum++)
   {
      const real_t scale = max(fabs(observed[datum]),
                               options.data_error_floor * maximum_amplitude);
      standard_deviation[datum] = options.data_relative_error * scale;
      observed[datum] += options.synthetic_noise * scale * normal(generator);
   }
}

void VerifyJacobian(CSEMProblem &problem, const CSEMOptions &options,
                    const vector<real_t> &frequencies,
                    const Vector &model, const SolverSettings &settings)
{
   Vector prediction;
   DenseMatrix analytic;
   SolveSurvey(problem, frequencies, model, prediction, &analytic, settings);

   DenseMatrix finite_difference(prediction.Size(), model.Size());
   for (int parameter = 0; parameter < model.Size(); parameter++)
   {
      Vector plus_model(model), minus_model(model);
      plus_model[parameter] += options.jacobian_fd_step;
      minus_model[parameter] -= options.jacobian_fd_step;
      Vector plus_response, minus_response;
      SolveSurvey(problem, frequencies, plus_model, plus_response, nullptr,
                  settings);
      SolveSurvey(problem, frequencies, minus_model, minus_response, nullptr,
                  settings);
      for (int datum = 0; datum < prediction.Size(); datum++)
      {
         finite_difference(datum, parameter) =
            (plus_response[datum] - minus_response[datum]) /
            (2.0 * options.jacobian_fd_step);
      }
   }

   if (problem.Rank() != 0) { return; }
   const string file_name =
      string(options.output_prefix) + "_jacobian_check.csv";
   ofstream output(file_name);
   MFEM_VERIFY(output, "Unable to open Jacobian verification output.");
   output << "frequency_hz,receiver,part,parameter,analytic,"
             "finite_difference,absolute_error,relative_error\n";
   output << setprecision(16);
   real_t maximum_relative_error = 0.0;
   real_t squared_error = 0.0;
   real_t squared_reference = 0.0;
   const int block_size = problem.DataSize();
   const char *parts[2] = {"real", "imag"};
   const char *parameters[2] = {"log_sigma_background",
                                "log_sigma_target"};
   for (int datum = 0; datum < prediction.Size(); datum++)
   {
      const int frequency_index = datum / block_size;
      const int local_datum = datum % block_size;
      const int receiver = local_datum / 2;
      const int part = local_datum % 2;
      for (int parameter = 0; parameter < model.Size(); parameter++)
      {
         const real_t exact = analytic(datum, parameter);
         const real_t approximation = finite_difference(datum, parameter);
         const real_t absolute_error = fabs(exact - approximation);
         const real_t relative_error = absolute_error /
            max(max(fabs(exact), fabs(approximation)), real_t(1.0e-30));
         maximum_relative_error = max(maximum_relative_error, relative_error);
         squared_error += absolute_error * absolute_error;
         squared_reference += approximation * approximation;
         output << frequencies[frequency_index] << ',' << receiver << ','
                << parts[part] << ',' << parameters[parameter] << ','
                << exact << ',' << approximation << ',' << absolute_error
                << ',' << relative_error << '\n';
      }
   }
   cout << "Jacobian check: relative Frobenius error="
        << sqrt(squared_error /
                max(squared_reference, numeric_limits<real_t>::min()))
        << ", maximum entrywise relative error=" << maximum_relative_error
        << " (details: " << file_name << ")\n";
}

void RunInversion(CSEMProblem &problem, const CSEMOptions &options,
                  const vector<real_t> &frequencies,
                  const SolverSettings &forward_settings,
                  const SolverSettings &inversion_settings)
{
   Vector observed, standard_deviation;
   if (string(options.observed_file).empty())
   {
      BuildSyntheticData(problem, options, frequencies, forward_settings,
                         observed,
                         standard_deviation);
      WriteData(string(options.output_prefix) + "_synthetic.csv", problem,
                frequencies, observed, standard_deviation);
   }
   else
   {
      ReadObservedData(options.observed_file,
                       (problem.DataSize() / 2) * frequencies.size(), observed,
                       standard_deviation);
   }

   Vector model(2), reference_model(2), lower(2), upper(2);
   model[0] = log(options.initial_background_conductivity);
   model[1] = log(options.initial_target_conductivity);
   reference_model = model;
   lower = log(options.minimum_conductivity);
   upper = log(options.maximum_conductivity);

   ofstream history;
   if (problem.Rank() == 0)
   {
      history.open(string(options.output_prefix) + "_inversion.csv");
      MFEM_VERIFY(history, "Unable to open inversion history output.");
      history << "iteration,objective,sigma_background,sigma_target,"
                 "gradient_norm,step_norm,line_search\n";
      history << setprecision(16);
   }

   for (int iteration = 0; iteration < options.inversion_iterations;
        iteration++)
   {
      Vector prediction;
      DenseMatrix jacobian;
      SolveSurvey(problem, frequencies, model, prediction, &jacobian,
                  inversion_settings);
      const real_t objective =
         Objective(prediction, observed, standard_deviation, model,
                   reference_model, options.regularization);

      Vector weighted_residual(observed.Size());
      for (int datum = 0; datum < observed.Size(); datum++)
      {
         const real_t variance =
            standard_deviation[datum] * standard_deviation[datum];
         weighted_residual[datum] =
            (prediction[datum] - observed[datum]) / variance;
      }
      Vector gradient(2);
      jacobian.MultTranspose(weighted_residual, gradient);
      for (int parameter = 0; parameter < 2; parameter++)
      {
         gradient[parameter] += options.regularization *
                                (model[parameter] -
                                 reference_model[parameter]);
      }
      const real_t gradient_norm = gradient.Norml2();

      GaussNewtonNormalOperator normal_operator(
         jacobian, standard_deviation,
         options.regularization + options.gauss_newton_damping);
      Vector step(2);
      step = 0.0;
      Vector negative_gradient(gradient);
      negative_gradient *= -1.0;
      CGSolver least_squares_cg;
      least_squares_cg.SetOperator(normal_operator);
      least_squares_cg.SetRelTol(1.0e-12);
      least_squares_cg.SetAbsTol(0.0);
      least_squares_cg.SetMaxIter(20);
      least_squares_cg.SetPrintLevel(0);
      least_squares_cg.Mult(negative_gradient, step);
      MFEM_VERIFY(least_squares_cg.GetConverged(),
                  "The Gauss-Newton least-squares solve did not converge.");
      const real_t step_norm = step.Norml2();

      real_t line_search = 1.0;
      bool accepted = false;
      Vector trial_model(2), trial_prediction;
      for (int trial = 0; trial < options.line_search_steps; trial++)
      {
         trial_model = model;
         trial_model.Add(line_search, step);
         for (int parameter = 0; parameter < 2; parameter++)
         {
            trial_model[parameter] =
               max(lower[parameter], min(upper[parameter],
                                          trial_model[parameter]));
         }
         SolveSurvey(problem, frequencies, trial_model, trial_prediction,
                     nullptr, inversion_settings);
         const real_t trial_objective =
            Objective(trial_prediction, observed, standard_deviation,
                      trial_model, reference_model, options.regularization);
         if (trial_objective < objective)
         {
            accepted = true;
            break;
         }
         line_search *= 0.5;
      }

      if (problem.Rank() == 0)
      {
         cout << "GN " << iteration << ": Phi=" << objective
              << ", sigma=[" << exp(model[0]) << ", " << exp(model[1])
              << "], |g|=" << gradient_norm << ", |dm|=" << step_norm
              << ", alpha=" << (accepted ? line_search : 0.0) << '\n';
         history << iteration << ',' << objective << ',' << exp(model[0])
                 << ',' << exp(model[1]) << ',' << gradient_norm << ','
                 << step_norm << ',' << (accepted ? line_search : 0.0)
                 << '\n';
      }

      if (gradient_norm <= options.inversion_gradient_tolerance ||
          step_norm <= options.inversion_step_tolerance || !accepted)
      {
         break;
      }
      model = trial_model;
   }

   Vector final_prediction;
   SolveSurvey(problem, frequencies, model, final_prediction, nullptr,
               forward_settings,
               string(options.output_prefix) + "_field");
   WriteData(string(options.output_prefix) + "_recovered.csv", problem,
             frequencies, final_prediction, standard_deviation);
   if (problem.Rank() == 0)
   {
      cout << "Recovered conductivities [S/m]: background=" << exp(model[0])
           << ", target=" << exp(model[1]) << '\n';
   }
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();

   CSEMOptions options;
   OptionsParser args(argc, argv);
   args.AddOption(&options.mesh_file, "-m", "--mesh",
                  "Physical-coordinate hexahedral mesh file.");
   args.AddOption(&options.solver, "-s", "--solver",
                  "Linear solver: presb, presb-mumps, or mumps.");
   args.AddOption(&options.order, "-o", "--order",
                  "Nedelec spectral element order.");
   args.AddOption(&options.refinements, "-r", "--refinements",
                  "Uniform mesh refinements.");
   args.AddOption(&options.nx, "-nx", "--elements-x",
                  "Generated mesh elements in x.");
   args.AddOption(&options.ny, "-ny", "--elements-y",
                  "Generated mesh elements in y.");
   args.AddOption(&options.nz, "-nz", "--elements-z",
                  "Generated mesh elements in z.");
   args.AddOption(&options.frequency, "-f", "--frequency",
                  "Source frequency in Hz.");
   args.AddOption(&options.frequencies, "-fs", "--frequencies",
                  "Comma-separated source frequencies in Hz; overrides -f.");
   args.AddOption(&options.background_conductivity, "-sb", "--sigma-background",
                  "True/background conductivity in S/m.");
   args.AddOption(&options.target_conductivity, "-st", "--sigma-target",
                  "True/target conductivity in S/m.");
   args.AddOption(&options.initial_background_conductivity, "-isb",
                  "--initial-sigma-background",
                  "Initial inversion background conductivity in S/m.");
   args.AddOption(&options.initial_target_conductivity, "-ist",
                  "--initial-sigma-target",
                  "Initial inversion target conductivity in S/m.");
   args.AddOption(&options.receivers, "-nr", "--receivers",
                  "Number of Gaussian receiver functionals.");
   args.AddOption(&options.receiver_file, "-rf", "--receiver-file",
                  "Receiver CSV with x,y,z[,azimuth,dip] rows.");
   args.AddOption(&options.receiver_field, "-rc", "--receiver-field",
                  "Receiver field quantity: E or H.");
   args.AddOption(&options.source_x, "-sx", "--source-x",
                  "Wire center x coordinate in metres.");
   args.AddOption(&options.source_file, "-sf", "--source-file",
                  "Transmitter polyline CSV with x,y,z vertex rows.");
   args.AddOption(&options.source_y, "-sy", "--source-y",
                  "Wire center y coordinate in metres.");
   args.AddOption(&options.source_z, "-sz", "--source-z",
                  "Wire center z coordinate in metres.");
   args.AddOption(&options.source_length, "-sl", "--source-length",
                  "Wire length in metres.");
   args.AddOption(&options.source_radius, "-sr", "--source-radius",
                  "Gaussian wire regularization radius in metres.");
   args.AddOption(&options.source_current, "-si", "--source-current",
                  "Source current in amperes.");
   args.AddOption(&options.source_azimuth, "-saz", "--source-azimuth",
                  "Wire azimuth in degrees, measured from +x toward +y.");
   args.AddOption(&options.source_dip, "-sdip", "--source-dip",
                  "Wire dip in degrees; positive points toward +z.");
   args.AddOption(&options.receiver_azimuth, "-raz",
                  "--receiver-azimuth",
                  "Electric receiver-axis azimuth in degrees.");
   args.AddOption(&options.receiver_dip, "-rdip", "--receiver-dip",
                  "Electric receiver-axis dip in degrees.");
   args.AddOption(&options.receiver_x, "-rx", "--receiver-x",
                  "Receiver-line center x coordinate in metres.");
   args.AddOption(&options.receiver_y, "-ry", "--receiver-y",
                  "Receiver-line center y coordinate in metres.");
   args.AddOption(&options.receiver_z, "-rz", "--receiver-z",
                  "Receiver-line center z coordinate in metres.");
   args.AddOption(&options.receiver_span, "-rsp", "--receiver-span",
                  "Receiver-line span in metres.");
   args.AddOption(&options.receiver_radius, "-rr", "--receiver-radius",
                  "Gaussian receiver averaging radius in metres.");
   args.AddOption(&options.receiver_line_azimuth, "-rlaz",
                  "--receiver-line-azimuth",
                  "Receiver-line azimuth in degrees.");
   args.AddOption(&options.observed_file, "-obs", "--observed",
                  "Observed CSV; empty selects a synthetic inversion.");
   args.AddOption(&options.output_prefix, "-p", "--output-prefix",
                  "Prefix for CSV and ParaView outputs.");
   args.AddOption(&options.forward_outer_tolerance, "-fot",
                  "--forward-outer-tolerance", "Forward outer tolerance.");
   args.AddOption(&options.forward_inner_tolerance, "-fit",
                  "--forward-inner-tolerance", "Forward inner tolerance.");
   args.AddOption(&options.inversion_outer_tolerance, "-iot",
                  "--inversion-outer-tolerance", "Inversion outer tolerance.");
   args.AddOption(&options.inversion_inner_tolerance, "-iit",
                  "--inversion-inner-tolerance", "Inversion inner tolerance.");
   args.AddOption(&options.inversion_iterations, "-ni",
                  "--inversion-iterations", "Gauss-Newton iterations.");
   args.AddOption(&options.regularization, "-reg", "--regularization",
                  "Log-conductivity Tikhonov weight.");
   args.AddOption(&options.synthetic_noise, "-noise", "--synthetic-noise",
                  "Synthetic Gaussian noise as a relative fraction.");
   args.AddOption(&options.jacobian_fd_step, "-jfd",
                  "--jacobian-fd-step",
                  "Centered finite-difference step in log conductivity.");
   args.AddOption(&options.verify_jacobian, "-vj", "--verify-jacobian",
                  "-no-vj", "--no-verify-jacobian",
                  "Compare analytic tangent sensitivities with finite "
                  "differences.");
   args.AddOption(&options.inversion, "-inv", "--inversion", "-no-inv",
                  "--no-inversion", "Run inversion instead of forward only.");
   args.AddOption(&options.visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization", "Save a ParaView field collection.");
   args.AddOption(&options.device, "-d", "--device", "Device configuration.");
   args.Parse();
   if (!args.Good())
   {
      if (Mpi::Root()) { args.PrintUsage(cout); }
      return 1;
   }
   if (Mpi::Root()) { args.PrintOptions(cout); }

   MFEM_VERIFY(options.order >= 1, "Element order must be positive.");
   MFEM_VERIFY(options.frequency > 0.0, "Frequency must be positive.");
   MFEM_VERIFY(options.jacobian_fd_step > 0.0,
               "The Jacobian finite-difference step must be positive.");
   MFEM_VERIFY(options.receivers > 0 ||
                  !string(options.receiver_file).empty(),
               "At least one receiver or a receiver file is required.");
   MFEM_VERIFY(string(options.receiver_field) == "E" ||
                  string(options.receiver_field) == "e" ||
                  string(options.receiver_field) == "H" ||
                  string(options.receiver_field) == "h",
               "Receiver field must be E or H.");
   MFEM_VERIFY((options.source_length > 0.0 ||
                   !string(options.source_file).empty()) &&
                  options.source_radius > 0.0 &&
                  options.receiver_radius > 0.0,
               "Source length (unless a source file is used) and the "
               "source/receiver radii must be positive.");
   MFEM_VERIFY(options.background_conductivity > 0.0 &&
                  options.target_conductivity > 0.0 &&
                  options.initial_background_conductivity > 0.0 &&
                  options.initial_target_conductivity > 0.0,
               "Conductivities must be positive.");

   Device device(options.device);
   if (Mpi::Root()) { device.Print(); }
   const vector<real_t> frequencies = ParseFrequencies(options);
   if (Mpi::Root())
   {
      for (const real_t frequency : frequencies)
      {
         const real_t skin_depth = PlaneWaveSkinDepth(options, frequency);
         cout << "Background plane-wave skin depth at " << frequency
              << " Hz: " << skin_depth << " m; paper guideline 5 delta = "
              << 5.0 * skin_depth << " m from the area of interest\n";
      }
   }
   CSEMProblem problem(MPI_COMM_WORLD, options);
   SolverSettings forward_settings{options.forward_outer_tolerance,
                                   options.forward_inner_tolerance, true};
   SolverSettings inversion_settings{options.inversion_outer_tolerance,
                                     options.inversion_inner_tolerance, false};

   if (options.inversion)
   {
      RunInversion(problem, options, frequencies, forward_settings,
                   inversion_settings);
   }
   else
   {
      Vector model(2), responses;
      model[0] = log(options.background_conductivity);
      model[1] = log(options.target_conductivity);
      SolveSurvey(problem, frequencies, model, responses, nullptr,
                  forward_settings,
                  options.visualization
                     ? string(options.output_prefix) + "_field"
                     : string());
      Vector standard_deviation(responses.Size());
      real_t maximum_amplitude = 0.0;
      for (int i = 0; i < responses.Size(); i++)
      {
         maximum_amplitude = max(maximum_amplitude, fabs(responses[i]));
      }
      for (int i = 0; i < responses.Size(); i++)
      {
         standard_deviation[i] = options.data_relative_error *
            max(fabs(responses[i]),
                options.data_error_floor * maximum_amplitude);
      }
      WriteData(string(options.output_prefix) + "_response.csv", problem,
                frequencies, responses, standard_deviation);
      if (options.verify_jacobian)
      {
         VerifyJacobian(problem, options, frequencies, model,
                        forward_settings);
      }
   }
   return 0;
}
