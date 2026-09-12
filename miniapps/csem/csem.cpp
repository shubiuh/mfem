// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

//                 MFEM Miniapp: CSEM (csem)
//
// This miniapp implements the total-electric-field controlled-source
// electromagnetic formulation, the Weiss PRESB inner-outer solver, and a
// anisotropic region/cell Gauss-Newton inversion.  It requires MPI and hypre.

#include "mfem.hpp"
#include "csem_solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <set>
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
   const char *material_file = "";
   const char *boundary_condition = "pec";
   const char *absorbing_boundary_attributes = "";
   const char *solver = "presb";
   const char *observed_file = "";
   const char *receiver_file = "";
   const char *receiver_field = "E";
   const char *receiver_model = "gaussian";
   const char *source_file = "";
   const char *source_catalog = "";
   const char *source_model = "gaussian";
   const char *model_parameterization = "two";
   const char *model_file = "";
   const char *initial_model_file = "";
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
   int source_quadrature_order = 8;
   int survey_groups = 1;
   int regularization_neighbors = 6;

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
   real_t background_vertical_resistivity = 1000.0;
   real_t target_vertical_resistivity = 20.0;
   real_t background_dip = 0.0;
   real_t target_dip = 0.0;
   real_t minimum_resistivity = 0.1;
   real_t maximum_resistivity = 1.0e6;
   real_t minimum_dip = -89.0;
   real_t maximum_dip = 89.0;

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
   real_t regularization_rh_scale = 1.0;
   real_t regularization_rv_scale = 1.0;
   real_t regularization_dip_scale = 1.0;
   real_t gauss_newton_damping = 1.0e-6;
   real_t gauss_newton_cg_tolerance = 1.0e-6;
   real_t inversion_gradient_tolerance = 1.0e-5;
   real_t inversion_step_tolerance = 1.0e-4;
   real_t jacobian_fd_step = 1.0e-4;
   real_t absorbing_conductivity = 0.0;
   real_t absorbing_relative_permittivity = 1.0;
   real_t absorbing_relative_permeability = 1.0;

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

vector<int> ParseIntegerList(const string &text)
{
   string normalized(text);
   replace(normalized.begin(), normalized.end(), ',', ' ');
   replace(normalized.begin(), normalized.end(), ';', ' ');
   istringstream input(normalized);
   vector<int> values;
   int value;
   while (input >> value) { values.push_back(value); }
   return values;
}

struct TaggedMaterial
{
   real_t rh = 1.0;
   real_t rv = 1.0;
   real_t dip = 0.0;
   real_t relative_permittivity = 1.0;
   bool invertible = false;
};

class TaggedMaterialCatalog
{
private:
   bool enabled = false;
   vector<TaggedMaterial> materials;
   vector<int> defined;

public:
   TaggedMaterialCatalog(const CSEMOptions &options, const ParMesh &mesh)
   {
      const string file_name(options.material_file);
      if (file_name.empty()) { return; }
      enabled = true;
      const int maximum_attribute = mesh.attributes.Size()
                                       ? mesh.attributes.Max() : 0;
      MFEM_VERIFY(maximum_attribute > 0,
                  "A tagged material table requires volume attributes.");
      materials.resize(maximum_attribute + 1);
      defined.assign(maximum_attribute + 1, 0);

      ifstream input(file_name);
      MFEM_VERIFY(input, "Unable to open tagged material file: " << file_name);
      string line;
      while (getline(input, line))
      {
         const size_t comment = line.find('#');
         if (comment != string::npos) { line.erase(comment); }
         replace(line.begin(), line.end(), ',', ' ');
         istringstream row(line);
         int attribute, invertible;
         real_t rh, rv, dip_degrees, relative_permittivity;
         if (!(row >> attribute >> rh >> rv >> dip_degrees >>
               relative_permittivity >> invertible))
         {
            continue;
         }
         MFEM_VERIFY(attribute > 0 && attribute <= maximum_attribute,
                     "Tagged material attribute is outside the mesh range.");
         MFEM_VERIFY(!defined[attribute],
                     "Duplicate tagged material attribute " << attribute);
         MFEM_VERIFY(rh > 0.0 && rv > 0.0 && relative_permittivity > 0.0,
                     "Tagged Rh, Rv, and relative permittivity must be "
                     "positive.");
         MFEM_VERIFY(invertible == 0 || invertible == 1,
                     "Tagged material invertible flag must be 0 or 1.");
         materials[attribute].rh = rh;
         materials[attribute].rv = rv;
         materials[attribute].dip = dip_degrees * M_PI / 180.0;
         materials[attribute].relative_permittivity = relative_permittivity;
         materials[attribute].invertible = invertible != 0;
         defined[attribute] = 1;
      }

      for (int i = 0; i < mesh.attributes.Size(); i++)
      {
         const int attribute = mesh.attributes[i];
         MFEM_VERIFY(attribute > 0 && attribute <= maximum_attribute &&
                        defined[attribute],
                     "Mesh volume attribute " << attribute
                     << " is missing from the tagged material table.");
      }
   }

   bool Enabled() const { return enabled; }

   const TaggedMaterial &Get(int attribute) const
   {
      MFEM_VERIFY(enabled && attribute > 0 &&
                     attribute < static_cast<int>(materials.size()) &&
                     defined[attribute],
                  "Undefined tagged material attribute " << attribute);
      return materials[attribute];
   }

   bool IsInvertible(int attribute) const
   {
      return Get(attribute).invertible;
   }
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

struct Transmitter
{
   vector<Vector> vertices;
   real_t current = 1.0;
   real_t radius = 1.0;
};

void AddStraightTransmitter(vector<Transmitter> &transmitters,
                            const Vector &center, real_t length,
                            real_t azimuth, real_t dip, real_t current,
                            real_t radius)
{
   MFEM_VERIFY(length > 0.0 && radius > 0.0,
               "Transmitter length and radius must be positive.");
   const Vector direction = UnitVectorFromAngles(azimuth, dip);
   Transmitter transmitter;
   transmitter.current = current;
   transmitter.radius = radius;
   Vector start(center), end(center);
   start.Add(-0.5 * length, direction);
   end.Add(0.5 * length, direction);
   transmitter.vertices.push_back(start);
   transmitter.vertices.push_back(end);
   transmitters.push_back(transmitter);
}

vector<Transmitter> LoadTransmitters(const CSEMOptions &options)
{
   vector<Transmitter> transmitters;
   if (!string(options.source_catalog).empty())
   {
      ifstream input(options.source_catalog);
      MFEM_VERIFY(input, "Unable to open transmitter catalog: "
                  << options.source_catalog);
      string line;
      while (getline(input, line))
      {
         replace(line.begin(), line.end(), ',', ' ');
         istringstream row(line);
         vector<real_t> values;
         real_t value;
         while (row >> value) { values.push_back(value); }
         if (values.empty()) { continue; }
         MFEM_VERIFY(values.size() == 7 || values.size() == 8,
                     "Each transmitter catalog row must contain "
                     "x,y,z,length,azimuth,dip,current[,radius].");
         Vector center(3);
         center[0] = values[0];
         center[1] = values[1];
         center[2] = values[2];
         AddStraightTransmitter(transmitters, center, values[3], values[4],
                                values[5], values[6],
                                values.size() == 8 ? values[7]
                                                   : options.source_radius);
      }
      MFEM_VERIFY(!transmitters.empty(),
                  "The transmitter catalog contains no stations.");
      return transmitters;
   }

   if (!string(options.source_file).empty())
   {
      ifstream input(options.source_file);
      MFEM_VERIFY(input, "Unable to open transmitter geometry file: "
                  << options.source_file);
      Transmitter transmitter;
      transmitter.current = options.source_current;
      transmitter.radius = options.source_radius;
      string line;
      while (getline(input, line))
      {
         replace(line.begin(), line.end(), ',', ' ');
         istringstream row(line);
         Vector vertex(3);
         if (row >> vertex[0] >> vertex[1] >> vertex[2])
         {
            transmitter.vertices.push_back(vertex);
         }
      }
      MFEM_VERIFY(transmitter.vertices.size() >= 2,
                  "A transmitter polyline needs at least two vertices.");
      transmitters.push_back(transmitter);
      return transmitters;
   }

   Vector center(3);
   center[0] = options.source_x;
   center[1] = options.source_y;
   center[2] = options.source_z;
   AddStraightTransmitter(transmitters, center, options.source_length,
                          options.source_azimuth, options.source_dip,
                          options.source_current, options.source_radius);
   return transmitters;
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

class CSEMConductivity : public MatrixCoefficient
{
private:
   const CSEMOptions &options;
   ParMesh &mesh;
   const TaggedMaterialCatalog *materials;
   real_t length_scale;
   int unit_count = 0;
   vector<int> unit_ids;
   Vector model;
   Vector true_model;
   Vector initial_model;

   int Unit(const Vector &scaled_x, int local_element) const
   {
      const string mode(options.model_parameterization);
      const int attribute = mesh.GetAttribute(local_element);
      if (materials->Enabled() && !materials->IsInvertible(attribute))
      {
         return -1;
      }
      if (mode == "cell")
      {
         const int global_element =
            static_cast<int>(mesh.GetGlobalElementNum(local_element));
         const auto found = lower_bound(unit_ids.begin(), unit_ids.end(),
                                        global_element);
         MFEM_VERIFY(found != unit_ids.end() && *found == global_element,
                     "Active cell is missing from the cell model.");
         return static_cast<int>(found - unit_ids.begin());
      }
      if (mode == "region")
      {
         const auto found = lower_bound(unit_ids.begin(), unit_ids.end(),
                                        attribute);
         MFEM_VERIFY(found != unit_ids.end() && *found == attribute,
                     "Element attribute is missing from the region model.");
         return static_cast<int>(found - unit_ids.begin());
      }
      return InTarget(scaled_x) ? 1 : 0;
   }

   void ReadProperties(const string &file_name, Vector &values) const
   {
      if (file_name.empty()) { return; }
      ifstream input(file_name);
      MFEM_VERIFY(input, "Unable to open resistivity model file: "
                  << file_name);
      string line;
      while (getline(input, line))
      {
         replace(line.begin(), line.end(), ',', ' ');
         istringstream row(line);
         long long id;
         real_t rh, rv, dip;
         if (!(row >> id >> rh >> rv >> dip)) { continue; }
         long long index = id;
          if (string(options.model_parameterization) != "two")
          {
            const auto found = lower_bound(unit_ids.begin(), unit_ids.end(),
                                           static_cast<int>(id));
            MFEM_VERIFY(found != unit_ids.end() && *found == id,
                        "Model property attribute is not present in the mesh.");
            index = found - unit_ids.begin();
         }
         MFEM_VERIFY(index >= 0 && index < unit_count,
                     "Model property id is outside the mesh parameter range.");
         MFEM_VERIFY(rh > 0.0 && rv > 0.0,
                     "Rh and Rv must be positive.");
         const int offset = 3 * static_cast<int>(index);
         values[offset] = log(rh);
         values[offset + 1] = log(rv);
         values[offset + 2] = dip * M_PI / 180.0;
      }
   }

   void Tensor(const Vector &scaled_x, int local_element,
               DenseMatrix &sigma) const
   {
      sigma.SetSize(3);
      sigma = 0.0;
      if (!materials->Enabled() && length_scale * scaled_x[2] < 0.0)
      {
         for (int d = 0; d < 3; d++) { sigma(d, d) = options.air_conductivity; }
         return;
      }
      const int unit = Unit(scaled_x, local_element);
      real_t sigma_h, sigma_v, dip;
      if (unit < 0)
      {
         const TaggedMaterial &material =
            materials->Get(mesh.GetAttribute(local_element));
         sigma_h = 1.0 / material.rh;
         sigma_v = 1.0 / material.rv;
         dip = material.dip;
      }
      else
      {
         sigma_h = exp(-model[3 * unit]);
         sigma_v = exp(-model[3 * unit + 1]);
         dip = model[3 * unit + 2];
      }
      Vector axis(3);
      axis[0] = sin(dip);
      axis[1] = 0.0;
      axis[2] = cos(dip);
      for (int i = 0; i < 3; i++)
      {
         sigma(i, i) = sigma_h;
         for (int j = 0; j < 3; j++)
         {
            sigma(i, j) += (sigma_v - sigma_h) * axis[i] * axis[j];
         }
      }
   }

public:
   CSEMConductivity(const CSEMOptions &options_, ParMesh &mesh_,
                    const TaggedMaterialCatalog &materials_,
                    real_t length_scale_)
      : MatrixCoefficient(3), options(options_), mesh(mesh_),
        materials(&materials_), length_scale(length_scale_)
   {
      const string mode(options.model_parameterization);
      vector<int> global_cell_attributes;
      if (mode == "cell")
      {
         const int global_elements = static_cast<int>(mesh.GetGlobalNE());
         vector<int> marker(global_elements, 0);
         global_cell_attributes.assign(global_elements, 0);
         for (int element = 0; element < mesh.GetNE(); element++)
         {
            const int global_element =
               static_cast<int>(mesh.GetGlobalElementNum(element));
            const int attribute = mesh.GetAttribute(element);
            global_cell_attributes[global_element] = attribute;
            if (!materials->Enabled() || materials->IsInvertible(attribute))
            {
               marker[global_element] = 1;
            }
         }
         MPI_Allreduce(MPI_IN_PLACE, marker.data(), marker.size(), MPI_INT,
                       MPI_MAX, mesh.GetComm());
         MPI_Allreduce(MPI_IN_PLACE, global_cell_attributes.data(),
                       global_cell_attributes.size(), MPI_INT, MPI_MAX,
                       mesh.GetComm());
         for (int element = 0; element < global_elements; element++)
         {
            if (marker[element]) { unit_ids.push_back(element); }
         }
         unit_count = static_cast<int>(unit_ids.size());
      }
      else if (mode == "region")
      {
         int local_max = mesh.attributes.Size() ? mesh.attributes.Max() : 0;
         int global_max = 0;
         MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX,
                       mesh.GetComm());
         vector<int> marker(global_max + 1, 0);
         for (int element = 0; element < mesh.GetNE(); element++)
         {
            const int attribute = mesh.GetAttribute(element);
            if (!materials->Enabled() || materials->IsInvertible(attribute))
            {
               marker[attribute] = 1;
            }
         }
         MPI_Allreduce(MPI_IN_PLACE, marker.data(), marker.size(), MPI_INT,
                       MPI_MAX, mesh.GetComm());
         for (int attribute = 1; attribute <= global_max; attribute++)
         {
            if (marker[attribute]) { unit_ids.push_back(attribute); }
         }
         unit_count = static_cast<int>(unit_ids.size());
      }
      else
      {
         MFEM_VERIFY(mode == "two", "Model parameterization must be two, "
                     "region, or cell.");
         MFEM_VERIFY(!materials->Enabled(),
                     "Tagged material assignment requires region or cell "
                     "parameterization.");
         unit_count = 2;
         unit_ids = {0, 1};
      }
      MFEM_VERIFY(unit_count > 0, "The conductivity model has no units.");
      true_model.SetSize(3 * unit_count);
      initial_model.SetSize(3 * unit_count);
      for (int unit = 0; unit < unit_count; unit++)
      {
         if (materials->Enabled())
         {
            const int attribute = mode == "region"
               ? unit_ids[unit] : global_cell_attributes[unit_ids[unit]];
            const TaggedMaterial &material = materials->Get(attribute);
            true_model[3 * unit] = log(material.rh);
            true_model[3 * unit + 1] = log(material.rv);
            true_model[3 * unit + 2] = material.dip;
            initial_model[3 * unit] = true_model[3 * unit];
            initial_model[3 * unit + 1] = true_model[3 * unit + 1];
            initial_model[3 * unit + 2] = true_model[3 * unit + 2];
            continue;
         }
         const bool target_unit = mode == "two" && unit == 1;
         const real_t sigma = target_unit ? options.target_conductivity
                                          : options.background_conductivity;
         const real_t initial_sigma = target_unit
            ? options.initial_target_conductivity
            : options.initial_background_conductivity;
         true_model[3 * unit] = log(1.0 / sigma);
         true_model[3 * unit + 1] = log(target_unit
            ? options.target_vertical_resistivity
            : options.background_vertical_resistivity);
         true_model[3 * unit + 2] = (target_unit ? options.target_dip
                                                 : options.background_dip) *
                                    M_PI / 180.0;
         initial_model[3 * unit] = log(1.0 / initial_sigma);
         initial_model[3 * unit + 1] = true_model[3 * unit + 1];
         initial_model[3 * unit + 2] = true_model[3 * unit + 2];
      }
      ReadProperties(options.model_file, true_model);
      ReadProperties(options.initial_model_file, initial_model);
      model = true_model;
   }

   void SetModel(const Vector &model_)
   {
      MFEM_VERIFY(model_.Size() == ParameterCount(),
                  "The CSEM model vector has the wrong size.");
      model = model_;
   }

   int ParameterCount() const { return 3 * unit_count; }
   int UnitCount() const { return unit_count; }
   int UnitId(int unit) const { return unit_ids[unit]; }
   const Vector &TrueModel() const { return true_model; }
   const Vector &InitialModel() const { return initial_model; }
   int UnitAt(const Vector &scaled_x, int local_element) const
   {
      return Unit(scaled_x, local_element);
   }

   void Properties(ElementTransformation &T, const IntegrationPoint &ip,
                   real_t &rh, real_t &rv, real_t &dip_degrees) const
   {
      Vector x;
      T.Transform(ip, x);
       if (!materials->Enabled() && length_scale * x[2] < 0.0)
       {
         rh = rv = 1.0 / options.air_conductivity;
         dip_degrees = 0.0;
         return;
       }
       const int unit = Unit(x, T.ElementNo);
       if (unit < 0)
       {
          const TaggedMaterial &material =
             materials->Get(mesh.GetAttribute(T.ElementNo));
          rh = material.rh;
          rv = material.rv;
          dip_degrees = material.dip * 180.0 / M_PI;
       }
       else
       {
          rh = exp(model[3 * unit]);
          rv = exp(model[3 * unit + 1]);
          dip_degrees = model[3 * unit + 2] * 180.0 / M_PI;
       }
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

   void Eval(DenseMatrix &sigma, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      Vector x;
      T.Transform(ip, x);
      Tensor(x, T.ElementNo, sigma);
   }

   void ParameterDerivative(const Vector &scaled_x, int local_element,
                             int parameter, DenseMatrix &derivative) const
   {
      derivative.SetSize(3);
      derivative = 0.0;
      const int point_unit = Unit(scaled_x, local_element);
      if ((!materials->Enabled() && length_scale * scaled_x[2] < 0.0) ||
          point_unit < 0 || parameter / 3 != point_unit) { return; }
      const int unit = parameter / 3;
      const int component = parameter % 3;
      const real_t sigma_h = exp(-model[3 * unit]);
      const real_t sigma_v = exp(-model[3 * unit + 1]);
      const real_t dip = model[3 * unit + 2];
      Vector axis(3), derivative_axis(3);
      axis[0] = sin(dip); axis[1] = 0.0; axis[2] = cos(dip);
      derivative_axis[0] = cos(dip);
      derivative_axis[1] = 0.0;
      derivative_axis[2] = -sin(dip);
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            if (component == 0)
            {
               derivative(i, j) = -sigma_h *
                  ((i == j ? 1.0 : 0.0) - axis[i] * axis[j]);
            }
            else if (component == 1)
            {
               derivative(i, j) = -sigma_v * axis[i] * axis[j];
            }
            else
            {
               derivative(i, j) = (sigma_v - sigma_h) *
                  (derivative_axis[i] * axis[j] +
                   axis[i] * derivative_axis[j]);
            }
         }
      }
   }

   void DirectionalDerivative(const Vector &scaled_x, int local_element,
                              const Vector &direction,
                              DenseMatrix &derivative) const
   {
      MFEM_VERIFY(direction.Size() == ParameterCount(),
                  "The conductivity direction has the wrong size.");
      derivative.SetSize(3);
      derivative = 0.0;
      if (!materials->Enabled() && length_scale * scaled_x[2] < 0.0)
      {
         return;
      }

      const int unit = Unit(scaled_x, local_element);
      if (unit < 0) { return; }
      DenseMatrix component(3);
      for (int property = 0; property < 3; property++)
      {
         const real_t weight = direction[3 * unit + property];
         if (weight == 0.0) { continue; }
         ParameterDerivative(scaled_x, local_element,
                             3 * unit + property, component);
         for (int i = 0; i < 3; i++)
         {
            for (int j = 0; j < 3; j++)
            {
               derivative(i, j) += weight * component(i, j);
            }
         }
      }
   }
};

class ConductivityPropertyCoefficient : public Coefficient
{
private:
   CSEMConductivity &conductivity;
   int property;

public:
   ConductivityPropertyCoefficient(CSEMConductivity &conductivity_,
                                   int property_)
      : conductivity(conductivity_), property(property_) { }

   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override
   {
      real_t rh, rv, dip;
      conductivity.Properties(T, ip, rh, rv, dip);
      if (property == 0) { return rh; }
      if (property == 1) { return rv; }
      if (property == 2) { return dip; }
      if (property == 3) { return 1.0 / rh; }
      return 1.0 / rv;
   }
};

class MaterialTagCoefficient : public Coefficient
{
private:
   ParMesh &mesh;
   const TaggedMaterialCatalog &materials;
   real_t length_scale;
   bool active_mask;

public:
   MaterialTagCoefficient(ParMesh &mesh_,
                          const TaggedMaterialCatalog &materials_,
                          real_t length_scale_, bool active_mask_)
      : mesh(mesh_), materials(materials_), length_scale(length_scale_),
        active_mask(active_mask_)
   {
   }

   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override
   {
      const int attribute = mesh.GetAttribute(T.ElementNo);
      if (!active_mask) { return static_cast<real_t>(attribute); }
      if (materials.Enabled())
      {
         return materials.IsInvertible(attribute) ? 1.0 : 0.0;
      }
      Vector x;
      T.Transform(ip, x);
      return length_scale * x[2] >= 0.0 ? 1.0 : 0.0;
   }
};

class RelativePermittivity : public Coefficient
{
private:
   const CSEMOptions &options;
   ParMesh &mesh;
   const TaggedMaterialCatalog *materials;
   real_t length_scale;

public:
   RelativePermittivity(const CSEMOptions &options_, ParMesh &mesh_,
                        const TaggedMaterialCatalog &materials_,
                        real_t length_scale_)
      : options(options_), mesh(mesh_), materials(&materials_),
        length_scale(length_scale_)
   {
   }

   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override
   {
      if (materials->Enabled())
      {
         return materials->Get(
            mesh.GetAttribute(T.ElementNo)).relative_permittivity;
      }
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

class ScaledMatrixCoefficient : public MatrixCoefficient
{
private:
   MatrixCoefficient &coefficient;
   real_t scale;

public:
   ScaledMatrixCoefficient(MatrixCoefficient &coefficient_, real_t scale_)
      : MatrixCoefficient(coefficient_.GetVDim()), coefficient(coefficient_),
        scale(scale_) { }

   void Eval(DenseMatrix &matrix, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      coefficient.Eval(matrix, T, ip);
      matrix *= scale;
   }
};

class ShiftedMassCoefficient : public MatrixCoefficient
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
      : MatrixCoefficient(3), conductivity(conductivity_),
        permittivity(permittivity_),
        conductivity_scale(conductivity_scale_),
        permittivity_scale(permittivity_scale_)
   {
   }

   void Eval(DenseMatrix &matrix, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      conductivity.Eval(matrix, T, ip);
      matrix *= conductivity_scale;
      const real_t dielectric =
         permittivity_scale * permittivity.Eval(T, ip);
      for (int d = 0; d < 3; d++) { matrix(d, d) -= dielectric; }
   }
};

class DirectionalDerivativeCoefficient : public MatrixCoefficient
{
private:
   CSEMConductivity &conductivity;
   const Vector &direction;
   real_t scale;

public:
   DirectionalDerivativeCoefficient(CSEMConductivity &conductivity_,
                                    const Vector &direction_, real_t scale_)
      : MatrixCoefficient(3), conductivity(conductivity_),
        direction(direction_), scale(scale_)
   {
   }

   void Eval(DenseMatrix &matrix, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      Vector x;
      T.Transform(ip, x);
      conductivity.DirectionalDerivative(x, T.ElementNo, direction, matrix);
      matrix *= scale;
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

   real_t length_scale;
   real_t right_hand_side_scale;
   real_t current;
   real_t radius;
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
   GaussianWireCurrent(const Transmitter &transmitter, real_t length_scale_,
                       real_t omega)
      : VectorCoefficient(3), length_scale(length_scale_),
        right_hand_side_scale(-omega * mu0 * length_scale_ * length_scale_),
        current(transmitter.current), radius(transmitter.radius)
   {
      for (int vertex = 0;
           vertex + 1 < static_cast<int>(transmitter.vertices.size());
           vertex++)
      {
         AddSegment(transmitter.vertices[vertex],
                    transmitter.vertices[vertex + 1]);
      }
   }

   void Eval(Vector &value, ElementTransformation &T,
             const IntegrationPoint &ip) override
   {
      Vector scaled_x;
      T.Transform(ip, scaled_x);
      const real_t radius2 = radius * radius;
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
                       radius) -
                   erf((axial_distance - segment.half_length) /
                       radius));
         const real_t profile =
            longitudinal * exp(-transverse2 / radius2) / (M_PI * radius2);
         value.Add(right_hand_side_scale * current * profile,
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
   const Operator &jacobian;
   const Vector &standard_deviation;
   const SparseMatrix &regularization;
   real_t regularization_weight;
   real_t damping;
   mutable Vector weighted_action;
   mutable Vector regularization_action;
   mutable Vector regularization_normal;

public:
   GaussNewtonNormalOperator(const Operator &jacobian_,
                             const Vector &standard_deviation_,
                             const SparseMatrix &regularization_,
                             real_t regularization_weight_, real_t damping_)
      : Operator(jacobian_.Width()), jacobian(jacobian_),
        standard_deviation(standard_deviation_),
        regularization(regularization_),
        regularization_weight(regularization_weight_), damping(damping_),
        weighted_action(jacobian_.Height()),
        regularization_action(regularization_.Height()),
        regularization_normal(jacobian_.Width())
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
      if (regularization.Height() > 0)
      {
         regularization.Mult(x, regularization_action);
         regularization.MultTranspose(regularization_action,
                                      regularization_normal);
         y.Add(regularization_weight, regularization_normal);
      }
      y.Add(damping, x);
   }
};

class CSEMProblem
{
private:
   MPI_Comm communicator;
   int rank;
   int world_rank;
   int survey_group;
   int survey_group_count;
   const CSEMOptions &options;
   real_t length_scale = 1.0;
   unique_ptr<ParMesh> mesh;
   unique_ptr<ND_FECollection> collection;
   unique_ptr<ParFiniteElementSpace> space;
   unique_ptr<ParComplexGridFunction> field;
   unique_ptr<TaggedMaterialCatalog> materials;
   unique_ptr<RT_FECollection> magnetic_collection;
   unique_ptr<ParFiniteElementSpace> magnetic_space;
   unique_ptr<ParDiscreteLinearOperator> curl;
   unique_ptr<CSEMConductivity> conductivity;
   unique_ptr<RelativePermittivity> permittivity;
   unique_ptr<SparseMatrix> regularization_operator;
   vector<Transmitter> transmitters;
   Array<int> essential_tdofs;
   Array<int> absorbing_boundary;
   vector<Vector> receiver_centers;
   vector<Vector> receiver_directions;
   vector<unique_ptr<GaussianReceiver>> receiver_coefficients;
   vector<unique_ptr<ParLinearForm>> receiver_functionals;
   Array<int> receiver_elements;
   Array<IntegrationPoint> receiver_integration_points;

   void BuildRegularization()
   {
      const int units = conductivity->UnitCount();
      const int parameters = conductivity->ParameterCount();
      set<pair<int, int>> edges;
      if (units > 1)
      {
         Vector centers(3 * units), counts(units);
         centers = 0.0;
         counts = 0.0;
         for (int element = 0; element < mesh->GetNE(); element++)
         {
            ElementTransformation *transformation =
               mesh->GetElementTransformation(element);
            const IntegrationPoint &center = Geometries.GetCenter(
               mesh->GetElementBaseGeometry(element));
            Vector x;
             transformation->Transform(center, x);
             const int unit = conductivity->UnitAt(x, element);
             if (unit < 0) { continue; }
             for (int d = 0; d < 3; d++)
            {
               centers[3 * unit + d] += length_scale * x[d];
            }
            counts[unit] += 1.0;
         }
         MPI_Allreduce(MPI_IN_PLACE, centers.GetData(), centers.Size(),
                       MPITypeMap<real_t>::mpi_type, MPI_SUM, communicator);
         MPI_Allreduce(MPI_IN_PLACE, counts.GetData(), counts.Size(),
                       MPITypeMap<real_t>::mpi_type, MPI_SUM, communicator);
         for (int unit = 0; unit < units; unit++)
         {
            if (counts[unit] > 0.0)
            {
               for (int d = 0; d < 3; d++)
               {
                  centers[3 * unit + d] /= counts[unit];
               }
            }
         }
         const int neighbors = min(options.regularization_neighbors,
                                   units - 1);
         for (int unit = 0; unit < units; unit++)
         {
            vector<pair<real_t, int>> distances;
            for (int other = 0; other < units; other++)
            {
               if (other == unit) { continue; }
               real_t distance2 = 0.0;
               for (int d = 0; d < 3; d++)
               {
                  const real_t difference = centers[3 * unit + d] -
                                             centers[3 * other + d];
                  distance2 += difference * difference;
               }
               distances.emplace_back(distance2, other);
            }
            sort(distances.begin(), distances.end());
            for (int neighbor = 0; neighbor < neighbors; neighbor++)
            {
               const int other = distances[neighbor].second;
               edges.insert(minmax(unit, other));
            }
         }
      }

      regularization_operator.reset(
         new SparseMatrix(3 * static_cast<int>(edges.size()), parameters));
      const real_t scales[3] = {options.regularization_rh_scale,
                                options.regularization_rv_scale,
                                options.regularization_dip_scale};
      int edge = 0;
      for (const auto &connection : edges)
      {
         for (int property = 0; property < 3; property++)
         {
            const int row = 3 * edge + property;
            regularization_operator->Add(
               row, 3 * connection.first + property, -scales[property]);
            regularization_operator->Add(
               row, 3 * connection.second + property, scales[property]);
         }
         edge++;
      }
      regularization_operator->Finalize();
   }

   bool MagneticReceivers() const
   {
      return string(options.receiver_field) == "H" ||
             string(options.receiver_field) == "h";
   }

   bool PointReceivers() const
   {
      return string(options.receiver_model) == "point";
   }

   bool AbsorbingBoundary() const
   {
      return string(options.boundary_condition) == "absorbing";
   }

   complex<real_t> AbsorbingBoundaryCoefficient(real_t frequency) const
   {
      const real_t omega = 2.0 * M_PI * frequency;
      const real_t wave_number = omega * sqrt(mu0 * epsilon0);
      const complex<real_t> imaginary_unit(0.0, 1.0);
      const complex<real_t> complex_permittivity(
         options.absorbing_relative_permittivity,
         -options.absorbing_conductivity / (omega * epsilon0));
      const complex<real_t> exterior_wave_number = wave_number *
         sqrt(options.absorbing_relative_permeability * complex_permittivity);
      return imaginary_unit * length_scale * exterior_wave_number /
             options.absorbing_relative_permeability;
   }

   void BuildBoundaryConditions()
   {
      const int maximum_attribute = mesh->bdr_attributes.Size()
                                       ? mesh->bdr_attributes.Max() : 0;
      Array<int> essential_boundary(maximum_attribute);
      essential_boundary = 1;
      absorbing_boundary.DeleteAll();
      if (AbsorbingBoundary())
      {
         MFEM_VERIFY(maximum_attribute > 0,
                     "Absorbing treatment requires boundary attributes.");
         absorbing_boundary.SetSize(maximum_attribute);
         absorbing_boundary = 0;
         const vector<int> requested = ParseIntegerList(
            options.absorbing_boundary_attributes);
         if (requested.empty())
         {
            absorbing_boundary = 1;
            essential_boundary = 0;
         }
         else
         {
            for (const int attribute : requested)
            {
               MFEM_VERIFY(attribute > 0 && attribute <= maximum_attribute &&
                              mesh->bdr_attributes.Find(attribute) >= 0,
                           "Absorbing boundary attribute " << attribute
                           << " is not present in the mesh.");
               absorbing_boundary[attribute - 1] = 1;
               essential_boundary[attribute - 1] = 0;
            }
         }
      }
      space->GetEssentialTrueDofs(essential_boundary, essential_tdofs);
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
         if (PointReceivers()) { return; }
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
      }
      else
      {
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

      if (PointReceivers())
      {
         DenseMatrix points(3, static_cast<int>(receiver_centers.size()));
         for (int receiver = 0;
              receiver < static_cast<int>(receiver_centers.size()); receiver++)
         {
            for (int d = 0; d < 3; d++)
            {
               points(d, receiver) = receiver_centers[receiver][d] /
                                     length_scale;
            }
         }
         const int found = mesh->FindPoints(points, receiver_elements,
                                            receiver_integration_points,
                                            false);
         MFEM_VERIFY(found == static_cast<int>(receiver_centers.size()),
                     "One or more point receivers lie outside the mesh.");
      }
   }

   void ExtractResponses(const ParComplexGridFunction &solution,
                          real_t frequency,
                          Vector &responses) const
   {
      responses.SetSize(2 * receiver_centers.size());
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
      if (PointReceivers())
      {
         responses = 0.0;
         Vector real_value(3), imaginary_value(3);
         for (int receiver = 0;
              receiver < static_cast<int>(receiver_centers.size()); receiver++)
         {
            const int element = receiver_elements[receiver];
            if (element < 0) { continue; }
            const ParGridFunction &real_field = MagneticReceivers()
               ? *magnetic_real : solution.real();
            const ParGridFunction &imaginary_field = MagneticReceivers()
               ? *magnetic_imaginary : solution.imag();
            real_field.GetVectorValue(
               element, receiver_integration_points[receiver], real_value);
            imaginary_field.GetVectorValue(
               element, receiver_integration_points[receiver],
               imaginary_value);
            responses[2 * receiver] =
               real_value * receiver_directions[receiver];
            responses[2 * receiver + 1] =
               imaginary_value * receiver_directions[receiver];
         }
         MPI_Allreduce(MPI_IN_PLACE, responses.GetData(), responses.Size(),
                       MPITypeMap<real_t>::mpi_type, MPI_SUM, communicator);
         return;
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

   void AssembleReceiverTranspose(const Vector &weights, real_t frequency,
                                  Vector &true_rhs) const
   {
      MFEM_VERIFY(!PointReceivers(),
                  "Point receivers are forward-regression probes and do not "
                  "define an adjoint transpose. Use Gaussian receivers for "
                  "Jacobian verification and inversion.");
      MFEM_VERIFY(weights.Size() == DataSize(),
                  "The adjoint receiver-weight vector has the wrong size.");
      Vector local_real(space->GetVSize()), local_imaginary(space->GetVSize());
      local_real = 0.0;
      local_imaginary = 0.0;

      if (MagneticReceivers())
      {
         Vector magnetic_real(magnetic_space->GetVSize());
         Vector magnetic_imaginary(magnetic_space->GetVSize());
         magnetic_real = 0.0;
         magnetic_imaginary = 0.0;
         const real_t magnetic_scale =
            1.0 / (2.0 * M_PI * frequency * mu0 *
                   options.relative_permeability * length_scale);
         for (int receiver = 0;
              receiver < static_cast<int>(receiver_functionals.size());
              receiver++)
         {
            // H_I = scale curl(E_R), H_R = -scale curl(E_I).
            magnetic_real.Add(magnetic_scale * weights[2 * receiver + 1],
                              *receiver_functionals[receiver]);
            magnetic_imaginary.Add(-magnetic_scale * weights[2 * receiver],
                                   *receiver_functionals[receiver]);
         }
         curl->MultTranspose(magnetic_real, local_real);
         curl->MultTranspose(magnetic_imaginary, local_imaginary);
      }
      else
      {
         for (int receiver = 0;
              receiver < static_cast<int>(receiver_functionals.size());
              receiver++)
         {
            local_real.Add(weights[2 * receiver],
                           *receiver_functionals[receiver]);
            local_imaginary.Add(weights[2 * receiver + 1],
                                *receiver_functionals[receiver]);
         }
      }

      const int block_size = space->GetTrueVSize();
      true_rhs.SetSize(2 * block_size);
      Vector true_real, true_imaginary;
      true_real.MakeRef(true_rhs, 0, block_size);
      true_imaginary.MakeRef(true_rhs, block_size, block_size);
      const Operator *prolongation = space->GetProlongationMatrix();
      prolongation->MultTranspose(local_real, true_real);
      prolongation->MultTranspose(local_imaginary, true_imaginary);
      for (int i = 0; i < essential_tdofs.Size(); i++)
      {
         true_real[essential_tdofs[i]] = 0.0;
         true_imaginary[essential_tdofs[i]] = 0.0;
      }
   }

   void AccumulateParameterGradient(
      const ParComplexGridFunction &adjoint_field,
      real_t scaled_conductivity, Vector &gradient) const
   {
      Vector local_gradient(conductivity->ParameterCount());
      local_gradient = 0.0;
      Vector electric_real(3), electric_imaginary(3);
      Vector adjoint_real(3), adjoint_imaginary(3);
      Vector derivative_electric_real(3), derivative_electric_imaginary(3);
      Vector x;
      DenseMatrix derivative(3);

      for (int element_number = 0; element_number < mesh->GetNE();
           element_number++)
      {
         ElementTransformation *transformation =
            mesh->GetElementTransformation(element_number);
         const FiniteElement *element = space->GetFE(element_number);
         const int integration_order =
            transformation->OrderW() + 2 * element->GetOrder();
         const IntegrationRule &rule =
            IntRules.Get(element->GetGeomType(), integration_order);
         for (int q = 0; q < rule.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = rule.IntPoint(q);
            transformation->SetIntPoint(&ip);
            field->real().GetVectorValue(element_number, ip, electric_real);
            field->imag().GetVectorValue(element_number, ip,
                                         electric_imaginary);
            adjoint_field.real().GetVectorValue(element_number, ip,
                                                adjoint_real);
            adjoint_field.imag().GetVectorValue(element_number, ip,
                                                adjoint_imaginary);
            transformation->Transform(ip, x);
            const int unit = conductivity->UnitAt(x, element_number);
            if (unit < 0) { continue; }
            const real_t weight = scaled_conductivity * ip.weight *
                                  transformation->Weight();
            for (int property = 0; property < 3; property++)
            {
               const int parameter = 3 * unit + property;
               conductivity->ParameterDerivative(x, element_number,
                                                  parameter, derivative);
               derivative.Mult(electric_imaginary,
                               derivative_electric_imaginary);
               derivative.Mult(electric_real, derivative_electric_real);
               local_gradient[parameter] += weight *
                  ((adjoint_real * derivative_electric_imaginary) -
                   (adjoint_imaginary * derivative_electric_real));
            }
         }
      }

      gradient.SetSize(local_gradient.Size());
      MPI_Allreduce(local_gradient.GetData(), gradient.GetData(),
                    local_gradient.Size(), MPITypeMap<real_t>::mpi_type,
                    MPI_SUM, communicator);
   }

public:
   CSEMProblem(MPI_Comm communicator_, const CSEMOptions &options_,
               int survey_group_, int survey_group_count_)
      : communicator(communicator_), survey_group(survey_group_),
        survey_group_count(survey_group_count_), options(options_)
   {
      MPI_Comm_rank(communicator, &rank);
      MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
      BuildMesh();
      materials.reset(new TaggedMaterialCatalog(options, *mesh));
      collection.reset(new ND_FECollection(
         options.order, 3, BasisType::GaussLobatto,
         BasisType::GaussLegendre));
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
      conductivity.reset(new CSEMConductivity(options, *mesh, *materials,
                                              length_scale));
      permittivity.reset(new RelativePermittivity(options, *mesh, *materials,
                                                  length_scale));
      BuildRegularization();
      transmitters = LoadTransmitters(options);

      BuildBoundaryConditions();
      BuildReceivers();

      if (rank == 0)
      {
         cout << "CSEM mesh: " << mesh->GetGlobalNE() << " hexahedra, "
              << space->GlobalTrueVSize() << " true H(curl) DOFs\n"
              << "Length scaling L = " << length_scale << " m\n"
              << "Model: " << conductivity->UnitCount() << " unit(s), "
              << conductivity->ParameterCount() << " parameter(s), "
              << regularization_operator->Height() / 3
              << " spatial edge(s)\n"
              << "Sensitivity: matrix-free tangent/adjoint actions; "
                 "no survey Jacobian is stored\n"
              << "Materials: "
              << (materials->Enabled()
                     ? "volume-attribute tagged (topography conforming)"
                     : "legacy z-split") << "\n"
              << "Boundary: " << options.boundary_condition;
         if (AbsorbingBoundary())
         {
            cout << " on " << absorbing_boundary.Sum()
                 << " boundary attribute(s)";
         }
         cout << '\n';
      }
   }

   int Rank() const { return rank; }
   int WorldRank() const { return world_rank; }
   int SurveyGroup() const { return survey_group; }
   int SurveyGroupCount() const { return survey_group_count; }
   int SourceCount() const { return static_cast<int>(transmitters.size()); }
   int ParameterCount() const { return conductivity->ParameterCount(); }
   const Vector &TrueModel() const { return conductivity->TrueModel(); }
   const Vector &InitialModel() const { return conductivity->InitialModel(); }
   const SparseMatrix &Regularization() const
   {
      return *regularization_operator;
   }
   string ModelMode() const { return options.model_parameterization; }
   int ModelUnitId(int unit) const
   {
      return conductivity->UnitId(unit);
   }
   int DataSize() const { return 2 * static_cast<int>(receiver_centers.size()); }
   const Vector &ReceiverCenter(int receiver) const
   {
      return receiver_centers[receiver];
   }
   const Vector &ReceiverDirection(int receiver) const
   {
      return receiver_directions[receiver];
   }

   void Solve(const Vector &log_model, real_t frequency, int source_index,
              Vector &responses, const SolverSettings &settings,
              const Vector *tangent_direction = nullptr,
              Vector *tangent_responses = nullptr,
              const Vector *adjoint_weights = nullptr,
              Vector *adjoint_gradient = nullptr)
   {
      MFEM_VERIFY(source_index >= 0 && source_index < SourceCount(),
                  "Transmitter index is out of range.");
      MFEM_VERIFY((tangent_direction == nullptr) ==
                     (tangent_responses == nullptr),
                  "Tangent direction and response must be supplied together.");
      MFEM_VERIFY((adjoint_weights == nullptr) ==
                     (adjoint_gradient == nullptr),
                  "Adjoint weights and gradient must be supplied together.");
      conductivity->SetModel(log_model);
      const real_t omega = 2.0 * M_PI * frequency;
      const real_t scaled_conductivity =
         omega * mu0 * length_scale * length_scale;
      const real_t wave_number = omega * sqrt(mu0 * epsilon0);
      const real_t scaled_permittivity =
         wave_number * wave_number * length_scale * length_scale;

      ConstantCoefficient inverse_relative_permeability(
         1.0 / options.relative_permeability);
      ScaledMatrixCoefficient conduction(*conductivity, scaled_conductivity);
      ScaledCoefficient dielectric(*permittivity, -scaled_permittivity);
      ShiftedMassCoefficient shifted_mass(
         *conductivity, *permittivity, scaled_conductivity,
         scaled_permittivity);
      unique_ptr<ConstantCoefficient> absorbing_real;
      unique_ptr<ConstantCoefficient> absorbing_imaginary;
      unique_ptr<ConstantCoefficient> absorbing_shifted;
      if (AbsorbingBoundary())
      {
         const complex<real_t> gamma =
            AbsorbingBoundaryCoefficient(frequency);
         absorbing_real.reset(new ConstantCoefficient(gamma.real()));
         absorbing_imaginary.reset(new ConstantCoefficient(gamma.imag()));
         absorbing_shifted.reset(new ConstantCoefficient(gamma.real() +
                                                          gamma.imag()));
         if (settings.verbose && rank == 0)
         {
            cout << "Absorbing boundary gamma: " << gamma.real() << " + i "
                 << gamma.imag() << '\n';
         }
      }

      ParComplexLinearForm right_hand_side(
         space.get(), ComplexOperator::HERMITIAN);
      unique_ptr<GaussianWireCurrent> gaussian_source;
      vector<unique_ptr<VectorDeltaCoefficient>> line_deltas;
      if (string(options.source_model) == "gaussian")
      {
         gaussian_source.reset(new GaussianWireCurrent(
            transmitters[source_index], length_scale, omega));
         right_hand_side.AddDomainIntegrator(
            nullptr, new VectorFEDomainLFIntegrator(*gaussian_source));
      }
      else
      {
         const Transmitter &transmitter = transmitters[source_index];
         const IntegrationRule &line_rule =
            IntRules.Get(Geometry::SEGMENT, options.source_quadrature_order);
         for (int segment = 0;
              segment + 1 < static_cast<int>(transmitter.vertices.size());
              segment++)
         {
            const Vector &start = transmitter.vertices[segment];
            const Vector &end = transmitter.vertices[segment + 1];
            Vector direction(end);
            direction -= start;
            const real_t physical_length = direction.Norml2();
            MFEM_VERIFY(physical_length > 0.0,
                        "Transmitter segments must have positive length.");
            direction /= physical_length;
            for (int q = 0; q < line_rule.GetNPoints(); q++)
            {
               const IntegrationPoint &point = line_rule.IntPoint(q);
               Vector center(start);
               center.Add(point.x, end);
               center.Add(-point.x, start);
               center /= length_scale;
               const real_t delta_scale =
                  -omega * mu0 * transmitter.current *
                  (physical_length / length_scale) * point.weight;
               line_deltas.emplace_back(new VectorDeltaCoefficient(
                  direction, center[0], center[1], center[2], delta_scale));
               right_hand_side.AddDomainIntegrator(
                  nullptr,
                  new VectorFEDomainLFIntegrator(*line_deltas.back()));
            }
         }
      }
      right_hand_side.Assemble();

      *field = complex<real_t>(0.0, 0.0);
      ParSesquilinearForm form(space.get(), ComplexOperator::HERMITIAN);
      form.AddDomainIntegrator(
         new CurlCurlIntegrator(inverse_relative_permeability), nullptr);
      form.AddDomainIntegrator(new VectorFEMassIntegrator(dielectric),
                               new VectorFEMassIntegrator(conduction));
      if (AbsorbingBoundary())
      {
         form.AddBoundaryIntegrator(
            new VectorFEMassIntegrator(*absorbing_real),
            new VectorFEMassIntegrator(*absorbing_imaginary),
            absorbing_boundary);
      }
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
         if (AbsorbingBoundary())
         {
            shifted_form->AddBoundaryIntegrator(
               new VectorFEMassIntegrator(*absorbing_shifted),
               absorbing_boundary);
         }
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

       if (tangent_direction)
       {
          MFEM_VERIFY(tangent_direction->Size() ==
                         conductivity->ParameterCount(),
                      "The tangent model direction has the wrong size.");
          const int block_size = space->GetTrueVSize();
          Vector electric_real, electric_imaginary;
          electric_real.MakeRef(X, 0, block_size);
          electric_imaginary.MakeRef(X, block_size, block_size);

          DirectionalDerivativeCoefficient derivative_coefficient(
             *conductivity, *tangent_direction, scaled_conductivity);
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
          ExtractResponses(tangent_field, frequency, *tangent_responses);
       }

       if (adjoint_weights)
       {
          Vector adjoint_rhs;
          AssembleReceiverTranspose(*adjoint_weights, frequency, adjoint_rhs);
          const int block_size = space->GetTrueVSize();
          Vector conjugated_rhs(adjoint_rhs);
          Vector conjugated_imaginary;
          conjugated_imaginary.MakeRef(conjugated_rhs, block_size, block_size);
          conjugated_imaginary *= -1.0;

          // The real block transpose represents the conjugated complex
          // system. Solve A y = conj(q), then z = conj(y) satisfies A^T z=q.
          Vector adjoint_solution;
          solve_linear_system(conjugated_rhs, adjoint_solution);
          Vector adjoint_real, adjoint_imaginary;
          adjoint_real.MakeRef(adjoint_solution, 0, block_size);
          adjoint_imaginary.MakeRef(adjoint_solution, block_size, block_size);
          adjoint_imaginary *= -1.0;

          ParComplexGridFunction adjoint_field(space.get());
          adjoint_field = complex<real_t>(0.0, 0.0);
          adjoint_field.real().SetFromTrueDofs(adjoint_real);
          adjoint_field.imag().SetFromTrueDofs(adjoint_imaginary);
          adjoint_field.Sync();
          AccumulateParameterGradient(adjoint_field, scaled_conductivity,
                                      *adjoint_gradient);
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
      ParGridFunction rh(&material_space), rv(&material_space);
      ParGridFunction sigma_h(&material_space), sigma_v(&material_space);
      ParGridFunction dip(&material_space);
      ParGridFunction epsilon_r(&material_space);
      ParGridFunction material_attribute(&material_space);
      ParGridFunction invertible_material(&material_space);
      ConductivityPropertyCoefficient rh_coefficient(*conductivity, 0);
      ConductivityPropertyCoefficient rv_coefficient(*conductivity, 1);
      ConductivityPropertyCoefficient dip_coefficient(*conductivity, 2);
      ConductivityPropertyCoefficient sigma_h_coefficient(*conductivity, 3);
      ConductivityPropertyCoefficient sigma_v_coefficient(*conductivity, 4);
      MaterialTagCoefficient material_attribute_coefficient(
         *mesh, *materials, length_scale, false);
      MaterialTagCoefficient invertible_material_coefficient(
         *mesh, *materials, length_scale, true);
      rh.ProjectCoefficient(rh_coefficient);
      rv.ProjectCoefficient(rv_coefficient);
      dip.ProjectCoefficient(dip_coefficient);
      sigma_h.ProjectCoefficient(sigma_h_coefficient);
      sigma_v.ProjectCoefficient(sigma_v_coefficient);
      epsilon_r.ProjectCoefficient(*permittivity);
      material_attribute.ProjectCoefficient(material_attribute_coefficient);
      invertible_material.ProjectCoefficient(invertible_material_coefficient);

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
      collection_output.RegisterField("Rh_ohm_m", &rh);
      collection_output.RegisterField("Rv_ohm_m", &rv);
      collection_output.RegisterField("sigma_h_S_per_m", &sigma_h);
      collection_output.RegisterField("sigma_v_S_per_m", &sigma_v);
      collection_output.RegisterField("dip_degrees", &dip);
      collection_output.RegisterField("relative_permittivity", &epsilon_r);
      collection_output.RegisterField("material_attribute",
                                      &material_attribute);
      collection_output.RegisterField("invertible_material",
                                      &invertible_material);
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
                  const SolverSettings &settings,
                  const string &field_prefix = "")
{
   const int block_size = problem.DataSize();
   const int frequency_count = static_cast<int>(frequencies.size());
   const int source_count = problem.SourceCount();
   responses.SetSize(block_size * source_count * frequency_count);
   responses = 0.0;

   for (int frequency_index = 0;
         frequency_index < frequency_count;
        frequency_index++)
   {
      const real_t frequency = frequencies[frequency_index];
      for (int source = 0; source < source_count; source++)
      {
         const int task = frequency_index * source_count + source;
         if (task % problem.SurveyGroupCount() != problem.SurveyGroup())
         {
            continue;
         }
         if (settings.verbose && problem.Rank() == 0)
         {
            cout << "Solving frequency " << frequency << " Hz ("
                 << frequency_index + 1 << '/' << frequencies.size()
                 << "), transmitter " << source + 1 << '/' << source_count
                 << '\n';
         }
         Vector source_responses;
         problem.Solve(log_model, frequency, source, source_responses,
                       settings);
         for (int datum = 0; datum < block_size; datum++)
         {
            const int global_datum =
               (frequency_index * source_count + source) * block_size + datum;
            responses[global_datum] = source_responses[datum];
         }
         if (!field_prefix.empty())
         {
            string source_prefix = field_prefix;
            if (source_count > 1)
            {
               ostringstream name;
               name << field_prefix << "_s" << setw(3) << setfill('0')
                    << source;
               source_prefix = name.str();
            }
            problem.SaveField(FrequencyFieldPrefix(
                                 source_prefix, frequency_index,
                                 frequency_count, frequency),
                              frequency);
         }
      }
   }

   // Only one rank per survey group contributes its already group-reduced
   // receiver values. Summation then assembles the disjoint task blocks.
   if (problem.Rank() != 0)
   {
      responses = 0.0;
   }
   MPI_Allreduce(MPI_IN_PLACE, responses.GetData(), responses.Size(),
                  MPITypeMap<real_t>::mpi_type, MPI_SUM, MPI_COMM_WORLD);
}

void ApplySurveyTangent(CSEMProblem &problem,
                        const vector<real_t> &frequencies,
                        const Vector &log_model, const Vector &direction,
                        Vector &tangent_responses,
                        const SolverSettings &settings)
{
   MFEM_VERIFY(direction.Size() == problem.ParameterCount(),
               "The survey tangent direction has the wrong size.");
   const int block_size = problem.DataSize();
   const int frequency_count = static_cast<int>(frequencies.size());
   const int source_count = problem.SourceCount();
   tangent_responses.SetSize(block_size * source_count * frequency_count);
   tangent_responses = 0.0;

   for (int frequency_index = 0; frequency_index < frequency_count;
        frequency_index++)
   {
      for (int source = 0; source < source_count; source++)
      {
         const int task = frequency_index * source_count + source;
         if (task % problem.SurveyGroupCount() != problem.SurveyGroup())
         {
            continue;
         }
         Vector source_responses, source_tangent;
         problem.Solve(log_model, frequencies[frequency_index], source,
                       source_responses, settings, &direction,
                       &source_tangent);
         for (int datum = 0; datum < block_size; datum++)
         {
            tangent_responses[task * block_size + datum] =
               source_tangent[datum];
         }
      }
   }

   if (problem.Rank() != 0) { tangent_responses = 0.0; }
   MPI_Allreduce(MPI_IN_PLACE, tangent_responses.GetData(),
                 tangent_responses.Size(), MPITypeMap<real_t>::mpi_type,
                 MPI_SUM, MPI_COMM_WORLD);
}

void ApplySurveyAdjoint(CSEMProblem &problem,
                        const vector<real_t> &frequencies,
                        const Vector &log_model, const Vector &data_weights,
                        Vector &model_gradient,
                        const SolverSettings &settings)
{
   const int block_size = problem.DataSize();
   const int frequency_count = static_cast<int>(frequencies.size());
   const int source_count = problem.SourceCount();
   MFEM_VERIFY(data_weights.Size() ==
                  block_size * source_count * frequency_count,
               "The survey adjoint data vector has the wrong size.");
   model_gradient.SetSize(problem.ParameterCount());
   model_gradient = 0.0;

   for (int frequency_index = 0; frequency_index < frequency_count;
        frequency_index++)
   {
      for (int source = 0; source < source_count; source++)
      {
         const int task = frequency_index * source_count + source;
         if (task % problem.SurveyGroupCount() != problem.SurveyGroup())
         {
            continue;
         }
         Vector source_weights(block_size);
         for (int datum = 0; datum < block_size; datum++)
         {
            source_weights[datum] =
               data_weights[task * block_size + datum];
         }
         Vector source_responses, source_gradient;
         problem.Solve(log_model, frequencies[frequency_index], source,
                       source_responses, settings, nullptr, nullptr,
                       &source_weights, &source_gradient);
         model_gradient += source_gradient;
      }
   }

   if (problem.Rank() != 0) { model_gradient = 0.0; }
   MPI_Allreduce(MPI_IN_PLACE, model_gradient.GetData(), model_gradient.Size(),
                 MPITypeMap<real_t>::mpi_type, MPI_SUM, MPI_COMM_WORLD);
}

class MatrixFreeSurveyJacobian : public Operator
{
private:
   CSEMProblem &problem;
   const vector<real_t> &frequencies;
   const Vector &model;
   const SolverSettings &settings;

public:
   MatrixFreeSurveyJacobian(CSEMProblem &problem_,
                            const vector<real_t> &frequencies_,
                            const Vector &model_,
                            const SolverSettings &settings_)
      : Operator(problem_.DataSize() * problem_.SourceCount() *
                    static_cast<int>(frequencies_.size()),
                 problem_.ParameterCount()),
        problem(problem_), frequencies(frequencies_), model(model_),
        settings(settings_)
   {
   }

   void Mult(const Vector &x, Vector &y) const override
   {
      ApplySurveyTangent(problem, frequencies, model, x, y, settings);
   }

   void MultTranspose(const Vector &x, Vector &y) const override
   {
      ApplySurveyAdjoint(problem, frequencies, model, x, y, settings);
   }
};

real_t Objective(const Vector &prediction, const Vector &observed,
                 const Vector &standard_deviation, const Vector &model,
                 const Vector &reference_model,
                 const SparseMatrix &regularization_operator,
                 real_t regularization)
{
   real_t value = 0.0;
   for (int datum = 0; datum < prediction.Size(); datum++)
   {
      const real_t residual =
         (prediction[datum] - observed[datum]) / standard_deviation[datum];
      value += residual * residual;
   }
   if (regularization_operator.Height() > 0)
   {
      Vector difference(model);
      difference -= reference_model;
      Vector regularized(regularization_operator.Height());
      regularization_operator.Mult(difference, regularized);
      value += regularization * (regularized * regularized);
   }
   return 0.5 * value;
}

real_t DataObjective(const Vector &prediction, const Vector &observed,
                     const Vector &standard_deviation)
{
   real_t value = 0.0;
   for (int datum = 0; datum < prediction.Size(); datum++)
   {
      const real_t residual =
         (prediction[datum] - observed[datum]) / standard_deviation[datum];
      value += residual * residual;
   }
   return 0.5 * value;
}

real_t SpatialObjective(const Vector &model, const Vector &reference_model,
                        const SparseMatrix &regularization_operator,
                        real_t regularization)
{
   if (regularization_operator.Height() == 0) { return 0.0; }
   Vector difference(model);
   difference -= reference_model;
   Vector regularized(regularization_operator.Height());
   regularization_operator.Mult(difference, regularized);
   return 0.5 * regularization * (regularized * regularized);
}

void WriteModel(const string &file_name, const CSEMProblem &problem,
                const Vector &model)
{
   if (problem.WorldRank() != 0) { return; }
   ofstream output(file_name);
   MFEM_VERIFY(output, "Unable to open model output file: " << file_name);
   output << "parameterization,id,Rh_ohm_m,Rv_ohm_m,dip_degrees\n";
   output << setprecision(16);
   for (int unit = 0; unit < problem.ParameterCount() / 3; unit++)
   {
      output << problem.ModelMode() << ',' << problem.ModelUnitId(unit) << ','
             << exp(model[3 * unit]) << ',' << exp(model[3 * unit + 1])
             << ',' << model[3 * unit + 2] * 180.0 / M_PI << '\n';
   }
}

void WriteData(const string &file_name, const CSEMProblem &problem,
               const vector<real_t> &frequencies, const Vector &data,
               const Vector &standard_deviation)
{
   if (problem.WorldRank() != 0) { return; }
   ofstream output(file_name);
   MFEM_VERIFY(output, "Unable to open output data file: " << file_name);
   output << "frequency_hz,source,receiver,x_m,y_m,z_m,direction_x,direction_y,"
             "direction_z,real,imag,"
             "stddev_real,stddev_imag\n";
   output << setprecision(16);
   const int frequency_block_size = problem.DataSize();
   const int source_count = problem.SourceCount();
   MFEM_VERIFY(data.Size() == frequency_block_size * source_count *
                              frequencies.size(),
               "Survey data size does not match the frequency list.");
   for (int frequency_index = 0;
        frequency_index < static_cast<int>(frequencies.size());
        frequency_index++)
   {
      for (int source = 0; source < source_count; source++)
      {
         for (int receiver = 0;
              receiver < frequency_block_size / 2; receiver++)
         {
            const int offset =
               (frequency_index * source_count + source) *
               frequency_block_size + 2 * receiver;
            const Vector &center = problem.ReceiverCenter(receiver);
            const Vector &direction = problem.ReceiverDirection(receiver);
            output << frequencies[frequency_index] << ',' << source << ','
                   << receiver << ',' << center[0] << ',' << center[1] << ','
                   << center[2] << ',' << direction[0] << ','
                   << direction[1] << ',' << direction[2] << ','
                   << data[offset] << ',' << data[offset + 1] << ','
                   << standard_deviation[offset] << ','
                   << standard_deviation[offset + 1] << '\n';
         }
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
      if (values.size() == 13)
      {
         data_values.push_back(values[9]);
         data_values.push_back(values[10]);
         error_values.push_back(values[11]);
         error_values.push_back(values[12]);
      }
      else if (values.size() == 12)
      {
         // Backward-compatible reader for the one-source oriented CSV.
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
               "transmitters times receivers.");
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
   Vector true_model(problem.TrueModel());
   SolveSurvey(problem, frequencies, true_model, observed, forward_settings);

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
   SolveSurvey(problem, frequencies, model, prediction, settings);
   const string file_name =
      string(options.output_prefix) + "_jacobian_check.csv";
   ofstream output;
   if (problem.WorldRank() == 0)
   {
      output.open(file_name);
      MFEM_VERIFY(output, "Unable to open Jacobian verification output.");
      output << "frequency_hz,source,receiver,part,parameter,analytic,"
                "finite_difference,absolute_error,relative_error\n";
      output << setprecision(16);
   }
   real_t maximum_relative_error = 0.0;
   real_t squared_error = 0.0;
   real_t squared_reference = 0.0;
   const int block_size = problem.DataSize();
   const int source_count = problem.SourceCount();
   const char *parts[2] = {"real", "imag"};
   Vector direction(model.Size());
   for (int parameter = 0; parameter < model.Size(); parameter++)
   {
      direction = 0.0;
      direction[parameter] = 1.0;
      Vector analytic_column;
      ApplySurveyTangent(problem, frequencies, model, direction,
                         analytic_column, settings);

      Vector plus_model(model), minus_model(model);
      plus_model[parameter] += options.jacobian_fd_step;
      minus_model[parameter] -= options.jacobian_fd_step;
      Vector plus_response, minus_response;
      SolveSurvey(problem, frequencies, plus_model, plus_response, settings);
      SolveSurvey(problem, frequencies, minus_model, minus_response, settings);
      for (int datum = 0; datum < prediction.Size(); datum++)
      {
         const real_t exact = analytic_column[datum];
         const real_t approximation =
            (plus_response[datum] - minus_response[datum]) /
            (2.0 * options.jacobian_fd_step);
         const real_t absolute_error = fabs(exact - approximation);
         const real_t relative_error = absolute_error /
            max(max(fabs(exact), fabs(approximation)), real_t(1.0e-30));
         maximum_relative_error = max(maximum_relative_error, relative_error);
         squared_error += absolute_error * absolute_error;
         squared_reference += approximation * approximation;
         if (problem.WorldRank() == 0)
         {
            const int frequency_index =
               datum / (source_count * block_size);
            const int frequency_datum =
               datum % (source_count * block_size);
            const int source = frequency_datum / block_size;
            const int local_datum = frequency_datum % block_size;
            const int receiver = local_datum / 2;
            const int part = local_datum % 2;
            const char *property = parameter % 3 == 0 ? "log_Rh" :
                                   parameter % 3 == 1 ? "log_Rv" :
                                                        "dip_rad";
            output << frequencies[frequency_index] << ',' << source << ','
                   << receiver << ',' << parts[part] << ','
                   << parameter / 3 << ':' << property << ','
                   << exact << ',' << approximation << ',' << absolute_error
                   << ',' << relative_error << '\n';
         }
      }
   }

   Vector dot_direction(model.Size()), dot_weights(prediction.Size());
   for (int i = 0; i < dot_direction.Size(); i++)
   {
      dot_direction[i] = sin(0.37 * (i + 1)) + 0.1 * cos(0.11 * (i + 1));
   }
   for (int i = 0; i < dot_weights.Size(); i++)
   {
      dot_weights[i] = cos(0.23 * (i + 1)) - 0.2 * sin(0.07 * (i + 1));
   }
   Vector tangent_action, adjoint_action;
   ApplySurveyTangent(problem, frequencies, model, dot_direction,
                      tangent_action, settings);
   ApplySurveyAdjoint(problem, frequencies, model, dot_weights,
                      adjoint_action, settings);
   const real_t tangent_dot = tangent_action * dot_weights;
   const real_t adjoint_dot = dot_direction * adjoint_action;
   const real_t adjoint_relative_error = fabs(tangent_dot - adjoint_dot) /
      max(max(fabs(tangent_dot), fabs(adjoint_dot)), real_t(1.0e-30));

   if (problem.WorldRank() == 0)
   {
      cout << "Matrix-free tangent check: relative Frobenius error="
           << sqrt(squared_error /
                   max(squared_reference, numeric_limits<real_t>::min()))
           << ", maximum entrywise relative error=" << maximum_relative_error
           << " (details: " << file_name << ")\n"
           << "Matrix-free adjoint dot test: <Jv,w>=" << tangent_dot
           << ", <v,J^T w>=" << adjoint_dot
           << ", relative error=" << adjoint_relative_error << '\n';
   }
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
                       (problem.DataSize() / 2) * problem.SourceCount() *
                       frequencies.size(), observed,
                       standard_deviation);
   }

   Vector model(problem.InitialModel()), reference_model(model);
   Vector lower(problem.ParameterCount()), upper(problem.ParameterCount());
   for (int parameter = 0; parameter < problem.ParameterCount(); parameter++)
   {
      if (parameter % 3 == 2)
      {
         lower[parameter] = options.minimum_dip * M_PI / 180.0;
         upper[parameter] = options.maximum_dip * M_PI / 180.0;
      }
      else
      {
         lower[parameter] = log(options.minimum_resistivity);
         upper[parameter] = log(options.maximum_resistivity);
      }
   }

   ofstream history;
   if (problem.WorldRank() == 0)
   {
      history.open(string(options.output_prefix) + "_inversion.csv");
      MFEM_VERIFY(history, "Unable to open inversion history output.");
      history << "iteration,objective,data_objective,spatial_objective,"
                 "gradient_norm,step_norm,line_search\n";
      history << setprecision(16);
   }

   for (int iteration = 0; iteration < options.inversion_iterations;
        iteration++)
   {
      Vector prediction;
      SolveSurvey(problem, frequencies, model, prediction,
                  inversion_settings);
      MatrixFreeSurveyJacobian jacobian(problem, frequencies, model,
                                        inversion_settings);
      const real_t objective =
         Objective(prediction, observed, standard_deviation, model,
                   reference_model, problem.Regularization(),
                   options.regularization);
      const real_t data_objective =
         DataObjective(prediction, observed, standard_deviation);
      const real_t spatial_objective =
         SpatialObjective(model, reference_model, problem.Regularization(),
                          options.regularization);

      Vector weighted_residual(observed.Size());
      for (int datum = 0; datum < observed.Size(); datum++)
      {
         const real_t variance =
            standard_deviation[datum] * standard_deviation[datum];
         weighted_residual[datum] =
            (prediction[datum] - observed[datum]) / variance;
      }
      Vector gradient(problem.ParameterCount());
      jacobian.MultTranspose(weighted_residual, gradient);
      if (problem.Regularization().Height() > 0)
      {
         Vector difference(model), regularized(problem.Regularization().Height());
         Vector regularization_gradient(problem.ParameterCount());
         difference -= reference_model;
         problem.Regularization().Mult(difference, regularized);
         problem.Regularization().MultTranspose(regularized,
                                                regularization_gradient);
         gradient.Add(options.regularization, regularization_gradient);
      }
      const real_t gradient_norm = gradient.Norml2();

      GaussNewtonNormalOperator normal_operator(
         jacobian, standard_deviation,
         problem.Regularization(), options.regularization,
         options.gauss_newton_damping);
      Vector step(problem.ParameterCount());
      step = 0.0;
      Vector negative_gradient(gradient);
      negative_gradient *= -1.0;
      CGSolver least_squares_cg;
      least_squares_cg.SetOperator(normal_operator);
      // Matrix-free J/J^T products inherit the accuracy of the PDE solves.
      // Requiring a substantially smaller normal residual causes false
      // non-convergence once the action error reaches its natural floor.
      least_squares_cg.SetRelTol(options.gauss_newton_cg_tolerance);
      least_squares_cg.SetAbsTol(0.0);
      least_squares_cg.SetMaxIter(max(50, 2 * problem.ParameterCount()));
      least_squares_cg.SetPrintLevel(0);
      least_squares_cg.Mult(negative_gradient, step);
      MFEM_VERIFY(least_squares_cg.GetConverged(),
                  "The Gauss-Newton least-squares solve did not converge.");
      const real_t step_norm = step.Norml2();

      real_t line_search = 1.0;
      bool accepted = false;
      Vector trial_model(problem.ParameterCount()), trial_prediction;
      for (int trial = 0; trial < options.line_search_steps; trial++)
      {
         trial_model = model;
         trial_model.Add(line_search, step);
         for (int parameter = 0; parameter < problem.ParameterCount(); parameter++)
         {
            trial_model[parameter] =
               max(lower[parameter], min(upper[parameter],
                                          trial_model[parameter]));
         }
         SolveSurvey(problem, frequencies, trial_model, trial_prediction,
                     inversion_settings);
         const real_t trial_objective =
            Objective(trial_prediction, observed, standard_deviation,
                      trial_model, reference_model, problem.Regularization(),
                      options.regularization);
         if (trial_objective < objective)
         {
            accepted = true;
            break;
         }
         line_search *= 0.5;
      }

      if (problem.WorldRank() == 0)
      {
         cout << "GN " << iteration << ": Phi=" << objective
              << " (data=" << data_objective
              << ", spatial=" << spatial_objective << ')'
              << ", first Rh/Rv/dip=[" << exp(model[0]) << ", "
              << exp(model[1]) << ", " << model[2] * 180.0 / M_PI
              << " deg], |g|=" << gradient_norm << ", |dm|=" << step_norm
              << ", alpha=" << (accepted ? line_search : 0.0) << '\n';
         history << iteration << ',' << objective << ',' << data_objective
                 << ',' << spatial_objective << ',' << gradient_norm << ','
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
   SolveSurvey(problem, frequencies, model, final_prediction,
               forward_settings,
               string(options.output_prefix) + "_field");
   WriteData(string(options.output_prefix) + "_recovered.csv", problem,
             frequencies, final_prediction, standard_deviation);
   WriteModel(string(options.output_prefix) + "_recovered_model.csv",
              problem, model);
   if (problem.WorldRank() == 0)
   {
      cout << "Recovered anisotropic model with "
           << problem.ParameterCount() / 3 << " unit(s).\n";
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
   args.AddOption(&options.material_file, "-mat", "--material-table",
                  "Volume-attribute material CSV for tagged/topographic "
                  "meshes.");
   args.AddOption(&options.boundary_condition, "-bc",
                  "--boundary-condition",
                  "Far-boundary treatment: pec or absorbing.");
   args.AddOption(&options.absorbing_boundary_attributes, "-ab",
                  "--absorbing-boundaries",
                  "Comma-separated absorbing boundary attributes; empty "
                  "selects all.");
   args.AddOption(&options.absorbing_conductivity, "-abcs",
                  "--absorbing-conductivity",
                  "Exterior conductivity for the absorbing condition.");
   args.AddOption(&options.absorbing_relative_permittivity, "-abce",
                  "--absorbing-relative-permittivity",
                  "Exterior relative permittivity for the absorber.");
   args.AddOption(&options.absorbing_relative_permeability, "-abcm",
                  "--absorbing-relative-permeability",
                  "Exterior relative permeability for the absorber.");
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
   args.AddOption(&options.model_parameterization, "-mp",
                  "--model-parameterization", "Model: two, region, or cell.");
   args.AddOption(&options.model_file, "-mf", "--model-file",
                  "True-model CSV: id,Rh,Rv,dip_degrees.");
   args.AddOption(&options.initial_model_file, "-imf",
                  "--initial-model-file",
                  "Initial-model CSV: id,Rh,Rv,dip_degrees.");
   args.AddOption(&options.background_vertical_resistivity, "-brv",
                  "--background-rv", "Background vertical resistivity.");
   args.AddOption(&options.target_vertical_resistivity, "-trv",
                  "--target-rv", "Target vertical resistivity.");
   args.AddOption(&options.background_dip, "-bdip", "--background-dip",
                  "Background symmetry-axis dip in degrees.");
   args.AddOption(&options.target_dip, "-tdip", "--target-dip",
                  "Target symmetry-axis dip in degrees.");
   args.AddOption(&options.minimum_resistivity, "-rmin",
                  "--minimum-resistivity", "Inversion Rh/Rv lower bound.");
   args.AddOption(&options.maximum_resistivity, "-rmax",
                  "--maximum-resistivity", "Inversion Rh/Rv upper bound.");
   args.AddOption(&options.minimum_dip, "-dmin", "--minimum-dip",
                  "Inversion dip lower bound in degrees.");
   args.AddOption(&options.maximum_dip, "-dmax", "--maximum-dip",
                  "Inversion dip upper bound in degrees.");
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
   args.AddOption(&options.receiver_model, "-rm", "--receiver-model",
                  "Receiver evaluation: gaussian or point.");
   args.AddOption(&options.source_x, "-sx", "--source-x",
                  "Wire center x coordinate in metres.");
   args.AddOption(&options.source_file, "-sf", "--source-file",
                  "Transmitter polyline CSV with x,y,z vertex rows.");
   args.AddOption(&options.source_catalog, "-sc", "--source-catalog",
                  "Straight-transmitter station catalog CSV.");
   args.AddOption(&options.source_model, "-sm", "--source-model",
                  "Source discretization: gaussian or line.");
   args.AddOption(&options.source_quadrature_order, "-sqo",
                  "--source-quadrature-order",
                  "Integration order for the line-delta source.");
   args.AddOption(&options.survey_groups, "-sg", "--survey-groups",
                  "MPI groups used for independent frequency/source tasks.");
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
                  "Spatial graph-regularization weight.");
   args.AddOption(&options.regularization_neighbors, "-rkn",
                  "--regularization-neighbors",
                  "Nearest model-unit neighbors in the spatial graph.");
   args.AddOption(&options.regularization_rh_scale, "-rrhs",
                  "--regularization-rh-scale", "Scale for log(Rh) edges.");
   args.AddOption(&options.regularization_rv_scale, "-rrvs",
                  "--regularization-rv-scale", "Scale for log(Rv) edges.");
   args.AddOption(&options.regularization_dip_scale, "-rds",
                  "--regularization-dip-scale", "Scale for dip edges.");
   args.AddOption(&options.gauss_newton_damping, "-gnd",
                  "--gauss-newton-damping",
                  "Diagonal damping added to the normal operator.");
   args.AddOption(&options.gauss_newton_cg_tolerance, "-gnct",
                  "--gauss-newton-cg-tolerance",
                  "Relative tolerance for the matrix-free normal solve.");
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

   if (Mpi::Root())
   {
      const filesystem::path output_path(options.output_prefix);
      if (output_path.has_parent_path())
      {
         error_code error;
         filesystem::create_directories(output_path.parent_path(), error);
         MFEM_VERIFY(!error, "Unable to create output directory: "
                     << output_path.parent_path().string() << " ("
                     << error.message() << ')');
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);

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
   MFEM_VERIFY(string(options.receiver_model) == "gaussian" ||
                  string(options.receiver_model) == "point",
               "Receiver model must be gaussian or point.");
   MFEM_VERIFY(string(options.receiver_model) != "point" ||
                  (!options.inversion && !options.verify_jacobian),
               "Point receivers are available for forward/reference runs; "
               "inversion and adjoint Jacobian checks require Gaussian "
               "receivers.");
   MFEM_VERIFY(string(options.source_file).empty() ||
                  string(options.source_catalog).empty(),
               "Use either --source-file or --source-catalog, not both.");
   MFEM_VERIFY(string(options.source_model) == "gaussian" ||
                  string(options.source_model) == "line",
               "Source model must be gaussian or line.");
   MFEM_VERIFY(options.source_quadrature_order >= 1,
               "Source quadrature order must be positive.");
   int world_rank = 0, world_size = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   MFEM_VERIFY(options.survey_groups >= 1 &&
                  options.survey_groups <= world_size,
               "Survey groups must be between one and the MPI world size.");
   MFEM_VERIFY(string(options.model_parameterization) != "cell" ||
                  world_size % options.survey_groups == 0,
               "Cell models require equal-size survey groups so global cell "
               "numbering remains identical in every group.");
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
   MFEM_VERIFY(string(options.model_parameterization) == "two" ||
                   string(options.model_parameterization) == "region" ||
                   string(options.model_parameterization) == "cell",
                "Model parameterization must be two, region, or cell.");
   MFEM_VERIFY(string(options.material_file).empty() ||
                  string(options.model_parameterization) == "region" ||
                  string(options.model_parameterization) == "cell",
               "Tagged material assignment requires region or cell model "
               "parameterization.");
   MFEM_VERIFY(string(options.boundary_condition) == "pec" ||
                  string(options.boundary_condition) == "absorbing",
               "Boundary condition must be pec or absorbing.");
   MFEM_VERIFY(options.absorbing_conductivity >= 0.0 &&
                  options.absorbing_relative_permittivity > 0.0 &&
                  options.absorbing_relative_permeability > 0.0,
               "Absorbing exterior conductivity must be nonnegative, and "
               "relative permittivity/permeability must be positive.");
   MFEM_VERIFY(options.minimum_resistivity > 0.0 &&
                  options.maximum_resistivity > options.minimum_resistivity,
               "Resistivity bounds must be positive and ordered.");
   MFEM_VERIFY(options.maximum_dip > options.minimum_dip,
               "Dip bounds must be ordered.");
   MFEM_VERIFY(options.regularization_neighbors >= 0 &&
                   options.regularization_rh_scale >= 0.0 &&
                   options.regularization_rv_scale >= 0.0 &&
                  options.regularization_dip_scale >= 0.0,
                "Regularization neighbors and component scales must be "
                "nonnegative.");
   MFEM_VERIFY(options.gauss_newton_damping >= 0.0 &&
                  options.gauss_newton_cg_tolerance > 0.0,
               "Gauss-Newton damping must be nonnegative and its CG "
               "tolerance must be positive.");

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
   const int survey_group =
      world_rank * options.survey_groups / world_size;
   MPI_Comm survey_communicator;
   MPI_Comm_split(MPI_COMM_WORLD, survey_group, world_rank,
                  &survey_communicator);
   if (Mpi::Root())
   {
      cout << "Survey scheduling: " << options.survey_groups
           << " MPI group(s) for "
           << frequencies.size() * LoadTransmitters(options).size()
           << " frequency/source task(s)\n";
   }
   CSEMProblem problem(survey_communicator, options, survey_group,
                       options.survey_groups);
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
      Vector model(problem.TrueModel()), responses;
      SolveSurvey(problem, frequencies, model, responses,
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
      WriteModel(string(options.output_prefix) + "_model.csv", problem,
                 model);
      if (options.verify_jacobian)
      {
         VerifyJacobian(problem, options, frequencies, model,
                        forward_settings);
      }
   }
   return 0;
}
