//                       MFEM UDAR 3D Forward Solver
//
// Compile with: make udar3p
//
// Sample runs:
//   mpirun -np 4 udar3p -m ../../data/inline-hex.mesh -rs 1 -rp 1
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
#include "udar3p_adapt.hpp"
#include "udar3p_frequency.hpp"
#include "udar3p_pml.hpp"
#include "udar3p_postprocess.hpp"
#include "udar3p_points.hpp"
#include "udar3p_solver.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
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
{ return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

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

/** @brief Parse a comma-separated list of nonnegative station indices.
    @param text Comma-separated integer indices, or an empty string.
    @return Sorted unique station indices; empty when @a text is empty.

    The caller checks indices against the configured survey station count. */
vector<int> ParseStationIndices(const char *text)
{
   const vector<real_t> values = ParseList(text);
   vector<int> indices;
   for (real_t value : values)
   {
      MFEM_VERIFY(isfinite(value) && value >= 0.0 && value == floor(value) &&
                      value <= numeric_limits<int>::max(),
                  "Station indices must be nonnegative integers.");
      indices.push_back(static_cast<int>(value));
   }
   sort(indices.begin(), indices.end());
   const auto duplicate = adjacent_find(indices.begin(), indices.end());
   MFEM_VERIFY(duplicate == indices.end(), "Station indices must be unique.");
   return indices;
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
         MFEM_VERIFY(fabs(tensor(i, j) - tensor(j, i)) <= symmetry_tolerance,
                     "Conductivity tensor must be symmetric.");
         const real_t symmetric_value = 0.5 * (tensor(i, j) + tensor(j, i));
         tensor(i, j) = tensor(j, i) = symmetric_value;
      }
   }

   DenseMatrix normalized(tensor);
   normalized *= 1.0 / tensor_scale;
   const real_t leading_minor_1 = normalized(0, 0);
   const real_t leading_minor_2 = normalized(0, 0) * normalized(1, 1) -
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

void BuildToolFrame(const Vector &start, const Vector &end, const Vector &up,
                    DenseMatrix &frame)
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
                             const Vector &slab_center_, real_t slab_thickness_,
                             real_t interface_width_,
                             const Vector &box_minimum_,
                             const Vector &box_maximum_)
       : MatrixCoefficient(3), model(model_), background(background_),
         target(target_), slab_normal(slab_normal_), slab_center(slab_center_),
         slab_thickness(slab_thickness_), interface_width(interface_width_),
         box_minimum(box_minimum_), box_maximum(box_maximum_)
   { Normalize3(slab_normal, "Slab normal"); }

   /** @brief Replace the tensors returned by Eval().
       @param background_ Symmetric positive-definite background tensor.
       @param target_ Symmetric positive-definite target tensor.

       The persistent server calls this between solves. The coefficient does
       not retain references to either input matrix. */
   void SetConductivities(const DenseMatrix &background_,
                          const DenseMatrix &target_)
   {
      background = background_;
      target = target_;
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
            target_weight =
                0.5 * (tanh((half_thickness + distance) / interface_width) +
                       tanh((half_thickness - distance) / interface_width));
         }
         else
         {
            target_weight = fabs(distance) <= 0.5 * slab_thickness ? 1.0 : 0.0;
         }
      }
      else if (model == 2)
      {
         target_weight = 1.0;
         if (interface_width > 0.0)
         {
            for (int d = 0; d < 3; d++)
            {
               target_weight *=
                   0.5 * (tanh((x[d] - box_minimum[d]) / interface_width) -
                          tanh((x[d] - box_maximum[d]) / interface_width));
            }
         }
         else
         {
            for (int d = 0; d < 3; d++)
            {
               target_weight *=
                   x[d] >= box_minimum[d] && x[d] <= box_maximum[d] ? 1.0 : 0.0;
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

/** @brief Scalar horizontal or vertical resistivity derived from conductivity.

    The horizontal value is the mean of the xx and yy entries of the inverse
    conductivity tensor; the vertical value is its zz entry. This convention
    gives Rh = 1/sigma_h and Rv = 1/sigma_v for a transversely isotropic model
    aligned with the global z axis. */
class ResistivityCoefficient : public Coefficient
{
 private:
   MatrixCoefficient &conductivity;
   bool horizontal;

 public:
   ResistivityCoefficient(MatrixCoefficient &conductivity_, bool horizontal_)
       : conductivity(conductivity_), horizontal(horizontal_)
   {
   }

   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
   {
      DenseMatrix resistivity;
      conductivity.Eval(resistivity, T, ip);
      resistivity.Invert();
      return horizontal ? 0.5 * (resistivity(0, 0) + resistivity(1, 1))
                        : resistivity(2, 2);
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
      const real_t current =
          moment_magnitude / (M_PI * loop_radius * loop_radius);
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
   const int candidate_owner =
       (local_element >= 0) ? Mpi::WorldRank() : no_owner;
   int owner = no_owner;
   MPI_Allreduce(&candidate_owner, &owner, 1, MPI_INT, MPI_MIN, mesh.GetComm());
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
         for (int d = 0; d < 3; d++) { value += moment[d] * curl_shape(i, d); }
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
      local_found = (T->TransformBack(point, ip, tolerance) ==
                     InverseElementTransformation::Inside);
   }
   int global_found = 0;
   MPI_Allreduce(&local_found, &global_found, 1, MPI_INT, MPI_MAX,
                 mesh.GetComm());
   return global_found != 0;
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

/** @brief Write a legacy VTK mesh with scalar cell fields for one evaluation.
    @param output Destination legacy VTK file.
    @param mesh Parallel mesh local to this MPI rank.
    @param rh Horizontal resistivity coefficient.
    @param rv Vertical resistivity coefficient.

    The model is deliberately written as legacy VTK rather than through
    ParaViewDataCollection because coefficient fields in that collection are
    always PointData. Rh and Rv are material properties and must be CellData.
    Each MPI rank writes its local mesh partition. */
void WriteModelCellDataVTK(const string &output, ParMesh &mesh, Coefficient &rh,
                           Coefficient &rv, const UDARCartesianPML *pml)
{
   ofstream stream(output);
   MFEM_VERIFY(stream.is_open(), "Unable to open VTK output " << output);
   stream.precision(16);

   const int refinement = 1;
   int point_count = 0;
   int cell_count = 0;
   int cell_size = 0;
   for (int element = 0; element < mesh.GetNE(); element++)
   {
      const Geometry::Type geometry = mesh.GetElementBaseGeometry(element);
      const int vertices = Geometries.GetVertices(geometry)->GetNPoints();
      RefinedGeometry *refined =
          GlobGeometryRefiner.Refine(geometry, refinement, 1);
      point_count += refined->RefPts.GetNPoints();
      const int cells = refined->RefGeoms.Size() / vertices;
      cell_count += cells;
      cell_size += cells * (vertices + 1);
   }

   stream << "# vtk DataFile Version 3.0\n"
          << "MFEM UDAR formation resistivity cell data\n"
          << "ASCII\n"
          << "DATASET UNSTRUCTURED_GRID\n"
          << "POINTS " << point_count << " double\n";
   for (int element = 0; element < mesh.GetNE(); element++)
   {
      const Geometry::Type geometry = mesh.GetElementBaseGeometry(element);
      RefinedGeometry *refined =
          GlobGeometryRefiner.Refine(geometry, refinement, 1);
      DenseMatrix points;
      mesh.GetElementTransformation(element)->Transform(refined->RefPts,
                                                        points);
      for (int point = 0; point < points.Width(); point++)
      {
         stream << points(0, point) << ' '
                << (points.Height() > 1 ? points(1, point) : 0.0) << ' '
                << (points.Height() > 2 ? points(2, point) : 0.0) << '\n';
      }
   }

   stream << "CELLS " << cell_count << ' ' << cell_size << '\n';
   int point_offset = 0;
   for (int element = 0; element < mesh.GetNE(); element++)
   {
      const Geometry::Type geometry = mesh.GetElementBaseGeometry(element);
      const int vertices = Geometries.GetVertices(geometry)->GetNPoints();
      RefinedGeometry *refined =
          GlobGeometryRefiner.Refine(geometry, refinement, 1);
      for (int point = 0; point < refined->RefGeoms.Size();)
      {
         stream << vertices;
         for (int vertex = 0; vertex < vertices; vertex++, point++)
         {
            stream << ' ' << point_offset + refined->RefGeoms[point];
         }
         stream << '\n';
      }
      point_offset += refined->RefPts.GetNPoints();
   }

   stream << "CELL_TYPES " << cell_count << '\n';
   for (int element = 0; element < mesh.GetNE(); element++)
   {
      const Geometry::Type geometry = mesh.GetElementBaseGeometry(element);
      const int vertices = Geometries.GetVertices(geometry)->GetNPoints();
      RefinedGeometry *refined =
          GlobGeometryRefiner.Refine(geometry, refinement, 1);
      for (int point = 0; point < refined->RefGeoms.Size(); point += vertices)
      {
         stream << VTKGeometry::Map[geometry] << '\n';
      }
   }

   stream << "CELL_DATA " << cell_count << '\n'
          << "SCALARS Rh_ohm_m double 1\n"
          << "LOOKUP_TABLE default\n";
   for (int element = 0; element < mesh.GetNE(); element++)
   {
      const Geometry::Type geometry = mesh.GetElementBaseGeometry(element);
      const int vertices = Geometries.GetVertices(geometry)->GetNPoints();
      RefinedGeometry *refined =
          GlobGeometryRefiner.Refine(geometry, refinement, 1);
      ElementTransformation *transformation =
          mesh.GetElementTransformation(element);
      const IntegrationPoint &center = Geometries.GetCenter(geometry);
      const real_t value = rh.Eval(*transformation, center);
      for (int point = 0; point < refined->RefGeoms.Size(); point += vertices)
      {
         stream << value << '\n';
      }
   }

   stream << "SCALARS Rv_ohm_m double 1\n"
          << "LOOKUP_TABLE default\n";
   for (int element = 0; element < mesh.GetNE(); element++)
   {
      const Geometry::Type geometry = mesh.GetElementBaseGeometry(element);
      const int vertices = Geometries.GetVertices(geometry)->GetNPoints();
      RefinedGeometry *refined =
          GlobGeometryRefiner.Refine(geometry, refinement, 1);
      ElementTransformation *transformation =
          mesh.GetElementTransformation(element);
      const IntegrationPoint &center = Geometries.GetCenter(geometry);
      const real_t value = rv.Eval(*transformation, center);
      for (int point = 0; point < refined->RefGeoms.Size(); point += vertices)
      {
         stream << value << '\n';
      }
   }

   stream << "SCALARS PML_region int 1\n"
          << "LOOKUP_TABLE default\n";
   for (int element = 0; element < mesh.GetNE(); element++)
   {
      const Geometry::Type geometry = mesh.GetElementBaseGeometry(element);
      const int vertices = Geometries.GetVertices(geometry)->GetNPoints();
      RefinedGeometry *refined =
          GlobGeometryRefiner.Refine(geometry, refinement, 1);
      Vector center;
      mesh.GetElementCenter(element, center);
      const int value = pml && pml->Contains(center) ? 1 : 0;
      for (int point = 0; point < refined->RefGeoms.Size(); point += vertices)
      {
         stream << value << '\n';
      }
   }

   stream << "SCALARS PML_profile double 1\n"
          << "LOOKUP_TABLE default\n";
   for (int element = 0; element < mesh.GetNE(); element++)
   {
      const Geometry::Type geometry = mesh.GetElementBaseGeometry(element);
      const int vertices = Geometries.GetVertices(geometry)->GetNPoints();
      RefinedGeometry *refined =
          GlobGeometryRefiner.Refine(geometry, refinement, 1);
      Vector center;
      mesh.GetElementCenter(element, center);
      const real_t value = pml ? pml->Profile(center) : 0.0;
      for (int point = 0; point < refined->RefGeoms.Size(); point += vertices)
      {
         stream << value << '\n';
      }
   }

   stream << "SCALARS PML_sigma_S_m double 1\n"
          << "LOOKUP_TABLE default\n";
   for (int element = 0; element < mesh.GetNE(); element++)
   {
      const Geometry::Type geometry = mesh.GetElementBaseGeometry(element);
      const int vertices = Geometries.GetVertices(geometry)->GetNPoints();
      RefinedGeometry *refined =
          GlobGeometryRefiner.Refine(geometry, refinement, 1);
      Vector center;
      mesh.GetElementCenter(element, center);
      const real_t value = pml ? pml->ArtificialConductivity(center) : 0.0;
      for (int point = 0; point < refined->RefGeoms.Size(); point += vertices)
      {
         stream << value << '\n';
      }
   }

   stream << "SCALARS cell_size double 1\n"
          << "LOOKUP_TABLE default\n";
   for (int element = 0; element < mesh.GetNE(); element++)
   {
      const Geometry::Type geometry = mesh.GetElementBaseGeometry(element);
      const int vertices = Geometries.GetVertices(geometry)->GetNPoints();
      RefinedGeometry *refined =
          GlobGeometryRefiner.Refine(geometry, refinement, 1);
      const real_t value = mesh.GetElementSize(element, 1);
      for (int point = 0; point < refined->RefGeoms.Size(); point += vertices)
      {
         stream << value << '\n';
      }
   }
}

/** @brief Write the formation resistivity fields for one forward evaluation.
    @param directory Root directory selected with @c -vtk.
    @param evaluation_name Stem of the CSV response file.
    @param mesh Parallel mesh containing the formation.
    @param sigma Physical formation conductivity tensor.
    @param pml PML geometry and grading, or nullptr when PML is disabled.

    The produced @c model.pvd references one legacy VTK file per MPI rank.
    Rh and Rv are expressed in ohm m, follow the global x/y horizontal and z
    vertical convention documented by ResistivityCoefficient, and are stored
    as VTK cell data. The PML mask, grading, artificial stretch conductivity,
    and cell size are also cell data for inspection in ParaView. */
void WriteVTKModel(const string &directory, const string &evaluation_name,
                   ParMesh &mesh, MatrixCoefficient &sigma,
                   const UDARCartesianPML *pml)
{
   if (directory.empty()) { return; }

   ResistivityCoefficient rh(sigma, true);
   ResistivityCoefficient rv(sigma, false);
   const filesystem::path model_directory =
       filesystem::path(directory) / evaluation_name / "model";
   filesystem::create_directories(model_directory);
   ostringstream filename;
   filename << "model_rank" << setw(6) << setfill('0') << Mpi::WorldRank()
            << ".vtk";
   WriteModelCellDataVTK((model_directory / filename.str()).string(), mesh, rh,
                         rv, pml);

   if (Mpi::WorldRank() == 0)
   {
      ofstream collection(model_directory / "model.pvd");
      MFEM_VERIFY(collection.is_open(), "Unable to open model PVD output.");
      collection << "<?xml version=\"1.0\"?>\n"
                 << "<VTKFile type=\"Collection\" version=\"1.0\" "
                 << "byte_order=\"LittleEndian\">\n<Collection>\n";
      for (int rank = 0; rank < Mpi::WorldSize(); rank++)
      {
         ostringstream piece;
         piece << "model_rank" << setw(6) << setfill('0') << rank << ".vtk";
         collection << "<DataSet timestep=\"0\" group=\"\" part=\"" << rank
                    << "\" file=\"" << piece.str() << "\"/>\n";
      }
      collection << "</Collection>\n</VTKFile>\n";
   }
}


} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   const int myid = Mpi::WorldRank();

   const char *mesh_file = "../../data/inline-hex.mesh";
   const char *frequency_text = "2000,12000,48000";
   const char *spacing_text = "";
   const char *station_indices_text = "";
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
   const char *server_directory = "";
   const char *vtk_directory = "";
   const char *device_config = "cpu";
   const char *iterative_preconditioner = "ilu";
   const char *ams_system_mode = "complex";
   int order = 1;
   int serial_refinements = 0;
   int parallel_refinements = 0;
   int stations = 7;
   int formation_model = 2;
   int boundary_condition = 1;
   int source_type = 0;
   int max_iterations = 2000;
   int krylov_dimension = 200;
   int octree_levels = 2;
   int pml_order = 2;
   int ams_cycle_type = 14;
   int ams_iterations = 2;
   int ams_smoothing_iterations = 1;
   int ams_print_level = 0;
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
   real_t octree_radius = -1.0;
   real_t octree_padding = 0.0;
   real_t pml_thickness = -1.0;
   real_t pml_conductivity = -1.0;
   bool use_mumps = false;
   bool reuse_mumps_ordering = true;
   bool use_octree = false;
   bool use_pml = false;
   bool use_ilu_fallback = true;
   const char *assembly_mode = "auto";
   bool profile_vtk = false;

   OptionsParser args(argc, argv);
   args.AddOption(&profile_vtk, "-vtk-profile", "--vtk-profile",
                  "-no-vtk-profile", "--no-vtk-profile",
                  "Measure host derived-field evaluation inside ParaView Save.");
   args.AddOption(&assembly_mode, "-assembly", "--assembly",
                  "Maxwell assembly: auto, full, or guarded partial.");
   args.AddOption(&mesh_file, "-m", "--mesh", "3D mesh file.");
   args.AddOption(&order, "-o", "--order", "Nedelec finite element order.");
   args.AddOption(&serial_refinements, "-rs", "--refinements-serial",
                  "Number of serial uniform refinements.");
   args.AddOption(&parallel_refinements, "-rp", "--refinements-parallel",
                  "Number of parallel uniform refinements.");
   args.AddOption(
       &use_octree, "-octree", "--octree", "-no-octree", "--no-octree",
       "Enable graded nonconforming refinement around the tool path.");
   args.AddOption(&octree_levels, "-octree-l", "--octree-levels",
                  "Number of locally refined octree levels near the tool.");
   args.AddOption(&octree_radius, "-octree-r", "--octree-radius",
                  "Radius of the finest transmitter-receiver refinement tube.");
   args.AddOption(
       &octree_padding, "-octree-pad", "--octree-padding",
       "Extra radial and axial refinement coverage around the tool.");
   args.AddOption(&use_pml, "-pml", "--pml", "-no-pml", "--no-pml",
                  "Enable a Cartesian coordinate-stretch PML at outer faces.");
   args.AddOption(
       &pml_thickness, "-pml-t", "--pml-thickness",
       "PML thickness in mesh units; default is 15 percent of span.");
   args.AddOption(&pml_conductivity, "-pml-s", "--pml-conductivity",
                  "PML maximum stretch conductivity in S/m.");
   args.AddOption(&pml_order, "-pml-o", "--pml-order",
                  "Polynomial grading order for the PML conductivity.");
   args.AddOption(&frequency_text, "-freq", "--frequencies",
                  "Comma-separated frequencies in Hz.");
   args.AddOption(
       &spacing_text, "-sp", "--spacings",
       "Comma-separated transmitter-receiver spacings in mesh units.");
   args.AddOption(&stations, "-ns", "--num-stations",
                  "Number of transmitter stations on the straight survey.");
   args.AddOption(
       &station_indices_text, "-si", "--station-indices",
       "Comma-separated survey stations to simulate; empty selects all.");
   args.AddOption(&survey_start_text, "-p0", "--survey-start",
                  "Survey start x,y,z; default is based on the mesh bounds.");
   args.AddOption(&survey_end_text, "-p1", "--survey-end",
                  "Survey end x,y,z; default is based on the mesh bounds.");
   args.AddOption(&measured_depth_origin, "-md0", "--measured-depth-origin",
                  "Measured depth assigned to the first station.");
   args.AddOption(
       &up_text, "-up", "--tool-up",
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
   args.AddOption(
       &sigma_background_text, "-sig-bg", "--sigma-background",
       "Background conductivity: 1, 3 diagonal, or 9 tensor values [S/m].");
   args.AddOption(
       &sigma_target_text, "-sig-target", "--sigma-target",
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
                  "--iterative-solver",
                  "Use MUMPS or iterative FGMRES with ILU/AMS.");
   args.AddOption(&reuse_mumps_ordering, "-mumps-reuse",
                  "--mumps-reuse-ordering", "-no-mumps-reuse",
                  "--no-mumps-reuse-ordering",
                  "Reuse MUMPS symbolic ordering for a fixed mesh.");
   args.AddOption(&iterative_preconditioner, "-pc", "--preconditioner",
                  "Iterative preconditioner: ilu or ams.");
   args.AddOption(&ams_system_mode, "-ams-mode", "--ams-system-mode",
                  "AMS system representation: complex or real.");
   args.AddOption(&ams_cycle_type, "-ams-cycle", "--ams-cycle-type",
                  "Hypre AMS cycle type; Palace-style default is 14.");
   args.AddOption(&ams_iterations, "-ams-it", "--ams-iterations",
                  "AMS cycles applied per preconditioner call.");
   args.AddOption(&ams_smoothing_iterations, "-ams-smooth",
                  "--ams-smoothing-iterations",
                  "AMS smoothing sweeps per cycle.");
   args.AddOption(&ams_print_level, "-ams-print", "--ams-print-level",
                  "Hypre AMS diagnostic print level.");
   args.AddOption(&use_ilu_fallback, "-ams-fallback", "--ams-ilu-fallback",
                  "-no-ams-fallback", "--no-ams-ilu-fallback",
                  "Retry a failed AMS solve with distributed ILU.");
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
   args.AddOption(
       &vtk_directory, "-vtk", "--vtk-directory",
       "ParaView VTK output directory; empty disables field output.");
   // Server mode receives rank-zero requests and writes one CSV per request.
   args.AddOption(
       &server_directory, "-server", "--server-directory",
       "Persistent inversion server output directory; empty disables it.");
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
   MFEM_VERIFY(pml_order >= 1, "PML polynomial order must be at least one.");
   MFEM_VERIFY(transmitter_moment != 0.0,
               "Transmitter moment must be nonzero.");
   MFEM_VERIFY(string(iterative_preconditioner) == "ilu" ||
                   string(iterative_preconditioner) == "ams",
               "Iterative preconditioner must be 'ilu' or 'ams'.");
   MFEM_VERIFY(string(ams_system_mode) == "complex" ||
                   string(ams_system_mode) == "real",
               "AMS system mode must be 'complex' or 'real'.");
   MFEM_VERIFY(ams_cycle_type >= 1 && ams_cycle_type <= 14,
               "AMS cycle type must be between 1 and 14.");
   MFEM_VERIFY(ams_iterations >= 1 && ams_smoothing_iterations >= 1,
               "AMS iteration counts must be positive.");

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

   Vector bounds_minimum(3), bounds_maximum(3);
   serial_mesh.GetBoundingBox(bounds_minimum, bounds_maximum);
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
   vector<int> station_indices = ParseStationIndices(station_indices_text);
   if (station_indices.empty())
   {
      station_indices.resize(stations);
      for (int station = 0; station < stations; station++)
      {
         station_indices[station] = station;
      }
   }
   for (int station : station_indices)
   {
      MFEM_VERIFY(station < stations, "Station index " << station
                                                       << " is outside [0, "
                                                       << stations - 1 << "].");
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
      survey_start[1] = survey_end[1] = bounds_minimum[1] + 0.47 * span[1];
      survey_start[2] = survey_end[2] = bounds_minimum[2] + 0.53 * span[2];
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
   if (use_octree)
   {
      MFEM_VERIFY(octree_levels >= 1,
                  "Octree refinement needs at least one level.");
      if (octree_radius <= 0.0)
      {
         octree_radius =
             max(real_t(0.04) * min_span,
                 real_t(0.5) * *max_element(spacings.begin(), spacings.end()));
      }
      MFEM_VERIFY(octree_radius > 0.0 && octree_padding >= 0.0,
                  "Octree radius must be positive and padding nonnegative.");
      serial_mesh.EnsureNCMesh(true);
   }

   ND_FECollection finite_elements(order, 3);
   auto *coarse_mesh = new ParMesh(MPI_COMM_WORLD, serial_mesh);
   serial_mesh.Clear();
   auto *coarse_space = new ParFiniteElementSpace(coarse_mesh, &finite_elements);
   ParFiniteElementSpaceHierarchy hierarchy(coarse_mesh, coarse_space,
                                            true, true);
   if (use_octree)
   {
      const int octree_passes =
          RefineUDARToolRegion(hierarchy, survey_start, survey_end, tool_axis,
                               *max_element(spacings.begin(), spacings.end()),
                               octree_radius, octree_padding, octree_levels);
      if (myid == 0)
      {
         cout << "Applied " << octree_passes
              << " graded nested refinement passes around the tool path.\n";
      }
   }

   for (int level = 0; level < parallel_refinements; level++)
   {
      hierarchy.AddUniformlyRefinedLevel();
   }
   ParFiniteElementSpace &space = hierarchy.GetFinestFESpace();
   ParMesh &mesh = *space.GetParMesh();
   for (int level = 0; level < hierarchy.GetNumLevels(); level++)
   {
      const auto unknowns = hierarchy.GetFESpaceAtLevel(level).GlobalTrueVSize();
      if (myid == 0)
      {
         cout << "H(curl) hierarchy level " << level << ": " << unknowns
              << " true unknowns\n";
      }
   }
   if (myid == 0) { cout << "UDAR unknowns: "; }
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
           << (source_type == 0 ? "point magnetic dipole"
                                : "finite current loop")
           << endl;
   }

   FormationSigmaCoefficient sigma(
       formation_model, sigma_background, sigma_target, slab_normal,
       slab_center, slab_thickness, interface_width, box_minimum, box_maximum);
   if (use_pml)
   {
      if (pml_thickness <= 0.0) { pml_thickness = 0.15 * min_span; }
      if (pml_conductivity <= 0.0)
      {
         pml_conductivity = 10.0 * sigma_background(0, 0);
      }
      MFEM_VERIFY(pml_thickness > 0.0 && pml_conductivity > 0.0,
                  "PML thickness and conductivity must be positive.");
   }
   unique_ptr<UDARCartesianPML> pml;
   if (use_pml)
   {
      pml.reset(new UDARCartesianPML(bounds_minimum, bounds_maximum,
                                     pml_thickness, pml_conductivity,
                                     pml_order));
      Vector first_receiver(survey_start);
      Vector last_receiver(survey_end);
      const real_t maximum_spacing =
          *max_element(spacings.begin(), spacings.end());
      first_receiver.Add(-maximum_spacing, tool_axis);
      last_receiver.Add(-maximum_spacing, tool_axis);
      MFEM_VERIFY(!pml->Contains(survey_start) && !pml->Contains(survey_end) &&
                      !pml->Contains(first_receiver) &&
                      !pml->Contains(last_receiver),
                  "The transmitter-receiver trajectory intersects the PML; "
                  "enlarge the mesh or reduce the PML thickness.");
   }

   if (vtk_directory[0] && myid == 0)
   {
      filesystem::create_directories(vtk_directory);
   }
   MPI_Barrier(mesh.GetComm());

   if (pml)
   {
      const vector<UDARPMLSideDiagnostics> pml_diagnostics =
          AnalyzeUDARPMLResolution(mesh, *pml);
      if (myid == 0)
      {
         filesystem::path diagnostics_path;
         if (vtk_directory[0])
         {
            diagnostics_path =
                filesystem::path(vtk_directory) / "pml_resolution.csv";
         }
         else if (server_directory[0])
         {
            diagnostics_path =
                filesystem::path(server_directory) / "pml_resolution.csv";
         }
         else
         {
            const filesystem::path response_path(output_file);
            diagnostics_path =
                response_path.parent_path() /
                (response_path.stem().string() + "_pml_resolution.csv");
         }

         cout << "PML mesh-resolution diagnostics:\n";
         for (const UDARPMLSideDiagnostics &entry : pml_diagnostics)
         {
            cout << "  " << entry.SideName() << ": "
                 << entry.conservative_layers << " conservative layer(s), "
                 << entry.intersected_elements << " intersected element(s), "
                 << "max profile jump " << entry.maximum_profile_jump
                 << ", max normal span " << entry.maximum_normal_span
                 << '\n';
            if (entry.conservative_layers < 4)
            {
               cout << "  Warning: " << entry.SideName()
                    << " PML has fewer than four conservative element layers; "
                    << "refine the outer mesh or increase PML thickness.\n";
            }
            if (entry.maximum_profile_jump > 0.35)
            {
               cout << "  Warning: " << entry.SideName()
                    << " PML grading changes by more than 0.35 in one element; "
                    << "refine the absorbing layer.\n";
            }
         }
         WriteUDARPMLResolutionCSV(diagnostics_path.string(), pml_diagnostics,
                                   *pml);
         cout << "Wrote PML diagnostics to " << diagnostics_path.string()
              << '\n';
      }
   }

   Array<int> essential_tdofs;
   if (boundary_condition == 0 && mesh.bdr_attributes.Size())
   {
      Array<int> essential_boundary(mesh.bdr_attributes.Max());
      essential_boundary = 1;
      space.GetEssentialTrueDofs(essential_boundary, essential_tdofs);
   }

   Vector trajectory(survey_end);
   trajectory -= survey_start;
   const real_t survey_length = trajectory.Norml2();
   const int true_size = space.GetTrueVSize();
   const ComplexOperator::Convention convention = ComplexOperator::HERMITIAN;
   unique_ptr<RT_FECollection> curl_finite_elements;
   unique_ptr<ParFiniteElementSpace> curl_space;
   unique_ptr<ParDiscreteLinearOperator> curl;
   if (vtk_directory[0])
   {
      curl_finite_elements.reset(new RT_FECollection(order - 1, 3));
      curl_space.reset(
          new ParFiniteElementSpace(&mesh, curl_finite_elements.get()));
      curl.reset(new ParDiscreteLinearOperator(&space, curl_space.get()));
      curl->AddDomainInterpolator(new CurlInterpolator);
      curl->Assemble();
      curl->Finalize();
   }

   // Cache receiver locations once. The mesh and survey geometry remain fixed
   // across frequencies and persistent inversion-server evaluations.
   vector<Vector> station_transmitters(station_indices.size());
   vector<vector<Vector>> station_receivers(station_indices.size());
   vector<unique_ptr<UDARBatchedPointEvaluator>> receiver_evaluators;
   receiver_evaluators.reserve(station_indices.size());
   for (int station_position = 0;
        station_position < static_cast<int>(station_indices.size());
        station_position++)
   {
      const int station = station_indices[station_position];
      const real_t fraction =
          (stations == 1) ? 0.0 : real_t(station) / real_t(stations - 1);
      station_transmitters[station_position] = survey_start;
      station_transmitters[station_position].Add(fraction, trajectory);
      MFEM_VERIFY(PointIsInside(mesh, station_transmitters[station_position],
                                point_tolerance),
                  "Transmitter station " << station
                                         << " lies outside the mesh.");
      station_receivers[station_position].reserve(spacings.size());
      for (real_t spacing : spacings)
      {
         Vector receiver(station_transmitters[station_position]);
         receiver.Add(-spacing, tool_axis);
         station_receivers[station_position].push_back(receiver);
      }
      receiver_evaluators.emplace_back(make_unique<UDARBatchedPointEvaluator>(
          mesh, station_receivers[station_position], point_tolerance));
   }

   vector<unique_ptr<UDARFrequencySystem>> frequency_systems;
   frequency_systems.reserve(frequencies.size());
   for (real_t frequency : frequencies)
   {
      UDARFrequencySystemOptions frequency_options;
      frequency_options.frequency = frequency;
      frequency_options.epsilon0 = eps0;
      frequency_options.mu0 = mu0;
      frequency_options.relative_permittivity = eps_r;
      frequency_options.relative_permeability = mu_r;
      frequency_options.boundary_condition = boundary_condition;
      frequency_options.use_mumps = use_mumps;
      frequency_options.reuse_mumps_ordering = reuse_mumps_ordering;
      frequency_options.assembly = assembly_mode;
      frequency_options.iterative.preconditioner = iterative_preconditioner;
      frequency_options.iterative.system_mode = ams_system_mode;
      frequency_options.iterative.maximum_iterations = max_iterations;
      frequency_options.iterative.krylov_dimension = krylov_dimension;
      frequency_options.iterative.relative_tolerance = relative_tolerance;
      frequency_options.iterative.absolute_tolerance = absolute_tolerance;
      frequency_options.iterative.ams_cycle_type = ams_cycle_type;
      frequency_options.iterative.ams_iterations = ams_iterations;
      frequency_options.iterative.ams_smoothing_iterations =
          ams_smoothing_iterations;
      frequency_options.iterative.ams_print_level = ams_print_level;
      frequency_options.iterative.ilu_fallback = use_ilu_fallback;
      frequency_systems.emplace_back(make_unique<UDARFrequencySystem>(
          space, sigma, pml.get(), essential_tdofs, frequency_options));
      if (myid == 0)
      {
         const auto &system = *frequency_systems.back();
         cout << "Maxwell assembly at " << frequency << " Hz: "
              << system.AssemblyMode();
         if (!system.AssemblyReason().empty())
         { cout << " (" << system.AssemblyReason() << ')'; }
         cout << '\n';
      }
   }

   /** @brief Solve and write one survey for the current conductivity tensors.
       @param current_output_file Rank-zero CSV destination.

       The enclosing scope retains mesh, finite-element space, tool geometry,
       and one UDARFrequencySystem per frequency. Material-dependent operators
       and iterative preconditioners are rebuilt; MUMPS objects retain reusable
       symbolic state for the fixed sparsity pattern. */
   auto RunForward = [&](const string &current_output_file)
   {
      const string evaluation_name =
          filesystem::path(current_output_file).stem().string();
      const long long global_elements = mesh.GetGlobalNE();
      WriteVTKModel(vtk_directory, evaluation_name, mesh, sigma, pml.get());
      ofstream output, timing_output;
      if (myid == 0)
      {
         output.open(current_output_file);
         MFEM_VERIFY(output.good(),
                     "Unable to open output file " << current_output_file);
         output << setprecision(16) << scientific;
         WriteCSVHeader(output);
         const filesystem::path response_path(current_output_file);
         const filesystem::path timing_path =
             response_path.parent_path() /
             (response_path.stem().string() + "_timing.csv");
         timing_output.open(timing_path);
         MFEM_VERIFY(timing_output.good(),
                     "Unable to open timing output " << timing_path.string());
         timing_output << setprecision(9) << scientific
                       << "frequency_hz,global_elements,complex_unknowns,"
                       << "solver,preconditioner,system_mode,assembly_mode,assembly_s,"
                       << "solver_setup_s,rhs_s,solve_s,field_output_s,"
                       << "sampling_s,vtk_save_s,derived_eval_s,"
                       << "krylov_iterations,ilu_fallbacks,vtk_profiled\n";
      }

      for (int frequency_index = 0;
           frequency_index < static_cast<int>(frequencies.size());
           frequency_index++)
      {
         UDARFrequencySystem &frequency_system =
             *frequency_systems[frequency_index];
         const real_t frequency = frequency_system.Frequency();
         const real_t omega = frequency_system.Omega();
         StopWatch phase_timer;
         phase_timer.Start();
         frequency_system.AssembleOperator(sigma_background(0, 0));
         phase_timer.Stop();
         const double assembly_seconds = phase_timer.RealTime();
         phase_timer.Restart();
         frequency_system.ConfigureSolver();
         phase_timer.Stop();
         const double solver_setup_seconds = phase_timer.RealTime();
         double rhs_seconds = 0.0;
         double solve_seconds = 0.0;
         double field_output_seconds = 0.0;
         double vtk_save_seconds = 0.0;
         double derived_eval_seconds = 0.0;
         double sampling_seconds = 0.0;
         int krylov_iterations = 0;
         int ilu_fallbacks = 0;

         if (myid == 0)
         {
            cout << "Solving " << frequency << " Hz with "
                 << frequency_system.SolverDescription()
                 << endl;
         }

         for (int station_position = 0;
              station_position < static_cast<int>(station_indices.size());
              station_position++)
         {
            const int station = station_indices[station_position];
            const real_t fraction =
                (stations == 1) ? 0.0 : real_t(station) / real_t(stations - 1);
            const Vector &transmitter = station_transmitters[station_position];
            vector<Vector> rhs(3), solution(3);
            phase_timer.Restart();
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
                              "Transmitter station "
                                  << station << " lies outside the mesh.");
               }
               else
               {
                  loop_current.reset(new FiniteLoopCurrentCoefficient(
                      transmitter, moment, loop_radius, wire_radius));
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
            phase_timer.Stop();
            rhs_seconds += phase_timer.RealTime();

            Array<const Vector *> rhs_array(3);
            Array<Vector *> solution_array(3);
            for (int axis = 0; axis < 3; axis++)
            {
               rhs_array[axis] = &rhs[axis];
               solution_array[axis] = &solution[axis];
            }
            phase_timer.Restart();
            const UDARFrequencySolveStats solve_stats =
                frequency_system.Solve(rhs_array, solution_array);
            krylov_iterations += solve_stats.krylov_iterations;
            ilu_fallbacks += solve_stats.ilu_fallbacks;
            if (solve_stats.ilu_fallbacks > 0 && myid == 0)
            {
               cout << "  AMS required distributed ILU recovery for "
                    << solve_stats.ilu_fallbacks
                    << " source solve(s) at station " << station << ".\n";
            }
            phase_timer.Stop();
            solve_seconds += phase_timer.RealTime();

            phase_timer.Restart();
            vector<unique_ptr<ParComplexGridFunction>> electric_fields(3);
            for (int axis = 0; axis < 3; axis++)
            {
               electric_fields[axis].reset(new ParComplexGridFunction(&space));
               Vector solution_real, solution_imaginary;
               solution_real.MakeRef(solution[axis], 0, true_size);
               solution_imaginary.MakeRef(solution[axis], true_size, true_size);
               electric_fields[axis]->real().SetFromTrueDofs(solution_real);
               electric_fields[axis]->imag().SetFromTrueDofs(
                   solution_imaginary);
            }

            vector<unique_ptr<ParGridFunction>> magnetic_real(3);
            vector<unique_ptr<ParGridFunction>> magnetic_imaginary(3);
            if (vtk_directory[0])
            {
               const real_t magnetic_scale = 1.0 / (omega * mu0 * mu_r);
               for (int axis = 0; axis < 3; axis++)
               {
                  magnetic_real[axis].reset(
                      new ParGridFunction(curl_space.get()));
                  magnetic_imaginary[axis].reset(
                      new ParGridFunction(curl_space.get()));
                  curl->Mult(electric_fields[axis]->imag(),
                             *magnetic_real[axis]);
                  curl->Mult(electric_fields[axis]->real(),
                             *magnetic_imaginary[axis]);
                  *magnetic_real[axis] *= -magnetic_scale;
                  *magnetic_imaginary[axis] *= magnetic_scale;
               }
               const auto vtk_stats = WriteUDARVTKFields(
                              vtk_directory, evaluation_name, frequency_index,
                              station, frequency, mesh, sigma, eps0 * eps_r,
                              mu0 * mu_r, electric_fields, magnetic_real,
                              magnetic_imaginary, profile_vtk);
               vtk_save_seconds += vtk_stats.save_seconds;
               derived_eval_seconds += vtk_stats.derived_evaluation_seconds;
            }
            phase_timer.Stop();
            field_output_seconds += phase_timer.RealTime();

            vector<const ParComplexGridFunction *> electric_field_views(3);
            for (int axis = 0; axis < 3; axis++)
            {
               electric_field_views[axis] = electric_fields[axis].get();
            }
            vector<DenseMatrix> receiver_field_real;
            vector<DenseMatrix> receiver_field_imaginary;
            phase_timer.Restart();
            receiver_evaluators[station_position]->EvaluateMagneticFields(
                electric_field_views, omega, mu0, mu_r, receiver_field_real,
                receiver_field_imaginary);
            phase_timer.Stop();
            sampling_seconds += phase_timer.RealTime();

            for (int spacing_index = 0;
                 spacing_index < static_cast<int>(spacings.size());
                 spacing_index++)
            {
               const real_t spacing = spacings[spacing_index];
               const Vector &receiver =
                   station_receivers[station_position][spacing_index];
               complex<real_t> response[3][3];
               const bool receiver_found =
                   receiver_evaluators[station_position]->Found(spacing_index);

               for (int transmitter_axis = 0; transmitter_axis < 3;
                    transmitter_axis++)
               {
                  Vector field_real(3), field_imaginary(3);
                  receiver_field_real[spacing_index].GetColumnReference(
                      transmitter_axis, field_real);
                  receiver_field_imaginary[spacing_index].GetColumnReference(
                      transmitter_axis, field_imaginary);

                  for (int receiver_axis = 0; receiver_axis < 3;
                       receiver_axis++)
                  {
                     Vector receiver_direction(3);
                     GetFrameAxis(tool_frame, receiver_axis,
                                  receiver_direction);
                     if (receiver_found)
                     {
                        response[receiver_axis][transmitter_axis] =
                            complex<real_t>(
                                Dot3(receiver_direction, field_real),
                                Dot3(receiver_direction, field_imaginary)) /
                            transmitter_moment;
                     }
                     else
                     {
                        const real_t nan = numeric_limits<real_t>::quiet_NaN();
                        response[receiver_axis][transmitter_axis] =
                            complex<real_t>(nan, nan);
                     }
                  }
               }

               if (!receiver_found && myid == 0)
               {
                  cerr << "Warning: receiver for station " << station
                       << ", spacing " << spacing
                       << " lies outside the mesh.\n";
               }

               if (myid == 0)
               {
                  output << station << ','
                         << measured_depth_origin + fraction * survey_length
                         << ',' << frequency << ',' << spacing;
                  for (int d = 0; d < 3; d++)
                  {
                     output << ',' << transmitter[d];
                  }
                  for (int d = 0; d < 3; d++) { output << ',' << receiver[d]; }
                  for (int receiver_axis = 0; receiver_axis < 3;
                       receiver_axis++)
                  {
                     for (int transmitter_axis = 0; transmitter_axis < 3;
                          transmitter_axis++)
                     {
                        output
                            << ','
                            << response[receiver_axis][transmitter_axis].real()
                            << ','
                            << response[receiver_axis][transmitter_axis].imag();
                     }
                  }
                  output << '\n';
               }
            }

            if (myid == 0)
            {
               cout << "  station " << station + 1 << '/' << stations
                    << " (selected " << station_position + 1 << '/'
                    << station_indices.size() << ")\n";
            }
         }

         const double local_timing[8] = {
             assembly_seconds, solver_setup_seconds, rhs_seconds,
             solve_seconds, field_output_seconds, sampling_seconds,
             vtk_save_seconds, derived_eval_seconds};
         double maximum_timing[8] = {};
         MPI_Reduce(local_timing, maximum_timing, 8, MPI_DOUBLE, MPI_MAX, 0,
                    mesh.GetComm());
         int global_krylov_iterations = 0;
         int maximum_fallbacks = 0;
         MPI_Reduce(&krylov_iterations, &global_krylov_iterations, 1, MPI_INT,
                    MPI_MAX, 0, mesh.GetComm());
         MPI_Reduce(&ilu_fallbacks, &maximum_fallbacks, 1, MPI_INT, MPI_MAX, 0,
                    mesh.GetComm());
         if (myid == 0)
         {
            const string solver_name = frequency_system.SolverDescription();
            timing_output << frequency << ',' << global_elements << ','
                          << space.GlobalTrueVSize() << ',' << solver_name
                          << ','
                          << (frequency_system.UsesMUMPS()
                                  ? "direct"
                                  : frequency_system.Preconditioner())
                          << ','
                          << (frequency_system.UsesMUMPS()
                                  ? "complex-or-real-block"
                                  : frequency_system.SystemMode())
                          << ',' << frequency_system.AssemblyMode();
            for (double value : maximum_timing)
            {
               timing_output << ',' << value;
            }
            timing_output << ',' << global_krylov_iterations << ','
                          << maximum_fallbacks << ','
                          << (profile_vtk && vtk_directory[0] ? 1 : 0) << '\n';
            cout << "  timing (max rank): assembly " << maximum_timing[0]
                 << " s, setup " << maximum_timing[1] << " s, RHS "
                 << maximum_timing[2] << " s, solve " << maximum_timing[3]
                 << " s, fields " << maximum_timing[4] << " s, sampling "
                 << maximum_timing[5] << " s\n";
         }
      }

      if (myid == 0)
      {
         cout << "Wrote full complex 3x3 H tensor responses to "
              << current_output_file << '\n';
      }
   };

   if (!server_directory[0])
   {
      RunForward(output_file);
      return 0;
   }

   if (myid == 0)
   {
      MFEM_VERIFY(
          filesystem::is_directory(server_directory),
          "Server output directory does not exist: " << server_directory);
      cout << "UDAR_SERVER_READY 1" << endl;
   }
   /** @brief Serve scalar isotropic conductivity requests from Python.

       Rank 0 accepts @c E request_id sigma_background sigma_target, with S/m
       conductivities, or @c Q. Valid requests are broadcast to all MPI ranks.
       The protocol supports only the two scalar conductivities; geometry and
       anisotropy remain fixed for the process lifetime. */
   while (true)
   {
      int opcode = 0;
      long long request = -1;
      real_t background_value = 0.0;
      real_t target_value = 0.0;
      if (myid == 0)
      {
         string line;
         if (!getline(cin, line)) { opcode = 2; }
         else
         {
            istringstream input(line);
            char command = 0;
            input >> command;
            string extra;
            if (command == 'Q' && !(input >> extra)) { opcode = 2; }
            else if (command == 'E' &&
                     (input >> request >> background_value >> target_value) &&
                     request >= 0 && !(input >> extra))
            {
               opcode = 1;
            }
            else
            {
               opcode = -1;
            }
         }
      }
      MPI_Bcast(&opcode, 1, MPI_INT, 0, mesh.GetComm());
      if (opcode == 2) { break; }
      if (opcode < 0)
      {
         if (myid == 0) { cout << "UDAR_SERVER_ERROR invalid_request" << endl; }
         continue;
      }
      MPI_Bcast(&request, 1, MPI_LONG_LONG, 0, mesh.GetComm());
      MPI_Bcast(&background_value, 1, MFEM_MPI_REAL_T, 0, mesh.GetComm());
      MPI_Bcast(&target_value, 1, MFEM_MPI_REAL_T, 0, mesh.GetComm());
      MFEM_VERIFY(isfinite(background_value) && background_value > 0.0 &&
                      isfinite(target_value) && target_value > 0.0,
                  "Server conductivities must be positive and finite.");
      sigma_background = 0.0;
      sigma_target = 0.0;
      for (int d = 0; d < 3; d++)
      {
         sigma_background(d, d) = background_value;
         sigma_target(d, d) = target_value;
      }
      sigma.SetConductivities(sigma_background, sigma_target);
      ostringstream filename;
      filename << "evaluation_" << setw(5) << setfill('0') << request << ".csv";
      const string path =
          (filesystem::path(server_directory) / filename.str()).string();
      RunForward(path);
      if (myid == 0) { cout << "UDAR_SERVER_DONE " << request << endl; }
   }
   return 0;
}
