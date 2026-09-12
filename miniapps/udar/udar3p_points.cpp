// Copyright (c) 2026, Lawrence Livermore National Security, LLC.
// See the top-level LICENSE file for details.

#include "udar3p_points.hpp"

namespace mfem
{

UDARBatchedPointEvaluator::UDARBatchedPointEvaluator(
    ParMesh &mesh_, const std::vector<Vector> &points, real_t tolerance)
   : mesh(mesh_), local_elements(points.size()), local_points(points.size()),
     global_matches(points.size(), 0)
{
   std::vector<int> local_matches(points.size(), 0);
   for (int element = 0; element < mesh.GetNE(); element++)
   {
      ElementTransformation *transformation =
          mesh.GetElementTransformation(element);
      for (int point = 0; point < static_cast<int>(points.size()); point++)
      {
         IntegrationPoint reference_point;
         if (transformation->TransformBack(points[point], reference_point,
                                           tolerance) ==
             InverseElementTransformation::Inside)
         {
            local_elements[point].push_back(element);
            local_points[point].push_back(reference_point);
            local_matches[point]++;
         }
      }
   }
   MPI_Allreduce(local_matches.data(), global_matches.data(),
                 static_cast<int>(points.size()), MPI_INT, MPI_SUM,
                 mesh.GetComm());
}

void UDARBatchedPointEvaluator::EvaluateMagneticFields(
    const std::vector<const ParComplexGridFunction *> &electric_fields,
    real_t omega, real_t mu0, real_t mu_r,
    std::vector<DenseMatrix> &field_real,
    std::vector<DenseMatrix> &field_imaginary) const
{
   MFEM_VERIFY(omega > 0.0 && mu0 > 0.0 && mu_r > 0.0,
               "Magnetic-field sampling requires positive physical values.");
   const int fields = static_cast<int>(electric_fields.size());
   const int entries = Size() * fields * 3;
   Vector local(2 * entries), global(2 * entries);
   local = 0.0;

   for (int point = 0; point < Size(); point++)
   {
      for (int match = 0;
           match < static_cast<int>(local_elements[point].size()); match++)
      {
         ElementTransformation *transformation =
             mesh.GetElementTransformation(local_elements[point][match]);
         transformation->SetIntPoint(&local_points[point][match]);
         for (int field = 0; field < fields; field++)
         {
            Vector curl_real(3), curl_imaginary(3);
            electric_fields[field]->real().GetCurl(*transformation, curl_real);
            electric_fields[field]->imag().GetCurl(*transformation,
                                                    curl_imaginary);
            for (int d = 0; d < 3; d++)
            {
               const int index = (point * fields + field) * 3 + d;
               local[index] += curl_real[d];
               local[index + entries] += curl_imaginary[d];
            }
         }
      }
   }

   MPI_Allreduce(local.GetData(), global.GetData(), 2 * entries,
                 MFEM_MPI_REAL_T, MPI_SUM, mesh.GetComm());
   const real_t magnetic_scale = 1.0 / (omega * mu0 * mu_r);
   field_real.resize(Size());
   field_imaginary.resize(Size());
   for (int point = 0; point < Size(); point++)
   {
      field_real[point].SetSize(3, fields);
      field_imaginary[point].SetSize(3, fields);
      if (!Found(point))
      {
         field_real[point] = 0.0;
         field_imaginary[point] = 0.0;
         continue;
      }
      for (int field = 0; field < fields; field++)
      {
         for (int d = 0; d < 3; d++)
         {
            const int index = (point * fields + field) * 3 + d;
            field_real[point](d, field) =
                -magnetic_scale * global[index + entries] /
                global_matches[point];
            field_imaginary[point](d, field) =
                magnetic_scale * global[index] / global_matches[point];
         }
      }
   }
}

} // namespace mfem
