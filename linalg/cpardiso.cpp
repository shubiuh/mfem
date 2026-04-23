// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "cpardiso.hpp"
#include "hypre.hpp"
#include <algorithm>
#include <vector>
#include <numeric>

#ifdef MFEM_USE_MPI
#ifdef MFEM_USE_MKL_CPARDISO

namespace mfem
{
CPardisoSolver::CPardisoSolver(MPI_Comm comm, int nrhs_)
{
   comm_ = MPI_Comm_c2f(comm);

   // Indicate that default parameters are changed
   iparm[0] = 1;
   // Use METIS for fill-in reordering
   iparm[1] = 2;
   // Do not write the solution into the x vector data
   iparm[5] = 0;
   // Maximum number of iterative refinement steps
   iparm[7] = 2;
   // Perturb the pivot elements with 1E-13
   iparm[9] = 13;
   // Use nonsymmetric permutation
   iparm[10] = 1;
   // Perform a check on the input data
   iparm[26] = 1;
#ifdef MFEM_USE_SINGLE
   // Single precision
   iparm[27] = 1;
#endif
   // 0-based indexing in CSR data structure
   iparm[34] = 1;
   // All inputs are distributed between MPI processes
   iparm[39] = 2;
   // Maximum number of numerical factorizations
   maxfct = 1;
   // Which factorization to use. This parameter is ignored and always assumed
   // to be equal to 1. See MKL documentation.
   mnum = 1;
   // Print statistical information in file
   msglvl = 0;
   // Initialize error flag
   error = 0;
   // Real nonsymmetric matrix
   mtype = MatType::REAL_NONSYMMETRIC;
   // Number of right hand sides
   nrhs = nrhs_;
};

void CPardisoSolver::SetOperator(const Operator &op)
{
   auto hypreParMat = dynamic_cast<const HypreParMatrix &>(op);

   MFEM_ASSERT(hypreParMat, "Must pass HypreParMatrix as Operator");

   auto parcsr_op = static_cast<hypre_ParCSRMatrix *>(
                       const_cast<HypreParMatrix &>(hypreParMat));

   hypreParMat.HostRead();
   hypre_CSRMatrix *csr_op = hypre_MergeDiagAndOffd(parcsr_op);
   hypreParMat.HypreRead();
#if MFEM_HYPRE_VERSION >= 21600
   hypre_CSRMatrixBigJtoJ(csr_op);
#endif

   m = parcsr_op->global_num_rows;
   first_row = parcsr_op->first_row_index;
   nnz_loc = csr_op->num_nonzeros;
   m_loc = csr_op->num_rows;

   height = m_loc;
   width = m_loc;

   real_t *csr_nzval = csr_op->data;
   int *csr_colind = csr_op->j;

   delete[] csr_rowptr;
   delete[] reordered_csr_colind;
   delete[] reordered_csr_nzval;
   csr_rowptr = new int[m_loc + 1];
   reordered_csr_colind = new int[nnz_loc];
   reordered_csr_nzval = new real_t[nnz_loc];

   for (int i = 0; i <= m_loc; i++)
   {
      csr_rowptr[i] = (csr_op->i)[i];
   }

   // CPardiso expects the column indices to be sorted for each row
   std::vector<int> permutation_idx(nnz_loc);
   std::iota(permutation_idx.begin(), permutation_idx.end(), 0);
   for (int i = 0; i < m_loc; i++)
   {
      std::sort(permutation_idx.begin() + csr_rowptr[i],
                permutation_idx.begin() + csr_rowptr[i + 1],
                [csr_colind](int i1, int i2)
      {
         return csr_colind[i1] < csr_colind[i2];
      });
   }

   for (int i = 0; i < nnz_loc; i++)
   {
      reordered_csr_colind[i] = csr_colind[permutation_idx[i]];
      reordered_csr_nzval[i] = csr_nzval[permutation_idx[i]];
   }

   hypre_CSRMatrixDestroy(csr_op);

   // iparm[40], the number of row in global matrix, rhs element and solution vector that
   // begins the input domain belonging to this MPI process

   // iparm[41], the number of row in global matrix, rhs element and solution vector that
   // ends the input domain belonging to this MPI process
   if (m_loc == 0 && first_row == 0)
   {
      // Workaround for the issue https://github.com/mfem/mfem/issues/4634
      iparm[40] = 1;
      iparm[41] = first_row;
   }
   else
   {
      iparm[40] = first_row;
      iparm[41] = first_row + m_loc - 1;
   }

   // Analyze inputs
   phase = 11;
   cluster_sparse_solver(pt,
                         &maxfct,
                         &mnum,
                         &mtype,
                         &phase,
                         &m,
                         reordered_csr_nzval,
                         csr_rowptr,
                         reordered_csr_colind,
                         &idum,
                         &nrhs,
                         iparm,
                         &msglvl,
                         &ddum,
                         &ddum,
                         &comm_,
                         &error);

   MFEM_ASSERT(error == 0, "CPardiso analyze input error");

   // Numerical factorization
   phase = 22;
   cluster_sparse_solver(pt,
                         &maxfct,
                         &mnum,
                         &mtype,
                         &phase,
                         &m,
                         reordered_csr_nzval,
                         csr_rowptr,
                         reordered_csr_colind,
                         &idum,
                         &nrhs,
                         iparm,
                         &msglvl,
                         &ddum,
                         &ddum,
                         &comm_,
                         &error);

   MFEM_ASSERT(error == 0, "CPardiso factorization input error");
}

void CPardisoSolver::Mult(const Vector &b, Vector &x) const
{
   // Solve
   phase = 33;
   cluster_sparse_solver(pt,
                         &maxfct,
                         &mnum,
                         &mtype,
                         &phase,
                         &m,
                         reordered_csr_nzval,
                         csr_rowptr,
                         reordered_csr_colind,
                         &idum,
                         &nrhs,
                         iparm,
                         &msglvl,
                         b.GetData(),
                         x.GetData(),
                         &comm_,
                         &error);

   MFEM_ASSERT(error == 0, "Pardiso solve error");
}

void CPardisoSolver::Mult(const DenseMatrix &B, DenseMatrix &X) const
{
   MFEM_ASSERT(B.Width() == X.Width(), "Incompatible matrix sizes");

   // Solve for multiple right-hand side
   phase = 33;
   cluster_sparse_solver(pt,
                         &maxfct,
                         &mnum,
                         &mtype,
                         &phase,
                         &m,
                         reordered_csr_nzval,
                         csr_rowptr,
                         reordered_csr_colind,
                         &idum,
                         &nrhs,
                         iparm,
                         &msglvl,
                         B.GetData(), // n x nrhs, column-major
                         X.GetData(), // solution n x nrhs
                         &comm_,
                         &error);

   MFEM_ASSERT(error == 0, "Pardiso solve error");
}

void CPardisoSolver::SetPrintLevel(int print_level)
{
   msglvl = print_level;
}

void CPardisoSolver::setRHSCount(int nrhs_)
{
   nrhs = nrhs_;
}

void CPardisoSolver::SetMatrixType(MatType mat_type)
{
   mtype = mat_type;
}

CPardisoSolver::~CPardisoSolver()
{
   // Release all internal memory
   phase = -1;
   cluster_sparse_solver(pt,
                         &maxfct,
                         &mnum,
                         &mtype,
                         &phase,
                         &m,
                         reordered_csr_nzval,
                         csr_rowptr,
                         reordered_csr_colind,
                         &idum,
                         &nrhs,
                         iparm,
                         &msglvl,
                         &ddum,
                         &ddum,
                         &comm_,
                         &error);

   MFEM_ASSERT(error == 0, "CPardiso free error");

   delete[] csr_rowptr;
   delete[] reordered_csr_colind;
   delete[] reordered_csr_nzval;
}

// ============================================================================
// ComplexCPardisoSolver implementation
// ============================================================================

ComplexCPardisoSolver::ComplexCPardisoSolver(MPI_Comm comm, int nrhs_)
{
   comm_ = MPI_Comm_c2f(comm);

   // Indicate that default parameters are changed
   iparm[0] = 1;
   // Use METIS for fill-in reordering
   iparm[1] = 2;
   // Do not write the solution into the x vector data
   iparm[5] = 0;
   // Maximum number of iterative refinement steps
   iparm[7] = 2;
   // Perturb the pivot elements with 1E-13
   iparm[9] = 13;
   // Use nonsymmetric permutation
   iparm[10] = 1;
   // Perform a check on the input data
   iparm[26] = 1;
#ifdef MFEM_USE_SINGLE
   // Single precision
   iparm[27] = 1;
#endif
   // 0-based indexing in CSR data structure
   iparm[34] = 1;
   // All inputs are distributed between MPI processes
   iparm[39] = 2;
   // Maximum number of numerical factorizations
   maxfct = 1;
   // Which factorization to use
   mnum = 1;
   // Print statistical information in file
   msglvl = 0;
   // Initialize error flag
   error = 0;
   // Complex nonsymmetric matrix (default)
   mtype = static_cast<int>(MatType::COMPLEX_NONSYMMETRIC);
   // Number of right hand sides
   nrhs = nrhs_;
   // Zero the dummy complex scalar
   ddum_cpx.real = 0;
   ddum_cpx.imag = 0;
}

void ComplexCPardisoSolver::SetOperator(const Operator &op)
{
   const auto *chpm = dynamic_cast<const ComplexHypreParMatrix *>(&op);
   MFEM_VERIFY(chpm,
               "ComplexCPardisoSolver::SetOperator requires a ComplexHypreParMatrix");
   SetOperator(*chpm);
}

void ComplexCPardisoSolver::SetOperator(const ComplexHypreParMatrix &op)
{
   const HypreParMatrix &A_r = op.real();
   const HypreParMatrix &A_i = op.imag();

   // Merge diagonal and off-diagonal blocks for A_r
   auto *parcsr_r = static_cast<hypre_ParCSRMatrix *>(
                       const_cast<HypreParMatrix &>(A_r));
   const_cast<HypreParMatrix &>(A_r).HostRead();
   hypre_CSRMatrix *csr_r = hypre_MergeDiagAndOffd(parcsr_r);
   const_cast<HypreParMatrix &>(A_r).HypreRead();
#if MFEM_HYPRE_VERSION >= 21600
   hypre_CSRMatrixBigJtoJ(csr_r);
#endif

   // Merge diagonal and off-diagonal blocks for A_i
   auto *parcsr_i = static_cast<hypre_ParCSRMatrix *>(
                       const_cast<HypreParMatrix &>(A_i));
   const_cast<HypreParMatrix &>(A_i).HostRead();
   hypre_CSRMatrix *csr_i = hypre_MergeDiagAndOffd(parcsr_i);
   const_cast<HypreParMatrix &>(A_i).HypreRead();
#if MFEM_HYPRE_VERSION >= 21600
   hypre_CSRMatrixBigJtoJ(csr_i);
#endif

   m         = static_cast<int>(parcsr_r->global_num_rows);
   first_row = static_cast<int>(parcsr_r->first_row_index);
   m_loc     = csr_r->num_rows;
   height    = m_loc;
   width     = m_loc;

   const int nnz_r = csr_r->num_nonzeros;
   const int nnz_i = csr_i->num_nonzeros;

   int *col_r = csr_r->j;
   int *col_i = csr_i->j;

   // Build sort permutations so column indices within each row are ascending.
   std::vector<int> perm_r(nnz_r), perm_i(nnz_i);
   std::iota(perm_r.begin(), perm_r.end(), 0);
   std::iota(perm_i.begin(), perm_i.end(), 0);
   for (int i = 0; i < m_loc; i++)
   {
      std::sort(perm_r.begin() + csr_r->i[i], perm_r.begin() + csr_r->i[i + 1],
                [col_r](int a, int b) { return col_r[a] < col_r[b]; });
      std::sort(perm_i.begin() + csr_i->i[i], perm_i.begin() + csr_i->i[i + 1],
                [col_i](int a, int b) { return col_i[a] < col_i[b]; });
   }

   // For symmetric/Hermitian types only the lower triangular part is stored.
   const bool lower_tri = (mtype != static_cast<int>(MatType::COMPLEX_NONSYMMETRIC));

   // First pass: count nonzeros in the merged complex CSR.
   nnz_loc = 0;
   for (int i = 0; i < m_loc; i++)
   {
      const int grow = first_row + i; // 0-based global row
      int jr = csr_r->i[i], jrend = csr_r->i[i + 1];
      int ji = csr_i->i[i], jiend = csr_i->i[i + 1];
      while (jr < jrend || ji < jiend)
      {
         int cr = (jr < jrend) ? col_r[perm_r[jr]] : INT_MAX;
         int ci = (ji < jiend) ? col_i[perm_i[ji]] : INT_MAX;
         int col = (cr < ci) ? cr : ci;
         if (cr == col) { jr++; }
         if (ci == col) { ji++; }
         if (lower_tri && col > grow) { continue; }
         nnz_loc++;
      }
   }

   // Allocate merged complex CSR storage.
   delete[] csr_rowptr;
   delete[] reordered_csr_colind;
   delete[] reordered_csr_nzval;
   csr_rowptr           = new int[m_loc + 1];
   reordered_csr_colind = new int[nnz_loc];
   reordered_csr_nzval  = new CPardisoComplexScalar[nnz_loc];

   // Second pass: fill merged complex CSR.
   csr_rowptr[0] = 0;
   int k = 0;
   for (int i = 0; i < m_loc; i++)
   {
      const int grow = first_row + i;
      int jr = csr_r->i[i], jrend = csr_r->i[i + 1];
      int ji = csr_i->i[i], jiend = csr_i->i[i + 1];
      while (jr < jrend || ji < jiend)
      {
         int cr = (jr < jrend) ? col_r[perm_r[jr]] : INT_MAX;
         int ci = (ji < jiend) ? col_i[perm_i[ji]] : INT_MAX;
         int col = (cr < ci) ? cr : ci;
         real_t vr = 0.0, vi = 0.0;
         if (cr == col) { vr = csr_r->data[perm_r[jr]]; jr++; }
         if (ci == col) { vi = csr_i->data[perm_i[ji]]; ji++; }
         if (lower_tri && col > grow) { continue; }
         reordered_csr_colind[k]       = col;
         reordered_csr_nzval[k].real   =
            static_cast<decltype(reordered_csr_nzval[k].real)>(vr);
         reordered_csr_nzval[k].imag   =
            static_cast<decltype(reordered_csr_nzval[k].imag)>(vi);
         k++;
      }
      csr_rowptr[i + 1] = k;
   }

   hypre_CSRMatrixDestroy(csr_r);
   hypre_CSRMatrixDestroy(csr_i);

   // Set the distributed range for this MPI rank.
   // iparm[40]: first row (0-based), iparm[41]: last row (0-based, inclusive).
   if (m_loc == 0 && first_row == 0)
   {
      // Workaround: same as CPardisoSolver for empty local share
      iparm[40] = 1;
      iparm[41] = first_row;
   }
   else
   {
      iparm[40] = first_row;
      iparm[41] = first_row + m_loc - 1;
   }

   // Symbolic analysis (phase 11)
   phase = 11;
   cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &m,
                         reordered_csr_nzval, csr_rowptr, reordered_csr_colind,
                         &idum, &nrhs, iparm, &msglvl,
                         &ddum_cpx, &ddum_cpx, &comm_, &error);
   MFEM_ASSERT(error == 0, "ComplexCPardiso analysis error");

   // Numerical factorization (phase 22)
   phase = 22;
   cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &m,
                         reordered_csr_nzval, csr_rowptr, reordered_csr_colind,
                         &idum, &nrhs, iparm, &msglvl,
                         &ddum_cpx, &ddum_cpx, &comm_, &error);
   MFEM_ASSERT(error == 0, "ComplexCPardiso factorization error");
}

void ComplexCPardisoSolver::Mult(const Vector &b_r, const Vector &b_i,
                                 Vector &x_r, Vector &x_i) const
{
   MFEM_ASSERT(b_r.Size() == m_loc && b_i.Size() == m_loc,
               "ComplexCPardisoSolver::Mult: RHS size mismatch");

   // Pack interleaved complex RHS buffer
   std::vector<CPardisoComplexScalar> b_cpx(m_loc), x_cpx(m_loc);
   for (int i = 0; i < m_loc; i++)
   {
      b_cpx[i].real = static_cast<decltype(b_cpx[i].real)>(b_r[i]);
      b_cpx[i].imag = static_cast<decltype(b_cpx[i].imag)>(b_i[i]);
   }

   phase = 33;
   cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &m,
                         reordered_csr_nzval, csr_rowptr, reordered_csr_colind,
                         &idum, &nrhs, iparm, &msglvl,
                         b_cpx.data(), x_cpx.data(), &comm_, &error);
   MFEM_ASSERT(error == 0, "ComplexCPardiso solve error");

   // Unpack interleaved complex solution
   x_r.SetSize(m_loc);
   x_i.SetSize(m_loc);
   for (int i = 0; i < m_loc; i++)
   {
      x_r[i] = static_cast<real_t>(x_cpx[i].real);
      x_i[i] = static_cast<real_t>(x_cpx[i].imag);
   }
}

void ComplexCPardisoSolver::Mult(const Vector &b, Vector &x) const
{
   MFEM_ASSERT(b.Size() == 2 * m_loc,
               "ComplexCPardisoSolver::Mult: b must have size 2*m_loc");
   x.SetSize(2 * m_loc);

   // View [re; im] sub-vectors without copying.
   const Vector b_r(const_cast<real_t *>(b.GetData()),         m_loc);
   const Vector b_i(const_cast<real_t *>(b.GetData()) + m_loc, m_loc);
   Vector x_r(x.GetData(),         m_loc);
   Vector x_i(x.GetData() + m_loc, m_loc);

   Mult(b_r, b_i, x_r, x_i);
}

void ComplexCPardisoSolver::SetPrintLevel(int print_lvl)
{
   msglvl = print_lvl;
}

void ComplexCPardisoSolver::setRHSCount(int nrhs_)
{
   nrhs = nrhs_;
}

void ComplexCPardisoSolver::SetMatrixType(MatType mat_type)
{
   mtype = static_cast<int>(mat_type);
}

ComplexCPardisoSolver::~ComplexCPardisoSolver()
{
   // Release all internal solver memory
   phase = -1;
   cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &m,
                         reordered_csr_nzval, csr_rowptr, reordered_csr_colind,
                         &idum, &nrhs, iparm, &msglvl,
                         &ddum_cpx, &ddum_cpx, &comm_, &error);
   MFEM_ASSERT(error == 0, "ComplexCPardiso free error");

   delete[] csr_rowptr;
   delete[] reordered_csr_colind;
   delete[] reordered_csr_nzval;
}

} // namespace mfem

#endif // MFEM_USE_MKL_CPARDISO
#endif // MFEM_USE_MPI
