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

#include "../config/config.hpp"
#include "../general/communication.hpp"

#ifdef MFEM_USE_MUMPS
#ifdef MFEM_USE_MPI

#include "mumps.hpp"

#include <algorithm>

#if MFEM_MUMPS_VERSION >= 530
#ifdef MUMPS_INTSIZE64
#error "Full 64-bit MUMPS is not yet supported"
#endif
#else
#ifdef INTSIZE64
#error "Full 64-bit MUMPS is not yet supported"
#endif
#endif

// Macro s.t. indices match MUMPS documentation
#define MUMPS_ICNTL(I) icntl[(I) -1]
#define MUMPS_CNTL(I) cntl[(I) -1]
#define MUMPS_INFO(I) info[(I) -1]
#define MUMPS_INFOG(I) infog[(I) -1]

namespace mfem
{

MUMPSSolver::MUMPSSolver(MPI_Comm comm_)
{
   Init(comm_);
}

MUMPSSolver::MUMPSSolver(const Operator &op)
{
   auto APtr = dynamic_cast<const HypreParMatrix *>(&op);
   MFEM_VERIFY(APtr, "Not a compatible matrix type");
   Init(APtr->GetComm());
   SetOperator(op);
}

void MUMPSSolver::Init(MPI_Comm comm_)
{
   id = nullptr;
   comm = comm_;
   MPI_Comm_size(comm, &numProcs);
   MPI_Comm_rank(comm, &myid);

   mat_type = MatType::UNSYMMETRIC;
   print_level = 0;
   reorder_method = ReorderingStrategy::AUTOMATIC;
   reorder_reuse = false;
#if MFEM_MUMPS_VERSION >= 510
   blr_mode = 0;
   blr_tol = 0.0;
   blr_compression_type = 0;
#endif
#if MFEM_MUMPS_VERSION >= 550
   blr_cb_compression = 0;
#endif
   mem_relaxation = 20;
   num_threads = 0;
   pivot_threshold = -1.0;
   out_of_core = 0;

#if MFEM_MUMPS_VERSION >= 530
   irhs_loc = nullptr;
   rhs_loc = nullptr;
   isol_loc = nullptr;
   sol_loc = nullptr;
#else
   recv_counts = nullptr;
   displs = nullptr;
   rhs_glob = nullptr;
#endif
}

MUMPSSolver::~MUMPSSolver()
{
#if MFEM_MUMPS_VERSION >= 530
   delete [] irhs_loc;
   delete [] rhs_loc;
   delete [] isol_loc;
   delete [] sol_loc;
#else
   delete [] recv_counts;
   delete [] displs;
   delete [] rhs_glob;
#endif
   if (id)
   {
      id->job = -2;
#ifdef MFEM_USE_SINGLE
      smumps_c(id);
#else
      dmumps_c(id);
#endif
      delete id;
   }
}

void MUMPSSolver::SetOperator(const Operator &op)
{
   auto APtr = dynamic_cast<const HypreParMatrix *>(&op);
   MFEM_VERIFY(APtr, "Not a compatible matrix type");

   height = op.Height();
   width = op.Width();

   auto parcsr_op = (hypre_ParCSRMatrix *)const_cast<HypreParMatrix &>(*APtr);
   APtr->HostRead();
   hypre_CSRMatrix *csr_op = hypre_MergeDiagAndOffd(parcsr_op);
   APtr->HypreRead();
   HYPRE_Int       *Iptr   = csr_op->i;
#if MFEM_HYPRE_VERSION >= 21600
   HYPRE_BigInt    *Jptr   = csr_op->big_j;
#else
   HYPRE_Int       *Jptr   = csr_op->j;
#endif

   int n_loc = internal::to_int(csr_op->num_rows);
   row_start = internal::to_int(parcsr_op->first_row_index);

   MUMPS_INT8 nnz = 0, k = 0;
   if (mat_type)
   {
      // Count nnz in case of symmetric mode
      for (int i = 0; i < n_loc; i++)
      {
         for (HYPRE_Int j = Iptr[i]; j < Iptr[i + 1]; j++)
         {
            int ii = row_start + i + 1;
#if MFEM_HYPRE_VERSION >= 21600
            HYPRE_BigInt jj = Jptr[k] + 1;
#else
            HYPRE_Int jj = Jptr[k] + 1;
#endif
            if (ii >= jj) { nnz++; }
            k++;
         }
      }
   }
   else
   {
      nnz = csr_op->num_nonzeros;
   }
   int *I = new int[nnz];
   int *J = new int[nnz];

   // Fill in I and J arrays for
   // COO format in 1-based indexing
   k = 0;
   real_t *data;
   if (mat_type)
   {
      MUMPS_INT8 l = 0;
      data = new real_t[nnz];
      for (int i = 0; i < n_loc; i++)
      {
         for (HYPRE_Int j = Iptr[i]; j < Iptr[i + 1]; j++)
         {
            int ii = row_start + i + 1;
#if MFEM_HYPRE_VERSION >= 21600
            HYPRE_BigInt jj = Jptr[k] + 1;
#else
            HYPRE_Int jj = Jptr[k] + 1;
#endif
            if (ii >= jj)
            {
               I[l] = ii;
               J[l] = internal::to_int(jj);
               data[l++] = csr_op->data[k];
            }
            k++;
         }
      }
   }
   else
   {
      for (int i = 0; i < n_loc; i++)
      {
         for (HYPRE_Int j = Iptr[i]; j < Iptr[i + 1]; j++)
         {
            I[k] = row_start + i + 1;
            J[k] = internal::to_int(Jptr[k] + 1);
            k++;
         }
      }
      data = csr_op->data;
   }

   // New MUMPS object or reuse the one from a previous matrix
   if (!id || !reorder_reuse)
   {
      if (id)
      {
         id->job = -2;
#ifdef MFEM_USE_SINGLE
         smumps_c(id);
#else
         dmumps_c(id);
#endif
         delete id;
      }
#ifdef MFEM_USE_SINGLE
      id = new SMUMPS_STRUC_C();
#else
      id = new DMUMPS_STRUC_C();
#endif
      id->sym = mat_type;

      // C to Fortran communicator
      id->comm_fortran = (MUMPS_INT)MPI_Comm_c2f(comm);

      // Host is involved in computation
      id->par = 1;

      // MUMPS init
      id->job = -1;
#ifdef MFEM_USE_SINGLE
      smumps_c(id);
#else
      dmumps_c(id);
#endif

      // Set MUMPS default parameters
      SetParameters();

      id->n = internal::to_int(parcsr_op->global_num_rows);
      id->nnz_loc = nnz;
      id->irn_loc = I;
      id->jcn_loc = J;
      id->a_loc = data;

      // MUMPS analysis
      id->job = 1;
#ifdef MFEM_USE_SINGLE
      smumps_c(id);
#else
      dmumps_c(id);
#endif
   }
   else
   {
      id->irn_loc = I;
      id->jcn_loc = J;
      id->a_loc = data;
   }

   // MUMPS factorization
   id->job = 2;
   {
      const int mem_relax_lim = 200;
      while (true)
      {
#ifdef MFEM_USE_SINGLE
         smumps_c(id);
#else
         dmumps_c(id);
#endif
         if (id->MUMPS_INFOG(1) < 0)
         {
            if (id->MUMPS_INFOG(1) == -8 || id->MUMPS_INFOG(1) == -9)
            {
               id->MUMPS_ICNTL(14) += 20;
               MFEM_VERIFY(id->MUMPS_ICNTL(14) <= mem_relax_lim,
                           "Memory relaxation limit reached for MUMPS factorization");
               if (myid == 0 && print_level > 0)
               {
                  mfem::out << "Re-running MUMPS factorization with memory relaxation "
                            << id->MUMPS_ICNTL(14) << '\n';
               }
            }
            else
            {
               MFEM_ABORT("Error during MUMPS numerical factorization");
            }
         }
         else { break; }
      }
   }

   hypre_CSRMatrixDestroy(csr_op);
   delete [] I;
   delete [] J;
   if (mat_type) { delete [] data; }

   id->nrhs = -1;  // Set up solution storage on first call to Mult
#if MFEM_MUMPS_VERSION >= 530
   delete [] irhs_loc;
   delete [] isol_loc;
   id->nloc_rhs = n_loc;
   id->lrhs_loc = n_loc;
   id->lsol_loc = id->MUMPS_INFO(23);
   irhs_loc = new int[id->lrhs_loc];
   isol_loc = new int[id->lsol_loc];
   for (int i = 0; i < n_loc; i++)
   {
      irhs_loc[i] = row_start + i + 1;
   }
   id->irhs_loc = irhs_loc;
   id->isol_loc = isol_loc;

   row_starts.SetSize(numProcs);
   MPI_Allgather(&row_start, 1, MPI_INT, row_starts, 1, MPI_INT, comm);
#else
   id->lrhs = id->n;
   if (myid == 0)
   {
      delete [] recv_counts;
      delete [] displs;
      recv_counts = new int[numProcs];
      displs = new int[numProcs];
   }
   MPI_Gather(&n_loc, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, comm);
   if (myid == 0)
   {
      displs[0] = 0;
      int s = 0;
      for (int k = 0; k < numProcs-1; k++)
      {
         s += recv_counts[k];
         displs[k+1] = s;
      }
   }
#endif
}

void MUMPSSolver::InitRhsSol(int nrhs) const
{
   if (id->nrhs != nrhs)
   {
#if MFEM_MUMPS_VERSION >= 530
      delete [] rhs_loc;
      delete [] sol_loc;
      rhs_loc = (nrhs > 1) ? new real_t[nrhs * id->lrhs_loc] : nullptr;
      sol_loc = new real_t[nrhs * id->lsol_loc];
      id->rhs_loc = rhs_loc;
      id->sol_loc = sol_loc;
#else
      if (myid == 0)
      {
         delete [] rhs_glob;
         rhs_glob = new real_t[nrhs * id->lrhs];
         id->rhs = rhs_glob;
      }
#endif
   }
   id->nrhs = nrhs;
}

void MUMPSSolver::Mult(const Vector &x, Vector &y) const
{
   Array<const Vector *> X(1);
   Array<Vector *> Y(1);
   X[0] = &x;
   Y[0] = &y;
   ArrayMult(X, Y);
}

void MUMPSSolver::ArrayMult(const Array<const Vector *> &X,
                            Array<Vector *> &Y) const
{
   MFEM_ASSERT(X.Size() == Y.Size(),
               "Number of columns mismatch in MUMPSSolver::Mult!");
   InitRhsSol(X.Size());
#if MFEM_MUMPS_VERSION >= 530
   if (id->nrhs == 1)
   {
      MFEM_ASSERT(X.Size() == 1 && X[0], "Missing Vector in MUMPSSolver::Mult!");
      X[0]->HostRead();
      id->rhs_loc = X[0]->GetData();
   }
   else
   {
      for (int i = 0; i < id->nrhs; i++)
      {
         MFEM_ASSERT(X[i], "Missing Vector in MUMPSSolver::Mult!");
         X[i]->HostRead();
         std::copy(X[i]->GetData(), X[i]->GetData() + X[i]->Size(),
                   id->rhs_loc + i * id->lrhs_loc);
      }
   }

   // MUMPS solve
   id->job = 3;
#ifdef MFEM_USE_SINGLE
   smumps_c(id);
#else
   dmumps_c(id);
#endif

   RedistributeSol(id->isol_loc, id->sol_loc, id->lsol_loc, Y);
#else
   for (int i = 0; i < id->nrhs; i++)
   {
      MFEM_ASSERT(X[i], "Missing Vector in MUMPSSolver::Mult!");
      X[i]->HostRead();
      MPI_Gatherv(X[i]->GetData(), X[i]->Size(), MPITypeMap<real_t>::mpi_type,
                  id->rhs + i * id->lrhs, recv_counts, displs, MPITypeMap<real_t>::mpi_type, 0,
                  comm);
   }

   // MUMPS solve
   id->job = 3;
#ifdef MFEM_USE_SINGLE
   smumps_c(id);
#else
   dmumps_c(id);
#endif

   for (int i = 0; i < id->nrhs; i++)
   {
      MFEM_ASSERT(Y[i], "Missing Vector in MUMPSSolver::Mult!");
      Y[i]->HostWrite();
      MPI_Scatterv(id->rhs + i * id->lrhs, recv_counts, displs,
                   MPITypeMap<real_t>::mpi_type,
                   Y[i]->GetData(), Y[i]->Size(), MPITypeMap<real_t>::mpi_type, 0, comm);
   }
#endif
}

void MUMPSSolver::MultTranspose(const Vector &x, Vector &y) const
{
   // Set flag for transpose solve
   id->MUMPS_ICNTL(9) = 0;
   Mult(x, y);

   // Reset the flag
   id->MUMPS_ICNTL(9) = 1;
}

void MUMPSSolver::ArrayMultTranspose(const Array<const Vector *> &X,
                                     Array<Vector *> &Y) const
{
   // Set flag for transpose solve
   id->MUMPS_ICNTL(9) = 0;
   ArrayMult(X, Y);

   // Reset the flag
   id->MUMPS_ICNTL(9) = 1;
}

void MUMPSSolver::SetPrintLevel(int print_lvl)
{
   print_level = print_lvl;
}

void MUMPSSolver::SetMatrixSymType(MatType mtype)
{
   mat_type = mtype;
}

void MUMPSSolver::SetReorderingStrategy(ReorderingStrategy method)
{
   reorder_method = method;
}

void MUMPSSolver::SetReorderingReuse(bool reuse)
{
   reorder_reuse = reuse;
}

#if MFEM_MUMPS_VERSION >= 510
void MUMPSSolver::SetBLRMode(int mode)
{
   blr_mode = mode;
}

void MUMPSSolver::SetBLRTol(double tol)
{
   blr_tol = tol;
}

void MUMPSSolver::SetBLRCompressionType(int type)
{
   blr_compression_type = type;
}
#endif

#if MFEM_MUMPS_VERSION >= 550
void MUMPSSolver::SetBLRCBCompression(int mode)
{
   blr_cb_compression = mode;
}
#endif

void MUMPSSolver::SetMemRelaxation(int pct)
{
   mem_relaxation = pct;
}

void MUMPSSolver::SetNumThreads(int n)
{
   num_threads = n;
}

void MUMPSSolver::SetPivotThreshold(double tol)
{
   pivot_threshold = tol;
}

void MUMPSSolver::SetOutOfCore(int mode)
{
   out_of_core = mode;
}

void MUMPSSolver::SetParameters()
{
   // Output stream for error messages
   id->MUMPS_ICNTL(1) = 6;
   // Output stream for diagnosting printing local to each proc
   id->MUMPS_ICNTL(2) = 0;
   // Output stream for global info
   id->MUMPS_ICNTL(3) = 6;
   // Level of error printing
   id->MUMPS_ICNTL(4) = print_level;
   // Input matrix format (assembled)
   id->MUMPS_ICNTL(5) = 0;
   // Use A or A^T
   id->MUMPS_ICNTL(9) = 1;
   // Iterative refinement (disabled)
   id->MUMPS_ICNTL(10) = 0;
   // Error analysis-statistics (disabled)
   id->MUMPS_ICNTL(11) = 0;
   // Use of ScaLAPACK (Parallel factorization on root)
   id->MUMPS_ICNTL(13) = 0;
   // Percentage increase of estimated workspace
   id->MUMPS_ICNTL(14) = mem_relaxation;
   // Number of OpenMP threads
   id->MUMPS_ICNTL(16) = num_threads;
   // Matrix input format (distributed)
   id->MUMPS_ICNTL(18) = 3;
   // Schur complement (no Schur complement matrix returned)
   id->MUMPS_ICNTL(19) = 0;
#if MFEM_MUMPS_VERSION >= 530
   // Distributed RHS
   id->MUMPS_ICNTL(20) = 10;
   // Distributed Sol
   id->MUMPS_ICNTL(21) = 1;
#else
   // Centralized RHS
   id->MUMPS_ICNTL(20) = 0;
   // Centralized Sol
   id->MUMPS_ICNTL(21) = 0;
#endif
   // Out of core factorization and solve
   id->MUMPS_ICNTL(22) = out_of_core;
   // Max size of working memory (default = based on estimates)
   id->MUMPS_ICNTL(23) = 0;
   // Configure reordering
   switch (reorder_method)
   {
      case ReorderingStrategy::AUTOMATIC:
         id->MUMPS_ICNTL(28) = 0;
         id->MUMPS_ICNTL(7) = 7;
         id->MUMPS_ICNTL(29) = 0;
         break;
      case ReorderingStrategy::AMD:
         id->MUMPS_ICNTL(28) = 1;
         id->MUMPS_ICNTL(7) = 0;
         break;
      case ReorderingStrategy::AMF:
         id->MUMPS_ICNTL(28) = 1;
         id->MUMPS_ICNTL(7) = 2;
         break;
      case ReorderingStrategy::PORD:
         id->MUMPS_ICNTL(28) = 1;
         id->MUMPS_ICNTL(7) = 4;
         break;
      case ReorderingStrategy::METIS:
         id->MUMPS_ICNTL(28) = 1;
         id->MUMPS_ICNTL(7) = 5;
         break;
      case ReorderingStrategy::PARMETIS:
         id->MUMPS_ICNTL(28) = 2;
         id->MUMPS_ICNTL(29) = 2;
         break;
      case ReorderingStrategy::SCOTCH:
         id->MUMPS_ICNTL(28) = 1;
         id->MUMPS_ICNTL(7) = 3;
         break;
      case ReorderingStrategy::PTSCOTCH:
         id->MUMPS_ICNTL(28) = 2;
         id->MUMPS_ICNTL(29) = 1;
         break;
      default:
         break; // This should be unreachable
   }
   // Pivot threshold
   if (pivot_threshold >= 0.0)
   {
      id->MUMPS_CNTL(1) = pivot_threshold;
   }
   // BLR options
#if MFEM_MUMPS_VERSION >= 510
   if (blr_mode > 0)
   {
      id->MUMPS_ICNTL(35) = blr_mode;
      id->MUMPS_ICNTL(36) = blr_compression_type;
      if (blr_tol > 0.0) { id->MUMPS_CNTL(7) = blr_tol; }
#if MFEM_MUMPS_VERSION >= 550
      id->MUMPS_ICNTL(37) = blr_cb_compression;
#endif
   }
#endif
}

#if MFEM_MUMPS_VERSION >= 530
int MUMPSSolver::GetRowRank(int i, const Array<int> &row_starts_) const
{
   if (row_starts_.Size() == 1)
   {
      return 0;
   }
   auto up = std::upper_bound(row_starts_.begin(), row_starts_.end(), i);
   return std::distance(row_starts_.begin(), up) - 1;
}

void MUMPSSolver::RedistributeSol(const int *rmap, const real_t *x,
                                  const int lx_loc, Array<Vector *> &Y) const
{
   int *send_count = new int[numProcs]();
   for (int i = 0; i < lx_loc; i++)
   {
      int j = rmap[i] - 1;
      int row_rank = GetRowRank(j, row_starts);
      if (myid == row_rank) { continue; }
      send_count[row_rank]++;
   }

   int *recv_count = new int[numProcs];
   MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

   int *send_displ = new int[numProcs]; send_displ[0] = 0;
   int *recv_displ = new int[numProcs]; recv_displ[0] = 0;
   int sbuff_size = send_count[numProcs-1];
   int rbuff_size = recv_count[numProcs-1];
   for (int k = 0; k < numProcs - 1; k++)
   {
      send_displ[k + 1] = send_displ[k] + send_count[k];
      recv_displ[k + 1] = recv_displ[k] + recv_count[k];
      sbuff_size += send_count[k];
      rbuff_size += recv_count[k];
   }

   int *sendbuf_index = new int[sbuff_size];
   real_t *sendbuf_values = new real_t[sbuff_size];
   int *recvbuf_index = new int[rbuff_size];
   real_t *recvbuf_values = new real_t[rbuff_size];
   int *soffs = new int[numProcs]();

   for (int i = 0; i < lx_loc; i++)
   {
      int j = rmap[i] - 1;
      int row_rank = GetRowRank(j, row_starts);
      if (myid != row_rank)
      {
         int k = send_displ[row_rank] + soffs[row_rank];
         sendbuf_index[k] = j;
         soffs[row_rank]++;
      }
   }

   MPI_Alltoallv(sendbuf_index, send_count, send_displ, MPI_INT,
                 recvbuf_index, recv_count, recv_displ, MPI_INT, comm);

   for (int rhs = 0; rhs < Y.Size(); rhs++)
   {
      MFEM_ASSERT(Y[rhs], "Missing Vector in MUMPSSolver::Mult!");
      Y[rhs]->HostWrite();

      std::fill(soffs, soffs + numProcs, 0);
      for (int i = 0; i < lx_loc; i++)
      {
         int j = rmap[i] - 1;
         int row_rank = GetRowRank(j, row_starts);
         if (myid == row_rank)
         {
            int local_index = j - row_start;
            (*Y[rhs])(local_index) = x[rhs * lx_loc + i];
         }
         else
         {
            int k = send_displ[row_rank] + soffs[row_rank];
            sendbuf_values[k] = x[rhs * lx_loc + i];
            soffs[row_rank]++;
         }
      }

      MPI_Alltoallv(sendbuf_values, send_count, send_displ,
                    MPITypeMap<real_t>::mpi_type,
                    recvbuf_values, recv_count, recv_displ, MPITypeMap<real_t>::mpi_type, comm);

      // Unpack recv buffer
      for (int i = 0; i < rbuff_size; i++)
      {
         int local_index = recvbuf_index[i] - row_start;
         (*Y[rhs])(local_index) = recvbuf_values[i];
      }
   }

   delete [] recvbuf_values;
   delete [] recvbuf_index;
   delete [] soffs;
   delete [] sendbuf_values;
   delete [] sendbuf_index;
   delete [] recv_displ;
   delete [] send_displ;
   delete [] recv_count;
   delete [] send_count;
}
#endif

// ============================================================================
// ComplexMUMPSSolver implementation
// ============================================================================
#ifdef MFEM_USE_COMPLEX_MUMPS
#if MFEM_MUMPS_VERSION < 550
#error ComplexMUMPSSolver requires MUMPS >= 5.5.0
#endif

// ---------------------------------------------------------------------------
// Portable helpers for MUMPS complex scalar types.
//
// ZMUMPS_COMPLEX / CMUMPS_COMPLEX is either a C struct {double/float r, i;}
// or the C99 type "double/float _Complex" depending on the MUMPS build.
// Both representations store two consecutive values of the base floating-point
// type, so reinterpret_cast between them and std::complex<T> is safe.
// ---------------------------------------------------------------------------

#ifdef MFEM_USE_SINGLE
// ---------- single-precision helpers ----------------------------------------
static inline void CMumpsSet(CMUMPS_COMPLEX &z, float re, float im)
{
   float *p = reinterpret_cast<float *>(&z);
   p[0] = re;  p[1] = im;
}
static inline float CMumpsReal(const CMUMPS_COMPLEX &z)
{ return reinterpret_cast<const float *>(&z)[0]; }
static inline float CMumpsImag(const CMUMPS_COMPLEX &z)
{ return reinterpret_cast<const float *>(&z)[1]; }
#define ComplexMumpsSet(z,re,im) CMumpsSet(z, static_cast<float>(re), static_cast<float>(im))
#define ComplexMumpsReal(z)      CMumpsReal(z)
#define ComplexMumpsImag(z)      CMumpsImag(z)
#define COMPLEX_MUMPS_CALL       cmumps_c
#else
// ---------- double-precision helpers ----------------------------------------
static inline void ZMumpsSet(ZMUMPS_COMPLEX &z, double re, double im)
{
   double *p = reinterpret_cast<double *>(&z);
   p[0] = re;  p[1] = im;
}
static inline double ZMumpsReal(const ZMUMPS_COMPLEX &z)
{ return reinterpret_cast<const double *>(&z)[0]; }
static inline double ZMumpsImag(const ZMUMPS_COMPLEX &z)
{ return reinterpret_cast<const double *>(&z)[1]; }
#define ComplexMumpsSet(z,re,im) ZMumpsSet(z, re, im)
#define ComplexMumpsReal(z)      ZMumpsReal(z)
#define ComplexMumpsImag(z)      ZMumpsImag(z)
#define COMPLEX_MUMPS_CALL       zmumps_c
#endif

void ComplexMUMPSSolver::Init(MPI_Comm comm_)
{
   id = nullptr;
   comm = comm_;
   MPI_Comm_size(comm, &numProcs);
   MPI_Comm_rank(comm, &myid);

   mat_type      = MatType::UNSYMMETRIC;
   print_level   = 0;
   reorder_method = ReorderingStrategy::AUTOMATIC;
   mem_relaxation    = 20;
   num_threads       = 0;
   pivot_threshold   = -1.0;
   out_of_core       = 0;
   reorder_reuse     = false;
   analysis_done_    = false;
   blr_mode = 0;
   blr_tol = 0.0;
   blr_compression_type = 0;
   blr_cb_compression = 0;

   I_coo = nullptr;  J_coo = nullptr;  data_coo = nullptr;
   nnz_loc_total = 0;
   row_start = 0;  m_loc = 0;  n_global = 0;

   irhs_loc = nullptr;  isol_loc  = nullptr;
   rhs_loc_buf = nullptr;  sol_loc_buf = nullptr;
   lrhs_loc = 0;  lsol_loc = 0;
   nrhs_cur_ = 0;
}

ComplexMUMPSSolver::ComplexMUMPSSolver(MPI_Comm comm_)
{
   Init(comm_);
}

ComplexMUMPSSolver::~ComplexMUMPSSolver()
{
   delete [] I_coo;
   delete [] J_coo;
   delete [] data_coo;

   delete [] irhs_loc;
   delete [] isol_loc;
   delete [] rhs_loc_buf;
   delete [] sol_loc_buf;

   if (id)
   {
      id->job = -2;
      COMPLEX_MUMPS_CALL(id);
      delete id;
   }
}

void ComplexMUMPSSolver::SetPrintLevel(int print_lvl)
{
   print_level = print_lvl;
}

void ComplexMUMPSSolver::SetMatrixSymType(MatType mtype)
{
   mat_type = mtype;
}

void ComplexMUMPSSolver::SetReorderingStrategy(ReorderingStrategy method)
{
   reorder_method = method;
}

void ComplexMUMPSSolver::SetBLRMode(int mode)
{
   blr_mode = mode;
}

void ComplexMUMPSSolver::SetBLRTol(double tol)
{
   blr_tol = tol;
}

void ComplexMUMPSSolver::SetBLRCompressionType(int type)
{
   blr_compression_type = type;
}

void ComplexMUMPSSolver::SetBLRCBCompression(int mode)
{
   blr_cb_compression = mode;
}

void ComplexMUMPSSolver::SetMemRelaxation(int pct)
{
   mem_relaxation = pct;
}

void ComplexMUMPSSolver::SetNumThreads(int n)
{
   num_threads = n;
}

void ComplexMUMPSSolver::SetPivotThreshold(double tol)
{
   pivot_threshold = tol;
}

void ComplexMUMPSSolver::SetOutOfCore(int mode)
{
   out_of_core = mode;
}

void ComplexMUMPSSolver::SetReorderingReuse(bool reuse)
{
   reorder_reuse = reuse;
}

void ComplexMUMPSSolver::SetParameters()
{
   id->MUMPS_ICNTL(1) = 6;    // error stream
   id->MUMPS_ICNTL(2) = 0;    // diagnostic stream (per proc)
   id->MUMPS_ICNTL(3) = 6;    // global info stream
   id->MUMPS_ICNTL(4) = print_level;
   id->MUMPS_ICNTL(5) = 0;    // assembled input (duplicates are summed)
   id->MUMPS_ICNTL(9) = 1;    // solve A x = b
   id->MUMPS_ICNTL(10) = 0;   // no iterative refinement
   id->MUMPS_ICNTL(11) = 0;   // no error statistics
   id->MUMPS_ICNTL(13) = 0;   // ScaLAPACK on root
   id->MUMPS_ICNTL(14) = mem_relaxation;  // workspace % above estimate
   id->MUMPS_ICNTL(16) = num_threads;     // OpenMP threads per rank
   id->MUMPS_ICNTL(18) = 3;   // distributed assembled input
   id->MUMPS_ICNTL(19) = 0;   // no Schur complement

   id->MUMPS_ICNTL(20) = 10;  // distributed RHS
   id->MUMPS_ICNTL(21) = 1;   // distributed solution
   id->MUMPS_ICNTL(22) = out_of_core;  // 0=in-core, 1=out-of-core
   id->MUMPS_ICNTL(23) = 0;   // default working memory

   if (pivot_threshold >= 0.0) { id->MUMPS_CNTL(1) = pivot_threshold; }

   // Reordering
   switch (reorder_method)
   {
      case ReorderingStrategy::AUTOMATIC:
         id->MUMPS_ICNTL(28) = 0;
         id->MUMPS_ICNTL(7)  = 7;
         id->MUMPS_ICNTL(29) = 0;
         break;
      case ReorderingStrategy::AMD:
         id->MUMPS_ICNTL(28) = 1;
         id->MUMPS_ICNTL(7)  = 0;
         break;
      case ReorderingStrategy::AMF:
         id->MUMPS_ICNTL(28) = 1;
         id->MUMPS_ICNTL(7)  = 2;
         break;
      case ReorderingStrategy::PORD:
         id->MUMPS_ICNTL(28) = 1;
         id->MUMPS_ICNTL(7)  = 4;
         break;
      case ReorderingStrategy::METIS:
         id->MUMPS_ICNTL(28) = 1;
         id->MUMPS_ICNTL(7)  = 5;
         break;
      case ReorderingStrategy::PARMETIS:
         id->MUMPS_ICNTL(28) = 2;
         id->MUMPS_ICNTL(29) = 2;
         break;
      case ReorderingStrategy::SCOTCH:
         id->MUMPS_ICNTL(28) = 1;
         id->MUMPS_ICNTL(7)  = 3;
         break;
      case ReorderingStrategy::PTSCOTCH:
         id->MUMPS_ICNTL(28) = 2;
         id->MUMPS_ICNTL(29) = 1;
         break;
      default:
         break;
   }

   if (blr_mode > 0)
   {
      id->MUMPS_ICNTL(35) = blr_mode;
      id->MUMPS_ICNTL(36) = blr_compression_type;
      id->MUMPS_ICNTL(37) = blr_cb_compression;
      if (blr_tol > 0.0) { id->MUMPS_CNTL(7) = blr_tol; }
   }
}

void ComplexMUMPSSolver::SetOperator(const Operator &op)
{
   const auto *chpm = dynamic_cast<const ComplexHypreParMatrix *>(&op);
   MFEM_VERIFY(chpm,
               "ComplexMUMPSSolver::SetOperator requires a ComplexHypreParMatrix.");
   SetOperator(*chpm);
}

void ComplexMUMPSSolver::SetOperator(const ComplexHypreParMatrix &op)
{
   height = op.Height();
   width  = op.Width();

   const HypreParMatrix &A_r = op.real();
   const HypreParMatrix &A_i = op.imag();

   // ---- Merge local CSR from A_r ----------------------------------------
   auto *parcsr_r =
      static_cast<hypre_ParCSRMatrix *>(
         const_cast<HypreParMatrix &>(A_r));
   const_cast<HypreParMatrix &>(A_r).HostRead();
   hypre_CSRMatrix *csr_r = hypre_MergeDiagAndOffd(parcsr_r);
   const_cast<HypreParMatrix &>(A_r).HypreRead();
#if MFEM_HYPRE_VERSION >= 21600
   hypre_CSRMatrixBigJtoJ(csr_r);
#endif

   // ---- Merge local CSR from A_i ----------------------------------------
   auto *parcsr_i =
      static_cast<hypre_ParCSRMatrix *>(
         const_cast<HypreParMatrix &>(A_i));
   const_cast<HypreParMatrix &>(A_i).HostRead();
   hypre_CSRMatrix *csr_i = hypre_MergeDiagAndOffd(parcsr_i);
   const_cast<HypreParMatrix &>(A_i).HypreRead();
#if MFEM_HYPRE_VERSION >= 21600
   hypre_CSRMatrixBigJtoJ(csr_i);
#endif

   m_loc      = static_cast<int>(csr_r->num_rows);
   row_start  = static_cast<int>(parcsr_r->first_row_index);
   n_global   = static_cast<int>(parcsr_r->global_num_rows);

   const int nnz_r = csr_r->num_nonzeros;
   const int nnz_i = csr_i->num_nonzeros;

   // ---- For symmetric modes keep only the lower triangle ----------------
   // Lower-triangle filter: keep entry (row, col) where global_row >= col.
   int nnz_r_use = 0, nnz_i_use = 0;
   if (mat_type != UNSYMMETRIC)
   {
      for (int i = 0; i < m_loc; i++)
      {
         const int grow = row_start + i + 1; // 1-based
         for (int j = csr_r->i[i]; j < csr_r->i[i + 1]; j++)
         {
            if (grow >= static_cast<int>(csr_r->j[j]) + 1) { nnz_r_use++; }
         }
         for (int j = csr_i->i[i]; j < csr_i->i[i + 1]; j++)
         {
            if (grow >= static_cast<int>(csr_i->j[j]) + 1) { nnz_i_use++; }
         }
      }
   }
   else
   {
      nnz_r_use = nnz_r;
      nnz_i_use = nnz_i;
   }

   nnz_loc_total = nnz_r_use + nnz_i_use;

   delete [] I_coo;
   delete [] J_coo;
   delete [] data_coo;
   I_coo    = new int[nnz_loc_total];
   J_coo    = new int[nnz_loc_total];
   data_coo = new ComplexMumpsScalar[nnz_loc_total];

   // ---- Fill COO: real entries first, then imaginary entries ------------
   // MUMPS sums duplicate (row,col) pairs, yielding (D_r, 0) + (0, D_i)
   // = (D_r, D_i) for coincident sparsity entries.
   int k = 0;
   for (int i = 0; i < m_loc; i++)
   {
      const int grow = row_start + i + 1; // 1-based
      for (int j = csr_r->i[i]; j < csr_r->i[i + 1]; j++)
      {
         const int gcol = static_cast<int>(csr_r->j[j]) + 1; // 1-based
         if (mat_type != UNSYMMETRIC && grow < gcol) { continue; }
         I_coo[k] = grow;
         J_coo[k] = gcol;
         ComplexMumpsSet(data_coo[k], csr_r->data[j], 0.0);
         k++;
      }
   }
   for (int i = 0; i < m_loc; i++)
   {
      const int grow = row_start + i + 1;
      for (int j = csr_i->i[i]; j < csr_i->i[i + 1]; j++)
      {
         const int gcol = static_cast<int>(csr_i->j[j]) + 1;
         if (mat_type != UNSYMMETRIC && grow < gcol) { continue; }
         I_coo[k] = grow;
         J_coo[k] = gcol;
         ComplexMumpsSet(data_coo[k], 0.0, csr_i->data[j]);
         k++;
      }
   }

   hypre_CSRMatrixDestroy(csr_r);
   hypre_CSRMatrixDestroy(csr_i);

   // ---- Initialise or reset the MUMPS object ----------------------------
   if (!id || !reorder_reuse)
   {
      if (id)
      {
         id->job = -2;
         COMPLEX_MUMPS_CALL(id);
         delete id;
      }
      id = new ComplexMumpsStruc();
      id->sym = static_cast<int>(mat_type);
      id->comm_fortran = static_cast<MUMPS_INT>(MPI_Comm_c2f(comm));
      id->par = 1;   // host participates in computation

      id->job = -1;  // initialise
      COMPLEX_MUMPS_CALL(id);

      SetParameters();

      id->n       = n_global;
      id->nnz_loc = nnz_loc_total;
      id->irn_loc = I_coo;
      id->jcn_loc = J_coo;
      id->a_loc   = data_coo;

      // ---- Symbolic analysis (phase 1) -------------------------------------
      id->job = 1;
      COMPLEX_MUMPS_CALL(id);
      analysis_done_ = true;
   }
   else
   {
      // Reuse existing symbolic analysis: just update data pointers.
      id->nnz_loc = nnz_loc_total;
      id->irn_loc = I_coo;
      id->jcn_loc = J_coo;
      id->a_loc   = data_coo;
   }

   // ---- Numeric factorisation (phase 2) with memory-relaxation retry ----
   id->job = 2;
   {
      const int mem_relax_lim = 200;
      while (true)
      {
         COMPLEX_MUMPS_CALL(id);
         if (id->MUMPS_INFOG(1) < 0)
         {
            if (id->MUMPS_INFOG(1) == -8 || id->MUMPS_INFOG(1) == -9)
            {
               id->MUMPS_ICNTL(14) += 20;
               MFEM_VERIFY(id->MUMPS_ICNTL(14) <= mem_relax_lim,
                           "ComplexMUMPS: memory relaxation limit reached.");
               if (myid == 0 && print_level > 0)
               {
                  mfem::out << "ComplexMUMPS: re-factorising with memory "
                            << "relaxation " << id->MUMPS_ICNTL(14) << '\n';
               }
            }
            else
            {
               MFEM_ABORT("ComplexMUMPS: error during numeric factorisation.");
            }
         }
         else { break; }
      }
   }

   // ---- Allocate RHS / solution index arrays ---------------------------
   delete [] irhs_loc;
   delete [] isol_loc;
   delete [] rhs_loc_buf;
   delete [] sol_loc_buf;
   rhs_loc_buf = nullptr;
   sol_loc_buf = nullptr;
   nrhs_cur_ = 0;

   lrhs_loc = m_loc;
   lsol_loc = id->MUMPS_INFO(23);

   irhs_loc = new int[lrhs_loc];
   isol_loc = new int[lsol_loc];

   for (int i = 0; i < m_loc; i++) { irhs_loc[i] = row_start + i + 1; }

   id->nloc_rhs = m_loc;
   id->lrhs_loc = lrhs_loc;
   id->lsol_loc = lsol_loc;
   id->irhs_loc = irhs_loc;
   id->isol_loc = isol_loc;

   row_starts.SetSize(numProcs);
   MPI_Allgather(&row_start, 1, MPI_INT, row_starts, 1, MPI_INT, comm);
}

void ComplexMUMPSSolver::InitRhsSol(int nrhs) const
{
   if (nrhs != nrhs_cur_)
   {
      delete [] rhs_loc_buf;
      delete [] sol_loc_buf;
      rhs_loc_buf = new ComplexMumpsScalar[nrhs * lrhs_loc];
      sol_loc_buf = new ComplexMumpsScalar[nrhs * lsol_loc];
      nrhs_cur_   = nrhs;
   }
   id->nrhs = nrhs;
}

void ComplexMUMPSSolver::Mult(const Vector &b_r, const Vector &b_i,
                              Vector &x_r, Vector &x_i) const
{
   MFEM_VERIFY(id, "ComplexMUMPSSolver: SetOperator has not been called.");
   MFEM_VERIFY(b_r.Size() == m_loc && b_i.Size() == m_loc,
               "ComplexMUMPSSolver::Mult: RHS size mismatch.");

   x_r.SetSize(m_loc);
   x_i.SetSize(m_loc);

   InitRhsSol(1);

   // ---- Pack local complex RHS ------------------------------------------
   for (int i = 0; i < m_loc; i++)
   {
      ComplexMumpsSet(rhs_loc_buf[i], b_r[i], b_i[i]);
   }
   id->rhs_loc = rhs_loc_buf;
   id->sol_loc = sol_loc_buf;

   id->job = 3;
   COMPLEX_MUMPS_CALL(id);

   MFEM_VERIFY(id->MUMPS_INFOG(1) >= 0,
               "ComplexMUMPS: error during solve, INFOG(1)=" <<
               id->MUMPS_INFOG(1));

   RedistributeSol(isol_loc, sol_loc_buf, lsol_loc, x_r, x_i);
}

void ComplexMUMPSSolver::Mult(const Vector &b, Vector &x) const
{
   MFEM_VERIFY(id, "ComplexMUMPSSolver: SetOperator has not been called.");
   const int n = b.Size() / 2;
   MFEM_VERIFY(b.Size() == 2 * n,
               "ComplexMUMPSSolver::Mult: b must have even size (2*N).");

   // View the first / second halves as real / imaginary parts.
   Vector b_r(const_cast<real_t *>(b.GetData()),     n);
   Vector b_i(const_cast<real_t *>(b.GetData()) + n, n);

   x.SetSize(2 * n);
   Vector x_r(x.GetData(),     n);
   Vector x_i(x.GetData() + n, n);

   Mult(b_r, b_i, x_r, x_i);
}

void ComplexMUMPSSolver::ArrayMult(const Array<const Vector *> &X,
                                   Array<Vector *> &Y) const
{
   MFEM_VERIFY(id, "ComplexMUMPSSolver: SetOperator has not been called.");
   MFEM_ASSERT(X.Size() == Y.Size(),
               "Number of columns mismatch in ComplexMUMPSSolver::ArrayMult!");

   const int nrhs = X.Size();
   InitRhsSol(nrhs);

   // ---- Pack: each X[k] has layout [real(0..m-1) ; imag(0..m-1)] --------
   for (int k = 0; k < nrhs; k++)
   {
      MFEM_ASSERT(X[k] && X[k]->Size() == 2 * m_loc,
                  "ComplexMUMPSSolver::ArrayMult: RHS size mismatch.");
      X[k]->HostRead();
      const real_t *b_r = X[k]->GetData();
      const real_t *b_i = X[k]->GetData() + m_loc;
      for (int i = 0; i < m_loc; i++)
      {
         ComplexMumpsSet(rhs_loc_buf[k * lrhs_loc + i], b_r[i], b_i[i]);
      }
   }

   id->rhs_loc = rhs_loc_buf;
   id->sol_loc = sol_loc_buf;
   id->job = 3;
   COMPLEX_MUMPS_CALL(id);

   MFEM_VERIFY(id->MUMPS_INFOG(1) >= 0,
               "ComplexMUMPS: error during solve, INFOG(1)=" <<
               id->MUMPS_INFOG(1));

   // ---- Unpack: solution column k starts at sol_loc_buf + k*lsol_loc ----
   for (int k = 0; k < nrhs; k++)
   {
      Y[k]->SetSize(2 * m_loc);
      Y[k]->HostWrite();
      Vector x_r(Y[k]->GetData(),          m_loc);
      Vector x_i(Y[k]->GetData() + m_loc,  m_loc);
      RedistributeSol(isol_loc, sol_loc_buf + k * lsol_loc, lsol_loc, x_r, x_i);
   }
}

int ComplexMUMPSSolver::GetRowRank(int i, const Array<int> &row_starts_) const
{
   if (row_starts_.Size() == 1) { return 0; }
   auto up = std::upper_bound(row_starts_.begin(), row_starts_.end(), i);
   return static_cast<int>(std::distance(row_starts_.begin(), up)) - 1;
}

void ComplexMUMPSSolver::RedistributeSol(const int *rmap,
                                         const ComplexMumpsScalar *x,
                                         const int lx_loc,
                                         Vector &yr, Vector &yi) const
{
   // Identify how many entries each rank needs to send to every other rank.
   int *send_count = new int[numProcs]();
   for (int i = 0; i < lx_loc; i++)
   {
      const int j         = rmap[i] - 1;
      const int row_rank  = GetRowRank(j, row_starts);
      if (myid != row_rank) { send_count[row_rank]++; }
   }

   int *recv_count = new int[numProcs];
   MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

   int *send_displ = new int[numProcs]; send_displ[0] = 0;
   int *recv_displ = new int[numProcs]; recv_displ[0] = 0;
   int sbuff_size = send_count[numProcs - 1];
   int rbuff_size = recv_count[numProcs - 1];
   for (int k = 0; k < numProcs - 1; k++)
   {
      send_displ[k + 1] = send_displ[k] + send_count[k];
      recv_displ[k + 1] = recv_displ[k] + recv_count[k];
      sbuff_size += send_count[k];
      rbuff_size += recv_count[k];
   }

   int *sendbuf_index        = new int[sbuff_size];
   ComplexMumpsScalar *sendbuf_val = new ComplexMumpsScalar[sbuff_size];
   int *recvbuf_index        = new int[rbuff_size];
   ComplexMumpsScalar *recvbuf_val = new ComplexMumpsScalar[rbuff_size];
   int *soffs = new int[numProcs]();

   // Build index send buffer.
   for (int i = 0; i < lx_loc; i++)
   {
      const int j        = rmap[i] - 1;
      const int row_rank = GetRowRank(j, row_starts);
      if (myid != row_rank)
      {
         const int k = send_displ[row_rank] + soffs[row_rank];
         sendbuf_index[k] = j;
         soffs[row_rank]++;
      }
   }
   MPI_Alltoallv(sendbuf_index, send_count, send_displ, MPI_INT,
                 recvbuf_index, recv_count, recv_displ, MPI_INT, comm);

   // Build value send buffer (one pass per RHS; here we have exactly one).
   yr.HostWrite();
   yi.HostWrite();
   std::fill(soffs, soffs + numProcs, 0);
   for (int i = 0; i < lx_loc; i++)
   {
      const int j        = rmap[i] - 1;
      const int row_rank = GetRowRank(j, row_starts);
      if (myid == row_rank)
      {
         const int local_index = j - row_start;
         yr[local_index] = ComplexMumpsReal(x[i]);
         yi[local_index] = ComplexMumpsImag(x[i]);
      }
      else
      {
         const int k = send_displ[row_rank] + soffs[row_rank];
         sendbuf_val[k] = x[i];
         soffs[row_rank]++;
      }
   }

   // Exchange complex values.  We send pairs of doubles – use a byte-based
   // MPI type to avoid creating a custom MPI datatype.
   const int bytes_per_elem = static_cast<int>(sizeof(ComplexMumpsScalar));
   static_assert(sizeof(ComplexMumpsScalar) == 2 * sizeof(double),
                 "ComplexMumpsScalar must be two consecutive doubles.");

   // Adjust counts/displs to byte-pairs (factor 2 doubles = 16 bytes each).
   // We use MPI_DOUBLE with doubled counts for portability.
   int *send_count2 = new int[numProcs];
   int *send_displ2 = new int[numProcs];
   int *recv_count2 = new int[numProcs];
   int *recv_displ2 = new int[numProcs];
   for (int k = 0; k < numProcs; k++)
   {
      send_count2[k] = send_count[k] * 2;
      send_displ2[k] = send_displ[k] * 2;
      recv_count2[k] = recv_count[k] * 2;
      recv_displ2[k] = recv_displ[k] * 2;
   }
   MPI_Alltoallv(reinterpret_cast<double *>(sendbuf_val),
                 send_count2, send_displ2, MPI_DOUBLE,
                 reinterpret_cast<double *>(recvbuf_val),
                 recv_count2, recv_displ2, MPI_DOUBLE, comm);

   // Unpack received values.
   for (int i = 0; i < rbuff_size; i++)
   {
      const int local_index = recvbuf_index[i] - row_start;
      yr[local_index] = ComplexMumpsReal(recvbuf_val[i]);
      yi[local_index] = ComplexMumpsImag(recvbuf_val[i]);
   }

   delete [] recvbuf_val;
   delete [] recvbuf_index;
   delete [] soffs;
   delete [] sendbuf_val;
   delete [] sendbuf_index;
   delete [] send_count2; delete [] send_displ2;
   delete [] recv_count2; delete [] recv_displ2;
   delete [] recv_displ;
   delete [] send_displ;
   delete [] recv_count;
   delete [] send_count;
}

#endif // MFEM_USE_COMPLEX_MUMPS

} // namespace mfem

#endif // MFEM_USE_MPI
#endif // MFEM_USE_MUMPS
