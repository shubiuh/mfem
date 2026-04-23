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

#ifndef MFEM_CPARDISO
#define MFEM_CPARDISO

#include "../config/config.hpp"

#ifdef MFEM_USE_MPI
#ifdef MFEM_USE_MKL_CPARDISO

#include "mkl_cluster_sparse_solver.h"
#include "operator.hpp"
#include "densemat.hpp"
#include "complex_operator.hpp"

namespace mfem
{
/**
 * @brief MKL Parallel Direct Sparse Solver for Clusters
 *
 * Interface to MKL CPardiso: the MPI-enabled Intel MKL version of Pardiso
 */
class CPardisoSolver : public Solver
{
public:
   enum MatType
   {
      REAL_STRUCTURE_SYMMETRIC = 1,
      REAL_SYMMETRIC_POSITIVE_DEFINITE = 2,
      REAL_SYMMETRIC_INDEFINITE = -2,
      REAL_NONSYMMETRIC = 11,

      COMPLEX_STRUCTURALLY_SYMMETRIC = 3,      ///< Complex structurally symmetric
      COMPLEX_HERMITIAN_POSITIVE_DEFINITE = 4, ///< Complex Hermitian positive definite
      COMPLEX_HERMITIAN_INDEFINITE = -4,       ///< Complex Hermitian indefinite
      COMPLEX_SYMMETRIC = 6,                   ///< Complex symmetric
      COMPLEX_NONSYMMETRIC = 13                ///< Complex nonsymmetric
   };

   /**
    * @brief Construct a new CPardisoSolver object
    *
    * @param comm MPI Communicator
    * @param nrhs_ Number of right hand sides, default is 1
    */
   CPardisoSolver(MPI_Comm comm, int nrhs_ = 1);

   /**
    * @brief Set the Operator object and perform factorization
    *
    * @a op needs to be of type HypreParMatrix. The contents are copied and
    * reordered in an internal CSR structure.
    *
    * @param op Operator to use in factorization and solve
    */
   void SetOperator(const Operator &op) override;

   /**
    * @brief Solve
    *
    * @param b RHS vector
    * @param x Solution vector
    */
   void Mult(const Vector &b, Vector &x) const override;

   /**
    * @brief Solve for multiple right-hand sides
    *
    * @param B RHS matrix
    * @param X Solution matrix
    */
   void Mult(const DenseMatrix &B, DenseMatrix &X) const;

   /**
    * @brief Set the print level for MKL CPardiso
    *
    * Prints statistics after the factorization and after each solve.
    *
    * @param print_lvl Print level
    */
   void SetPrintLevel(int print_lvl);

   /**
    * @brief Set the number of right-hand sides
    *
    * @param nrhs_ Number of right-hand sides
    */
   void setRHSCount(int nrhs_);

   /**
    * @brief Set the matrix type
    *
    * The matrix type supported is either real and symmetric or real and
    * non-symmetric.
    *
    * @param mat_type Matrix type
    */
   void SetMatrixType(MatType mat_type);

   ~CPardisoSolver();

private:
   MPI_Fint comm_;

   // Global number of rows
   int m;

   // First row index of the global matrix on the local MPI rank
   int first_row;

   // Local number of nonzero entries
   int nnz_loc;

   // Local number of rows, obtained from a ParCSR matrix
   int m_loc;

   // CSR data structure for the copy data of the local CSR matrix
   int *csr_rowptr = nullptr;
   real_t *reordered_csr_nzval = nullptr;
   int *reordered_csr_colind = nullptr;

   // Internal solver memory pointer pt,
   // 32-bit: int pt[64]
   // 64-bit: long int pt[64] or void *pt[64] should be OK on both architectures
   mutable void *pt[64] = {0};

   // Solver control parameters, detailed description can be found in the
   // constructor.
   mutable int iparm[64] = {0};
   mutable int maxfct, mnum, msglvl, phase, error;
   int mtype;
   int nrhs;

   // Dummy variables
   mutable int idum;
   mutable real_t ddum;
};
/// MKL complex scalar type (single or double precision).
#ifdef MFEM_USE_SINGLE
using CPardisoComplexScalar = MKL_Complex8;
#else
using CPardisoComplexScalar = MKL_Complex16;
#endif

/**
 * @brief MKL Cluster Sparse Solver for complex-valued distributed systems.
 *
 * Solves  (A_r + i A_i)(x_r + i x_i) = (b_r + i b_i)  where the operator is
 * provided as a ComplexHypreParMatrix.  Real and imaginary CSR blocks are
 * merged into a single complex CSR and passed natively to cluster_sparse_solver
 * using a complex mtype (no real 2x2 block expansion).
 *
 * Supported matrix types:
 *   - COMPLEX_NONSYMMETRIC (default): full matrix assembled.
 *   - COMPLEX_SYMMETRIC: lower triangular part only (A = A^T).
 *   - COMPLEX_HERMITIAN_POSITIVE_DEFINITE: lower triangular (A = A^H, SPD).
 *   - COMPLEX_HERMITIAN_INDEFINITE: lower triangular (A = A^H, indefinite).
 */
class ComplexCPardisoSolver : public Solver
{
public:
   enum MatType
   {
      COMPLEX_HERMITIAN_POSITIVE_DEFINITE = 4,  ///< Complex Hermitian positive definite
      COMPLEX_HERMITIAN_INDEFINITE        = -4, ///< Complex Hermitian indefinite
      COMPLEX_SYMMETRIC                   = 6,  ///< Complex symmetric (A = A^T)
      COMPLEX_NONSYMMETRIC                = 13  ///< Complex nonsymmetric (default)
   };

   /**
    * @brief Construct a new ComplexCPardisoSolver.
    *
    * @param comm  MPI communicator
    * @param nrhs_ Number of right-hand sides (default 1)
    */
   ComplexCPardisoSolver(MPI_Comm comm, int nrhs_ = 1);

   /**
    * @brief Set the operator and perform symbolic + numeric factorization.
    *
    * @a op must be of type ComplexHypreParMatrix.
    */
   void SetOperator(const Operator &op) override;

   /// Overload accepting a ComplexHypreParMatrix directly.
   void SetOperator(const ComplexHypreParMatrix &op);

   /**
    * @brief Solve  A (x_r + i x_i) = (b_r + i b_i)  with explicit real/imag
    *        vectors.  Each vector has local size m_loc.
    */
   void Mult(const Vector &b_r, const Vector &b_i,
             Vector &x_r, Vector &x_i) const;

   /**
    * @brief Solve using combined [real ; imag] block vectors of size 2*m_loc.
    *
    * Layout: b = [b_r_0 ... b_r_{m-1}  b_i_0 ... b_i_{m-1}].
    */
   void Mult(const Vector &b, Vector &x) const override;

   /**
    * @brief Set the print level for MKL CPardiso.
    *
    * @param print_lvl  Print level (0 = silent, 1 = statistics).
    */
   void SetPrintLevel(int print_lvl);

   /**
    * @brief Set the number of right-hand sides.
    *
    * @param nrhs_ Number of right-hand sides.
    */
   void setRHSCount(int nrhs_);

   /**
    * @brief Set the matrix type.
    *
    * Must be called before SetOperator.
    *
    * @param mat_type  One of the ComplexCPardisoSolver::MatType values.
    */
   void SetMatrixType(MatType mat_type);

   ~ComplexCPardisoSolver();

private:
   MPI_Fint comm_;

   /// Global number of rows.
   int m;

   /// First row index of the global matrix on this MPI rank (0-based).
   int first_row;

   /// Number of nonzeros in the local merged complex CSR.
   int nnz_loc;

   /// Number of local rows.
   int m_loc;

   /// CSR data for the local complex matrix.
   int                    *csr_rowptr           = nullptr;
   CPardisoComplexScalar  *reordered_csr_nzval  = nullptr;
   int                    *reordered_csr_colind  = nullptr;

   // Internal solver memory pointer (64 entries, see MKL docs).
   mutable void *pt[64] = {0};

   // Solver control parameters.
   mutable int iparm[64] = {0};
   mutable int maxfct, mnum, msglvl, phase, error;
   int mtype;
   int nrhs;

   // Dummy variables used for unused b/x arguments during analysis/factorization.
   mutable int                   idum;
   mutable CPardisoComplexScalar ddum_cpx;
};

} // namespace mfem

#endif
#endif // MFEM_USE_MKL_CPARDISO
#endif // MFEM_USE_MPI
