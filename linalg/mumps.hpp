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

#ifndef MFEM_MUMPS
#define MFEM_MUMPS

#include "../config/config.hpp"

#ifdef MFEM_USE_MUMPS
#ifdef MFEM_USE_MPI

#include "operator.hpp"
#include "hypre.hpp"
#include <mpi.h>

#ifdef MFEM_USE_SINGLE
#include "smumps_c.h"
#else
#include "dmumps_c.h"
#endif

#ifdef MFEM_USE_COMPLEX_MUMPS
#ifdef MFEM_USE_SINGLE
#include "cmumps_c.h"
#else
#include "zmumps_c.h"
#endif
#include "complex_operator.hpp"
#endif // MFEM_USE_COMPLEX_MUMPS

namespace mfem
{

/**
 * @brief MUMPS: A Parallel Sparse Direct Solver
 *
 * Interface for the distributed MUMPS solver
 */
class MUMPSSolver : public Solver
{
public:
   /// Specify the type of matrix we are applying the solver to
   enum MatType
   {
      /// General sparse matrix, no symmetry is assumed
      UNSYMMETRIC = 0,
      /// A sparse symmetric positive definite matrix
      SYMMETRIC_POSITIVE_DEFINITE = 1,
      /// A sparse symmetric matrix that is not necessarily positive definite
      SYMMETRIC_INDEFINITE = 2
   };

   /// Specify the reordering strategy for the MUMPS solver
   enum ReorderingStrategy
   {
      /// Let MUMPS automatically decide the reording strategy
      AUTOMATIC = 0,
      /// Approximate Minimum Degree with auto quasi-dense row detection is used
      AMD,
      /// Approximate Minimum Fill method will be used
      AMF,
      /// The PORD library will be used
      PORD,
      /// The METIS library will be used
      METIS,
      /// The ParMETIS library will be used
      PARMETIS,
      /// The Scotch library will be used
      SCOTCH,
      /// The PTScotch library will be used
      PTSCOTCH
   };

   /**
    * @brief Constructor with MPI_Comm parameter.
    */
   MUMPSSolver(MPI_Comm comm_);

   /**
    * @brief Constructor with a HypreParMatrix Operator.
    */
   MUMPSSolver(const Operator &op);

   /**
    * @brief Set the Operator and perform factorization
    *
    * @a op needs to be of type HypreParMatrix.
    *
    * @param op Operator used in factorization and solve
    */
   void SetOperator(const Operator &op);

   /**
    * @brief Solve $ y = Op^{-1} x $
    *
    * @param x RHS vector
    * @param y Solution vector
    */
   void Mult(const Vector &x, Vector &y) const;

   /**
    * @brief Solve $ Y_i = Op^{-T} X_i $
    *
    * @param X Array of RHS vectors
    * @param Y Array of Solution vectors
    */
   void ArrayMult(const Array<const Vector *> &X, Array<Vector *> &Y) const;

   /**
    * @brief Transpose Solve $ y = Op^{-T} x $
    *
    * @param x RHS vector
    * @param y Solution vector
    */
   void MultTranspose(const Vector &x, Vector &y) const;

   /**
    * @brief Transpose Solve $ Y_i = Op^{-T} X_i $
    *
    * @param X Array of RHS vectors
    * @param Y Array of Solution vectors
    */
   void ArrayMultTranspose(const Array<const Vector *> &X,
                           Array<Vector *> &Y) const;

   /**
    * @brief Set the error print level for MUMPS
    *
    * Supported values are:
    * - 0:  No output printed
    * - 1:  Only errors printed
    * - 2:  Errors, warnings, and main stats printed
    * - 3:  Errors, warning, main stats, and terse diagnostics printed
    * - 4:  Errors, warning, main stats, diagnostics, and input/output printed
    *
    * @param print_lvl Print level, default is 2
    *
    * @note This method has to be called before SetOperator
    */
   void SetPrintLevel(int print_lvl);

   /**
    * @brief Set the matrix type
    *
    * Supported matrix types: MUMPSSolver::UNSYMMETRIC,
    * MUMPSSolver::SYMMETRIC_POSITIVE_DEFINITE,
    * and MUMPSSolver::SYMMETRIC_INDEFINITE
    *
    * @param mtype Matrix type
    *
    * @note This method has to be called before SetOperator
    */
   void SetMatrixSymType(MatType mtype);

   /**
    * @brief Set the reordering strategy
    *
    * Supported reorderings are: MUMPSSolver::AUTOMATIC,
    * MUMPSSolver::AMD, MUMPSSolver::AMF, MUMPSSolver::PORD,
    * MUMPSSolver::METIS, MUMPSSolver::PARMETIS,
    * MUMPSSolver::SCOTCH, and MUMPSSolver::PTSCOTCH
    *
    * @param method Reordering method
    *
    * @note This method has to be called before SetOperator
    */
   void SetReorderingStrategy(ReorderingStrategy method);

   /**
    * @brief Set the flag controlling reuse of the symbolic factorization
    * for multiple operators
    *
    * @param reuse Flag to reuse symbolic factorization
    *
    * @note This method has to be called before repeated calls to SetOperator
    */
   void SetReorderingReuse(bool reuse);

   /**
    * @brief Set the BLR (Block Low-Rank) factorization mode (ICNTL(35)).
    *
    * 0 = standard (no BLR, default), 1 = auto, 2 = BLR in factorization and
    * solution, 3 = BLR in factorization only.
    *
    * @note This method has to be called before SetOperator
    */
#if MFEM_MUMPS_VERSION >= 510
   void SetBLRMode(int mode);

   /**
    * @brief Set the tolerance for block low-rank (BLR) approximate
    * factorization (CNTL(7)).
    *
    * @param tol Tolerance; only used when SetBLRMode > 0.
    *
    * @note This method has to be called before SetOperator
    */
   void SetBLRTol(double tol);

   /**
    * @brief Set the BLR factorization variant (ICNTL(36)).
    *
    * 0 = UFSC (default), 1 = UCFS / low-rank update accumulation.
    *
    * @note This method has to be called before SetOperator
    */
   void SetBLRCompressionType(int type);
#endif

#if MFEM_MUMPS_VERSION >= 550
   /**
    * @brief Enable compression of BLR contribution blocks (ICNTL(37)).
    *
    * 0 = not compressed (default), 1 = compressed.
    *
    * @note This method has to be called before SetOperator
    */
   void SetBLRCBCompression(int mode);
#endif

   /**
    * @brief Set MUMPS workspace relaxation percentage above its estimate
    * (ICNTL(14)).
    *
    * Default is 20.  Increase (up to ~200) if MUMPS aborts with error -8/-9.
    */
   void SetMemRelaxation(int pct);

   /**
    * @brief Set number of OpenMP threads per MPI rank (ICNTL(16)).
    *
    * 0 (default) lets the OS / OpenMP runtime decide.
    */
   void SetNumThreads(int n);

   /**
    * @brief Set numerical pivoting threshold (CNTL(1)).
    *
    * Pass a negative value to keep the MUMPS default (0.01 for unsymmetric).
    */
   void SetPivotThreshold(double tol);

   /**
    * @brief Enable out-of-core factorization (ICNTL(22)).
    *
    * 0 = in-core (default), 1 = out-of-core.
    */
   void SetOutOfCore(int mode);

   // Destructor
   ~MUMPSSolver();

private:
   // MPI communicator
   MPI_Comm comm;

   // Number of procs
   int numProcs;

   // MPI rank
   int myid;

   // Parameter controlling the matrix type
   MatType mat_type;

   // Parameter controlling the printing level
   int print_level;

   // Parameter controlling the reordering strategy
   ReorderingStrategy reorder_method;

   // Parameter controlling whether or not to reuse the symbolic factorization
   // for multiple calls to SetOperator
   bool reorder_reuse;

#if MFEM_MUMPS_VERSION >= 510
   int    blr_mode;             ///< ICNTL(35): 0=off, 1=auto, 2=fact+sol, 3=fact only
   double blr_tol;              ///< CNTL(7):  BLR approximation tolerance
   int    blr_compression_type; ///< ICNTL(36): 0=UFSC, 1=UCFS/LUA
#endif
#if MFEM_MUMPS_VERSION >= 550
   int    blr_cb_compression;   ///< ICNTL(37): 0=CB not compressed, 1=CB compressed
#endif

   int    mem_relaxation;  ///< ICNTL(14): workspace % above estimate. Default 20.
   int    num_threads;     ///< ICNTL(16): OpenMP threads per rank. Default 0.
   double pivot_threshold; ///< CNTL(1):  pivot threshold (<0 = MUMPS default).
   int    out_of_core;     ///< ICNTL(22): 0=in-core, 1=out-of-core.

   // Local row offsets
   int row_start;

   // MUMPS object
#ifdef MFEM_USE_SINGLE
   SMUMPS_STRUC_C *id;
#else
   DMUMPS_STRUC_C *id;
#endif

   /// Method for initialization
   void Init(MPI_Comm comm_);

   /// Method for setting MUMPS internal parameters
   void SetParameters();

   /// Method for configuring storage for distributed/centralized RHS and
   /// solution
   void InitRhsSol(int nrhs) const;

#if MFEM_MUMPS_VERSION >= 530
   // Row offests array on all procs
   Array<int> row_starts;

   // Row maps and storage for distributed RHS and solution
   int *irhs_loc, *isol_loc;
   mutable real_t *rhs_loc, *sol_loc;

   // These two methods are needed to distribute the local solution
   // vectors returned by MUMPS to the original MFEM parallel partition
   int GetRowRank(int i, const Array<int> &row_starts_) const;
   void RedistributeSol(const int *rmap, const real_t *x, const int lx_loc,
                        Array<Vector *> &Y) const;
#else
   // Arrays needed for MPI_Gatherv and MPI_Scatterv
   int *recv_counts, *displs;
   mutable real_t *rhs_glob;
#endif
}; // mfem::MUMPSSolver class


#ifdef MFEM_USE_COMPLEX_MUMPS
#if MFEM_MUMPS_VERSION < 550
#error ComplexMUMPSSolver requires MUMPS >= 5.5.0
#endif
// ---------------------------------------------------------------------------
// Portable helpers for native complex MUMPS (ZMUMPS/CMUMPS).
// MUMPS defines ZMUMPS_COMPLEX as either a struct {double r, i;} or as the
// C99 type "double _Complex" depending on how the library was compiled.
// Both representations store two consecutive doubles, so reinterpret_cast
// between them and std::complex<double> is safe in all cases.
// ---------------------------------------------------------------------------

#ifdef MFEM_USE_SINGLE
/// Alias for the MUMPS single-precision complex scalar type.
using ComplexMumpsScalar = CMUMPS_COMPLEX;
#define ComplexMumpsStruc CMUMPS_STRUC_C
#else
/// Alias for the MUMPS double-precision complex scalar type.
using ComplexMumpsScalar = ZMUMPS_COMPLEX;
#define ComplexMumpsStruc ZMUMPS_STRUC_C
#endif

/// @brief MPI-distributed complex sparse direct solver via native ZMUMPS/CMUMPS.
///
/// Solves the complex system \f$ (A_r + i A_i)\,(x_r + i x_i) = b_r + i b_i \f$
/// where the operator is given as a ComplexHypreParMatrix.  The two
/// HypreParMatrix blocks are combined into a single distributed complex COO
/// array and passed directly to the ZMUMPS/CMUMPS driver; no real block-2x2
/// expansion is performed.
///
/// Supported matrix symmetry types:
///   - UNSYMMETRIC (default): whole matrix is assembled on every rank.
///   - COMPLEX_SYMMETRIC: only the lower-triangular part is used (A = A^T,
///     NOT Hermitian).
///   - COMPLEX_HERMITIAN: only the lower-triangular part is used (A = A^H).
///
/// Both symmetry modes expect the lower triangle to be provided by A_r/A_i.
class ComplexMUMPSSolver : public Solver
{
public:
   /// Matrix symmetry type fed to ZMUMPS (stored in id->sym).
   enum MatType
   {
      UNSYMMETRIC        = 0, ///< General complex unsymmetric
      COMPLEX_SYMMETRIC  = 2, ///< Complex symmetric  (A = A^T)
      COMPLEX_HERMITIAN  = 4  ///< Complex Hermitian  (A = A^H)
   };

   /// Fill-reducing reordering strategy (same options as MUMPSSolver).
   enum ReorderingStrategy
   {
      AUTOMATIC = 0, ///< MUMPS chooses automatically
      AMD,           ///< Approximate Minimum Degree
      AMF,           ///< Approximate Minimum Fill
      PORD,          ///< PORD library
      METIS,         ///< METIS library
      PARMETIS,      ///< ParMETIS library
      SCOTCH,        ///< SCOTCH library
      PTSCOTCH       ///< PTScotch library
   };

   /// Construct with an MPI communicator.
   explicit ComplexMUMPSSolver(MPI_Comm comm_);

   ~ComplexMUMPSSolver();

   /// @brief Set the operator and perform symbolic analysis + numeric
   ///        factorization.  @a op must be a ComplexHypreParMatrix.
   void SetOperator(const Operator &op) override;

   /// Overload that accepts a ComplexHypreParMatrix directly.
   /// Both real() and imag() blocks must be valid HypreParMatrix objects.
   void SetOperator(const ComplexHypreParMatrix &op);

   /// @brief Solve  A (x_r + i x_i) = (b_r + i b_i)  with explicit
   ///        real and imaginary RHS / solution vectors.
   void Mult(const Vector &b_r, const Vector &b_i,
             Vector &x_r, Vector &x_i) const;

   /// @brief Solve using 2N block vectors in the layout  [real ; imag].
   ///
   /// b and x must each have size 2*N where N is the global problem size.
   void Mult(const Vector &b, Vector &x) const override;

   /// @brief Solve multiple RHS in a single MUMPS solve call.
   ///
   /// Each element of @a X and @a Y is a 2*N block vector with layout
   /// [real ; imag].  All solves share the same factorization.
   void ArrayMult(const Array<const Vector *> &X,
                  Array<Vector *> &Y) const override;

   /// Set MUMPS diagnostic verbosity (0 = silent, 2 = normal, 4 = full).
   /// Must be called before SetOperator.
   void SetPrintLevel(int print_lvl);

   /// Set matrix symmetry type.  Must be called before SetOperator.
   void SetMatrixSymType(MatType mtype);

   /// Set fill-reducing reordering strategy.  Must be called before
   /// SetOperator.
   void SetReorderingStrategy(ReorderingStrategy method);

   /// Control the BLR (Block Low-Rank) factorization mode (ICNTL(35)).
   /// 0 = standard multifrontal, no BLR (default).
   /// 1 = BLR enabled; software selects factorization/solution variant.
   /// 2 = BLR during both factorization and solution phases (best memory gain).
   /// 3 = BLR during factorization only; solution is full-rank.
   /// Must be set to 1, 2, or 3 *before* the first SetOperator call so that
   /// the BLR preprocessing is performed during the analysis phase.
   void SetBLRMode(int mode);

   /// Set BLR approximation tolerance (CNTL(7)).  Only used when BLR is
   /// active (SetBLRMode > 0).  Smaller values give higher accuracy but
   /// less compression.  Default: 0.0 (MUMPS picks its own default).
   void SetBLRTol(double tol);

   /// Set BLR factorization variant (ICNTL(36)); only effective when BLR is
   /// active.  0 = UFSC (default), 1 = UCFS / low-rank update accumulation.
   void SetBLRCompressionType(int type);

   /// Enable compression of BLR contribution blocks (ICNTL(37)).  Only
   /// effective when BLR is active.  0 = not compressed (default),
   /// 1 = compressed (reduces memory; recommended together with ICNTL(36)=1).
   void SetBLRCBCompression(int mode);

   /// Set MUMPS workspace percentage above its internal estimate (ICNTL(14)).
   /// Default is 20.  Increase (up to ~200) if MUMPS aborts with error -8/-9.
   void SetMemRelaxation(int pct);

   /// Set number of OpenMP threads per MPI rank (ICNTL(16)).
   /// 0 (default) lets the OS / OpenMP runtime decide.
   void SetNumThreads(int n);

   /// Set numerical pivoting threshold (CNTL(1)).
   /// Pass a negative value to keep the MUMPS default (0.01 for unsymmetric).
   void SetPivotThreshold(double tol);

   /// Enable out-of-core factorization (ICNTL(22)).
   /// 0 = in-core (default), 1 = out-of-core.
   void SetOutOfCore(int mode);

   /// When true, the symbolic analysis (phase 1) from the first SetOperator
   /// call is reused in subsequent calls (only the numeric factorization is
   /// repeated).  Useful for AMR loops where the sparsity pattern changes
   /// little between refinements.  Default is false.
   void SetReorderingReuse(bool reuse);

private:
   // ---- MPI state -------------------------------------------------------
   MPI_Comm comm;
   int numProcs, myid;

   // ---- Solver configuration --------------------------------------------
   int print_level;
   MatType mat_type;
   ReorderingStrategy reorder_method;
   int  mem_relaxation;       ///< ICNTL(14): workspace % above estimate.
   int  num_threads;          ///< ICNTL(16): OpenMP threads per rank.
   double pivot_threshold;    ///< CNTL(1):  pivot threshold (<0 = MUMPS default).
   int  out_of_core;          ///< ICNTL(22): 0=in-core, 1=out-of-core.
   bool reorder_reuse;        ///< Skip analysis phase after first factorization.
   bool analysis_done_;       ///< True after the first successful phase 1.
   int    blr_mode;             ///< ICNTL(35): 0=off, 1=auto, 2=fact+sol, 3=fact only
   double blr_tol;              ///< CNTL(7):  BLR approximation tolerance
   int    blr_compression_type; ///< ICNTL(36): 0=UFSC, 1=UCFS/LUA
   int    blr_cb_compression;   ///< ICNTL(37): 0=CB not compressed, 1=CB compressed

   // ---- Matrix partition info -------------------------------------------
   int row_start; ///< First global row on this rank (0-indexed).
   int m_loc;     ///< Number of local rows.
   int n_global;  ///< Global matrix order.

   // ---- MUMPS internal object -------------------------------------------
   ComplexMumpsStruc *id;

   // ---- Distributed complex COO storage ---------------------------------
   // Separate real (I_r,J_r,D_r) and imaginary (I_i,J_i,D_i) COO arrays
   // are concatenated: first nnz_r entries come from A_r, the next nnz_i
   // from A_i.  MUMPS sums duplicate (row,col) entries, yielding the
   // desired complex value  D_r[k] + i*D_i[k].
   int *I_coo, *J_coo;
   ComplexMumpsScalar *data_coo;
   int nnz_loc_total; ///< Total COO entries on this rank (nnz_r + nnz_i).

   // ---- Helpers ---------------------------------------------------------
   void Init(MPI_Comm comm_);
   void SetParameters();
   void InitRhsSol(int nrhs) const; ///< Reallocate RHS/sol bufs when nrhs changes.

   // Distributed RHS / solution
   Array<int> row_starts;         ///< First row index on each rank.
   int        lrhs_loc, lsol_loc; ///< Sizes of local RHS / solution arrays.
   int *irhs_loc, *isol_loc;      ///< Global index maps.
   mutable ComplexMumpsScalar *rhs_loc_buf, *sol_loc_buf;
   mutable int nrhs_cur_;         ///< nrhs for which rhs/sol bufs are allocated.

   int  GetRowRank(int i, const Array<int> &row_starts_) const;
   void RedistributeSol(const int *rmap, const ComplexMumpsScalar *x,
                        int lx_loc, Vector &yr, Vector &yi) const;
}; // mfem::ComplexMUMPSSolver

#endif // MFEM_USE_COMPLEX_MUMPS

} // namespace mfem

#endif // MFEM_USE_MPI
#endif // MFEM_USE_MUMPS
#endif // MFEM_MUMPS
