# Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

# Defines the following variables:
#   - MKL_CPARDISO_FOUND
#   - MKL_CPARDISO_LIBRARIES
#   - MKL_CPARDISO_INCLUDE_DIRS

if(NOT MKL_MPI_WRAPPER_LIB)
  message(FATAL_ERROR "MKL CPardiso enabled but no MKL MPI Wrapper lib specified")
endif()

if(NOT MKL_LIBRARY_DIR)
  message(WARNING "Using default MKL library path. Double check the variable MKL_LIBRARY_DIR")
  set(MKL_LIBRARY_DIR "lib")
endif()

# Pre-cache the include dir to MKL_CPARDISO_DIR/include so the find_path probe
# inside mfem_find_component uses the local oneMKL headers.
if(MKL_CPARDISO_DIR AND NOT MKL_CPARDISO_INCLUDE_DIR)
  set(MKL_CPARDISO_INCLUDE_DIR "${MKL_CPARDISO_DIR}/include" CACHE PATH
    "MKL CPardiso include directory" FORCE)
  message(STATUS "MKL CPardiso: pre-cached include dir at ${MKL_CPARDISO_DIR}/include")
endif()

include(MfemCmakeUtilities)
# LibSuffixes must be a relative PATH_SUFFIX appended to HINTS (MKL_CPARDISO_DIR).
# Use "lib" so CMake searches MKL_CPARDISO_DIR/lib/ for all MKL shared libraries.
mfem_find_package(MKL_CPARDISO MKL_CPARDISO
    MKL_CPARDISO_DIR "include" mkl_cluster_sparse_solver.h "lib" mkl_core
  "Paths to headers required by MKL CPardiso." "Libraries required by MKL CPARDISO."
  ADD_COMPONENT MKL_LP64 "include" "" "lib" mkl_gf_lp64
  ADD_COMPONENT MKL_SEQUENTIAL "include" "" "lib" mkl_sequential
  ADD_COMPONENT MKL_MPI_WRAPPER "include" "" "lib" ${MKL_MPI_WRAPPER_LIB})
