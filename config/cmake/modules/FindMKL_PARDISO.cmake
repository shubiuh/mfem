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
#   - MKL_PARDISO_FOUND
#   - MKL_PARDISO_LIBRARIES
#   - MKL_PARDISO_INCLUDE_DIRS

if(NOT MKL_LIBRARY_DIR)
  message(WARNING "Using default MKL library path. Double check the variable MKL_LIBRARY_DIR")
  set(MKL_LIBRARY_DIR "lib/intel64")
  message(STATUS "MKL_LIBRARY_DIR set to ${MKL_LIBRARY_DIR}")
else()
  message(STATUS "MKL_LIBRARY_DIR set to ${MKL_LIBRARY_DIR}")
endif()

if(NOT MKL_COMPILER_DIR)
  message(WARNING "Using default MKL compiler path. Double check the variable MKL_COMPILER_DIR")
  set(MKL_COMPILER_DIR "compiler/latest/lib")
  message(STATUS "MKL_COMPILER_DIR set to ${MKL_COMPILER_DIR}")
else()
  message(STATUS "MKL_COMPILER_DIR set to ${MKL_COMPILER_DIR}")
endif()

if(WIN32)
  set(PLATFORM_PTHREAD "")
  set(PLATFORM_M "")
  set(PLATFORM_DL "")
else()
  set(PLATFORM_PTHREAD pthread)
  set(PLATFORM_M m)
  set(PLATFORM_DL dl)
endif()

include(MfemCmakeUtilities)

if(WIN32) # Windows
  mfem_find_package(MKL_PARDISO MKL_PARDISO
    MKL_PARDISO_DIR "include" mkl_pardiso.h ${MKL_LIBRARY_DIR} mkl_core_dll
    "Paths to headers required by MKL Pardiso." "Libraries required by MKL PARDISO."

    # LP64 interface
    ADD_COMPONENT MKL_LP64 "include" "" ${MKL_LIBRARY_DIR} mkl_intel_lp64_dll

    # Threaded OpenMP version (remove sequential)
    ADD_COMPONENT MKL_OMP "include" "" ${MKL_LIBRARY_DIR} mkl_intel_thread_dll

    # MKL core
    ADD_COMPONENT MKL_CORE "include" "" ${MKL_LIBRARY_DIR} mkl_core_dll

    # OpenMP runtime
    ADD_COMPONENT MKL_OMP_RUNTIME "include" "" ${MKL_COMPILER_DIR} libiomp5md
  )
else() # Linux, macOS, etc.
  mfem_find_package(MKL_PARDISO MKL_PARDISO
    MKL_PARDISO_DIR "include" mkl_pardiso.h ${MKL_LIBRARY_DIR} mkl_core
    "Paths to headers required by MKL Pardiso." "Libraries required by MKL PARDISO."
    ADD_COMPONENT MKL_LP64 "include" "" ${MKL_LIBRARY_DIR} mkl_intel_lp64
    # ADD_COMPONENT MKL_SEQUENTIAL "include" "" ${MKL_LIBRARY_DIR} mkl_sequential
    ADD_COMPONENT MKL_OMP "include" "" ${MKL_LIBRARY_DIR} mkl_intel_thread
    # OpenMP runtime
    ADD_COMPONENT MKL_OMP_RUNTIME "include" "" ${MKL_COMPILER_DIR} iomp5
    ADD_COMPONENT MKL_PTHREAD "include" "" ${MKL_LIBRARY_DIR} ${PLATFORM_PTHREAD}
    ADD_COMPONENT MKL_M "include" "" ${MKL_LIBRARY_DIR} ${PLATFORM_M}
    ADD_COMPONENT MKL_DL "include" "" ${MKL_LIBRARY_DIR} ${PLATFORM_DL}
  )
endif()
