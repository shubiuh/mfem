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
  if(WIN32)
    set(MKL_LIBRARY_DIR "lib/intel64")
  elseif(UNIX)
    set(MKL_LIBRARY_DIR "lib")
  endif()
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

# Pre-cache the gomp library path so mfem_find_component's NO_DEFAULT_PATH
# find_library probe can locate it. Use gcc explicitly (CMAKE_C_COMPILER may
# be mpicc which does not support -print-file-name).
if(NOT MKL_PARDISO_MKL_OMP_RUNTIME_LIBRARY)
  find_program(_gcc_exe NAMES gcc)
  if(_gcc_exe)
    execute_process(
      COMMAND ${_gcc_exe} -print-file-name=libgomp.so
      OUTPUT_VARIABLE _mkl_gomp_path
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET)
    if(_mkl_gomp_path AND NOT _mkl_gomp_path STREQUAL "libgomp.so")
      set(MKL_PARDISO_MKL_OMP_RUNTIME_LIBRARY "${_mkl_gomp_path}" CACHE FILEPATH
        "gomp library for MKL Pardiso OMP runtime" FORCE)
      message(STATUS "MKL Pardiso: pre-cached gomp at ${_mkl_gomp_path}")
    endif()
  endif()
endif()

# Pre-cache the include dir to MKL_PARDISO_DIR/include so the find_path probe
# inside mfem_find_component uses the local oneMKL headers and not a stale
# system path that resolves to include/mkl instead.
if(MKL_PARDISO_DIR AND NOT MKL_PARDISO_INCLUDE_DIR)
  set(MKL_PARDISO_INCLUDE_DIR "${MKL_PARDISO_DIR}/include" CACHE PATH
    "MKL Pardiso include directory" FORCE)
  message(STATUS "MKL Pardiso: pre-cached include dir at ${MKL_PARDISO_DIR}/include")
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
  # LibSuffixes must be a relative PATH_SUFFIX appended to HINTS (MKL_PARDISO_DIR).
  # Use "lib" so CMake searches MKL_PARDISO_DIR/lib/ for MKL shared libraries.
  # gomp/pthread/m/dl are system libs: use "" suffix so the fallback system
  # search finds them (or the pre-cached path is used directly).
  mfem_find_package(MKL_PARDISO MKL_PARDISO
    MKL_PARDISO_DIR "include" mkl_pardiso.h "lib" mkl_core
    "Paths to headers required by MKL Pardiso." "Libraries required by MKL PARDISO."
    ADD_COMPONENT MKL_LP64 "include" "" "lib" mkl_gf_lp64
    ADD_COMPONENT MKL_OMP "include" "" "lib" mkl_gnu_thread
    # OpenMP runtime: gomp lives in the GCC lib dir, pre-cached above.
    ADD_COMPONENT MKL_OMP_RUNTIME "include" "" "" gomp
    ADD_COMPONENT MKL_PTHREAD "include" "" "" ${PLATFORM_PTHREAD}
    ADD_COMPONENT MKL_M "include" "" "" ${PLATFORM_M}
    ADD_COMPONENT MKL_DL "include" "" "" ${PLATFORM_DL}
  )
  set(MKL_PARDISO_INCLUDE_DIR ${MKL_PARDISO_DIR}/include)
  set(MKL_PARDISO_INCLUDE_DIRS ${MKL_PARDISO_INCLUDE_DIR})
endif()
