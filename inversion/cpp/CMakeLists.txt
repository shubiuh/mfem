cmake_minimum_required(VERSION 3.23)
project(udar_inversion LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(UDAR_BUILD_ROL "Build the Trilinos/ROL inversion backend" ON)
option(UDAR_BUILD_TAO "Build the PETSc/TAO inversion backend" OFF)

# Build only the ROL components used by this adapter. No changes to MFEM or
# dependency downloads occur during configuration. Match standalone ROL's
# compatibility layer, with static support objects inside our exported DLL.
if(UDAR_BUILD_ROL)
  set(ROL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/src/trilinos-rol/packages/rol"
      CACHE PATH "Trilinos ROL source directory (tested: release 17.2.1)")
  set(PUGIXML_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/src/pugixml-rol"
      CACHE PATH "pugixml 1.14 source directory")
  if(NOT EXISTS "${ROL_SOURCE_DIR}/src/algorithm/TypeB/ROL_TypeB_QuasiNewtonAlgorithm.hpp")
    message(FATAL_ERROR "ROL source missing. Follow inversion/README.md dependency setup.")
  endif()
  if(NOT EXISTS "${PUGIXML_SOURCE_DIR}/src/pugixml.cpp")
    message(FATAL_ERROR "pugixml source missing. Follow inversion/README.md dependency setup.")
  endif()

  set(UDAR_BLAS_LIBRARIES "" CACHE STRING "BLAS/LAPACK libraries; optional on Unix, required for oneMKL on Windows")
  if(NOT UDAR_BLAS_LIBRARIES)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    set(UDAR_BLAS_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  endif()
  # oneMKL on Windows exports uppercase Fortran symbols; Unix BLAS uses name_.
  if(WIN32)
    set(F77_BLAS_MANGLE "(name,NAME) NAME")
  else()
    set(F77_BLAS_MANGLE "(name,NAME) name ## _")
  endif()
  configure_file("${ROL_SOURCE_DIR}/cmake/ROL_config.h.in" ROL_config.h @ONLY)
  set(rol_src "${ROL_SOURCE_DIR}/src")
  set(rol_includes "${CMAKE_CURRENT_BINARY_DIR}")
  foreach(part algorithm elementwise function mol oed sol status step utils vector zoo)
    file(GLOB_RECURSE headers "${rol_src}/${part}/*.hpp")
    foreach(header IN LISTS headers)
      get_filename_component(dir "${header}" DIRECTORY)
      list(APPEND rol_includes "${dir}")
    endforeach()
  endforeach()
  list(REMOVE_DUPLICATES rol_includes)
  foreach(part teuchos/blas teuchos/la simple/lapack simple/mpi
               simple/parameterlist std/shared_ptr teuchos/stacktrace teuchos-lite)
    list(APPEND rol_includes "${rol_src}/compatibility/${part}")
  endforeach()

  add_library(udar_rol SHARED udar_rol.cpp
    "${rol_src}/compatibility/teuchos-lite/Teuchos_BLAS.cpp"
    "${rol_src}/compatibility/teuchos-lite/Teuchos_CompObject.cpp"
    "${rol_src}/compatibility/teuchos-lite/Teuchos_TestForException.cpp"
    "${PUGIXML_SOURCE_DIR}/src/pugixml.cpp")
  target_include_directories(udar_rol PRIVATE ${rol_includes} "${PUGIXML_SOURCE_DIR}/src")
  target_link_libraries(udar_rol PRIVATE ${UDAR_BLAS_LIBRARIES})
  if(MSVC)
    target_compile_options(udar_rol PRIVATE /bigobj /EHsc /Zc:__cplusplus /wd4819)
    target_compile_definitions(udar_rol PRIVATE NOMINMAX)
  endif()
  set_target_properties(udar_rol PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
endif()

if(UDAR_BUILD_TAO)
  set(PETSC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/install/petsc-oneapi"
      CACHE PATH "Installed PETSc prefix containing include/ and lib/")
  find_path(PETSC_INCLUDE_DIR petsctao.h PATHS "${PETSC_DIR}/include" NO_DEFAULT_PATH REQUIRED)
  set(PETSC_LIBRARIES)
  foreach(component tao ts snes ksp dm mat vec sys)
    find_library(PETSC_${component}_LIBRARY NAMES "petsc${component}"
      "libpetsc${component}" PATHS "${PETSC_DIR}/lib" NO_DEFAULT_PATH REQUIRED)
    list(APPEND PETSC_LIBRARIES "${PETSC_${component}_LIBRARY}")
  endforeach()
  file(READ "${PETSC_INCLUDE_DIR}/petscconf.h" petsc_configuration)
  if(petsc_configuration MATCHES "PETSC_USE_COMPLEX" OR
     NOT petsc_configuration MATCHES "PETSC_USE_REAL_DOUBLE" OR
     petsc_configuration MATCHES "PETSC_USE_64BIT_INDICES" OR
     petsc_configuration MATCHES "PETSC_HAVE_64BIT_BLAS_INDICES" OR
     petsc_configuration MATCHES "PETSC_USE_DEBUG")
    message(FATAL_ERROR
      "udar_tao requires Release PETSc with real doubles, 32-bit PetscInt, and LP64 BLAS")
  endif()
  set(UDAR_PETSC_EXTERNAL_LIBRARIES "" CACHE STRING
      "PETSc external link libraries; inferred for the Windows oneAPI stack")
  if(NOT UDAR_PETSC_EXTERNAL_LIBRARIES AND WIN32)
    set(UDAR_PETSC_EXTERNAL_LIBRARIES
      "C:/OPT/Intel/oneAPI/mpi/latest/lib/impi.lib;C:/OPT/Intel/oneAPI/mkl/latest/lib/mkl_intel_lp64_dll.lib;C:/OPT/Intel/oneAPI/mkl/latest/lib/mkl_sequential_dll.lib;C:/OPT/Intel/oneAPI/mkl/latest/lib/mkl_core_dll.lib;gdi32.lib;user32.lib;advapi32.lib;kernel32.lib;ws2_32.lib")
  endif()
  add_library(udar_tao SHARED udar_tao.cpp)
  target_include_directories(udar_tao PRIVATE "${PETSC_INCLUDE_DIR}")
  target_link_libraries(udar_tao PRIVATE ${PETSC_LIBRARIES}
    ${UDAR_PETSC_EXTERNAL_LIBRARIES})
  if(MSVC)
    set_property(TARGET udar_tao PROPERTY MSVC_RUNTIME_LIBRARY MultiThreadedDLL)
    target_compile_options(udar_tao PRIVATE /EHsc /Zc:__cplusplus
      /Zc:preprocessor /wd4819)
    target_compile_definitions(udar_tao PRIVATE NOMINMAX)
  endif()
  set_target_properties(udar_tao PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
endif()
