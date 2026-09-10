cmake_minimum_required(VERSION 3.23)
project(udar_inversion LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Build only the ROL components used by this adapter. No changes to MFEM or
# dependency downloads occur during configuration. Match standalone ROL's
# compatibility layer, with static support objects inside our exported DLL.
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
