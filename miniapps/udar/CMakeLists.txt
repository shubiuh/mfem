# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.

if (MFEM_USE_MPI)
  add_mfem_miniapp(udar3p
    MAIN udar3p.cpp
    EXTRA_SOURCES
      udar3p_adapt.cpp
      udar3p_frequency.cpp
      udar3p_pml.cpp
      udar3p_postprocess.cpp
      udar3p_points.cpp
      udar3p_solver.cpp
    EXTRA_HEADERS
      udar3p_adapt.hpp
      udar3p_frequency.hpp
      udar3p_pml.hpp
      udar3p_postprocess.hpp
      udar3p_points.hpp
      udar3p_solver.hpp
    LIBRARIES mfem)

  # Explicit target: exercises retained ND transfers without a Maxwell solve.
  add_executable(udar_hierarchy_test EXCLUDE_FROM_ALL
    tests/hierarchy.cpp udar3p_adapt.cpp)
  target_link_libraries(udar_hierarchy_test PRIVATE mfem)

  add_executable(udar_assembly_test EXCLUDE_FROM_ALL
    tests/assembly.cpp udar3p_frequency.cpp udar3p_solver.cpp udar3p_pml.cpp)
  target_link_libraries(udar_assembly_test PRIVATE mfem)

  add_executable(udar_output_benchmark EXCLUDE_FROM_ALL
    tests/output_benchmark.cpp udar3p_postprocess.cpp)
  target_link_libraries(udar_output_benchmark PRIVATE mfem)

  file(COPY README.md PALACE_UPGRADE_PLAN.md OUTPUT_PERFORMANCE.md plot_udar_tensor.py
       run_udar3p_oneapi.cmd
       DESTINATION ${PROJECT_BINARY_DIR}/miniapps/udar)
endif()
