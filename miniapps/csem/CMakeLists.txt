# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.

if (MFEM_USE_MPI)
  add_mfem_miniapp(csem
    MAIN csem.cpp
    csem_solver.cpp
    EXTRA_HEADERS csem_solver.hpp
    LIBRARIES mfem)

  file(COPY README.md PAPER_METHODS.md sample_receivers.csv
            sample_transmitter.csv sample_transmitters.csv
            sample_region_model.csv sample_region_initial.csv
            sample_cell_model_2x2x2.csv
            sample_topography.mesh sample_topography_materials.csv
            sample_topography_initial.csv paper_scenario2_reference.csv
            csem_convergence.py generate_layered_reference.py
       DESTINATION ${PROJECT_BINARY_DIR}/miniapps/csem)

  if (MFEM_ENABLE_TESTING AND MFEM_USE_MUMPS)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if (Python3_Interpreter_FOUND)
      add_test(NAME csem_layered_convergence
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/csem_convergence.py
          --executable $<TARGET_FILE:csem>
          --reference ${CMAKE_CURRENT_SOURCE_DIR}/paper_scenario2_reference.csv
          --output-dir ${PROJECT_BINARY_DIR}/miniapps/csem/convergence-test
          --suite paper --orders 1,2 --h-levels 0,1 --h-order 1
          --components Ex,Hy,Hz --solver mumps
          --p-tolerance-percent 100 --h-tolerance-percent 50)
      set_tests_properties(csem_layered_convergence
        PROPERTIES LABELS "CSEM;regression" TIMEOUT 600)
    endif()
  endif()
endif()
