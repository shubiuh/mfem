# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.

if (MFEM_USE_MPI)
  add_mfem_miniapp(udar3p
    MAIN udar3p.cpp
    EXTRA_SOURCES
      udar3p_adapt.cpp
      udar3p_pml.cpp
      udar3p_points.cpp
      udar3p_solver.cpp
    EXTRA_HEADERS
      udar3p_adapt.hpp
      udar3p_pml.hpp
      udar3p_points.hpp
      udar3p_solver.hpp
    LIBRARIES mfem)

  file(COPY README.md PALACE_UPGRADE_PLAN.md plot_udar_tensor.py
       run_udar3p_oneapi.cmd
       DESTINATION ${PROJECT_BINARY_DIR}/miniapps/udar)
endif()
