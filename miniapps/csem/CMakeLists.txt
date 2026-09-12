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
            sample_transmitter.csv
       DESTINATION ${PROJECT_BINARY_DIR}/miniapps/csem)
endif()
