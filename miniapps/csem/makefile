# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.

MFEM_DIR ?= ../..
MFEM_BUILD_DIR ?= ../..
MFEM_INSTALL_DIR ?= ../../mfem
SRC = $(if $(MFEM_DIR:../..=),$(MFEM_DIR)/miniapps/csem/,)
CONFIG_MK = $(or $(wildcard $(MFEM_BUILD_DIR)/config/config.mk),\
   $(wildcard $(MFEM_INSTALL_DIR)/share/mfem/config.mk))

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

.SUFFIXES:
.SUFFIXES: .o .cpp .mk
.PHONY: all clean clean-build clean-exec
.PRECIOUS: %.o

all: csem

csem: $(SRC)csem.cpp csem_solver.o $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) -c $(<) -o csem.o
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ csem.o csem_solver.o $(MFEM_LIBS)

csem_solver.o: $(SRC)csem_solver.cpp $(SRC)csem_solver.hpp $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) -c $(<) -o $(@)

MFEM_TESTS = csem
include $(MFEM_TEST_MK)

RUN_MPI = $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $(MFEM_MPI_NP)
csem-test-par: csem
	@$(call mfem-test,$<, $(RUN_MPI), CSEM miniapp,\
	-nx 3 -ny 3 -nz 3 -o 1 -nr 3 -no-vis)

$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: clean-build clean-exec

clean-build:
	rm -f *.o *~ csem
	rm -rf *.dSYM *.TVD.*breakpoints

clean-exec:
	@rm -rf csem*_field csem*.csv
