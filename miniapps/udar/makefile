# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.

MFEM_DIR ?= ../..
MFEM_BUILD_DIR ?= ../..
MFEM_INSTALL_DIR ?= ../../mfem
SRC = $(if $(MFEM_DIR:../..=),$(MFEM_DIR)/miniapps/udar/,)
CONFIG_MK = $(or $(wildcard $(MFEM_BUILD_DIR)/config/config.mk),\
   $(wildcard $(MFEM_INSTALL_DIR)/share/mfem/config.mk))

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

UDAR_SOURCES = udar3p.cpp udar3p_adapt.cpp udar3p_frequency.cpp \
   udar3p_pml.cpp \
   udar3p_postprocess.cpp \
   udar3p_points.cpp udar3p_solver.cpp
UDAR_OBJECTS = $(UDAR_SOURCES:.cpp=.o)

.PHONY: all clean
.SUFFIXES:
.SUFFIXES: .o .cpp

all: udar3p

udar3p: $(UDAR_OBJECTS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(UDAR_OBJECTS) $(MFEM_LIBS)

%.o: $(SRC)%.cpp $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) -c $< -o $@

udar3p.o: $(SRC)udar3p_adapt.hpp $(SRC)udar3p_frequency.hpp \
   $(SRC)udar3p_pml.hpp \
   $(SRC)udar3p_postprocess.hpp \
   $(SRC)udar3p_points.hpp $(SRC)udar3p_solver.hpp
udar3p_adapt.o: $(SRC)udar3p_adapt.hpp
udar3p_frequency.o: $(SRC)udar3p_frequency.hpp $(SRC)udar3p_pml.hpp \
   $(SRC)udar3p_solver.hpp
udar3p_pml.o: $(SRC)udar3p_pml.hpp
udar3p_postprocess.o: $(SRC)udar3p_postprocess.hpp
udar3p_points.o: $(SRC)udar3p_points.hpp
udar3p_solver.o: $(SRC)udar3p_solver.hpp $(SRC)udar3p_pml.hpp

$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean:
	rm -f udar3p $(UDAR_OBJECTS) *~
