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

# Variables corresponding to defines in config.hpp (YES, NO, or value)
MFEM_VERSION           = 40800
MFEM_VERSION_STRING    = 4.8
MFEM_SOURCE_DIR        = /tmp/shubin/spack-stage/spack-stage-mfem-4.8.1-3smr4cewczdcrdyr4yzh4ire5wyyl5ey/spack-src
MFEM_INSTALL_DIR       = /home/shubin/spack/opt/spack/linux-skylake/mfem-4.8.1-3smr4cewczdcrdyr4yzh4ire5wyyl5ey
MFEM_GIT_STRING        =
MFEM_USE_MPI           = YES
MFEM_USE_METIS         = YES
MFEM_USE_METIS_5       = YES
MFEM_USE_DOUBLE        = YES
MFEM_USE_SINGLE        = NO
MFEM_DEBUG             = YES
MFEM_USE_EXCEPTIONS    = YES
MFEM_USE_ZLIB          = YES
MFEM_USE_LIBUNWIND     = NO
MFEM_USE_LAPACK        = NO
MFEM_THREAD_SAFE       = YES
MFEM_USE_LEGACY_OPENMP = YES
MFEM_USE_OPENMP        = YES
MFEM_USE_MEMALLOC      = YES
MFEM_TIMER_TYPE        = 2
MFEM_USE_SUNDIALS      = YES
MFEM_USE_SUITESPARSE   = YES
MFEM_USE_SUPERLU       = YES
MFEM_USE_SUPERLU5      = NO
MFEM_USE_MUMPS         = YES
MFEM_USE_STRUMPACK     = NO
MFEM_USE_GINKGO        = NO
MFEM_USE_AMGX          = NO
MFEM_USE_MAGMA         = NO
MFEM_USE_GNUTLS        = NO
MFEM_USE_NETCDF        = YES
MFEM_USE_PETSC         = NO
MFEM_USE_SLEPC         = NO
MFEM_USE_MPFR          = NO
MFEM_USE_SIDRE         = NO
MFEM_USE_FMS           = YES
MFEM_USE_CONDUIT       = NO
MFEM_USE_PUMI          = NO
MFEM_USE_HIOP          = YES
MFEM_USE_GSLIB         = YES
MFEM_USE_CUDA          = NO
MFEM_USE_HIP           = NO
MFEM_USE_RAJA          = NO
MFEM_USE_OCCA          = NO
MFEM_USE_CEED          = YES
MFEM_USE_CALIPER       = NO
MFEM_USE_UMPIRE        = NO
MFEM_USE_SIMD          = NO
MFEM_USE_ADIOS2        = NO
MFEM_USE_MKL_CPARDISO  = YES
MFEM_USE_MKL_PARDISO   = YES
MFEM_USE_MOONOLITH     = NO
MFEM_USE_ADFORWARD     = NO
MFEM_USE_CODIPACK      = NO
MFEM_USE_BENCHMARK     = NO
MFEM_USE_PARELAG       = NO
MFEM_USE_TRIBOL        = NO
MFEM_USE_ENZYME        = NO

# Compiler, compile options, and link options
MFEM_CXX       = /home/shubin/spack/opt/spack/linux-skylake/intel-oneapi-mpi-2021.16.0-5hheospttsdenb5pjlaicjvblrmphqpb/mpi/2021.16/bin/mpiicpx
MFEM_HOST_CXX  = /home/shubin/spack/opt/spack/linux-skylake/intel-oneapi-mpi-2021.16.0-5hheospttsdenb5pjlaicjvblrmphqpb/mpi/2021.16/bin/mpiicpx
MFEM_CPPFLAGS  =
MFEM_CXXFLAGS  = -g -O0 -std=c++17
MFEM_TPLFLAGS  = -I/home/shubin/spack/opt/spack/linux-skylake/hypre-2.33.0-nydtorgebuplgvdb5nblnihoejmcqdvj/include -I/home/shubin/spack/opt/spack/linux-skylake/superlu-dist-9.1.0-gv62ouswnzooinkogesuef3xgmj76ywg/include -I/home/shubin/spack/opt/spack/linux-skylake/parmetis-4.0.3-7cnpoavovuyu2evmrkpkr3eriadkptxy/include -I/home/shubin/spack/opt/spack/linux-skylake/mumps-5.8.1-d7ovatxpbv42lv4vi4ncov6eaulkl4hw/include -fiopenmp -I/home/shubin/spack/opt/spack/linux-skylake/metis-5.1.0-ptgixmjg5q4yc5sddocxqh2gbzyb5bml/include -I/home/shubin/spack/opt/spack/linux-skylake/libfms-0.2.0-n2a6ozcomeunptvip4qm6ndfn5a5kkin/include -I/home/shubin/spack/opt/spack/linux-skylake/sundials-7.4.0-rfhy5gxpt5kmpshjkpk4h2h7iwlush4k/include -I/home/shubin/spack/opt/spack/linux-skylake/suite-sparse-7.8.3-bm7ilyqey3pbnmar433gx3z5q5cq326e/include -I/home/shubin/spack/opt/spack/linux-skylake/netcdf-c-4.9.3-ibthin4b6ckllnkwptsjnratnjgpfuud/include -I/home/shubin/spack/opt/spack/linux-skylake/hiop-1.1.1-2ludq2rxksdnsskz6skzz7m7zhdfmodz/include -I/home/shubin/spack/opt/spack/linux-skylake/openblas-0.3.29-jlt2hcnipwgblbn5znf7ppag5m2ztxsw/include -I/home/shubin/spack/opt/spack/linux-skylake/gslib-1.0.9-gp5tfc4hizh5zikv7gdekvqgm7c2wt5y/include -I/home/shubin/spack/opt/spack/linux-skylake/libceed-0.12.0-xk2hogkkykvpi7fepuewkfgk42iucnuv/include -I/home/shubin/spack/opt/spack/linux-skylake/intel-oneapi-mkl-2024.2.2-25577eemsh26sbnl2unimld7joi4sf6z/mkl/latest/include -I/home/shubin/spack/opt/spack/linux-skylake/intel-oneapi-mkl-2024.2.2-25577eemsh26sbnl2unimld7joi4sf6z/mkl/latest/include -fiopenmp -fiopenmp -I/home/shubin/spack/opt/spack/linux-skylake/zlib-ng-2.2.4-cqnc2q6cuzvwthu7dbl6pemaaqwypv3l/include
MFEM_INCFLAGS  = -I$(MFEM_INC_DIR) $(MFEM_TPLFLAGS)
MFEM_PICFLAG   =
MFEM_FLAGS     = $(MFEM_CPPFLAGS) $(MFEM_CXXFLAGS) $(MFEM_INCFLAGS)
MFEM_EXT_LIBS  = -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/hypre-2.33.0-nydtorgebuplgvdb5nblnihoejmcqdvj/lib -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/openblas-0.3.29-jlt2hcnipwgblbn5znf7ppag5m2ztxsw/lib -L/home/shubin/spack/opt/spack/linux-skylake/hypre-2.33.0-nydtorgebuplgvdb5nblnihoejmcqdvj/lib -L/home/shubin/spack/opt/spack/linux-skylake/openblas-0.3.29-jlt2hcnipwgblbn5znf7ppag5m2ztxsw/lib -lHYPRE -lopenblas -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/superlu-dist-9.1.0-gv62ouswnzooinkogesuef3xgmj76ywg/lib -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/parmetis-4.0.3-7cnpoavovuyu2evmrkpkr3eriadkptxy/lib -L/home/shubin/spack/opt/spack/linux-skylake/superlu-dist-9.1.0-gv62ouswnzooinkogesuef3xgmj76ywg/lib -L/home/shubin/spack/opt/spack/linux-skylake/parmetis-4.0.3-7cnpoavovuyu2evmrkpkr3eriadkptxy/lib -lsuperlu_dist -lparmetis -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/openblas-0.3.29-jlt2hcnipwgblbn5znf7ppag5m2ztxsw/lib -L/home/shubin/spack/opt/spack/linux-skylake/openblas-0.3.29-jlt2hcnipwgblbn5znf7ppag5m2ztxsw/lib -lopenblas -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/mumps-5.8.1-d7ovatxpbv42lv4vi4ncov6eaulkl4hw/lib -L/home/shubin/spack/opt/spack/linux-skylake/mumps-5.8.1-d7ovatxpbv42lv4vi4ncov6eaulkl4hw/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/metis-5.1.0-ptgixmjg5q4yc5sddocxqh2gbzyb5bml/lib -L/home/shubin/spack/opt/spack/linux-skylake/metis-5.1.0-ptgixmjg5q4yc5sddocxqh2gbzyb5bml/lib -lmetis -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/libfms-0.2.0-n2a6ozcomeunptvip4qm6ndfn5a5kkin/lib -L/home/shubin/spack/opt/spack/linux-skylake/libfms-0.2.0-n2a6ozcomeunptvip4qm6ndfn5a5kkin/lib -lfms -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/sundials-7.4.0-rfhy5gxpt5kmpshjkpk4h2h7iwlush4k/lib -L/home/shubin/spack/opt/spack/linux-skylake/sundials-7.4.0-rfhy5gxpt5kmpshjkpk4h2h7iwlush4k/lib -lsundials_arkode -lsundials_cvodes -lsundials_nvecserial -lsundials_kinsol -lsundials_nvecparallel -lsundials_nvecmpiplusx -lsundials_core -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/suite-sparse-7.8.3-bm7ilyqey3pbnmar433gx3z5q5cq326e/lib -L/home/shubin/spack/opt/spack/linux-skylake/suite-sparse-7.8.3-bm7ilyqey3pbnmar433gx3z5q5cq326e/lib -lklu -lbtf -lumfpack -lcholmod -lcolamd -lamd -lcamd -lccolamd -lsuitesparseconfig -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/netcdf-c-4.9.3-ibthin4b6ckllnkwptsjnratnjgpfuud/lib -L/home/shubin/spack/opt/spack/linux-skylake/netcdf-c-4.9.3-ibthin4b6ckllnkwptsjnratnjgpfuud/lib -lnetcdf -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/hiop-1.1.1-2ludq2rxksdnsskz6skzz7m7zhdfmodz/lib -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/openblas-0.3.29-jlt2hcnipwgblbn5znf7ppag5m2ztxsw/lib -L/home/shubin/spack/opt/spack/linux-skylake/hiop-1.1.1-2ludq2rxksdnsskz6skzz7m7zhdfmodz/lib -L/home/shubin/spack/opt/spack/linux-skylake/openblas-0.3.29-jlt2hcnipwgblbn5znf7ppag5m2ztxsw/lib -lhiop -lopenblas -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/gslib-1.0.9-gp5tfc4hizh5zikv7gdekvqgm7c2wt5y/lib -L/home/shubin/spack/opt/spack/linux-skylake/gslib-1.0.9-gp5tfc4hizh5zikv7gdekvqgm7c2wt5y/lib -lgs -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/libceed-0.12.0-xk2hogkkykvpi7fepuewkfgk42iucnuv/lib -L/home/shubin/spack/opt/spack/linux-skylake/libceed-0.12.0-xk2hogkkykvpi7fepuewkfgk42iucnuv/lib -lceed -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/intel-oneapi-mkl-2024.2.2-25577eemsh26sbnl2unimld7joi4sf6z/mkl/latest/lib/intel64 -L/home/shubin/spack/opt/spack/linux-skylake/intel-oneapi-mkl-2024.2.2-25577eemsh26sbnl2unimld7joi4sf6z/mkl/latest/lib/intel64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/intel-oneapi-mkl-2024.2.2-25577eemsh26sbnl2unimld7joi4sf6z/mkl/latest/lib/intel64 -L/home/shubin/spack/opt/spack/linux-skylake/intel-oneapi-mkl-2024.2.2-25577eemsh26sbnl2unimld7joi4sf6z/mkl/latest/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L/home/shubin/spack/opt/spack/linux-skylake/intel-oneapi-mkl-2024.2.2-25577eemsh26sbnl2unimld7joi4sf6z/compiler/latest/lib -liomp5 -lpthread -lm -ldl -lrt -Wl,-rpath,/home/shubin/spack/opt/spack/linux-skylake/zlib-ng-2.2.4-cqnc2q6cuzvwthu7dbl6pemaaqwypv3l/lib -L/home/shubin/spack/opt/spack/linux-skylake/zlib-ng-2.2.4-cqnc2q6cuzvwthu7dbl6pemaaqwypv3l/lib -lz
MFEM_LIBS      =  -L$(MFEM_LIB_DIR) -lmfem $(MFEM_EXT_LIBS)
MFEM_LIB_FILE  = $(MFEM_LIB_DIR)/libmfem.a
MFEM_STATIC    = YES
MFEM_SHARED    = NO
MFEM_BUILD_TAG = Linux Dell7690 x86_64
MFEM_PREFIX    = /home/shubin/spack/opt/spack/linux-skylake/mfem-4.8.1-3smr4cewczdcrdyr4yzh4ire5wyyl5ey
MFEM_INC_DIR   = /home/shubin/spack/opt/spack/linux-skylake/mfem-4.8.1-3smr4cewczdcrdyr4yzh4ire5wyyl5ey/include
MFEM_LIB_DIR   = /home/shubin/spack/opt/spack/linux-skylake/mfem-4.8.1-3smr4cewczdcrdyr4yzh4ire5wyyl5ey/lib

# Location of test.mk
MFEM_TEST_MK = /home/shubin/spack/opt/spack/linux-skylake/mfem-4.8.1-3smr4cewczdcrdyr4yzh4ire5wyyl5ey/share/mfem/test.mk

# Command used to launch MPI jobs
MFEM_MPIEXEC    = mpirun
MFEM_MPIEXEC_NP = -np
MFEM_MPI_NP     = 4

# The NVCC compiler cannot link with -x=cu
MFEM_LINK_FLAGS := $(filter-out -x=cu -xhip, $(MFEM_FLAGS))

# Optional extra configuration
