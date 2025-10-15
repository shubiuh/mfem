# How to Use Spack-Installed MFEM

This guide explains how to use MFEM (Modular Finite Element Methods) libraries that have been installed via Spack to build and run your own MFEM examples and applications.

## Overview

Spack is a package manager that builds and installs multiple versions and configurations of software packages. In your workspace, there are two different MFEM installations:

1. **CUDA-enabled build** (LLVM 18.1.3): `/data/shubin/opt/spack/store/llvm-18.1.3/mfem-4.8.1-3lrmot6fydgnecja7lr2626ecsgc52n5`
2. **CPU-only build** (GCC 15.2.0): `/data/shubin/opt/spack/store/gcc-15.2.0/mfem-4.8.1-boqsa3oms4wytmhnbzwrtfxqw2gt75bm`

## Installation Features

### CUDA-enabled MFEM (LLVM build)
- **Compiler**: NVCC with MPICXX backend
- **GPU Support**: CUDA 12.9.0 with SM_70 architecture
- **MPI**: Enabled with MPICH 4.3.1
- **Dependencies**: HYPRE, MUMPS, METIS, PETSc, SuiteSparse, Intel MKL, NetCDF, MPFR, GSLib, FMS
- **OpenMP**: Enabled
- **Features**: CUDA kernels, GPU acceleration

### CPU-only MFEM (GCC build)
- **Compiler**: GCC 15.2.0 with MPICXX
- **MPI**: Enabled with MPICH
- **Dependencies**: HYPRE, MUMPS, METIS, PETSc, SuiteSparse, NetCDF, MPFR, GSLib, FMS
- **OpenMP**: Disabled
- **Features**: CPU-only computations

## Directory Structure

Each Spack-installed MFEM provides:
```
<MFEM_INSTALL_DIR>/
├── include/                    # Header files
│   ├── mfem.hpp               # Main MFEM header
│   ├── mfem-performance.hpp   # Performance header
│   └── mfem/                  # All MFEM headers organized by module
├── lib/                       # Library files
│   └── libmfem.so.4.8        # Shared library
└── share/mfem/               # Configuration and examples
    ├── config.mk             # Makefile configuration
    ├── test.mk               # Testing utilities
    ├── examples/             # Pre-compiled examples
    ├── miniapps/             # Pre-compiled miniapps
    └── data/                 # Mesh files and test data
```

## Building New MFEM Applications

### Method 1: Using the Provided Makefile Template

The easiest way is to copy and modify the makefile from the installed examples:

1. **Copy the template makefile**:
   ```bash
   # For CUDA build
   cp /data/shubin/opt/spack/store/llvm-18.1.3/mfem-4.8.1-3lrmot6fydgnecja7lr2626ecsgc52n5/share/mfem/examples/makefile ./

   # For CPU build  
   cp /data/shubin/opt/spack/store/gcc-15.2.0/mfem-4.8.1-boqsa3oms4wytmhnbzwrtfxqw2gt75bm/share/mfem/examples/makefile ./
   ```

2. **Modify the MFEM_INSTALL_DIR variable** in the makefile:
   ```makefile
   # For CUDA build
   MFEM_INSTALL_DIR = /data/shubin/opt/spack/store/llvm-18.1.3/mfem-4.8.1-3lrmot6fydgnecja7lr2626ecsgc52n5

   # For CPU build
   MFEM_INSTALL_DIR = /data/shubin/opt/spack/store/gcc-15.2.0/mfem-4.8.1-boqsa3oms4wytmhnbzwrtfxqw2gt75bm
   ```

3. **Build your program**:
   ```bash
   # Build a specific example
   make ex25
   
   # Build all examples
   make all
   
   # Clean build artifacts
   make clean
   ```

### Method 2: Manual Compilation

For more control, you can compile manually using the configuration variables:

1. **Set the MFEM installation path**:
   ```bash
   # For CUDA build
   export MFEM_DIR="/data/shubin/opt/spack/store/llvm-18.1.3/mfem-4.8.1-3lrmot6fydgnecja7lr2626ecsgc52n5"
   
   # For CPU build
   export MFEM_DIR="/data/shubin/opt/spack/store/gcc-15.2.0/mfem-4.8.1-boqsa3oms4wytmhnbzwrtfxqw2gt75bm"
   ```

2. **Load the configuration**:
   ```bash
   source $MFEM_DIR/share/mfem/config.mk
   ```

3. **Compile your program**:
   ```bash
   # For CUDA build
   $MFEM_CXX $MFEM_FLAGS your_program.cpp -o your_program $MFEM_LIBS
   
   # Example:
   /data/shubin/opt/spack/store/none-none/cuda-12.9.0-xxvugy5yp7rx7hcg3xc5qthbgzr3ahgi/bin/nvcc \
     -O3 -std=c++17 -x=cu --expt-extended-lambda -arch=sm_70 \
     -I$MFEM_DIR/include [other flags] \
     your_program.cpp -o your_program \
     -L$MFEM_DIR/lib -lmfem [other libraries]
   ```

### Method 3: Using CMake

**Note**: Spack-installed MFEM typically doesn't include CMake configuration files, so we need to manually specify paths and use the MPI compilers.

Create a `CMakeLists.txt` file:

```cmake
cmake_minimum_required(VERSION 3.16)

# Set compilers before project() - use MPI compilers
set(CMAKE_C_COMPILER "/data/shubin/opt/spack/store/llvm-18.1.3/mpich-4.3.1-ytw5t77wuontb3t4eusjri7wlq73xciw/bin/mpicc")
set(CMAKE_CXX_COMPILER "/data/shubin/opt/spack/store/llvm-18.1.3/mpich-4.3.1-ytw5t77wuontb3t4eusjri7wlq73xciw/bin/mpicxx")

project(MyMFEMProject)

# Set MFEM installation directory (using CPU build)
set(MFEM_DIR "/data/shubin/opt/spack/store/gcc-15.2.0/mfem-4.8.1-boqsa3oms4wytmhnbzwrtfxqw2gt75bm" CACHE PATH "MFEM installation directory")

# Check if MFEM installation exists
if(NOT EXISTS ${MFEM_DIR}/include/mfem.hpp)
    message(FATAL_ERROR "MFEM not found at ${MFEM_DIR}. Please check the MFEM_DIR path.")
endif()

# Add your executable (change the source file as needed)
add_executable(my_program my_program.cpp)

# Set C++ standard
set_property(TARGET my_program PROPERTY CXX_STANDARD 17)

# Set include directories - MFEM and all its dependencies
target_include_directories(my_program PRIVATE 
    ${MFEM_DIR}/include
    # External library includes (from config.mk MFEM_TPLFLAGS)
    /data/shubin/opt/spack/store/gcc-15.2.0/hypre-2.33.0-cd6txaynkragykxtkdel6x53frvnphyh/include
    /data/shubin/opt/spack/store/gcc-15.2.0/mumps-5.8.1-o6ot7mab7oy5mghksb2afvt4r3nbqotl/include
    /data/shubin/opt/spack/store/llvm-18.1.3/metis-5.1.0-6auhyeng75yqkadngsfj66hshxzstlew/include
    /data/shubin/opt/spack/store/llvm-18.1.3/libfms-0.2.0-puu6cm6q3uru6k33dhgfybmxhfyjkhma/include
    /data/shubin/opt/spack/store/gcc-15.2.0/suite-sparse-7.8.3-n3tkghfxjiaoneiffre2ifhqdgqxcya6/include
    /data/shubin/opt/spack/store/llvm-18.1.3/netcdf-c-4.9.3-bajkwx26gn6ekx2cwlhgustbe37zej4z/include
    /data/shubin/opt/spack/store/gcc-15.2.0/petsc-3.24.0-pl2etcmcpjl3xkqalhtoex4omra5l5sl/include
    /data/shubin/opt/spack/store/llvm-18.1.3/mpfr-4.2.1-gapwxxowe6bta2gqpwtpf5i3qua4bunj/include
    /data/shubin/opt/spack/store/llvm-18.1.3/gslib-1.0.9-ynzipi2fjfkeyn4i5lzj43dgyknnzk2w/include
    /data/shubin/opt/spack/store/gcc-15.2.0/intel-oneapi-mkl-2024.2.2-qe2jd2xh7fwecityp3jdh4z6wdp2p3j7/mkl/latest/include
    /data/shubin/opt/spack/store/llvm-18.1.3/zlib-ng-2.2.4-yu364ogrhkzpy5z6mglcrulto3yfjnda/include
)

# Link with MFEM library
target_link_libraries(my_program ${MFEM_DIR}/lib/libmfem.so)
```

**Key Points for CMake with Spack-installed MFEM:**

1. **MPI Compilers**: You must set the MPI compilers before the `project()` command
2. **Include Directories**: All external library include paths must be explicitly listed
3. **Standard**: Set C++17 standard to match MFEM's requirements
4. **Library Linking**: The MFEM shared library contains all dependencies, so you typically only need to link with `libmfem.so`

Build with:
```bash
mkdir build && cd build
cmake ..
make
```

## Running MFEM Programs

### Sequential Programs
```bash
./ex25 -m ../data/beam-tri.mesh
```

### Parallel Programs (MPI)
```bash
# Run on 4 processes
mpirun -np 4 ./ex25p -m ../data/beam-tri.mesh

# The installed config.mk sets default values:
# MFEM_MPIEXEC = mpirun
# MFEM_MPIEXEC_NP = -np  
# MFEM_MPI_NP = 4
```

### GPU Programs (CUDA build only)
```bash
# Run with GPU acceleration
./ex25 -d cuda -m ../data/beam-tri.mesh

# For parallel GPU runs
mpirun -np 2 ./ex25p -d cuda -m ../data/beam-tri.mesh
```

## Available Features by Build

### CUDA Build Features
- ✅ MPI (MPICH 4.3.1)
- ✅ CUDA 12.9.0 (SM_70)
- ✅ OpenMP
- ✅ HYPRE
- ✅ MUMPS
- ✅ METIS
- ✅ PETSc 3.24.0
- ✅ SuiteSparse
- ✅ Intel MKL (PARDISO)
- ✅ NetCDF
- ✅ MPFR
- ✅ GSLib
- ✅ FMS
- ✅ ZLIB

### CPU Build Features  
- ✅ MPI (MPICH)
- ❌ CUDA
- ❌ OpenMP
- ✅ HYPRE
- ✅ MUMPS
- ✅ METIS
- ✅ PETSc
- ✅ SuiteSparse
- ✅ NetCDF
- ✅ MPFR
- ✅ GSLib
- ✅ FMS
- ✅ ZLIB

## Accessing Example Programs and Data

### Pre-compiled Examples
Both installations include pre-compiled examples in `share/mfem/examples/`:
```bash
ls $MFEM_DIR/share/mfem/examples/ex*
```

### Example Source Code
Copy example source files to your working directory:
```bash
cp $MFEM_DIR/share/mfem/examples/ex25.cpp ./
cp $MFEM_DIR/share/mfem/examples/ex25p.cpp ./
```

### Mesh Data Files
Access mesh files for testing:
```bash
ls $MFEM_DIR/share/mfem/data/
# Files include: beam-tri.mesh, beam-quad.mesh, star.mesh, etc.
```

## Working with Your Own Code

1. **Include MFEM headers**:
   ```cpp
   #include "mfem.hpp"
   using namespace mfem;
   ```

2. **Initialize MPI (for parallel programs)**:
   ```cpp
   MPI_Init(&argc, &argv);
   // Your MFEM code here
   MPI_Finalize();
   ```

3. **Use MFEM device support (CUDA build)**:
   ```cpp
   Device device("cuda"); // or "cpu", "debug"
   // GPU-accelerated computations
   ```

## Troubleshooting

### Common Issues

1. **Library not found errors**:
   - Ensure `LD_LIBRARY_PATH` includes the MFEM lib directory
   - Use the full RPATH settings from config.mk

2. **Compilation errors**:
   - Check that you're using the correct compiler (nvcc for CUDA build, mpicxx for CPU build)
   - Verify all required flags from `MFEM_FLAGS` are included

3. **Runtime errors with MPI**:
   - Ensure MPI is properly initialized
   - Check that the number of processes matches your mesh partitioning

4. **CUDA errors (CUDA build only)**:
   - Verify GPU compatibility with SM_70 architecture
   - Check CUDA runtime installation

### Getting Help

- **MFEM Documentation**: https://mfem.org
- **MFEM Examples**: Study the provided examples in `share/mfem/examples/`
- **Configuration Details**: Check `share/mfem/config.mk` for all compile flags and libraries

## Example Workflow

Here's a complete example of building and running a custom MFEM program:

```bash
# 1. Set up working directory
mkdir my_mfem_project
cd my_mfem_project

# 2. Copy example source and makefile
export MFEM_DIR="/data/shubin/opt/spack/store/llvm-18.1.3/mfem-4.8.1-3lrmot6fydgnecja7lr2626ecsgc52n5"
cp $MFEM_DIR/share/mfem/examples/ex25.cpp ./
cp $MFEM_DIR/share/mfem/examples/makefile ./

# 3. Edit makefile to set correct MFEM_INSTALL_DIR
sed -i "s|MFEM_INSTALL_DIR.*|MFEM_INSTALL_DIR = $MFEM_DIR|" makefile

# 4. Build the program
make ex25

# 5. Run the program
./ex25 -m $MFEM_DIR/share/mfem/data/beam-tri.mesh

# 6. For parallel execution
make ex25p
mpirun -np 4 ./ex25p -m $MFEM_DIR/share/mfem/data/beam-tri.mesh
```

This workflow demonstrates the complete process from setup to execution using Spack-installed MFEM libraries.
