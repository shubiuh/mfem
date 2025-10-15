#!/bin/bash
source ../spack-build-env.txt
export LD_LIBRARY_PATH="/data/shubin/opt/spack/store/none-none/cuda-12.9.0-xxvugy5yp7rx7hcg3xc5qthbgzr3ahgi/lib64:$LD_LIBRARY_PATH"
make ex25p

mpirun -np 3 ./ex25p
mpirun -np 3 ./ex25p -pa -d cuda

# make clean 

exit 0


### Instructions to build and run MFEM examples using Spack-installed MFEM (NOT RUNNING THIS SCRIPT)##########################################################
# example workflow
# 1. Set up working directory
mkdir my_mfem_project
cd my_mfem_project

# 2. Copy example source and makefile
export MFEM_DIR="/data/shubin/opt/spack/store/gcc-15.2.0/mfem-4.8.1-boqsa3oms4wytmhnbzwrtfxqw2gt75bm" # CPU version
export MFEM_DIR="/data/shubin/opt/spack/store/llvm-18.1.3/mfem-4.8.1-3lrmot6fydgnecja7lr2626ecsgc52n5" # GPU version
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
##########################################################################################################################################################
##########################################################################################################################################################