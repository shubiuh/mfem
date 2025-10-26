## MFEM Docker Image build instructions
if MUMPS cannot be linked due to MPI::MPI_C, replace it with
/home/pymfem/opt/spack/opt/spack/linux-skylake/intel-oneapi-mkl-2024.2.2-kthrgl5xg6j54ajsoj6ihva2buqazgvt/compiler/2024.2/lib/libiomp5.so
in the CMakeCache.txt file inside the build directory.

To build the MFEM Docker image, use the following command in the terminal from the `docker/pymfem_build` directory:

```bash
docker build -t pymfem .
```

### mfem repo and build
The MFEM repository is cloned from https://github.com/shubiuh/mfem.git and checeked out to v48_dev branch.
Edit the `config/user.cmake` file to enable AmgX, CUDA, MPI support, etc and the paths to the libraries.
Then build MFEM with the following commands:
```bash
spack load cuda@12.8.1
spack load intel-oneapi-mkl
spack load intel-oneapi-mpi
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/docker_ws/thirdparty/mfem_libs/mfem_pcuda ..
make -j16
```

IF AMGX is enabled, make sure to set the correct path to MPI (use intel-oneapi-mpi and intel-oneapi-mkl) for all libs and cuda@12.8.1 for cuda related libs.