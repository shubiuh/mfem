## MFEM Docker Image
if MUMPS cannot be linked due to MPI::MPI_C, replace it with
/home/pymfem/opt/spack/opt/spack/linux-skylake/intel-oneapi-mkl-2024.2.2-kthrgl5xg6j54ajsoj6ihva2buqazgvt/compiler/2024.2/lib/libiomp5.so
in the CMakeCache.txt file inside the build directory.

To build the MFEM Docker image, use the following command in the terminal from the `docker/pymfem_build` directory:

```bash
docker build -t pymfem .
```