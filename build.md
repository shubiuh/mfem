1. use vcpkg to install openblas, lapack, suitesparse, metis. MUMPS is not available on vcpkg.
2. Or use conda to install mumps, suitesparse, metis, parmetis



Then use the cmake command to build mfem:

```powershell
cmake -G "Visual Studio 17 2022" `
  -T "Intel C++ Compiler 2025" `
  -DCMAKE_BUILD_TYPE=Release `
  -DBUILD_SHARED_LIBS=OFF `
  -DMFEM_PRECISION=double `
  -DMFEM_USE_OPENMP=NO `
  -DMFEM_THREAD_SAFE=YES `
  -DMFEM_USE_MEMALLOC=NO `
  -DMFEM_USE_GSLIB=ON `
  -DGSLIB_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\gslib" `
  -DMFEM_USE_MPI=ON `
  -DCMAKE_CXX_COMPILER="C:\OPT\Intel\oneAPI\mpi\2021.12\bin\mpiexec.exe" `
  -DMFEM_USE_LAPACK=OFF `
  -DBLA_VENDOR=OpenBLAS `
  -DBLAS_LIBRARIES="C:\vcpkg\packages\openblas_x64-windows\lib\openblas.lib" `
  -DLAPACK_LIBRARIES="C:\OPT\Intel\oneAPI\mkl\2025.2\lib\mkl_lapack95_lp64.lib" `
  -DMFEM_USE_METIS=ON `
  -DHYPRE_DIR="..\..\..\contrib\opt\hypre" `
  -DMFEM_USE_METIS_5=ON `
  -DMETIS_DIR="..\..\..\contrib\opt\metis\Library" `
  -DMFEM_USE_ZLIB=ON `
  -DZLIB_LIBRARY="C:\vcpkg\packages\zlib_x64-windows\lib\zlib.lib" `
  -DZLIB_INCLUDE_DIR="C:\vcpkg\packages\zlib_x64-windows\include" `
  -DParMETIS_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\ParMETIS" `
  -Dscalapack_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\scalapack_mkl" `
  -DMFEM_USE_MUMPS=OFF `
  -DMUMPS_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\mumps-5.6.1" `
  -DMFEM_USE_SUITESPARSE=ON `
  -DSuiteSparse_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\suitesparse_7.8.3_vcpkg" `
  -DMFEM_USE_NETCDF=ON `
  -DNETCDF_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\netCDF4.9.2" `
  -DHDF5_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\hdf5" `
  -DMFEM_USE_MKL_PARDISO=ON `
  -DMKL_PARDISO_DIR="C:\OPT\Intel/oneAPI/mkl/2025.2" `
  -DMFEM_USE_MKL_CPARDISO=ON `
  -DMKL_CPARDISO_DIR="C:\OPT\Intel/oneAPI/mkl/2025.2" `
  -DMKL_MPI_WRAPPER_LIB="mkl_blacs_intelmpi_lp64" `
  -DCMAKE_INSTALL_PREFIX="..\..\..\contrib\opt\build_trash_xps8940" `
  -DMFEM_ENABLE_TESTING=OFF `
  -DMFEM_USE_CUDA=OFF `
  ..

```


In cmd, it is
Release mode (if use Pardiso OMP)
```cmd
cmake -G "Visual Studio 17 2022" -T "Intel C++ Compiler 2025" -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DMFEM_PRECISION=double -DMFEM_USE_OPENMP=NO -DMFEM_THREAD_SAFE=YES -DMFEM_USE_MEMALLOC=NO -DMFEM_USE_GSLIB=ON -DGSLIB_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\gslib" -DMFEM_USE_MPI=ON -DCMAKE_CXX_COMPILER="C:\OPT\Intel\oneAPI\mpi\2021.12\bin\mpiexec.exe" -DMFEM_USE_LAPACK=OFF -DBLA_VENDOR=OpenBLAS -DBLAS_LIBRARIES="C:\vcpkg\packages\openblas_x64-windows\lib\openblas.lib" -DLAPACK_LIBRARIES="C:\OPT\Intel\oneAPI\mkl\2025.2\lib\mkl_lapack95_lp64.lib" -DMFEM_USE_METIS=ON -DHYPRE_DIR="..\..\..\contrib\opt\hypre" -DMFEM_USE_METIS_5=ON -DMETIS_DIR="..\..\..\contrib\opt\metis\Library" -DMFEM_USE_ZLIB=ON -DZLIB_LIBRARY="C:\vcpkg\packages\zlib_x64-windows\lib\zlib.lib" -DZLIB_INCLUDE_DIR="C:\vcpkg\packages\zlib_x64-windows\include" -DParMETIS_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\ParMETIS" -Dscalapack_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\scalapack_mkl" -DMFEM_USE_MUMPS=OFF -DMUMPS_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\mumps-5.6.1" -DMFEM_USE_SUITESPARSE=ON -DSuiteSparse_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\suitesparse_7.8.3_vcpkg" -DMFEM_USE_NETCDF=ON -DNETCDF_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\netCDF4.9.2" -DHDF5_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\hdf5" -DMFEM_USE_MKL_PARDISO=ON -DMKL_PARDISO_DIR="C:\OPT\Intel/oneAPI/mkl/2025.2" -DMFEM_USE_MKL_CPARDISO=OFF -DMKL_CPARDISO_DIR="C:\OPT\Intel/oneAPI/mkl/2025.2" -DMKL_MPI_WRAPPER_LIB="mkl_blacs_intelmpi_lp64" -DCMAKE_INSTALL_PREFIX="..\..\..\contrib\opt\build_trash_xps8940" -DMFEM_ENABLE_TESTING=OFF -DMFEM_USE_CUDA=OFF ..
```

If use Cluster Pardiso
```
cmake -G "Visual Studio 17 2022" -T "Intel C++ Compiler 2025" -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DMFEM_PRECISION=double -DMFEM_USE_OPENMP=YES -DMFEM_THREAD_SAFE=YES -DMFEM_USE_MEMALLOC=YES -DMFEM_USE_GSLIB=ON -DGSLIB_DIR="D:\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\FMCSEM\contrib\opt\gslib" -DMFEM_USE_MPI=ON -DCMAKE_CXX_COMPILER="C:\Program Files (x86)\Intel\oneAPI\mpi\2021.14\bin\mpiexec.exe" -DMFEM_USE_LAPACK=ON -DBLA_VENDOR=Intel10_64lp -DBLAS_LIBRARIES="C:\Program Files (x86)\Intel\oneAPI\mkl\2025.0\lib\mkl_blas95_lp64.lib" -DLAPACK_LIBRARIES="C:\Program Files (x86)\Intel\oneAPI\mkl\2025.0\lib\mkl_lapack95_lp64.lib" -DMFEM_USE_METIS=ON -DHYPRE_DIR="..\..\..\contrib\opt\hypre" -DMFEM_USE_METIS_5=ON -DMETIS_DIR="..\..\..\contrib\opt\metis\Library" -DMFEM_USE_ZLIB=ON -DZLIB_LIBRARY="C:\Users\shubi\Documents\vcpkg\packages\zlib_x64-windows\lib\zlib.lib" -DZLIB_INCLUDE_DIR="C:\Users\shubi\Documents\vcpkg\packages\zlib_x64-windows\include" -DParMETIS_DIR="D:\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\ParMETIS" -Dscalapack_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\scalapack_mkl" -DMFEM_USE_MUMPS=OFF -DMUMPS_DIR="D:\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\mumps-5.6.1" -DMFEM_USE_SUITESPARSE=ON -DSuiteSparse_DIR="D:\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\SuiteSparse5.4.0" -DMFEM_USE_NETCDF=ON -DNETCDF_DIR="D:\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\netCDF4.9.2" -DHDF5_DIR="D:\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\hdf5" -DMFEM_USE_MKL_CPARDISO=ON -DMKL_CPARDISO_DIR="C:/Program Files (x86)/Intel/oneAPI/mkl/latest" -DMKL_CPARDISO_MKL_MPI_WRAPPER_LIBRARY="C:\Program Files (x86)\Intel\oneAPI\mkl\latest\lib\mkl_blacs_msmpi_lp64.lib" -DCMAKE_INSTALL_PREFIX="..\..\..\contrib\opt\build_lib_lapack" -DMFEM_ENABLE_TESTING=OFF ..
```
Debug mode (if use Pardiso OMP)
```cmd
cmake -G "Visual Studio 17 2022" -T "Intel C++ Compiler 2025" -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=OFF -DMFEM_PRECISION=double -DMFEM_USE_OPENMP=NO -DMFEM_THREAD_SAFE=YES -DMFEM_USE_MEMALLOC=NO -DMFEM_USE_GSLIB=ON -DGSLIB_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\gslib" -DMFEM_USE_MPI=ON -DCMAKE_CXX_COMPILER="C:\OPT\Intel\oneAPI\mpi\2021.12\bin\mpiexec.exe" -DMFEM_USE_LAPACK=OFF -DBLA_VENDOR=OpenBLAS -DBLAS_LIBRARIES="C:\vcpkg\packages\openblas_x64-windows\lib\openblas.lib" -DLAPACK_LIBRARIES="C:\OPT\Intel\oneAPI\mkl\2025.2\lib\mkl_lapack95_lp64.lib" -DMFEM_USE_METIS=ON -DHYPRE_DIR="..\..\..\contrib\opt\hypre" -DMFEM_USE_METIS_5=ON -DMETIS_DIR="..\..\..\contrib\opt\metis\Library" -DMFEM_USE_ZLIB=ON -DZLIB_LIBRARY="C:\vcpkg\packages\zlib_x64-windows\lib\zlib.lib" -DZLIB_INCLUDE_DIR="C:\vcpkg\packages\zlib_x64-windows\include" -DParMETIS_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\ParMETIS" -Dscalapack_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\scalapack_mkl" -DMFEM_USE_MUMPS=OFF -DMUMPS_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\mumps-5.6.1" -DMFEM_USE_SUITESPARSE=ON -DSuiteSparse_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\suitesparse_7.8.3_vcpkg" -DMFEM_USE_NETCDF=ON -DNETCDF_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\netCDF4.9.2" -DHDF5_DIR="C:\Users\shubi\Dropbox\Cyentech_Projects\06_projects\05_WFD_logging\01_coding\thirdParty\hdf5" -DMFEM_USE_MKL_PARDISO=ON -DMKL_PARDISO_DIR="C:\OPT\Intel/oneAPI/mkl/2025.2" -DMFEM_USE_MKL_CPARDISO=OFF -DMKL_CPARDISO_DIR="C:\OPT\Intel/oneAPI/mkl/2025.2" -DMKL_MPI_WRAPPER_LIB="mkl_blacs_intelmpi_lp64" -DCMAKE_INSTALL_PREFIX="..\..\..\contrib\opt\build_trash_xps8940" -DMFEM_ENABLE_TESTING=OFF -DMFEM_USE_CUDA=OFF ..
```