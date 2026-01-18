# on cot42, only build in the main spack no environment is set up
source /export/home/mnle8/01_szeng/oneapi.sh
# copy my personal mfem package file
cp package.py /export/home/mnle8/.spack/package_repos/fncqgg4/repos/spack_repo/builtin/packages/mfem/package.py # on cot-42

# below are the built commands for mfem (ON cot42)
# CPU +MPI version (MFEM v4.8.2)
spack -v install --keep-stage mfem+mpi

# GPU +MPI version (MFEM v4.8.2)
spack -v install --keep-stage mfem+cuda+mpi+mumps+netcdf~occa~miniapps cuda_arch=70 cxxstd=17


# CPU + MPI version (v4.8.1)
spack install mfem@4.8.1~amgx~conduit~cuda~debug+examples~exceptions+fms~ginkgo~gnutls+gslib~hiop~lapack~legacy-openmp~libceed~libunwind+metis~miniapps~mklcpardiso+mklpardiso+mpfr+mpi+mumps+netcdf~occa~openmp+petsc~pumi~raja~rocm+shared~slepc~static~strumpack+suite-sparse~sundials~superlu-dist~threadsafe~umpire+zlib build_system=generic cxxstd=17 precision=double timer=auto ^mpfr@3.1.6%clang@19.1.7 ^gslib+shared %clang@19.1.7

# GPU + MPI version (v4.8.1)
spack install mfem@4.8.1~amgx~conduit+cuda~debug+examples+exceptions+fms~ginkgo~gnutls+gslib~hiop~lapack~legacy-openmp+libceed~libunwind+metis~miniapps~mklcpardiso+mklpardiso+mpfr+mpi+mumps+netcdf~occa~openmp+petsc~pumi~raja~rocm+shared~slepc~static~strumpack+suite-sparse~sundials~superlu-dist~threadsafe~umpire+zlib build_system=generic cuda_arch=70 cxxstd=17 precision=double timer=auto ^mpfr@3.1.6%clang@19.1.7 ^gslib+shared