source /home/shubin/spack/share/spack/setup-env.sh
spack load intel-oneapi-mpi/5hhe
spack load intel-oneapi-compilers
make -j
make -j examples