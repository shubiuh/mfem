unset PYTHONHOME
unset PYTHONPATH
spack load intel-oneapi-compilers
spack load intel-oneapi-mkl
spack load intel-oneapi-mpi
source /home/shubin/spack/opt/spack/linux-skylake/intel-oneapi-compilers-2025.2.0-u5v2zoeq5ywjo2texzzrbn45g7dy4ech/setvars.sh
# source /home/shubin/spack/opt/spack/linux-skylake/intel-oneapi-compilers-2025.2.0-u5v2zoeq5ywjo2texzzrbn45g7dy4ech/debugger/latest/env/vars.sh
# export PYTHONHOME=/home/shubin/spack/opt/spack/linux-skylake/intel-oneapi-compilers-2025.2.0-u5v2zoeq5ywjo2texzzrbn45g7dy4ech/debugger/2025.2/opt/debugger/lib/
# export PYTHONPATH=$ONEAPI_DEBUGGER_ROOT/lib/python3.13

make ex1