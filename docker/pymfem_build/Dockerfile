# syntax=docker/dockerfile:1.7
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# --- Base OS, locale, timezone ------------------------------------------------
RUN apt-get update -y && apt-get install -y --no-install-recommends \
      tzdata locales ca-certificates \
   && ln -fs /usr/share/zoneinfo/Europe/Madrid /etc/localtime \
   && dpkg-reconfigure -f noninteractive tzdata \
   && sed -i -e 's/# es_ES.UTF-8 UTF-8/es_ES.UTF-8 UTF-8/' /etc/locale.gen \
   && dpkg-reconfigure --frontend=noninteractive locales \
   && update-locale LANG=es_ES.UTF-8 \
   && rm -rf /var/lib/apt/lists/*
ENV LANG=es_ES.UTF-8

# --- Build deps (pyenv + toolchain + OpenMPI) ---------------------------------
RUN apt-get update -y && apt-get install -y --no-install-recommends \
      build-essential make cmake git curl xz-utils \
      libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev \
      libffi-dev liblzma-dev tk-dev libedit-dev \
      openmpi-bin libopenmpi-dev \
      swig \
      python3 python3.12-venv \
   && rm -rf /var/lib/apt/lists/*

# --- Non-root user ------------------------------------------------------------
RUN useradd -m -s /bin/bash pymfem
USER pymfem
WORKDIR /home/pymfem
ENV HOME=/home/pymfem

# --- pyenv install ------------------------------------------------------------
RUN curl -fsSL https://pyenv.run | bash
ENV PYENV_ROOT="$HOME/.pyenv"
ENV PATH="$PYENV_ROOT/bin:$PYENV_ROOT/shims:$PATH"

# Use bash -lc for subsequent RUN
SHELL ["/bin/bash","-lc"]

# Install Python 3.10.17 with pyenv and basic tools
RUN echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc && \
    echo 'export PATH="$PYENV_ROOT/bin:$PYENV_ROOT/shims:$PATH"' >> ~/.bashrc && \
    echo 'eval "$(pyenv init -)"' >> ~/.bashrc && \
    pyenv update && \
    pyenv install -s 3.10.17 && \
    pyenv global 3.10.17 && \
    python --version && \
    python -m ensurepip --upgrade || true && \
    python -m pip install --upgrade pip setuptools wheel

# --- Python packages (single layer) -------------------------------------------
# (mpi4py will use the system OpenMPI we installed above)
RUN python -m pip install -U \
      numpy==1.26.4 \
      scipy \
      numba \
      numba-scipy \
      ipykernel \
      jupyter \
      mpi4py \
      glvis

# --- Install Spack -------------------------------------------------------------
RUN mkdir -p /home/pymfem/opt && \
    git clone -c feature.manyFiles=true https://github.com/spack/spack.git /home/pymfem/opt/spack && \
    echo 'export SPACK_ROOT="/home/pymfem/opt/spack"' >> ~/.bashrc && \
    echo 'source $SPACK_ROOT/share/spack/setup-env.sh' >> ~/.bashrc

# --- Add custom MFEM package.py ------------------------------------------------
# (make sure package.py exists next to this Dockerfile)
ADD package.py /home/pymfem/package.py
RUN mkdir -p /home/pymfem/.spack/package_repos/fncqgg4/repos/spack_repo/builtin/packages/mfem && \
    cp /home/pymfem/package.py /home/pymfem/.spack/package_repos/fncqgg4/repos/spack_repo/builtin/packages/mfem/package.py

# --- Install Intel OneAPI components with Spack -------------------------------
RUN source /home/pymfem/opt/spack/share/spack/setup-env.sh && \
    spack install intel-oneapi-compilers@2025.2.1 && \
    spack install intel-oneapi-mkl@2024.2.2 && \
    spack install intel-oneapi-mpi@2021.16.1

# --- Install MFEM with Spack ---------------------------------------------------
RUN source /home/pymfem/opt/spack/share/spack/setup-env.sh && \
    spack install mfem@4.8.1~amgx~conduit+cuda~debug+examples+exceptions+fms~ginkgo~gnutls+gslib~hiop~lapack~legacy-openmp+libceed~libunwind+metis~miniapps~mklcpardiso+mklpardiso+mpfr+mpi+mumps+netcdf~occa~openmp+petsc~pumi~raja~rocm+shared~slepc~static~strumpack+suite-sparse~sundials~superlu-dist~threadsafe~umpire+zlib build_system=generic cuda_arch=70 cxxstd=17 precision=double timer=auto

# --- Build & install PyMFEM (parallel) ----------------------------------------
RUN git clone https://github.com/mfem/PyMFEM.git && \
    cd PyMFEM && \
    python setup.py install --with-parallel

# --- Add your script -----------------------------------------------------------
# (make sure script.sh exists next to this Dockerfile)
ADD script.sh /home/pymfem/script.sh
RUN chmod +x /home/pymfem/script.sh

# Default shell
SHELL ["/bin/bash","-lc"]


