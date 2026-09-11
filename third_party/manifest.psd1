@{
   # Change these only after rebuilding the complete stack, including MFEM.
   HypreRepository = 'https://github.com/hypre-space/hypre.git'
   HypreTag        = 'v2.33.0'
   MumpsRepository = 'https://github.com/scivision/mumps-superbuild.git'
   MumpsTag        = 'main'
   MetisArchive = 'https://mfem.github.io/tpls/metis-5.1.0.tar.gz'
   MetisSha256  = '76faebe03f6c963127dbb73c13eab58c9a3faeae48779f049066a21c087c5db2'

   # Optional standalone C++ inversion core (inversion/cpp).
   RolRepository  = 'https://github.com/trilinos/Trilinos.git'
   RolTag         = 'trilinos-release-17-2-1'
   RolCommit      = 'a0de2fe38d6b68f88fc4557abb003ff379561933'
   PugiRepository = 'https://github.com/zeux/pugixml.git'
   PugiTag        = 'v1.14'
   PugiCommit     = 'db78afc2b7d8f043b4bc6b185635d949ea2ed2a8'

   # Optional PETSc/TAO inversion core (inversion/cpp/udar_tao.cpp).
   PetscRepository = 'https://gitlab.com/petsc/petsc.git'
   PetscTag        = 'v3.25.5'
   PetscCommit     = 'a0b507cc1cc9eac739c34d6e71521f4854e84e36'
   Win32feRepository = 'https://bitbucket.org/petsc/win32fe.git'
   Win32feCommit     = 'f2ff43464a6f25d7ccd55fe164585feee2c8bd76'
}
