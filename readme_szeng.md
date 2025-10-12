## Using MFEM with the Windows Subsystem for Linux

## Create MSVC project
use vslinux package to create visual studio project and filters.
See https://github.com/dealii/dealii/wiki/Windows for details.

## Note
Currently the gdb-oneapi cannot work in vs2022.
If built using intel-oneapi-compilers and/or intel-oneapi-mpi, use vscode debugger and set tasks.json and launch.json to debug.
Will test gcc gdb.