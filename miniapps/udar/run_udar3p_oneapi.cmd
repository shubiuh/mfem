@echo off
setlocal

if defined ONEAPI_ROOT (
  set "UDAR_ONEAPI_SETUP=%ONEAPI_ROOT%setvars.bat"
) else (
  set "UDAR_ONEAPI_SETUP=C:\OPT\Intel\oneAPI\setvars.bat"
)

if not exist "%UDAR_ONEAPI_SETUP%" (
  echo oneAPI setup script was not found: "%UDAR_ONEAPI_SETUP%" 1>&2
  echo Set ONEAPI_ROOT to the Intel oneAPI installation directory. 1>&2
  exit /b 2
)

call "%UDAR_ONEAPI_SETUP%" --force
if errorlevel 1 exit /b %errorlevel%

if not defined UDAR_MPI_NP set "UDAR_MPI_NP=1"
set "UDAR_REPOSITORY=%~dp0..\.."
set "UDAR_EXECUTABLE=%UDAR_REPOSITORY%\build\oneapi-mpi-mumps\miniapps\udar\udar3p.exe"
if not exist "%UDAR_EXECUTABLE%" (
  echo UDAR executable was not found: "%UDAR_EXECUTABLE%" 1>&2
  exit /b 3
)

rem -localroot launches the local rank directly and inherits the initialized
rem oneAPI PATH, avoiding impi.dll failures from a stale Hydra service.
mpiexec -localroot -genv PATH "%PATH%" -n %UDAR_MPI_NP% ^
  "%UDAR_EXECUTABLE%" %*
exit /b %errorlevel%
