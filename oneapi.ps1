$Env:CC="icx"
$Env:FC="ifx"
$Env:CXX="icx"

cmd.exe /c '"%PROGRAMFILES(X86)%\Intel\oneAPI\setvars.bat" intel64 && pwsh'
cmd.exe /c '"C:\OPT\Intel\oneAPI\setvars.bat" intel64 && pwsh'

# Continue with the rest of your PowerShell script
Write-Host "OneAPI environment initialized. Continuing with the script..."