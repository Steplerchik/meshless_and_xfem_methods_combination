@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="usr-vaio" (taskkill /f /pid 7728)
if /i "%LOCALHOST%"=="usr-vaio" (taskkill /f /pid 4052)

del /F cleanup-ansys-usr-vaio-4052.bat
