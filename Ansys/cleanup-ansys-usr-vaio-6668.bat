@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="usr-vaio" (taskkill /f /pid 5232)
if /i "%LOCALHOST%"=="usr-vaio" (taskkill /f /pid 6668)

del /F cleanup-ansys-usr-vaio-6668.bat
