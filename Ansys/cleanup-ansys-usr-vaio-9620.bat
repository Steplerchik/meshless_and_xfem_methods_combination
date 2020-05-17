@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="usr-vaio" (taskkill /f /pid 15960)
if /i "%LOCALHOST%"=="usr-vaio" (taskkill /f /pid 9620)

del /F cleanup-ansys-usr-vaio-9620.bat
