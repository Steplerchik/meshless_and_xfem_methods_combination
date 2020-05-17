@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="usr-vaio" (taskkill /f /pid 9028)
if /i "%LOCALHOST%"=="usr-vaio" (taskkill /f /pid 10056)

del /F cleanup-ansys-usr-vaio-10056.bat
