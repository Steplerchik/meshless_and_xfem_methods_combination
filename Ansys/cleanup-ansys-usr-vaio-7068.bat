@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="usr-vaio" (taskkill /f /pid 17708)
if /i "%LOCALHOST%"=="usr-vaio" (taskkill /f /pid 7068)

del /F cleanup-ansys-usr-vaio-7068.bat
