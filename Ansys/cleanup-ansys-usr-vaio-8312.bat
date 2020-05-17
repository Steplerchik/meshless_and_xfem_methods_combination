@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="usr-vaio" (taskkill /f /pid 9440)
if /i "%LOCALHOST%"=="usr-vaio" (taskkill /f /pid 8312)

del /F cleanup-ansys-usr-vaio-8312.bat
