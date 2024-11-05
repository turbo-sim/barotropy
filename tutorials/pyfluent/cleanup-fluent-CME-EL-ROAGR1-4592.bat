echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v242\fluent/ntbin/win64/winkill.exe"

start "tell.exe" /B "C:\PROGRA~1\ANSYSI~1\v242\fluent\ntbin\win64\tell.exe" CME-EL-ROAGR1 58494 CLEANUP_EXITING
timeout /t 1
"C:\PROGRA~1\ANSYSI~1\v242\fluent\ntbin\win64\kill.exe" tell.exe
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 28688) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 6176) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 19724) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 28448) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 28396) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 27312) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 27824) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 27620) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 3192) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 27488) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 4592) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 19612)
del "C:\Users\roagr\OneDrive - Danmarks Tekniske Universitet\Roberto\Projects\barotropy\tutorials\pyfluent\cleanup-fluent-CME-EL-ROAGR1-4592.bat"
