echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v231\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\v231\fluent\ntbin\win64\tell.exe" CME-EL-ROAGR1 50756 CLEANUP_EXITING
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 15972) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 20872) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 20616) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 20512) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 21000) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 18580) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 20868) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 17572) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 7856) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 20684) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 16560) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 19800) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 20040) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 17432)
del "C:\Users\roagr\OneDrive - Danmarks Tekniske Universitet\Roberto\Projects\barotropic-model\validation\Elliot_1979\mixture_model\cleanup-fluent-CME-EL-ROAGR1-20040.bat"
