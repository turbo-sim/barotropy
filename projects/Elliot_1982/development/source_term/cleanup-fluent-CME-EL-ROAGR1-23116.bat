echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v231\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\v231\fluent\ntbin\win64\tell.exe" CME-EL-ROAGR1 64964 CLEANUP_EXITING
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 1876) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 23116) 
if /i "%LOCALHOST%"=="CME-EL-ROAGR1" (%KILL_CMD% 12904)
del "C:\Users\roagr\OneDrive - Danmarks Tekniske Universitet\Roberto\Projects\barotropic-model\validation\Elliot_1979\source_term\cleanup-fluent-CME-EL-ROAGR1-23116.bat"
