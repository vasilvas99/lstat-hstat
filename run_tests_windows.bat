@echo off
setlocal enableextensions enabledelayedexpansion

echo ============Running with input:============
echo input.txt
echo ===========================================

echo Compiling fortran...
gfortran.exe -O3 fortran/lstat_hstat_test.f90 -o test_fort.exe
echo Done
echo ===========================================

echo Python output:
echo ===========================================
python hstat.py
echo ===========================================
echo Fortran output:
echo ===========================================
test_fort.exe
echo ===========================================
echo Cleaning up
del test_fort.exe
pause
endlocal
