#!/bin/sh
set -o errexit

echo "============Running with input:============"
echo input.txt
echo "==========================================="

echo "Compiling fortran..."
gfortran -O3 lstat_hstat_test.f90 -o test_fort.exe
chmod +x $PWD/test_fort.exe
echo "Done"
echo "==========================================="

echo "Python output:"
echo "==========================================="
python3  hstat.py
echo "==========================================="
echo "Fortran output:"
echo "==========================================="
$PWD/test_fort.exe
echo "==========================================="
echo "Cleaning up"
rm $PWD/test_fort.exe
