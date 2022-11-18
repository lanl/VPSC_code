@REM This file should be named compileVPSC8.bat
gfortran -fdefault-real-8 -std=legacy vpsc8dim.for vpsc8main.for vpsc8sub.for -o vpsc8.exe
@del *.mod
