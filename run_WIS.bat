@echo off

g++ -std=c++11 -O2 WIS.cpp -o WIS.exe

(
echo %1
echo %2
echo %3
) | WIS.exe

pause
