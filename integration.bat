@ECHO OFF 
:: This batch file run QR_SearchingAlgorithm and creates a graphic based on output data

ECHO ============================
ECHO Run QR_SearchingAlgorithm
ECHO ============================ 
cd Debug

@echo off
SET /A "index = 1"
SET /A "count = 1"
:while
if %index% leq %count% (
   @ECHO OFF 
   ::echo The value of index is %index%
   SET /A "index = index + 1"
   QR_SearchingAlgorithm.exe %1 %2 %3 %4 %5
   goto :while
)

ECHO ============================
ECHO Operation successfully ended
ECHO ============================
 
::cd ../python_script
::ECHO Run graphics
::python graphics.py
ECHO ============================
ECHO Operation successfully ended
ECHO ============================

