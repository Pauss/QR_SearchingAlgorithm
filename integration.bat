@ECHO OFF 
:: This batch file run QR_SearchingAlgorithm and creates a graphic based on output data

ECHO ============================
ECHO Run QR_SearchingAlgorithm
ECHO ============================ 
cd Debug 
QR_SearchingAlgorithm.exe
ECHO ============================
ECHO Operation successfully ended
ECHO ============================
cd ../python_script
ECHO Run graphics
python graphics.py
ECHO ============================
ECHO Operation successfully ended
ECHO ============================

