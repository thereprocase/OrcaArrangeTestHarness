@echo off
cd /d F:\Claude\OrcaTestHarness
set PATH=F:\Claude\OrcaSlicer\build\OrcaSlicer;F:\Claude\OrcaSlicer\deps\build\dep_OCCT-prefix\src\dep_OCCT-build\bin\occt;%PATH%

REM Cap TBB threads inside the arranger to 4 cores
set TBB_NUM_THREADS=4

build\arrange_stress.exe %*
