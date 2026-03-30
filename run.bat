@echo off
cd /d F:\Claude\OrcaTestHarness
set PATH=F:\Claude\OrcaSlicer\build\OrcaSlicer;F:\Claude\OrcaSlicer\deps\build\dep_OCCT-prefix\src\dep_OCCT-build\bin\occt;%PATH%
build\arrange_stress.exe %*
