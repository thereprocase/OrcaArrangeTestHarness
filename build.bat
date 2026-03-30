@echo off
set PATH=F:\Claude\strawberry-perl\perl\bin;F:\Claude\strawberry-perl\c\bin;%PATH%
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64
set CMAKE_POLICY_VERSION_MINIMUM=3.5
cd /d F:\Claude\OrcaTestHarness
if not exist build mkdir build
cd build
cmake .. -G "Ninja" -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
