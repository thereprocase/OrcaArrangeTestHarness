@echo off
REM Build lock — shared with OrcaSlicer builds
set LOCKFILE=F:\Claude\.build_lock
if exist "%LOCKFILE%" (
    echo ERROR: Another build is already running. Lockfile: %LOCKFILE%
    exit /b 1
)
echo %DATE% %TIME% > "%LOCKFILE%"

set PATH=F:\Claude\strawberry-perl\perl\bin;F:\Claude\strawberry-perl\c\bin;%PATH%
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64
set CMAKE_POLICY_VERSION_MINIMUM=3.5
cd /d F:\Claude\OrcaTestHarness
if not exist build mkdir build
cd build
cmake .. -G "Ninja" -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release -j 4

del "%LOCKFILE%" 2>NUL
