cmake_minimum_required(VERSION 3.13)
project(OrcaArrangeTests CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Paths to OrcaSlicer source and built deps
set(ORCA_ROOT "F:/Claude/OrcaSlicer")
set(ORCA_BUILD "${ORCA_ROOT}/build")
set(ORCA_DEPS "${ORCA_ROOT}/deps/build/OrcaSlicer_dep/usr/local")

# Include paths
include_directories(
    ${ORCA_ROOT}/src
    ${ORCA_ROOT}/src/libslic3r
    ${ORCA_ROOT}/deps_src
    ${ORCA_ROOT}/deps_src/Clipper2/Clipper2Lib/include
    ${ORCA_ROOT}/deps_src/Eigen
    ${ORCA_DEPS}/include
    ${ORCA_DEPS}/include/boost-1_84
    ${ORCA_DEPS}/include/occt
    ${ORCA_BUILD}/src/libslic3r  # libslic3r_version.h
    ${ORCA_BUILD}/src            # other generated headers
)

# Boost setup
set(BOOST_ROOT ${ORCA_DEPS})
set(Boost_USE_STATIC_LIBS ON)
add_definitions(-DBOOST_ALL_NO_LIB)

# Link directories
link_directories(
    ${ORCA_BUILD}/src/libslic3r/Release
    ${ORCA_DEPS}/lib
)

# The test harness executable
add_executable(arrange_stress src/arrange_stress.cpp)

target_link_libraries(arrange_stress
    libslic3r
    tbb12
    # Boost
    libboost_log-vc144-mt-x64-1_84
    libboost_log_setup-vc144-mt-x64-1_84
    libboost_filesystem-vc144-mt-x64-1_84
    libboost_thread-vc144-mt-x64-1_84
    libboost_nowide-vc144-mt-x64-1_84
    libboost_locale-vc144-mt-x64-1_84
    libboost_iostreams-vc144-mt-x64-1_84
    # CGAL deps
    libgmp-10
    libmpfr-4
    # Crypto/network
    libcrypto
    libssl
    libcurl
    # Compression/image
    zlib
    libpng16_static
    expat
    # OCCT (OpenCASCADE) — libslic3r links STEP support
    ${ORCA_ROOT}/deps/build/dep_OCCT-prefix/src/dep_OCCT-build/lib/occt/TKernel.lib
    ${ORCA_ROOT}/deps/build/dep_OCCT-prefix/src/dep_OCCT-build/lib/occt/TKCDF.lib
    ${ORCA_ROOT}/deps/build/dep_OCCT-prefix/src/dep_OCCT-build/lib/occt/TKCAF.lib
    ${ORCA_ROOT}/deps/build/dep_OCCT-prefix/src/dep_OCCT-build/lib/occt/TKLCAF.lib
    tbbmalloc
    # Windows
    Ws2_32
    Crypt32
    Wldap32
    Normaliz
)

target_compile_definitions(arrange_stress PRIVATE
    SLIC3R_LOGLEVEL=3
    _USE_MATH_DEFINES
    NOMINMAX
    WIN32_LEAN_AND_MEAN
    _CRT_SECURE_NO_WARNINGS
)
