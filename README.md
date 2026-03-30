# BitmapArranger Stress Test Harness

Standalone test harness for the concave bitmap arrangement algorithm in [OrcaSlicer](https://github.com/SoftFever/OrcaSlicer). Compiles outside the OrcaSlicer build tree, linking against a pre-built `libslic3r.lib`. No OrcaSlicer source code is redistributed — this project only contains test scenarios and build configuration.

## What it does

Exercises `BitmapArranger::arrange()` with synthetic part geometries (rectangles, L-shapes, crescents, T-shapes, U-channels, circles) in realistic combinations: tight packing, overflow, mixed sizes, edge cases, and stress loads. Reports placement counts, plate usage, timing, and overlap detection.

## Prerequisites

1. **OrcaSlicer source and build** — You need a built copy of the `feature/concave-arrange` branch:
   ```
   git clone https://github.com/thereprocase/OrcaSlicer.git
   cd OrcaSlicer
   git checkout feature/concave-arrange
   ```
   Follow OrcaSlicer's build instructions to compile deps and the slicer. On Windows:
   ```
   build_release_vs.bat
   ```
   The harness needs these artifacts from the OrcaSlicer build:
   - `build/src/libslic3r/Release/libslic3r.lib`
   - `deps/build/OrcaSlicer_dep/usr/local/` (Boost, TBB, GMP, MPFR, OpenSSL, zlib, etc.)
   - `deps/build/dep_OCCT-prefix/src/dep_OCCT-build/lib/occt/` (OpenCASCADE import libs)
   - Source headers in `src/` and `deps_src/`

2. **Build tools:**
   - CMake 3.13+ (`pip install cmake` works)
   - Ninja (`pip install ninja` works)
   - MSVC (Visual Studio 2022 Build Tools or full VS)

## Setup

### 1. Clone this repo

```
git clone <this-repo-url>
cd OrcaTestHarness
```

### 2. Set the OrcaSlicer path

Edit `CMakeLists.txt` line 8 — set `ORCA_ROOT` to your OrcaSlicer source directory:

```cmake
set(ORCA_ROOT "F:/your/path/to/OrcaSlicer")
```

### 3. Build

**Windows (CMD or PowerShell):**
```
build.bat
```

Or manually:
```
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64
mkdir build && cd build
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
```

### 4. Run

```
run.bat --all
```

Or run individual scenarios:
```
run.bat --scenario stress_50_mixed
```

## Usage

```
arrange_stress [--all] [--scenario NAME] [--verbose]
```

- `--all` — Run all scenarios (default if no args)
- `--scenario NAME` — Run one named scenario
- `--verbose` — Extra output (not yet implemented)

Results are logged to `logs/stress_YYYYMMDD_HHMMSS.log`.

## Scenarios

| Category | Scenario | What it tests |
|----------|----------|---------------|
| Sanity | `single_rect`, `empty_input`, `single_oversized` | Basic placement, empty input, too-large item |
| Packing | `6_rects_one_plate`, `tight_pack_25_squares`, `tight_pack_26_squares` | Density and overflow boundary |
| Concave | `L_shapes_interlock`, `crescents_nest`, `T_shapes`, `U_channels` | Concave shape nesting advantage |
| Competition | `large_blocks_small_bed`, `ordering_matters`, `barely_fits`, `barely_doesnt_fit` | Size conflicts, BLF ordering, boundary cases |
| Multi-plate | `overflow_3_plates` | Overflow to multiple plates |
| Mixed | `mixed_sizes` | Heterogeneous part sizes |
| Stress | `stress_100_rects`, `stress_50_mixed` | Volume and timing |
| Edge | `zero_area_item`, `very_thin_item`, `no_spacing` | Degenerate input handling |

## Adding scenarios

Add a new entry in `build_scenarios()` in `src/arrange_stress.cpp`:

```cpp
scenarios["my_test"] = [=]() {
    std::vector<ArrangePolygon> items;
    items.push_back(make_ap(make_rect(50, 50)));
    // ... add more items
    return run_scenario("my_test", items, bed, params, true, "description");
};
```

Rebuild and run: `build.bat && run.bat --scenario my_test`

## Known limitations

- **No concave_triangles** — Test shapes use ExPolygon (Clipper rasterization path), not the concave triangle path used in real OrcaSlicer. Tests that expect concave interlocking (L-shapes on small beds) may underperform.
- **Overlap checker is approximate** — Uses vertex sampling, not full bitmap intersection. May produce false positives on shapes with few vertices (circles, crescents).
- **Windows only** — Build scripts assume MSVC. CMakeLists.txt is cross-platform in principle but only tested on Windows.

## License

This test harness contains no OrcaSlicer source code. It links against `libslic3r.lib` at build time, which is licensed under AGPLv3. See [OrcaSlicer's license](https://github.com/SoftFever/OrcaSlicer/blob/main/LICENSE.txt) for details.

Test harness code (`src/arrange_stress.cpp`, `CMakeLists.txt`, build scripts) is provided under the same AGPLv3 license for compatibility.
