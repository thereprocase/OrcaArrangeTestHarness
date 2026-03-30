# Claude Code Guide — OrcaTestHarness

## What this is
Standalone stress test harness for BitmapArranger. Links against a pre-built
OrcaSlicer `libslic3r.lib`. Does NOT contain or redistribute OrcaSlicer source.

## Build
```
build.bat         # Windows: sets up MSVC env, runs cmake + ninja
```
Requires: ORCA_ROOT in CMakeLists.txt pointing to a built OrcaSlicer checkout.

## Run
```
run.bat --all             # all scenarios
run.bat --scenario NAME   # one scenario
```

## Adding tests
Add entries in `build_scenarios()` in `src/arrange_stress.cpp`.
Shape builders: `make_rect`, `make_circle`, `make_L`, `make_crescent`, `make_T`, `make_U`.
Use `run_scenario(name, items, bed, params, expect_all_placed, notes)`.

## Key constraint
Tests use ExPolygon shapes (Clipper path), NOT concave_triangles (the real
concave path). To test true concave behavior, populate ArrangePolygon.concave_triangles
with projected triangle vertices.

## Dependencies (not redistributed)
- OrcaSlicer source + built libs (AGPLv3)
- MSVC 2022 Build Tools
- CMake 3.13+, Ninja
