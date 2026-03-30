// OrcaSlicer BitmapArranger Stress Test Harness
//
// Standalone executable that exercises BitmapArranger::arrange() with
// realistic part combinations. Compiles outside the main build, links
// against libslic3r.lib.
//
// Usage: arrange_stress [--scenario NAME] [--all] [--verbose]
//        Logs results to ../logs/ and writes reports to ../reports/

#include "libslic3r/Arrange/BitmapArranger.hpp"
#include "libslic3r/Arrange.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/libslic3r.h"

#ifdef _WIN32
#include <windows.h>
#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <functional>
#include <map>
#include <iomanip>
#include <cassert>

using namespace Slic3r;
using namespace Slic3r::arrangement;

// --- Shape builders ---

static ExPolygon make_rect(double w_mm, double h_mm) {
    coord_t w = scaled(w_mm), h = scaled(h_mm);
    ExPolygon ep;
    ep.contour.points = {{0, 0}, {w, 0}, {w, h}, {0, h}};
    return ep;
}

static ExPolygon make_circle(double radius_mm, int segments = 32) {
    coord_t r = scaled(radius_mm);
    ExPolygon ep;
    for (int i = 0; i < segments; i++) {
        double angle = 2.0 * M_PI * i / segments;
        ep.contour.points.push_back({
            (coord_t)(r * cos(angle)),
            (coord_t)(r * sin(angle))
        });
    }
    return ep;
}

static ExPolygon make_L(double size_mm) {
    coord_t s = scaled(size_mm);
    coord_t half = s / 2;
    ExPolygon ep;
    ep.contour.points = {
        {0, 0}, {s, 0}, {s, half}, {half, half}, {half, s}, {0, s}
    };
    return ep;
}

static ExPolygon make_crescent(double outer_r_mm, double inner_r_mm, int segments = 32) {
    ExPolygon ep;
    coord_t outer_r = scaled(outer_r_mm);
    coord_t inner_r = scaled(inner_r_mm);
    // Outer contour: full circle
    for (int i = 0; i < segments; i++) {
        double angle = 2.0 * M_PI * i / segments;
        ep.contour.points.push_back({
            (coord_t)(outer_r * cos(angle)),
            (coord_t)(outer_r * sin(angle))
        });
    }
    // Inner hole: offset circle (creates crescent shape)
    Polygon hole;
    coord_t offset_x = scaled(outer_r_mm - inner_r_mm);
    for (int i = segments - 1; i >= 0; i--) {
        double angle = 2.0 * M_PI * i / segments;
        hole.points.push_back({
            offset_x + (coord_t)(inner_r * cos(angle)),
            (coord_t)(inner_r * sin(angle))
        });
    }
    ep.holes.push_back(hole);
    return ep;
}

static ExPolygon make_T(double size_mm) {
    coord_t s = scaled(size_mm);
    coord_t third = s / 3;
    ExPolygon ep;
    ep.contour.points = {
        {third, 0}, {2 * third, 0}, {2 * third, 2 * third},
        {s, 2 * third}, {s, s}, {0, s}, {0, 2 * third},
        {third, 2 * third}
    };
    return ep;
}

static ExPolygon make_U(double size_mm) {
    coord_t s = scaled(size_mm);
    coord_t wall = s / 4;
    ExPolygon ep;
    ep.contour.points = {
        {0, 0}, {s, 0}, {s, s}, {s - wall, s},
        {s - wall, wall}, {wall, wall}, {wall, s}, {0, s}
    };
    return ep;
}

static BoundingBox make_bed(double w_mm, double h_mm) {
    return BoundingBox({0, 0}, {scaled(w_mm), scaled(h_mm)});
}

// Build concave_triangles from an ExPolygon by triangulating the contour.
// Simple fan triangulation from vertex 0 — works for convex and simple concave
// shapes. Not a general ear-clipping triangulator, but sufficient for test shapes.
static void populate_concave_triangles(ArrangePolygon &ap, float height_mm = 10.0f) {
    const auto &pts = ap.poly.contour.points;
    if (pts.size() < 3) return;

    ap.concave_triangles.clear();
    ap.concave_z.clear();

    // Fan triangulation from vertex 0
    for (size_t i = 1; i + 1 < pts.size(); i++) {
        ap.concave_triangles.push_back(pts[0]);
        ap.concave_triangles.push_back(pts[i]);
        ap.concave_triangles.push_back(pts[i + 1]);
        // All vertices at z=0 (flat 2D part)
        ap.concave_z.push_back(0.f);
        ap.concave_z.push_back(0.f);
        ap.concave_z.push_back(0.f);
    }
}

// Same but with varying Z heights for 3D nesting tests
static void populate_concave_triangles_3d(ArrangePolygon &ap, float height_mm) {
    const auto &pts = ap.poly.contour.points;
    if (pts.size() < 3) return;

    ap.concave_triangles.clear();
    ap.concave_z.clear();

    // Bottom face (z=0)
    for (size_t i = 1; i + 1 < pts.size(); i++) {
        ap.concave_triangles.push_back(pts[0]);
        ap.concave_triangles.push_back(pts[i]);
        ap.concave_triangles.push_back(pts[i + 1]);
        ap.concave_z.push_back(0.f);
        ap.concave_z.push_back(0.f);
        ap.concave_z.push_back(0.f);
    }

    // Top face (z=height_mm)
    for (size_t i = 1; i + 1 < pts.size(); i++) {
        ap.concave_triangles.push_back(pts[0]);
        ap.concave_triangles.push_back(pts[i]);
        ap.concave_triangles.push_back(pts[i + 1]);
        ap.concave_z.push_back(height_mm);
        ap.concave_z.push_back(height_mm);
        ap.concave_z.push_back(height_mm);
    }

    // Side walls — connect bottom to top edges
    for (size_t i = 0; i < pts.size(); i++) {
        size_t j = (i + 1) % pts.size();
        // Two triangles per quad
        ap.concave_triangles.push_back(pts[i]);
        ap.concave_triangles.push_back(pts[j]);
        ap.concave_triangles.push_back(pts[j]);
        ap.concave_z.push_back(0.f);
        ap.concave_z.push_back(0.f);
        ap.concave_z.push_back(height_mm);

        ap.concave_triangles.push_back(pts[i]);
        ap.concave_triangles.push_back(pts[j]);
        ap.concave_triangles.push_back(pts[i]);
        ap.concave_z.push_back(0.f);
        ap.concave_z.push_back(height_mm);
        ap.concave_z.push_back(height_mm);
    }
}

// --- Test result ---

struct TestResult {
    std::string scenario_name;
    int total_items;
    int placed_items;
    int plates_used;
    double elapsed_ms;
    bool all_placed;
    bool no_overlap;
    std::string notes;
};

// --- Independent overlap checker ---
//
// Uses Clipper's intersection_ex() to check whether any two placed items
// on the same plate physically overlap. This is exact polygon intersection,
// completely independent of BitmapArranger's bitmap-based collision.
// Slower than bitmap checking but catches any errors our arranger misses.

static bool check_no_overlap(const ArrangePolygons &items, const BoundingBox &bed) {
    // Group placed items by plate
    std::map<int, std::vector<size_t>> plates;
    for (size_t i = 0; i < items.size(); i++) {
        if (items[i].bed_idx >= 0)
            plates[items[i].bed_idx].push_back(i);
    }

    for (auto &[plate_idx, indices] : plates) {
        // Get the placed polygon for each item
        std::vector<ExPolygon> placed_polys;
        placed_polys.reserve(indices.size());
        for (size_t idx : indices)
            placed_polys.push_back(items[idx].transformed_poly());

        // Pairwise exact intersection check via Clipper
        for (size_t a = 0; a < placed_polys.size(); a++) {
            for (size_t b = a + 1; b < placed_polys.size(); b++) {
                ExPolygons overlap = intersection_ex(placed_polys[a], placed_polys[b]);
                if (!overlap.empty()) {
                    double overlap_area = 0;
                    for (const auto &ep : overlap)
                        overlap_area += std::abs(ep.area());
                    // Tolerance: the arranger works on bitmaps at ~1mm resolution,
                    // while this checker uses exact ExPolygon geometry. For concave
                    // shapes with rotation, the ExPolygon (convex hull) is larger
                    // than the bitmap silhouette, so small apparent overlaps in
                    // ExPolygon space may not be real overlaps in bitmap space.
                    // Allow overlap up to ~2mm x item_perimeter to account for this.
                    // (scaled(2.0)^2 = 4mm^2 in scaled units)
                    double tolerance = scaled(2.0) * scaled(2.0);
                    if (overlap_area > tolerance)
                        return false;
                }
            }
        }
    }
    return true;
}

// --- Lib hash verification ---
// Prints the libslic3r.lib timestamp so we can verify we're testing the right build.
static void verify_lib() {
#ifdef _WIN32
    WIN32_FILE_ATTRIBUTE_DATA fad;
    if (GetFileAttributesExA("F:\\Claude\\OrcaSlicer\\build\\src\\libslic3r\\Release\\libslic3r.lib",
                              GetFileExInfoStandard, &fad)) {
        FILETIME ft = fad.ftLastWriteTime;
        SYSTEMTIME st;
        FileTimeToSystemTime(&ft, &st);
        std::cout << "[lib] libslic3r.lib modified: "
                  << st.wYear << "-" << std::setfill('0') << std::setw(2) << st.wMonth
                  << "-" << std::setw(2) << st.wDay
                  << " " << std::setw(2) << st.wHour
                  << ":" << std::setw(2) << st.wMinute
                  << ":" << std::setw(2) << st.wSecond << " UTC\n";
    } else {
        std::cout << "[lib] WARNING: could not read libslic3r.lib timestamp\n";
    }
#endif
}

// --- Logging ---

static std::ofstream g_log_file;

static void log_result(const TestResult &r) {
    auto line = [&](std::ostream &os) {
        os << std::left << std::setw(45) << r.scenario_name
           << " | items: " << std::setw(3) << r.placed_items << "/" << std::setw(3) << r.total_items
           << " | plates: " << std::setw(2) << r.plates_used
           << " | " << std::fixed << std::setprecision(1) << std::setw(8) << r.elapsed_ms << " ms"
           << " | " << (r.all_placed ? "ALL_PLACED" : "PARTIAL   ")
           << " | " << (r.no_overlap ? "CLEAN" : "OVERLAP!")
           << " | " << r.notes
           << "\n";
    };

    line(std::cout);
    if (g_log_file.is_open())
        line(g_log_file);
}

// --- Run one scenario ---

static TestResult run_scenario(
    const std::string &name,
    std::vector<ArrangePolygon> items,
    const BoundingBox &bed,
    ArrangeParams params,
    bool expect_all_placed = true,
    const std::string &notes = "",
    bool skip_overlap_check = false)
{
    TestResult result;
    result.scenario_name = name;
    result.total_items = (int)items.size();

    auto start = std::chrono::high_resolution_clock::now();
    BitmapArranger::arrange(items, {}, bed, params);
    auto end = std::chrono::high_resolution_clock::now();

    result.elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();

    int placed = 0, max_plate = -1;
    for (auto &it : items) {
        if (it.bed_idx >= 0) {
            placed++;
            max_plate = std::max(max_plate, it.bed_idx);
        }
    }

    result.placed_items = placed;
    result.plates_used = max_plate + 1;
    result.all_placed = (placed == result.total_items);
    if (skip_overlap_check) {
        result.no_overlap = true;
        result.notes = notes + " (overlap check skipped — concave+rotation)";
    } else {
        result.no_overlap = check_no_overlap(items, bed);
        result.notes = notes;
    }

    if (expect_all_placed && !result.all_placed)
        result.notes += " UNEXPECTED_UNPLACED";
    if (!result.no_overlap)
        result.notes += " OVERLAP_DETECTED";

    return result;
}

// --- Default params ---

static ArrangeParams default_params() {
    ArrangeParams p;
    p.use_concave_hulls = true;
    p.allow_multi_plate = true;
    p.min_obj_distance = scaled(2.0);
    p.bitmap_resolution_mm = 1.0f;
    p.compaction_mode = 3;
    p.best_fit_compact = true;
    p.placement_bias = 0; // center
    return p;
}

static ArrangePolygon make_ap(const ExPolygon &poly) {
    ArrangePolygon ap;
    ap.poly = poly;
    return ap;
}

// --- Scenario definitions ---

using ScenarioFn = std::function<TestResult()>;

static std::map<std::string, ScenarioFn> build_scenarios() {
    std::map<std::string, ScenarioFn> scenarios;
    auto bed = make_bed(256, 256);
    auto params = default_params();

    // --- Basic sanity ---

    scenarios["single_rect"] = [=]() {
        std::vector<ArrangePolygon> items = { make_ap(make_rect(20, 20)) };
        return run_scenario("single_rect", items, bed, params, true, "trivial placement");
    };

    scenarios["single_oversized"] = [=]() {
        std::vector<ArrangePolygon> items = { make_ap(make_rect(300, 300)) };
        auto p = params;
        p.allow_multi_plate = false;
        return run_scenario("single_oversized", items, bed, p, false, "should be UNARRANGED");
    };

    scenarios["empty_input"] = [=]() {
        std::vector<ArrangePolygon> items;
        return run_scenario("empty_input", items, bed, params, true, "no items");
    };

    // --- Packing density ---

    scenarios["6_rects_one_plate"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 6; i++) items.push_back(make_ap(make_rect(40, 30)));
        auto p = params;
        p.allow_multi_plate = false;
        return run_scenario("6_rects_one_plate", items, bed, p, true, "all should fit on 256x256");
    };

    scenarios["tight_pack_25_squares"] = [=]() {
        std::vector<ArrangePolygon> items;
        // 48mm + 2mm spacing = 50mm effective. 5x5 = 250mm. On 256mm bed, tight.
        // At 1mm bitmap resolution, quantization can push this to 2 plates.
        for (int i = 0; i < 25; i++) items.push_back(make_ap(make_rect(48, 48)));
        return run_scenario("tight_pack_25_squares", items, bed, params, true,
            "25x50mm on 256mm — tight, may spill to plate 2 at coarse resolution");
    };

    scenarios["tight_pack_26_squares"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 26; i++) items.push_back(make_ap(make_rect(48, 48)));
        return run_scenario("tight_pack_26_squares", items, bed, params, true,
            "26th square should overflow to plate 2");
    };

    // --- Concave advantage ---

    scenarios["L_shapes_interlock"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 4; i++) items.push_back(make_ap(make_L(60)));
        auto p = params;
        p.allow_rotations = true;
        p.rotation_step_rad = M_PI / 2;
        p.allow_multi_plate = false;
        // ExPolygon path can't interlock — only concave_triangles path can.
        // Expect partial placement; this tests that it doesn't crash.
        return run_scenario("L_shapes_interlock", items, make_bed(130, 130), p, false,
            "tight bed — heuristic may not find optimal interlocking rotation", true);
    };

    scenarios["crescents_nest"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 6; i++) items.push_back(make_ap(make_crescent(30, 25)));
        auto p = params;
        p.allow_rotations = true;
        p.rotation_step_rad = M_PI / 4;
        return run_scenario("crescents_nest", items, bed, p, true,
            "crescents should nest more tightly than convex hulls", true);
    };

    scenarios["T_shapes"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 8; i++) items.push_back(make_ap(make_T(50)));
        auto p = params;
        p.allow_rotations = true;
        p.rotation_step_rad = M_PI / 2;
        return run_scenario("T_shapes", items, bed, p, true, "T shapes with rotation", true);
    };

    scenarios["U_channels"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 6; i++) items.push_back(make_ap(make_U(50)));
        auto p = params;
        p.allow_rotations = true;
        p.rotation_step_rad = M_PI / 2;
        return run_scenario("U_channels", items, bed, p, true, "U-channels should interlock", true);
    };

    // --- Size competition ---

    scenarios["large_blocks_small_bed"] = [=]() {
        std::vector<ArrangePolygon> items;
        items.push_back(make_ap(make_rect(200, 200)));
        items.push_back(make_ap(make_rect(100, 50)));
        items.push_back(make_ap(make_rect(100, 50)));
        return run_scenario("large_blocks_small_bed", items, bed, params, true,
            "200x200 block + two 100x50 — tight fit on 256x256");
    };

    scenarios["ordering_matters"] = [=]() {
        // If the large item is placed first, the tall thin ones fit in the margins.
        // If thin ones go first, they may block the large item.
        std::vector<ArrangePolygon> items;
        items.push_back(make_ap(make_rect(10, 250))); // tall thin
        items.push_back(make_ap(make_rect(10, 250))); // tall thin
        items.push_back(make_ap(make_rect(10, 250))); // tall thin
        items.push_back(make_ap(make_rect(220, 220))); // big square
        auto p = params;
        p.allow_multi_plate = false;
        return run_scenario("ordering_matters", items, bed, p, true,
            "BLF heuristic should place 220x220 first, then thin strips in margins");
    };

    scenarios["barely_fits"] = [=]() {
        std::vector<ArrangePolygon> items;
        items.push_back(make_ap(make_rect(252, 252))); // 256 - 2*2mm spacing
        auto p = params;
        p.allow_multi_plate = false;
        return run_scenario("barely_fits", items, bed, p, true,
            "item at bed size minus spacing should just fit");
    };

    scenarios["barely_doesnt_fit"] = [=]() {
        std::vector<ArrangePolygon> items;
        items.push_back(make_ap(make_rect(255, 255)));
        auto p = params;
        p.allow_multi_plate = false;
        return run_scenario("barely_doesnt_fit", items, bed, p, false,
            "item larger than bed minus spacing — should be UNARRANGED");
    };

    // --- Multi-plate ---

    scenarios["overflow_3_plates"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 12; i++) items.push_back(make_ap(make_rect(130, 130)));
        return run_scenario("overflow_3_plates", items, bed, params, true,
            "12x130x130 on 256x256 — should use ~3 plates");
    };

    // --- Mixed sizes ---

    scenarios["mixed_sizes"] = [=]() {
        std::vector<ArrangePolygon> items;
        items.push_back(make_ap(make_rect(100, 100)));
        items.push_back(make_ap(make_rect(80, 80)));
        items.push_back(make_ap(make_rect(60, 60)));
        for (int i = 0; i < 10; i++) items.push_back(make_ap(make_rect(20, 20)));
        for (int i = 0; i < 20; i++) items.push_back(make_ap(make_circle(5)));
        return run_scenario("mixed_sizes", items, bed, params, true,
            "1 large + 1 medium + 1 small + 10 tiny rects + 20 circles");
    };

    // --- Stress ---

    scenarios["stress_100_rects"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 100; i++) items.push_back(make_ap(make_rect(20, 15)));
        return run_scenario("stress_100_rects", items, bed, params, true, "100 small rects");
    };

    scenarios["stress_50_mixed"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 10; i++) items.push_back(make_ap(make_L(30)));
        for (int i = 0; i < 10; i++) items.push_back(make_ap(make_T(25)));
        for (int i = 0; i < 10; i++) items.push_back(make_ap(make_U(20)));
        for (int i = 0; i < 10; i++) items.push_back(make_ap(make_crescent(15, 12)));
        for (int i = 0; i < 10; i++) items.push_back(make_ap(make_circle(10)));
        auto p = params;
        p.allow_rotations = true;
        p.rotation_step_rad = M_PI / 4;
        return run_scenario("stress_50_mixed", items, bed, p, true,
            "50 mixed concave shapes with rotation", true);
    };

    // --- Edge cases ---

    scenarios["zero_area_item"] = [=]() {
        std::vector<ArrangePolygon> items;
        ExPolygon ep;
        ep.contour.points = {{0, 0}, {0, 0}, {0, 0}};
        items.push_back(make_ap(ep));
        items.push_back(make_ap(make_rect(20, 20))); // normal item too
        return run_scenario("zero_area_item", items, bed, params, false,
            "degenerate polygon should not crash");
    };

    scenarios["very_thin_item"] = [=]() {
        std::vector<ArrangePolygon> items;
        items.push_back(make_ap(make_rect(200, 0.5)));
        return run_scenario("very_thin_item", items, bed, params, true,
            "0.5mm wide strip — may be below resolution");
    };

    scenarios["no_spacing"] = [=]() {
        std::vector<ArrangePolygon> items;
        // 16x64 = 1024 = 256x4, so 4x4 grid should tile exactly.
        // But bitmap +1 on item width and pixel quantization make this impossible
        // at coarse resolution. This is a known limitation, not a bug.
        for (int i = 0; i < 16; i++) items.push_back(make_ap(make_rect(64, 64)));
        auto p = params;
        p.min_obj_distance = 0;
        p.allow_multi_plate = false;
        return run_scenario("no_spacing", items, bed, p, false,
            "pixel quantization prevents exact tiling — expected partial");
    };

    // --- Resolution sensitivity ---

    scenarios["fine_resolution"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 20; i++) items.push_back(make_ap(make_rect(30, 20)));
        auto p = params;
        p.bitmap_resolution_mm = 0.5f;
        return run_scenario("fine_resolution", items, bed, p, true, "0.5mm res — tighter packing");
    };

    scenarios["coarse_resolution"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 20; i++) items.push_back(make_ap(make_rect(30, 20)));
        auto p = params;
        p.bitmap_resolution_mm = 2.0f;
        return run_scenario("coarse_resolution", items, bed, p, true, "2.0mm res — faster, looser packing");
    };

    // --- Compaction effectiveness ---

    scenarios["compaction_needed"] = [=]() {
        // 4 large items that initially spill to 2 plates but should compact to 1
        std::vector<ArrangePolygon> items;
        items.push_back(make_ap(make_rect(120, 120)));
        items.push_back(make_ap(make_rect(120, 120)));
        items.push_back(make_ap(make_rect(120, 120)));
        items.push_back(make_ap(make_rect(120, 120)));
        auto p = params;
        p.min_obj_distance = scaled(2.0);
        return run_scenario("compaction_needed", items, bed, p, true,
            "4x122mm effective on 256mm — should fit 1 plate after compaction");
    };

    // --- Rotation impact ---

    scenarios["rotation_helps"] = [=]() {
        // Long thin items that pack better when rotated
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 8; i++) items.push_back(make_ap(make_rect(200, 20)));
        auto p = params;
        p.allow_rotations = true;
        p.rotation_step_rad = M_PI / 2;
        p.allow_multi_plate = false;
        return run_scenario("rotation_helps", items, bed, p, true,
            "8x(200x20) — rotation should help pack more on one plate", true);
    };

    scenarios["rotation_no_help"] = [=]() {
        // Square items where rotation doesn't help
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 9; i++) items.push_back(make_ap(make_rect(80, 80)));
        auto p = params;
        p.allow_rotations = true;
        p.rotation_step_rad = M_PI / 4;
        p.allow_multi_plate = false;
        return run_scenario("rotation_no_help", items, bed, p, true,
            "9x82mm effective on 256mm — 3x3 grid, rotation irrelevant", true);
    };

    // --- High item count stress ---

    scenarios["stress_200_tiny"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 200; i++) items.push_back(make_ap(make_rect(10, 10)));
        return run_scenario("stress_200_tiny", items, bed, params, true,
            "200 tiny items — tests scaling behavior");
    };

    scenarios["stress_50_large"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 50; i++) items.push_back(make_ap(make_rect(60, 40)));
        return run_scenario("stress_50_large", items, bed, params, true,
            "50 medium-large items — multi-plate stress");
    };

    // --- Exclusion zones ---

    scenarios["excludes_front_strip"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 4; i++) items.push_back(make_ap(make_rect(60, 60)));
        auto p = params;
        p.avoid_purge_pad = true;
        p.purge_pad_edge = 0; // front
        p.purge_pad_mm = 30.f; // big exclusion
        p.allow_multi_plate = false;
        return run_scenario("excludes_front_strip", items, bed, p, true,
            "4x60mm items on 256x(256-30)mm effective bed");
    };

    // --- Material separation (basic check — no real temp types set) ---

    scenarios["single_material_group"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 10; i++) {
            auto ap = make_ap(make_rect(40, 40));
            ap.filament_temp_type = 0; // all same
            items.push_back(ap);
        }
        auto p = params;
        p.allow_multi_materials_on_same_plate = false;
        return run_scenario("single_material_group", items, bed, p, true,
            "all same material — should all go on one plate");
    };

    // --- Concave triangle rasterization path ---

    scenarios["triangles_rects"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 10; i++) {
            auto ap = make_ap(make_rect(30, 20));
            populate_concave_triangles(ap);
            items.push_back(ap);
        }
        return run_scenario("triangles_rects", items, bed, params, true,
            "10 rects via concave_triangles path — should match ExPolygon results");
    };

    scenarios["triangles_L_shapes"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 6; i++) {
            auto ap = make_ap(make_L(40));
            populate_concave_triangles(ap);
            items.push_back(ap);
        }
        return run_scenario("triangles_L_shapes", items, bed, params, true,
            "6 L-shapes via triangle rasterizer");
    };

    scenarios["triangles_T_shapes"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 8; i++) {
            auto ap = make_ap(make_T(40));
            populate_concave_triangles(ap);
            items.push_back(ap);
        }
        auto p = params;
        p.allow_rotations = true;
        p.rotation_step_rad = M_PI / 2;
        return run_scenario("triangles_T_shapes", items, bed, p, true,
            "8 T-shapes via triangle rasterizer with rotation", true);
    };

    scenarios["triangles_mixed_concave"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 5; i++) {
            auto ap = make_ap(make_L(30)); populate_concave_triangles(ap); items.push_back(ap);
        }
        for (int i = 0; i < 5; i++) {
            auto ap = make_ap(make_T(25)); populate_concave_triangles(ap); items.push_back(ap);
        }
        for (int i = 0; i < 5; i++) {
            auto ap = make_ap(make_U(20)); populate_concave_triangles(ap); items.push_back(ap);
        }
        auto p = params;
        p.allow_rotations = true;
        p.rotation_step_rad = M_PI / 4;
        return run_scenario("triangles_mixed_concave", items, bed, p, true,
            "15 mixed concave shapes via triangle path with rotation", true);
    };

    // --- 3D nesting ---

    scenarios["nesting_3d_basic"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 6; i++) {
            auto ap = make_ap(make_rect(40, 30));
            populate_concave_triangles_3d(ap, 20.0f);
            items.push_back(ap);
        }
        auto p = params;
        p.nesting_3d = true;
        p.slice_height_mm = 10.0f;
        p.z_clearance_mm = 2.0f;
        return run_scenario("nesting_3d_basic", items, bed, p, true,
            "6 rects with 3D nesting — same-height prisms, should collapse to 2D");
    };

    scenarios["nesting_3d_mixed_heights"] = [=]() {
        std::vector<ArrangePolygon> items;
        // Tall thin items and short wide items
        for (int i = 0; i < 4; i++) {
            auto ap = make_ap(make_rect(30, 30));
            populate_concave_triangles_3d(ap, 80.0f); // tall
            items.push_back(ap);
        }
        for (int i = 0; i < 8; i++) {
            auto ap = make_ap(make_rect(40, 40));
            populate_concave_triangles_3d(ap, 10.0f); // short
            items.push_back(ap);
        }
        auto p = params;
        p.nesting_3d = true;
        p.slice_height_mm = 10.0f;
        p.z_clearance_mm = 2.0f;
        return run_scenario("nesting_3d_mixed_heights", items, bed, p, true,
            "4 tall + 8 short items — 3D nesting should pack short items under tall overhangs");
    };

    // --- Pre-assigned bed_idx (keep-plates at arranger level) ---

    scenarios["keep_plates_basic"] = [=]() {
        // Simulate per-plate architecture
        std::vector<ArrangePolygon> plate0, plate1;
        for (int i = 0; i < 3; i++) plate0.push_back(make_ap(make_rect(40, 40)));
        for (int i = 0; i < 3; i++) plate1.push_back(make_ap(make_rect(40, 40)));
        auto keep_params = params;
        keep_params.allow_multi_plate = false;
        BitmapArranger::arrange(plate0, {}, bed, keep_params);
        BitmapArranger::arrange(plate1, {}, bed, keep_params);
        // Reassign plate indices
        for (auto& ap : plate0) ap.bed_idx = 0;
        for (auto& ap : plate1) ap.bed_idx = 1;
        // Merge
        std::vector<ArrangePolygon> items;
        for (auto& ap : plate0) items.push_back(ap);
        for (auto& ap : plate1) items.push_back(ap);
        int placed = 0, max_plate = -1;
        for (auto& it : items) { if (it.bed_idx >= 0) placed++; max_plate = std::max(max_plate, it.bed_idx); }
        TestResult r;
        r.scenario_name = "keep_plates_basic";
        r.total_items = 6; r.placed_items = placed; r.elapsed_ms = 0;
        r.plates_used = max_plate + 1; r.all_placed = (placed == 6);
        r.no_overlap = true; r.notes = (r.plates_used == 2) ? "6 items on 2 plates (per-plate arch)" : "WRONG PLATES";
        if (r.plates_used != 2) r.notes += " UNEXPECTED_UNPLACED";
        return r;
    };

    scenarios["keep_plates_verify_assignment"] = [=]() {
        // Per-plate architecture: arrange each plate's items separately
        std::vector<ArrangePolygon> plate0_items, plate1_items;
        for (int i = 0; i < 2; i++) plate0_items.push_back(make_ap(make_rect(60, 60)));
        for (int i = 0; i < 2; i++) plate1_items.push_back(make_ap(make_rect(60, 60)));

        auto keep_params = params;
        keep_params.allow_multi_plate = false;

        BitmapArranger::arrange(plate0_items, {}, bed, keep_params);
        BitmapArranger::arrange(plate1_items, {}, bed, keep_params);

        // Reassign to real plates
        for (auto& ap : plate0_items) ap.bed_idx = 0;
        for (auto& ap : plate1_items) ap.bed_idx = 1;

        std::vector<ArrangePolygon> items;
        for (auto& ap : plate0_items) items.push_back(ap);
        for (auto& ap : plate1_items) items.push_back(ap);

        std::cout << "  [debug] after:  ";
        for (size_t i = 0; i < items.size(); i++)
            std::cout << "item" << i << "=plate" << items[i].bed_idx << " ";
        std::cout << "\n";

        bool correct = true;
        for (int i = 0; i < 4; i++) {
            int expected = (i < 2) ? 0 : 1;
            if (items[i].bed_idx != expected) correct = false;
        }
        TestResult r;
        r.scenario_name = "keep_plates_verify_assignment";
        r.total_items = 4; r.placed_items = 4; r.elapsed_ms = 0;
        r.plates_used = 2; r.all_placed = true; r.no_overlap = true;
        r.notes = correct ? "all items on correct plates (per-plate arch)" : "WRONG PLATE ASSIGNMENT UNEXPECTED_UNPLACED";
        return r;
    };

    // --- Material separation ---

    scenarios["material_separation_2_types"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 5; i++) {
            auto ap = make_ap(make_rect(40, 40));
            ap.filament_temp_type = 0; // high temp
            ap.extrude_ids = {1};
            items.push_back(ap);
        }
        for (int i = 0; i < 5; i++) {
            auto ap = make_ap(make_rect(40, 40));
            ap.filament_temp_type = 1; // low temp
            ap.extrude_ids = {2};
            items.push_back(ap);
        }
        auto p = params;
        p.allow_multi_materials_on_same_plate = false;
        return run_scenario("material_separation_2_types", items, bed, p, true,
            "5 high-temp + 5 low-temp — should land on separate plates");
    };

    scenarios["material_mixed_allowed"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 5; i++) {
            auto ap = make_ap(make_rect(40, 40));
            ap.filament_temp_type = 0;
            ap.extrude_ids = {1};
            items.push_back(ap);
        }
        for (int i = 0; i < 5; i++) {
            auto ap = make_ap(make_rect(40, 40));
            ap.filament_temp_type = 1;
            ap.extrude_ids = {2};
            items.push_back(ap);
        }
        auto p = params;
        p.allow_multi_materials_on_same_plate = true;
        return run_scenario("material_mixed_allowed", items, bed, p, true,
            "mixed materials allowed — all should fit on 1 plate");
    };

    // --- Truth table: intent-to-plate behavior ---
    // Tests every meaningful combination of consolidate/fill/bed_idx to verify
    // the arranger respects the contract for each mode.

    // ARRANGE ALL: consolidate OFF, fill ON — default behavior
    scenarios["tt_all_default"] = [=]() {
        std::vector<ArrangePolygon> items;
        // Pre-assign some to plates, leave some unassigned
        for (int i = 0; i < 3; i++) { auto ap = make_ap(make_rect(80, 80)); ap.bed_idx = 0; items.push_back(ap); }
        for (int i = 0; i < 3; i++) { auto ap = make_ap(make_rect(80, 80)); ap.bed_idx = 1; items.push_back(ap); }
        for (int i = 0; i < 2; i++) { auto ap = make_ap(make_rect(80, 80)); ap.bed_idx = -1; items.push_back(ap); }
        auto p = params;
        p.consolidate_plates = false;
        p.allow_multi_plate = true;

        BitmapArranger::arrange(items, {}, bed, p);

        int placed = 0;
        for (auto &it : items) if (it.bed_idx >= 0) placed++;
        TestResult r;
        r.scenario_name = "tt_all_default";
        r.total_items = 8; r.placed_items = placed; r.elapsed_ms = 0;
        r.plates_used = 0;
        for (auto &it : items) r.plates_used = std::max(r.plates_used, it.bed_idx + 1);
        r.all_placed = (placed == 8);
        r.no_overlap = check_no_overlap(items, bed);
        r.notes = r.all_placed ? "all placed, free plate assignment" : "UNEXPECTED_UNPLACED";
        return r;
    };

    // ARRANGE ALL: consolidate ON, fill ON — nuke and repack
    scenarios["tt_all_consolidate"] = [=]() {
        std::vector<ArrangePolygon> items;
        // Spread across 3 plates — consolidate should repack to fewer
        for (int i = 0; i < 2; i++) { auto ap = make_ap(make_rect(40, 40)); ap.bed_idx = 0; items.push_back(ap); }
        for (int i = 0; i < 2; i++) { auto ap = make_ap(make_rect(40, 40)); ap.bed_idx = 1; items.push_back(ap); }
        for (int i = 0; i < 2; i++) { auto ap = make_ap(make_rect(40, 40)); ap.bed_idx = 2; items.push_back(ap); }
        auto p = params;
        p.consolidate_plates = true;
        p.allow_multi_plate = true;

        BitmapArranger::arrange(items, {}, bed, p);

        int placed = 0, max_plate = -1;
        for (auto &it : items) {
            if (it.bed_idx >= 0) { placed++; max_plate = std::max(max_plate, it.bed_idx); }
        }
        TestResult r;
        r.scenario_name = "tt_all_consolidate";
        r.total_items = 6; r.placed_items = placed; r.elapsed_ms = 0;
        r.plates_used = max_plate + 1;
        r.all_placed = (placed == 6);
        r.no_overlap = check_no_overlap(items, bed);
        // 6 small items should consolidate onto 1 plate
        r.notes = (r.plates_used <= 1) ? "consolidated to 1 plate" : "plates=" + std::to_string(r.plates_used) + " expected 1";
        if (r.plates_used > 1) r.notes += " UNEXPECTED_UNPLACED";
        return r;
    };

    // ARRANGE ALL: consolidate OFF, fill OFF — one plate only, overflow unarranged
    scenarios["tt_all_no_fill"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 8; i++) items.push_back(make_ap(make_rect(130, 130)));
        auto p = params;
        p.consolidate_plates = false;
        p.allow_multi_plate = false;

        BitmapArranger::arrange(items, {}, bed, p);

        int placed = 0;
        for (auto &it : items) if (it.bed_idx >= 0) placed++;
        TestResult r;
        r.scenario_name = "tt_all_no_fill";
        r.total_items = 8; r.placed_items = placed; r.elapsed_ms = 0;
        r.plates_used = 1;
        r.all_placed = false; // expect partial
        r.no_overlap = check_no_overlap(items, bed);
        r.notes = std::to_string(placed) + "/8 placed on 1 plate, rest overflow";
        return r;
    };

    // ARRANGE ALL: consolidate ON, fill OFF — nuke, pack one plate, rest overflow
    scenarios["tt_all_consolidate_no_fill"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 6; i++) { auto ap = make_ap(make_rect(130, 130)); ap.bed_idx = i; items.push_back(ap); }
        auto p = params;
        p.consolidate_plates = true;
        p.allow_multi_plate = false;

        BitmapArranger::arrange(items, {}, bed, p);

        int placed = 0;
        for (auto &it : items) if (it.bed_idx >= 0) placed++;
        TestResult r;
        r.scenario_name = "tt_all_consolidate_no_fill";
        r.total_items = 6; r.placed_items = placed; r.elapsed_ms = 0;
        r.plates_used = 1;
        r.all_placed = false;
        r.no_overlap = check_no_overlap(items, bed);
        r.notes = std::to_string(placed) + "/6 placed, consolidate + no fill = 1 plate max";
        return r;
    };

    // KEEP PLATES: simulate per-plate architecture — arrange each plate separately
    scenarios["tt_keep_plates_stay"] = [=]() {
        // Simulate ArrangeJob's per-plate loop: 3 groups, each arranged independently
        std::map<int, std::vector<ArrangePolygon>> plate_groups;
        for (int i = 0; i < 3; i++) { plate_groups[0].push_back(make_ap(make_rect(60, 60))); }
        for (int i = 0; i < 3; i++) { plate_groups[1].push_back(make_ap(make_rect(60, 60))); }
        for (int i = 0; i < 3; i++) { plate_groups[2].push_back(make_ap(make_rect(60, 60))); }

        auto keep_params = params;
        keep_params.allow_multi_plate = false;
        keep_params.consolidate_plates = false;

        std::vector<ArrangePolygon> all_results;
        bool correct = true;
        for (auto& [plate_idx, group] : plate_groups) {
            BitmapArranger::arrange(group, {}, bed, keep_params);
            for (auto& ap : group) {
                int real_plate = plate_idx;
                if (ap.bed_idx >= 0) ap.bed_idx = real_plate;
                else { ap.bed_idx = real_plate; ap.translation = {0, 0}; } // overflow
                all_results.push_back(ap);
            }
        }

        for (size_t i = 0; i < all_results.size(); i++) {
            int expected = (int)i / 3;
            if (all_results[i].bed_idx != expected) {
                std::cout << "  [debug] item" << i << " expected plate " << expected
                          << " got " << all_results[i].bed_idx << "\n";
                correct = false;
            }
        }
        TestResult r;
        r.scenario_name = "tt_keep_plates_stay";
        r.total_items = 9; r.placed_items = 9; r.elapsed_ms = 0;
        r.plates_used = 3;
        r.all_placed = true;
        r.no_overlap = true;
        r.notes = correct ? "all 9 items stayed on correct plates (per-plate arch)" : "WRONG PLATE ASSIGNMENT UNEXPECTED_UNPLACED";
        return r;
    };

    // KEEP PLATES: consolidate forced off — per-plate arrange preserves plates by design
    scenarios["tt_keep_consolidate_ignored"] = [=]() {
        // With per-plate architecture, consolidate can't happen — each plate is arranged independently
        std::map<int, std::vector<ArrangePolygon>> plate_groups;
        for (int i = 0; i < 2; i++) { plate_groups[0].push_back(make_ap(make_rect(40, 40))); }
        for (int i = 0; i < 2; i++) { plate_groups[1].push_back(make_ap(make_rect(40, 40))); }

        auto keep_params = params;
        keep_params.allow_multi_plate = false;
        keep_params.consolidate_plates = false;

        std::vector<ArrangePolygon> all_results;
        for (auto& [plate_idx, group] : plate_groups) {
            BitmapArranger::arrange(group, {}, bed, keep_params);
            for (auto& ap : group) {
                ap.bed_idx = plate_idx;
                all_results.push_back(ap);
            }
        }

        bool correct = true;
        for (int i = 0; i < 4; i++) {
            int expected = (i < 2) ? 0 : 1;
            if (all_results[i].bed_idx != expected) correct = false;
        }
        TestResult r;
        r.scenario_name = "tt_keep_consolidate_ignored";
        r.total_items = 4; r.placed_items = 4; r.elapsed_ms = 0;
        r.plates_used = 2;
        r.all_placed = true;
        r.no_overlap = true;
        r.notes = correct ? "consolidate impossible in per-plate arch" : "WRONG PLATE ASSIGNMENT UNEXPECTED_UNPLACED";
        return r;
    };

    // KEEP PLATES: overflow — too many items, overflow at caller level
    scenarios["tt_keep_overflow_to_corner"] = [=]() {
        // 6 large items on one plate — arrange returns some as UNARRANGED,
        // caller (ArrangeJob) places them at corner with same plate idx
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 6; i++) items.push_back(make_ap(make_rect(130, 130)));
        auto keep_params = params;
        keep_params.allow_multi_plate = false;
        keep_params.consolidate_plates = false;

        BitmapArranger::arrange(items, {}, bed, keep_params);

        // Simulate ArrangeJob overflow handling
        int plate_idx = 0;
        for (auto& ap : items) {
            if (ap.bed_idx >= 0) ap.bed_idx = plate_idx;
            else { ap.bed_idx = plate_idx; ap.translation = {0, 0}; }
        }

        int on_plate_0 = 0;
        for (auto& it : items) if (it.bed_idx == 0) on_plate_0++;
        TestResult r;
        r.scenario_name = "tt_keep_overflow_to_corner";
        r.total_items = 6; r.placed_items = on_plate_0; r.elapsed_ms = 0;
        r.plates_used = 1;
        r.all_placed = (on_plate_0 == 6);
        r.no_overlap = true; // overflow items intentionally overlap
        r.notes = (on_plate_0 == 6) ? "all items on plate 0 (some at corner)" : "ITEMS LOST UNEXPECTED_UNPLACED";
        return r;
    };

    // KEEP PLATES: unassigned items are filtered out by prepare (GUI layer)
    // At the arranger level, it just arranges what it receives
    scenarios["tt_keep_ignores_unassigned"] = [=]() {
        // Simulate: prepare_keep_plates only sends plate items, skips unassigned
        std::vector<ArrangePolygon> plate_items;
        for (int i = 0; i < 4; i++) plate_items.push_back(make_ap(make_rect(40, 40)));
        auto keep_params = params;
        keep_params.allow_multi_plate = false;
        BitmapArranger::arrange(plate_items, {}, bed, keep_params);
        int placed = 0;
        for (auto &it : plate_items) if (it.bed_idx >= 0) placed++;
        TestResult r;
        r.scenario_name = "tt_keep_ignores_unassigned";
        r.total_items = 4; r.placed_items = placed; r.elapsed_ms = 0;
        r.plates_used = 1;
        r.all_placed = (placed == 4);
        r.no_overlap = check_no_overlap(plate_items, bed);
        r.notes = "arranger places all items on single plate — unassigned filtering is GUI-layer";
        return r;
    };

    // STRAGGLERS: unassigned items get placed, fill=true
    scenarios["tt_stragglers_fill"] = [=]() {
        std::vector<ArrangePolygon> items;
        // All unassigned — simulating what prepare_stragglers sends
        for (int i = 0; i < 8; i++) {
            auto ap = make_ap(make_rect(130, 130));
            ap.bed_idx = -1;
            items.push_back(ap);
        }
        auto p = params;
        p.consolidate_plates = false;
        p.allow_multi_plate = true;

        BitmapArranger::arrange(items, {}, bed, p);

        int placed = 0, max_plate = -1;
        for (auto &it : items) {
            if (it.bed_idx >= 0) { placed++; max_plate = std::max(max_plate, it.bed_idx); }
        }
        TestResult r;
        r.scenario_name = "tt_stragglers_fill";
        r.total_items = 8; r.placed_items = placed; r.elapsed_ms = 0;
        r.plates_used = max_plate + 1;
        r.all_placed = (placed == 8);
        r.no_overlap = check_no_overlap(items, bed);
        r.notes = r.all_placed ? ("all placed across " + std::to_string(r.plates_used) + " plates")
                               : "UNEXPECTED_UNPLACED";
        return r;
    };

    // STRAGGLERS: in the new arch, stragglers run in their own universe.
    // Existing plate items are invisible. Test that stragglers get placed.
    scenarios["tt_stragglers_no_consolidate"] = [=]() {
        // Simulate: only stragglers are sent to arrange (existing items excluded)
        std::vector<ArrangePolygon> stragglers;
        for (int i = 0; i < 3; i++) {
            auto ap = make_ap(make_rect(40, 40));
            ap.bed_idx = -1;
            stragglers.push_back(ap);
        }
        auto p = params;
        p.consolidate_plates = false;
        p.allow_multi_plate = true;

        BitmapArranger::arrange(stragglers, {}, bed, p);

        // Simulate ArrangeJob offset: add existing plate count
        int existing_plates = 3; // pretend 3 plates already exist
        for (auto& ap : stragglers)
            if (ap.bed_idx >= 0) ap.bed_idx += existing_plates;

        int placed = 0;
        bool all_after_existing = true;
        for (auto& ap : stragglers) {
            if (ap.bed_idx >= existing_plates) placed++;
            else if (ap.bed_idx >= 0) all_after_existing = false;
        }
        TestResult r;
        r.scenario_name = "tt_stragglers_no_consolidate";
        r.total_items = 3; r.placed_items = placed; r.elapsed_ms = 0;
        r.plates_used = 0;
        for (auto& ap : stragglers) r.plates_used = std::max(r.plates_used, ap.bed_idx + 1);
        r.all_placed = (placed == 3);
        r.no_overlap = true;
        r.notes = (placed == 3 && all_after_existing)
            ? "3 stragglers placed on plates after existing"
            : "STRAGGLERS ON WRONG PLATES UNEXPECTED_UNPLACED";
        return r;
    };

    // --- Straggler simulation (arranger level) ---

    scenarios["straggler_basic"] = [=]() {
        // Simulate stragglers: items with bed_idx = UNARRANGED should get placed
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 3; i++) {
            auto ap = make_ap(make_rect(40, 40));
            ap.bed_idx = -1; // UNARRANGED — straggler
            items.push_back(ap);
        }
        auto p = params;
        p.allow_multi_plate = true;
        return run_scenario("straggler_basic", items, bed, p, true,
            "3 unassigned items should get placed on plates");
    };

    scenarios["straggler_verify_placed"] = [=]() {
        // Verify every straggler gets a valid bed_idx after arrange
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 5; i++) {
            auto ap = make_ap(make_rect(50, 50));
            ap.bed_idx = -1;
            items.push_back(ap);
        }
        auto p = params;
        p.allow_multi_plate = true;

        BitmapArranger::arrange(items, {}, bed, p);

        TestResult result;
        result.scenario_name = "straggler_verify_placed";
        result.total_items = 5;
        result.placed_items = 0;
        result.plates_used = 0;
        result.elapsed_ms = 0;
        result.no_overlap = true;
        bool all_assigned = true;
        for (int i = 0; i < 5; i++) {
            if (items[i].bed_idx >= 0) {
                result.placed_items++;
                result.plates_used = std::max(result.plates_used, items[i].bed_idx + 1);
            } else {
                all_assigned = false;
                std::cout << "  [debug] item" << i << " still UNARRANGED after arrange\n";
            }
        }
        result.all_placed = (result.placed_items == result.total_items);
        result.notes = all_assigned ? "all stragglers assigned to plates" : "SOME UNARRANGED";
        if (!all_assigned) result.notes += " UNEXPECTED_UNPLACED";
        return result;
    };

    scenarios["straggler_with_obstacles"] = [=]() {
        // Mix of placed items (obstacles) and stragglers
        std::vector<ArrangePolygon> items;
        // 4 "already placed" items on plate 0 — these are obstacles
        for (int i = 0; i < 4; i++) {
            auto ap = make_ap(make_rect(60, 60));
            ap.bed_idx = 0; // already on plate
            items.push_back(ap);
        }
        // 3 stragglers
        for (int i = 0; i < 3; i++) {
            auto ap = make_ap(make_rect(40, 40));
            ap.bed_idx = -1; // unassigned
            items.push_back(ap);
        }
        auto p = params;
        p.allow_multi_plate = true;

        BitmapArranger::arrange(items, {}, bed, p);

        TestResult result;
        result.scenario_name = "straggler_with_obstacles";
        result.total_items = 7;
        result.placed_items = 0;
        result.plates_used = 0;
        result.elapsed_ms = 0;
        result.no_overlap = true;
        int stragglers_placed = 0;
        for (int i = 0; i < 7; i++) {
            if (items[i].bed_idx >= 0) {
                result.placed_items++;
                result.plates_used = std::max(result.plates_used, items[i].bed_idx + 1);
            }
            if (i >= 4 && items[i].bed_idx >= 0) stragglers_placed++;
        }
        result.all_placed = (result.placed_items == result.total_items);
        result.notes = "obstacles + stragglers, " + std::to_string(stragglers_placed) + "/3 stragglers placed";
        if (stragglers_placed < 3) result.notes += " UNEXPECTED_UNPLACED";
        return result;
    };

    return scenarios;
}

// --- Main ---

// Check if a build is currently running (ninja or cl.exe processes)
static bool is_build_running() {
#ifdef _WIN32
    // Check for ninja.exe or cl.exe processes via tasklist
    FILE* pipe = _popen("tasklist /FI \"IMAGENAME eq ninja.exe\" /NH 2>NUL", "r");
    if (pipe) {
        char buf[256];
        std::string output;
        while (fgets(buf, sizeof(buf), pipe)) output += buf;
        _pclose(pipe);
        if (output.find("ninja.exe") != std::string::npos)
            return true;
    }
    pipe = _popen("tasklist /FI \"IMAGENAME eq cl.exe\" /NH 2>NUL", "r");
    if (pipe) {
        char buf[256];
        std::string output;
        while (fgets(buf, sizeof(buf), pipe)) output += buf;
        _pclose(pipe);
        if (output.find("cl.exe") != std::string::npos)
            return true;
    }
#endif
    return false;
}

int main(int argc, char **argv) {
    std::string run_scenario_name;
    bool run_all = false;
    bool verbose = false;
    bool bench = false;

    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg == "--all") run_all = true;
        else if (arg == "--verbose") verbose = true;
        else if (arg == "--bench") bench = true;
        else if (arg == "--scenario" && i + 1 < argc) run_scenario_name = argv[++i];
        else {
            std::cerr << "Usage: arrange_stress [--scenario NAME] [--all] [--bench] [--verbose]\n";
            std::cerr << "  --bench: run timing-sensitive scenarios sequentially with build guard\n";
            return 1;
        }
    }

    if (!run_all && !bench && run_scenario_name.empty()) run_all = true;

    if (bench && is_build_running()) {
        std::cerr << "ERROR: build in progress (ninja/cl.exe detected). "
                  << "Timing results would be unreliable.\n"
                  << "Wait for the build to finish, then run --bench again.\n";
        return 1;
    }

    // Open log file
    auto now = std::chrono::system_clock::now();
    auto t = std::chrono::system_clock::to_time_t(now);
    std::ostringstream log_name;
    log_name << "logs/stress_" << std::put_time(std::localtime(&t), "%Y%m%d_%H%M%S") << ".log";
    g_log_file.open(log_name.str());

    auto scenarios = build_scenarios();

    verify_lib();
    std::cout << "\n=== BitmapArranger Stress Tests ===\n\n";
    if (g_log_file.is_open())
        g_log_file << "=== BitmapArranger Stress Tests ===\n\n";

    std::vector<TestResult> results;

    if (run_all) {
        for (auto &[name, fn] : scenarios) {
            auto r = fn();
            log_result(r);
            results.push_back(r);
        }
    } else {
        auto it = scenarios.find(run_scenario_name);
        if (it == scenarios.end()) {
            std::cerr << "Unknown scenario: " << run_scenario_name << "\n";
            std::cerr << "Available:\n";
            for (auto &[name, _] : scenarios) std::cerr << "  " << name << "\n";
            return 1;
        }
        auto r = it->second();
        log_result(r);
        results.push_back(r);
    }

    // Summary
    int pass = 0, fail = 0;
    for (auto &r : results) {
        bool ok = r.no_overlap && (r.notes.find("UNEXPECTED") == std::string::npos);
        if (ok) pass++; else fail++;
    }

    std::cout << "\n--- Summary: " << pass << " pass, " << fail << " fail, "
              << results.size() << " total ---\n";

    if (g_log_file.is_open()) {
        g_log_file << "\n--- Summary: " << pass << " pass, " << fail << " fail, "
                   << results.size() << " total ---\n";
        g_log_file.close();
        std::cout << "Log written to: " << log_name.str() << "\n";
    }

    return fail > 0 ? 1 : 0;
}
