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

// --- Overlap checker ---

static bool check_no_overlap(const ArrangePolygons &items, const BoundingBox &bed) {
    // Group items by plate
    std::map<int, std::vector<size_t>> plates;
    for (size_t i = 0; i < items.size(); i++) {
        if (items[i].bed_idx >= 0)
            plates[items[i].bed_idx].push_back(i);
    }

    // Per plate: rasterize each item and check for collisions
    for (auto &[plate_idx, indices] : plates) {
        coord_t res = scaled(1.0); // 1mm resolution for checking
        int bed_w = (int)std::ceil((double)(bed.max.x() - bed.min.x()) / res);
        int bed_h = (int)std::ceil((double)(bed.max.y() - bed.min.y()) / res);
        int bed_ww = (bed_w + 63) / 64;
        std::vector<uint64_t> occupied(bed_ww * bed_h, 0);

        for (size_t idx : indices) {
            ExPolygon placed = items[idx].transformed_poly();
            BoundingBox item_bb = get_extents(placed);

            for (const auto &pt : placed.contour.points) {
                int px = (int)((pt.x() - bed.min.x()) / res);
                int py = (int)((pt.y() - bed.min.y()) / res);
                if (px < 0 || px >= bed_w || py < 0 || py >= bed_h) continue;
                size_t word_idx = (size_t)py * bed_ww + px / 64;
                uint64_t bit = uint64_t(1) << (px % 64);
                if (occupied[word_idx] & bit)
                    return false; // overlap detected
                occupied[word_idx] |= bit;
            }
        }
    }
    return true;
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
    const std::string &notes = "")
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
    result.no_overlap = check_no_overlap(items, bed);
    result.notes = notes;

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
        for (int i = 0; i < 25; i++) items.push_back(make_ap(make_rect(48, 48)));
        return run_scenario("tight_pack_25_squares", items, bed, params, true,
            "25x(48+2mm gap)=50*5=250 should just fit");
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
        return run_scenario("L_shapes_interlock", items, make_bed(130, 130), p, true,
            "concave L-shapes should interlock on small bed");
    };

    scenarios["crescents_nest"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 6; i++) items.push_back(make_ap(make_crescent(30, 25)));
        auto p = params;
        p.allow_rotations = true;
        p.rotation_step_rad = M_PI / 4;
        return run_scenario("crescents_nest", items, bed, p, true,
            "crescents should nest more tightly than convex hulls");
    };

    scenarios["T_shapes"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 8; i++) items.push_back(make_ap(make_T(50)));
        auto p = params;
        p.allow_rotations = true;
        p.rotation_step_rad = M_PI / 2;
        return run_scenario("T_shapes", items, bed, p, true, "T shapes with rotation");
    };

    scenarios["U_channels"] = [=]() {
        std::vector<ArrangePolygon> items;
        for (int i = 0; i < 6; i++) items.push_back(make_ap(make_U(50)));
        auto p = params;
        p.allow_rotations = true;
        p.rotation_step_rad = M_PI / 2;
        return run_scenario("U_channels", items, bed, p, true, "U-channels should interlock");
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
            "50 mixed concave shapes with rotation");
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
        for (int i = 0; i < 16; i++) items.push_back(make_ap(make_rect(64, 64)));
        auto p = params;
        p.min_obj_distance = 0;
        p.allow_multi_plate = false;
        return run_scenario("no_spacing", items, bed, p, true,
            "16x64x64 = 256x256 exactly, zero spacing — perfect tile");
    };

    return scenarios;
}

// --- Main ---

int main(int argc, char **argv) {
    std::string run_scenario_name;
    bool run_all = false;
    bool verbose = false;

    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg == "--all") run_all = true;
        else if (arg == "--verbose") verbose = true;
        else if (arg == "--scenario" && i + 1 < argc) run_scenario_name = argv[++i];
        else {
            std::cerr << "Usage: arrange_stress [--scenario NAME] [--all] [--verbose]\n";
            return 1;
        }
    }

    if (!run_all && run_scenario_name.empty()) run_all = true;

    // Open log file
    auto now = std::chrono::system_clock::now();
    auto t = std::chrono::system_clock::to_time_t(now);
    std::ostringstream log_name;
    log_name << "logs/stress_" << std::put_time(std::localtime(&t), "%Y%m%d_%H%M%S") << ".log";
    g_log_file.open(log_name.str());

    auto scenarios = build_scenarios();

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
