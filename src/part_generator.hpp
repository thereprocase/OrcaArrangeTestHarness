// Part Generator — creates 1000 simulated parts across 4 categories
//
// Categories:
//   grid/       — rectangles, squares (mechanical housings, PCBs)
//   mechanical/ — L-shapes, T-shapes, U-channels, brackets, flanges
//   triangles/  — triangular profiles, wedges, trapezoids, pentagons
//   organic/    — blobs, donuts, crescents, rorschach, amoebas, kidney shapes
//
// Each part is an ExPolygon stored as a serialized point list.
// Parts range from tiny (~5mm) to large (~150mm, ~60% of 256mm bed).

#pragma once

#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/libslic3r.h"

#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <string>
#include <filesystem>
#include <iostream>

using namespace Slic3r;

// ---------------------------------------------------------------------------
// Shape primitives with randomized parameters
// ---------------------------------------------------------------------------

static ExPolygon gen_rect(std::mt19937& rng, double min_mm, double max_mm) {
    std::uniform_real_distribution<double> sz(min_mm, max_mm);
    double w = sz(rng), h = sz(rng);
    // Allow aspect ratios up to 4:1
    std::uniform_real_distribution<double> aspect(0.25, 1.0);
    h *= aspect(rng);
    coord_t cw = scaled(w), ch = scaled(h);
    ExPolygon ep;
    ep.contour.points = {{0,0},{cw,0},{cw,ch},{0,ch}};
    return ep;
}

static ExPolygon gen_square(std::mt19937& rng, double min_mm, double max_mm) {
    std::uniform_real_distribution<double> sz(min_mm, max_mm);
    double s = sz(rng);
    coord_t cs = scaled(s);
    ExPolygon ep;
    ep.contour.points = {{0,0},{cs,0},{cs,cs},{0,cs}};
    return ep;
}

static ExPolygon gen_L(std::mt19937& rng, double min_mm, double max_mm) {
    std::uniform_real_distribution<double> sz(min_mm, max_mm);
    std::uniform_real_distribution<double> ratio(0.3, 0.6);
    double s = sz(rng);
    double cut = s * ratio(rng);
    coord_t cs = scaled(s), cc = scaled(cut);
    ExPolygon ep;
    ep.contour.points = {{0,0},{cs,0},{cs,cc},{cc,cc},{cc,cs},{0,cs}};
    return ep;
}

static ExPolygon gen_T(std::mt19937& rng, double min_mm, double max_mm) {
    std::uniform_real_distribution<double> sz(min_mm, max_mm);
    std::uniform_real_distribution<double> ratio(0.25, 0.45);
    double s = sz(rng);
    double stem_w = s * ratio(rng);
    coord_t cs = scaled(s), cw = scaled(stem_w);
    coord_t offset = (cs - cw) / 2;
    ExPolygon ep;
    ep.contour.points = {
        {offset, 0}, {offset + cw, 0}, {offset + cw, cs/2},
        {cs, cs/2}, {cs, cs}, {0, cs}, {0, cs/2}, {offset, cs/2}
    };
    return ep;
}

static ExPolygon gen_U(std::mt19937& rng, double min_mm, double max_mm) {
    std::uniform_real_distribution<double> sz(min_mm, max_mm);
    std::uniform_real_distribution<double> ratio(0.2, 0.35);
    double s = sz(rng);
    double wall = s * ratio(rng);
    coord_t cs = scaled(s), cw = scaled(wall);
    ExPolygon ep;
    ep.contour.points = {
        {0,0},{cs,0},{cs,cs},{cs-cw,cs},{cs-cw,cw},{cw,cw},{cw,cs},{0,cs}
    };
    return ep;
}

static ExPolygon gen_bracket(std::mt19937& rng, double min_mm, double max_mm) {
    // L-shape with a flange (thickened base)
    std::uniform_real_distribution<double> sz(min_mm, max_mm);
    double s = sz(rng);
    coord_t cs = scaled(s);
    coord_t flange = cs / 5;
    coord_t arm = cs * 2 / 3;
    ExPolygon ep;
    ep.contour.points = {
        {0,0},{cs,0},{cs,flange},{flange*2,flange},
        {flange*2,arm},{flange,arm},{flange,flange},{0,flange}
    };
    return ep;
}

static ExPolygon gen_circle(std::mt19937& rng, double min_mm, double max_mm, int segs = 24) {
    std::uniform_real_distribution<double> sz(min_mm/2, max_mm/2);
    double r = sz(rng);
    coord_t cr = scaled(r);
    ExPolygon ep;
    for (int i = 0; i < segs; i++) {
        double a = 2.0 * M_PI * i / segs;
        ep.contour.points.push_back({(coord_t)(cr*cos(a)), (coord_t)(cr*sin(a))});
    }
    return ep;
}

static ExPolygon gen_donut(std::mt19937& rng, double min_mm, double max_mm, int segs = 24) {
    std::uniform_real_distribution<double> sz(min_mm/2, max_mm/2);
    std::uniform_real_distribution<double> hole_ratio(0.3, 0.7);
    double outer_r = sz(rng);
    double inner_r = outer_r * hole_ratio(rng);
    coord_t cor = scaled(outer_r), cir = scaled(inner_r);
    ExPolygon ep;
    // Outer contour (CCW)
    for (int i = 0; i < segs; i++) {
        double a = 2.0 * M_PI * i / segs;
        ep.contour.points.push_back({(coord_t)(cor*cos(a)), (coord_t)(cor*sin(a))});
    }
    // Inner hole (CW)
    Polygon hole;
    for (int i = segs - 1; i >= 0; i--) {
        double a = 2.0 * M_PI * i / segs;
        hole.points.push_back({(coord_t)(cir*cos(a)), (coord_t)(cir*sin(a))});
    }
    ep.holes.push_back(hole);
    return ep;
}

static ExPolygon gen_crescent(std::mt19937& rng, double min_mm, double max_mm, int segs = 24) {
    std::uniform_real_distribution<double> sz(min_mm/2, max_mm/2);
    std::uniform_real_distribution<double> offset_r(0.2, 0.5);
    double outer_r = sz(rng);
    double inner_r = outer_r * 0.8;
    double offset = outer_r * offset_r(rng);
    coord_t cor = scaled(outer_r), cir = scaled(inner_r), coff = scaled(offset);
    ExPolygon ep;
    for (int i = 0; i < segs; i++) {
        double a = 2.0 * M_PI * i / segs;
        ep.contour.points.push_back({(coord_t)(cor*cos(a)), (coord_t)(cor*sin(a))});
    }
    Polygon hole;
    for (int i = segs - 1; i >= 0; i--) {
        double a = 2.0 * M_PI * i / segs;
        hole.points.push_back({coff + (coord_t)(cir*cos(a)), (coord_t)(cir*sin(a))});
    }
    ep.holes.push_back(hole);
    return ep;
}

static ExPolygon gen_blob(std::mt19937& rng, double min_mm, double max_mm, int segs = 20) {
    // Perturbed circle — organic blob shape
    std::uniform_real_distribution<double> sz(min_mm/2, max_mm/2);
    std::uniform_real_distribution<double> perturb(0.7, 1.3);
    double base_r = sz(rng);
    ExPolygon ep;
    for (int i = 0; i < segs; i++) {
        double a = 2.0 * M_PI * i / segs;
        double r = base_r * perturb(rng);
        ep.contour.points.push_back({scaled(r*cos(a)), scaled(r*sin(a))});
    }
    return ep;
}

static ExPolygon gen_rorschach(std::mt19937& rng, double min_mm, double max_mm, int segs = 16) {
    // Symmetric blob — generate half, mirror
    std::uniform_real_distribution<double> sz(min_mm/2, max_mm/2);
    std::uniform_real_distribution<double> perturb(0.5, 1.5);
    double base_r = sz(rng);
    ExPolygon ep;
    // Right half (top to bottom)
    std::vector<Point> right_half;
    for (int i = 0; i <= segs; i++) {
        double a = -M_PI/2 + M_PI * i / segs; // -90° to +90°
        double r = base_r * perturb(rng);
        right_half.push_back({scaled(r*cos(a)), scaled(r*sin(a))});
    }
    // Left half is mirror of right
    for (auto& p : right_half) ep.contour.points.push_back(p);
    for (int i = (int)right_half.size() - 2; i >= 1; i--)
        ep.contour.points.push_back({-right_half[i].x(), right_half[i].y()});
    return ep;
}

static ExPolygon gen_kidney(std::mt19937& rng, double min_mm, double max_mm, int segs = 24) {
    // Kidney/bean shape — circle with concave indent on one side
    std::uniform_real_distribution<double> sz(min_mm/2, max_mm/2);
    std::uniform_real_distribution<double> indent(0.15, 0.35);
    double r = sz(rng);
    double dent = r * indent(rng);
    ExPolygon ep;
    for (int i = 0; i < segs; i++) {
        double a = 2.0 * M_PI * i / segs;
        double mod_r = r;
        // Indent on the right side
        double indent_strength = std::max(0.0, cos(a)) * dent;
        mod_r -= indent_strength;
        ep.contour.points.push_back({scaled(mod_r*cos(a)), scaled(mod_r*sin(a))});
    }
    return ep;
}

static ExPolygon gen_triangle(std::mt19937& rng, double min_mm, double max_mm) {
    std::uniform_real_distribution<double> sz(min_mm, max_mm);
    double base = sz(rng), height = sz(rng);
    std::uniform_real_distribution<double> skew(0.2, 0.8);
    double apex_x = base * skew(rng); // apex offset for non-isoceles
    ExPolygon ep;
    ep.contour.points = {{0,0},{scaled(base),0},{scaled(apex_x),scaled(height)}};
    return ep;
}

static ExPolygon gen_trapezoid(std::mt19937& rng, double min_mm, double max_mm) {
    std::uniform_real_distribution<double> sz(min_mm, max_mm);
    std::uniform_real_distribution<double> ratio(0.4, 0.8);
    double base = sz(rng), height = sz(rng);
    double top = base * ratio(rng);
    double offset = (base - top) / 2;
    ExPolygon ep;
    ep.contour.points = {
        {0,0},{scaled(base),0},
        {scaled(offset + top), scaled(height)},{scaled(offset), scaled(height)}
    };
    return ep;
}

static ExPolygon gen_pentagon(std::mt19937& rng, double min_mm, double max_mm) {
    std::uniform_real_distribution<double> sz(min_mm/2, max_mm/2);
    std::uniform_real_distribution<double> perturb(0.85, 1.15);
    double r = sz(rng);
    ExPolygon ep;
    for (int i = 0; i < 5; i++) {
        double a = 2.0 * M_PI * i / 5 - M_PI/2;
        double pr = r * perturb(rng);
        ep.contour.points.push_back({scaled(pr*cos(a)), scaled(pr*sin(a))});
    }
    return ep;
}

static ExPolygon gen_star(std::mt19937& rng, double min_mm, double max_mm, int points = 5) {
    std::uniform_real_distribution<double> sz(min_mm/2, max_mm/2);
    std::uniform_real_distribution<double> inner_ratio(0.35, 0.55);
    double outer = sz(rng);
    double inner = outer * inner_ratio(rng);
    ExPolygon ep;
    for (int i = 0; i < points * 2; i++) {
        double a = M_PI * i / points - M_PI/2;
        double r = (i % 2 == 0) ? outer : inner;
        ep.contour.points.push_back({scaled(r*cos(a)), scaled(r*sin(a))});
    }
    return ep;
}

static ExPolygon gen_amoeba(std::mt19937& rng, double min_mm, double max_mm) {
    // Multi-lobed organic shape using summed sinusoids
    std::uniform_real_distribution<double> sz(min_mm/2, max_mm/2);
    std::uniform_int_distribution<int> lobes(2, 5);
    std::uniform_real_distribution<double> amp(0.15, 0.4);
    std::uniform_real_distribution<double> phase(0, 2*M_PI);
    double base_r = sz(rng);
    int n_lobes = lobes(rng);
    double lobe_amp = base_r * amp(rng);
    double lobe_phase = phase(rng);
    int segs = 24;
    ExPolygon ep;
    for (int i = 0; i < segs; i++) {
        double a = 2.0 * M_PI * i / segs;
        double r = base_r + lobe_amp * sin(n_lobes * a + lobe_phase);
        r = std::max(r, base_r * 0.3); // prevent negative radii
        ep.contour.points.push_back({scaled(r*cos(a)), scaled(r*sin(a))});
    }
    return ep;
}

// ---------------------------------------------------------------------------
// Part library generation
// ---------------------------------------------------------------------------

struct PartEntry {
    std::string name;
    std::string category;
    ExPolygon poly;
    double area_mm2; // actual polygon area
};

static double poly_area_mm2(const ExPolygon& ep) {
    return std::abs(unscaled(unscaled(ep.area())));
}

// ---------------------------------------------------------------------------
// Serialization: save/load parts to/from simple text files
// Format: one part per file, first line = name, then contour points, then holes
// ---------------------------------------------------------------------------

static void save_part(const PartEntry& pe, const std::string& dir) {
    std::string filename = dir + "/" + pe.category;
    std::filesystem::create_directories(filename);
    // Sanitize name for filename
    std::string safe_name = pe.name;
    for (auto& c : safe_name) if (c == '/') c = '_';
    filename += "/" + safe_name + ".part";

    std::ofstream f(filename);
    if (!f) return;
    f << pe.name << "\n";
    f << pe.category << "\n";
    f << "contour " << pe.poly.contour.points.size() << "\n";
    for (auto& p : pe.poly.contour.points)
        f << p.x() << " " << p.y() << "\n";
    f << "holes " << pe.poly.holes.size() << "\n";
    for (auto& hole : pe.poly.holes) {
        f << "hole " << hole.points.size() << "\n";
        for (auto& p : hole.points)
            f << p.x() << " " << p.y() << "\n";
    }
}

static PartEntry load_part(const std::string& filepath) {
    PartEntry pe;
    std::ifstream f(filepath);
    if (!f) return pe;

    std::getline(f, pe.name);
    std::getline(f, pe.category);

    std::string line;
    // Contour
    std::getline(f, line);
    int n_contour = 0;
    sscanf(line.c_str(), "contour %d", &n_contour);
    for (int i = 0; i < n_contour; i++) {
        coord_t x, y;
        f >> x >> y;
        pe.poly.contour.points.push_back({x, y});
    }
    std::getline(f, line); // consume newline

    // Holes
    std::getline(f, line);
    int n_holes = 0;
    sscanf(line.c_str(), "holes %d", &n_holes);
    for (int h = 0; h < n_holes; h++) {
        std::getline(f, line);
        int n_pts = 0;
        sscanf(line.c_str(), "hole %d", &n_pts);
        Polygon hole;
        for (int i = 0; i < n_pts; i++) {
            coord_t x, y;
            f >> x >> y;
            hole.points.push_back({x, y});
        }
        pe.poly.holes.push_back(hole);
        std::getline(f, line); // consume newline
    }

    pe.poly.contour.make_counter_clockwise();
    for (auto& hole : pe.poly.holes) hole.make_clockwise();
    pe.area_mm2 = poly_area_mm2(pe.poly);
    return pe;
}

static std::vector<PartEntry> load_part_library(const std::string& dir) {
    std::vector<PartEntry> parts;
    for (auto& entry : std::filesystem::recursive_directory_iterator(dir)) {
        if (entry.path().extension() == ".part") {
            auto pe = load_part(entry.path().string());
            if (pe.poly.contour.points.size() >= 3 && pe.area_mm2 > 1.0)
                parts.push_back(std::move(pe));
        }
    }
    std::cout << "  Loaded " << parts.size() << " parts from " << dir << "\n";
    return parts;
}

static std::vector<PartEntry> generate_part_library(unsigned seed = 2026) {
    std::mt19937 rng(seed);
    std::vector<PartEntry> parts;

    // Size ranges: tiny(5-15), small(15-35), medium(35-70), large(70-150)
    struct SizeRange { double min_mm, max_mm; const char* label; };
    SizeRange sizes[] = {
        {5, 15, "tiny"}, {15, 35, "small"}, {35, 70, "med"}, {70, 150, "large"}
    };

    int id = 0;
    auto add = [&](const char* cat, const char* shape, const SizeRange& sz, ExPolygon ep) {
        PartEntry pe;
        pe.name = std::string(cat) + "/" + shape + "_" + sz.label + "_" + std::to_string(id++);
        pe.category = cat;
        pe.poly = std::move(ep);
        pe.area_mm2 = poly_area_mm2(pe.poly);
        parts.push_back(std::move(pe));
    };

    // GRID (250 parts): rectangles, squares
    for (auto& sz : sizes) {
        for (int i = 0; i < 35; i++) add("grid", "rect", sz, gen_rect(rng, sz.min_mm, sz.max_mm));
        for (int i = 0; i < 28; i++) add("grid", "square", sz, gen_square(rng, sz.min_mm, sz.max_mm));
    }

    // MECHANICAL (250 parts): L, T, U, brackets, circles
    for (auto& sz : sizes) {
        for (int i = 0; i < 12; i++) add("mechanical", "L", sz, gen_L(rng, sz.min_mm, sz.max_mm));
        for (int i = 0; i < 12; i++) add("mechanical", "T", sz, gen_T(rng, sz.min_mm, sz.max_mm));
        for (int i = 0; i < 10; i++) add("mechanical", "U", sz, gen_U(rng, sz.min_mm, sz.max_mm));
        for (int i = 0; i < 10; i++) add("mechanical", "bracket", sz, gen_bracket(rng, sz.min_mm, sz.max_mm));
        for (int i = 0; i < 18; i++) add("mechanical", "circle", sz, gen_circle(rng, sz.min_mm, sz.max_mm));
    }

    // TRIANGLES (250 parts): triangles, trapezoids, pentagons, stars
    for (auto& sz : sizes) {
        for (int i = 0; i < 20; i++) add("triangles", "tri", sz, gen_triangle(rng, sz.min_mm, sz.max_mm));
        for (int i = 0; i < 15; i++) add("triangles", "trap", sz, gen_trapezoid(rng, sz.min_mm, sz.max_mm));
        for (int i = 0; i < 15; i++) add("triangles", "pent", sz, gen_pentagon(rng, sz.min_mm, sz.max_mm));
        for (int i = 0; i < 12; i++) add("triangles", "star", sz, gen_star(rng, sz.min_mm, sz.max_mm));
    }

    // ORGANIC (250 parts): blobs, donuts, crescents, rorschach, kidney, amoeba
    for (auto& sz : sizes) {
        for (int i = 0; i < 10; i++) add("organic", "blob", sz, gen_blob(rng, sz.min_mm, sz.max_mm));
        for (int i = 0; i < 12; i++) add("organic", "donut", sz, gen_donut(rng, sz.min_mm, sz.max_mm));
        for (int i = 0; i < 10; i++) add("organic", "crescent", sz, gen_crescent(rng, sz.min_mm, sz.max_mm));
        for (int i = 0; i < 10; i++) add("organic", "rorschach", sz, gen_rorschach(rng, sz.min_mm, sz.max_mm));
        for (int i = 0; i < 10; i++) add("organic", "kidney", sz, gen_kidney(rng, sz.min_mm, sz.max_mm));
        for (int i = 0; i < 10; i++) add("organic", "amoeba", sz, gen_amoeba(rng, sz.min_mm, sz.max_mm));
    }

    std::cout << "  Generated " << parts.size() << " parts across 4 categories\n";

    // Stats
    std::map<std::string, int> cat_counts;
    std::map<std::string, double> cat_min_area, cat_max_area;
    for (auto& p : parts) {
        cat_counts[p.category]++;
        if (cat_min_area.find(p.category) == cat_min_area.end() || p.area_mm2 < cat_min_area[p.category])
            cat_min_area[p.category] = p.area_mm2;
        if (p.area_mm2 > cat_max_area[p.category])
            cat_max_area[p.category] = p.area_mm2;
    }
    for (auto& [cat, count] : cat_counts) {
        std::cout << "    " << cat << ": " << count << " parts, area range "
                  << (int)cat_min_area[cat] << "-" << (int)cat_max_area[cat] << " mm2\n";
    }

    return parts;
}

// Generate and save to disk (run once)
static void generate_and_save_library(const std::string& dir, unsigned seed = 2026) {
    auto parts = generate_part_library(seed);
    std::filesystem::create_directories(dir);
    for (auto& pe : parts) save_part(pe, dir);
    std::cout << "  Saved " << parts.size() << " parts to " << dir << "\n";
}

// Load or generate: checks if library exists on disk, generates if not
static std::vector<PartEntry> get_part_library(const std::string& dir) {
    if (std::filesystem::exists(dir) && !std::filesystem::is_empty(dir)) {
        return load_part_library(dir);
    }
    std::cout << "  Part library not found at " << dir << ", generating...\n";
    generate_and_save_library(dir);
    return load_part_library(dir);
}

// ---------------------------------------------------------------------------
// Random plate builder: pull N random parts from a category (or mixed)
// ---------------------------------------------------------------------------

static std::vector<size_t> pick_random_parts(
    const std::vector<PartEntry>& library,
    int count,
    const std::string& category,  // "" = any, "grid", "mechanical", "triangles", "organic", "mixed"
    std::mt19937& rng,
    double max_total_area_mm2 = 0) // 0 = no limit
{
    // Build index of eligible parts
    std::vector<size_t> eligible;
    if (category == "mixed" || category.empty()) {
        for (size_t i = 0; i < library.size(); i++) eligible.push_back(i);
    } else {
        for (size_t i = 0; i < library.size(); i++)
            if (library[i].category == category) eligible.push_back(i);
    }

    std::shuffle(eligible.begin(), eligible.end(), rng);

    std::vector<size_t> picked;
    double total_area = 0;
    for (size_t idx : eligible) {
        if ((int)picked.size() >= count) break;
        if (max_total_area_mm2 > 0 && total_area + library[idx].area_mm2 > max_total_area_mm2)
            continue;
        picked.push_back(idx);
        total_area += library[idx].area_mm2;
    }
    return picked;
}

// Convert picked parts into ArrangePolygons ready for arrangement
static ArrangePolygons parts_to_arrange_polys(
    const std::vector<PartEntry>& library,
    const std::vector<size_t>& indices)
{
    ArrangePolygons items;
    for (size_t idx : indices) {
        const auto& pe = library[idx];
        // Skip degenerate polygons
        if (pe.poly.contour.points.size() < 3) continue;
        if (pe.area_mm2 < 1.0) continue; // skip sub-1mm² parts

        ArrangePolygon ap;
        ap.poly = pe.poly;
        // Ensure correct orientation (CCW contour, CW holes)
        ap.poly.contour.make_counter_clockwise();
        for (auto& h : ap.poly.holes) h.make_clockwise();
        // Normalize: shift so bounding box starts at (0,0)
        // Many generated shapes (circles, blobs) are centered at origin with
        // negative coordinates. The arranger expects all-positive contours.
        BoundingBox bb = get_extents(ap.poly);
        if (bb.min.x() < 0 || bb.min.y() < 0) {
            coord_t dx = -bb.min.x();
            coord_t dy = -bb.min.y();
            ap.poly.translate(dx, dy);
        }
        ap.translation = {0, 0};
        ap.rotation = 0;
        ap.bed_idx = -1;
        ap.priority = 0;
        items.push_back(std::move(ap));
    }
    return items;
}
