// Gravity V2 — experimental packing optimizer
//
// Improvements over v1:
// 1. Binary search for collision boundary (10-50x fewer collision checks)
// 2. Random order permutations with best-of-N selection
// 3. Jiggle-and-settle to escape local minima
// 4. Multi-angle gravity sweep (not just toward target — try tangent offsets)
// 5. Pair swap optimization post-gravity
// 6. Stochastic tournament: run N random configs, keep best

#pragma once

#include "gravity.hpp"
#include <random>
#include <numeric>

// ---------------------------------------------------------------------------
// Binary search slide: find collision boundary in O(log n) steps instead of O(n)
// ---------------------------------------------------------------------------

static SlideResult slide_binary(
    ArrangePolygons& items,
    int item_idx,
    double nx, double ny,  // normalized direction
    const BoundingBox& bed,
    const GravityParams& gp)
{
    auto& ap = items[item_idx];
    int plate = ap.bed_idx;
    Vec2crd origin = ap.translation;
    double rotation = ap.rotation;

    SlideResult result;
    result.best_translation = origin;
    result.best_rotation = rotation;
    result.distance_moved_mm = 0;
    result.contacts = 0;
    result.moved = false;

    // First, find upper bound: exponential probe
    coord_t step = scaled(gp.step_mm);
    int lo = 0, hi = 0;
    for (int probe = 1; probe <= gp.max_steps; probe *= 2) {
        Vec2crd candidate = {
            origin.x() + (coord_t)(nx * step * probe),
            origin.y() + (coord_t)(ny * step * probe)
        };
        ExPolygon test = ap.poly;
        test.rotate(rotation);
        test.translate(candidate.x(), candidate.y());

        if (!within_bed(test, bed, gp.padding_mm)) { hi = probe; break; }

        ap.translation = candidate;
        bool hit = collides_with_others(test, gp.padding_mm, items, item_idx, plate);
        ap.translation = origin;

        if (hit) { hi = probe; break; }
        lo = probe;
        if (probe >= gp.max_steps) { hi = probe; break; }
    }

    if (lo == 0 && hi == 0) return result; // nowhere to go
    if (hi == 0) hi = lo; // never hit anything — use max

    // Binary search between lo and hi
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        Vec2crd candidate = {
            origin.x() + (coord_t)(nx * step * mid),
            origin.y() + (coord_t)(ny * step * mid)
        };
        ExPolygon test = ap.poly;
        test.rotate(rotation);
        test.translate(candidate.x(), candidate.y());

        bool out = !within_bed(test, bed, gp.padding_mm);
        bool hit = false;
        if (!out) {
            ap.translation = candidate;
            hit = collides_with_others(test, gp.padding_mm, items, item_idx, plate);
            ap.translation = origin;
        }

        if (out || hit) hi = mid;
        else lo = mid;
    }

    if (lo > 0) {
        result.best_translation = {
            origin.x() + (coord_t)(nx * step * lo),
            origin.y() + (coord_t)(ny * step * lo)
        };
        result.distance_moved_mm = lo * gp.step_mm;
        result.moved = true;
    }
    return result;
}

// ---------------------------------------------------------------------------
// Multi-angle slide: try the main direction + tangent offsets
// ---------------------------------------------------------------------------

static SlideResult slide_multi_angle(
    ArrangePolygons& items,
    int item_idx,
    Point target,
    const BoundingBox& bed,
    const GravityParams& gp,
    int n_angles = 7)  // main + 3 each side
{
    auto& ap = items[item_idx];
    BoundingBox cur_bb = get_extents(item_footprint(ap));
    Point cur_center = cur_bb.center();

    double dx = (double)(target.x() - cur_center.x());
    double dy = (double)(target.y() - cur_center.y());
    double dist = std::sqrt(dx * dx + dy * dy);
    if (dist < scaled(0.5)) return {ap.translation, ap.rotation, 0, 0, false};

    double base_angle = std::atan2(dy, dx);
    double spread = M_PI / 6.0; // ±30° fan

    SlideResult best = {ap.translation, ap.rotation, 0, 0, false};

    for (int ai = 0; ai < n_angles; ai++) {
        double angle = base_angle + spread * (2.0 * ai / (n_angles - 1) - 1.0);
        double nx = std::cos(angle);
        double ny = std::sin(angle);

        auto sr = slide_binary(items, item_idx, nx, ny, bed, gp);
        if (sr.distance_moved_mm > best.distance_moved_mm) {
            best = sr;
        }
    }
    return best;
}

// ---------------------------------------------------------------------------
// Compute cluster area for a plate (used as the optimization metric)
// ---------------------------------------------------------------------------

static double cluster_area_mm2(const ArrangePolygons& items, int plate_idx) {
    BoundingBox bb = plate_cluster_bb(items, plate_idx);
    if (!bb.defined) return 0;
    return unscaled(bb.size().x()) * unscaled(bb.size().y());
}

// Total cluster area across all plates
static double total_cluster_area(const ArrangePolygons& items) {
    int max_plate = -1;
    for (const auto& ap : items)
        if (ap.bed_idx > max_plate) max_plate = ap.bed_idx;
    double total = 0;
    for (int pi = 0; pi <= max_plate; pi++)
        total += cluster_area_mm2(items, pi);
    return total;
}

// ---------------------------------------------------------------------------
// Save/restore item positions for rollback
// ---------------------------------------------------------------------------

struct ItemSnapshot {
    Vec2crd translation;
    double rotation;
    int bed_idx;
};

static std::vector<ItemSnapshot> save_positions(const ArrangePolygons& items) {
    std::vector<ItemSnapshot> snap(items.size());
    for (size_t i = 0; i < items.size(); i++) {
        snap[i] = {items[i].translation, items[i].rotation, items[i].bed_idx};
    }
    return snap;
}

static void restore_positions(ArrangePolygons& items, const std::vector<ItemSnapshot>& snap) {
    for (size_t i = 0; i < items.size() && i < snap.size(); i++) {
        items[i].translation = snap[i].translation;
        items[i].rotation = snap[i].rotation;
        items[i].bed_idx = snap[i].bed_idx;
    }
}

// ---------------------------------------------------------------------------
// Gravity V2: one pass with given item order
// ---------------------------------------------------------------------------

static double gravity_pass_v2(
    ArrangePolygons& items,
    const std::vector<int>& order,  // item indices in processing order
    Point target,
    const BoundingBox& bed,
    const GravityParams& gp,
    int plate_idx)
{
    int moved = 0;
    for (int idx : order) {
        if (items[idx].bed_idx != plate_idx) continue;

        auto sr = slide_multi_angle(items, idx, target, bed, gp);
        if (sr.moved) {
            items[idx].translation = sr.best_translation;
            if (sr.best_rotation != items[idx].rotation)
                items[idx].rotation = sr.best_rotation;
            moved++;
        }
    }
    return cluster_area_mm2(items, plate_idx);
}

// ---------------------------------------------------------------------------
// Jiggle: randomly perturb items by a small amount
// ---------------------------------------------------------------------------

static void jiggle_items(
    ArrangePolygons& items,
    int plate_idx,
    const BoundingBox& bed,
    const GravityParams& gp,
    std::mt19937& rng,
    double jiggle_mm = 3.0)
{
    std::uniform_real_distribution<double> dist(-jiggle_mm, jiggle_mm);
    coord_t jig = scaled(jiggle_mm);

    for (auto& ap : items) {
        if (ap.bed_idx != plate_idx) continue;

        Vec2crd orig = ap.translation;
        Vec2crd candidate = {
            orig.x() + scaled(dist(rng)),
            orig.y() + scaled(dist(rng))
        };

        ExPolygon test = ap.poly;
        test.rotate(ap.rotation);
        test.translate(candidate.x(), candidate.y());

        if (!within_bed(test, bed, gp.padding_mm)) continue;

        // Quick collision check — accept jiggle only if no collision
        ap.translation = candidate;
        bool hit = false;
        for (size_t j = 0; j < items.size(); j++) {
            if (&items[j] == &ap) continue;
            if (items[j].bed_idx != plate_idx) continue;
            ExPolygon other = item_footprint(items[j]);
            ExPolygons padded_test = offset_ex(test, scaled(gp.padding_mm / 2.0));
            ExPolygons padded_other = offset_ex(other, scaled(gp.padding_mm / 2.0));
            for (const auto& pt : padded_test) {
                for (const auto& po : padded_other) {
                    ExPolygons overlap = intersection_ex(pt, po);
                    double area = 0;
                    for (const auto& o : overlap) area += std::abs(o.area());
                    if (area > scaled(0.1) * scaled(0.1)) { hit = true; break; }
                }
                if (hit) break;
            }
            if (hit) break;
        }

        if (hit) ap.translation = orig; // rollback
    }
}

// ---------------------------------------------------------------------------
// Pair swap: try swapping positions of item pairs to reduce cluster area
// ---------------------------------------------------------------------------

static int try_pair_swaps(
    ArrangePolygons& items,
    int plate_idx,
    const BoundingBox& bed,
    const GravityParams& gp,
    std::mt19937& rng,
    int max_attempts = 50)
{
    std::vector<int> plate_items;
    for (size_t i = 0; i < items.size(); i++)
        if (items[i].bed_idx == plate_idx) plate_items.push_back((int)i);

    if (plate_items.size() < 2) return 0;

    double best_area = cluster_area_mm2(items, plate_idx);
    int swaps_accepted = 0;

    for (int attempt = 0; attempt < max_attempts; attempt++) {
        // Pick two random items
        std::uniform_int_distribution<int> pick(0, (int)plate_items.size() - 1);
        int ai = pick(rng), bi = pick(rng);
        if (ai == bi) continue;
        int a = plate_items[ai], b = plate_items[bi];

        // Swap translations
        std::swap(items[a].translation, items[b].translation);
        std::swap(items[a].rotation, items[b].rotation);

        // Check validity
        auto vr = validate_arrangement(items, bed, gp.padding_mm, plate_idx);
        double new_area = cluster_area_mm2(items, plate_idx);

        if (vr.valid && new_area < best_area) {
            best_area = new_area;
            swaps_accepted++;
        } else {
            // Rollback
            std::swap(items[a].translation, items[b].translation);
            std::swap(items[a].rotation, items[b].rotation);
        }
    }
    return swaps_accepted;
}

// ---------------------------------------------------------------------------
// MAIN: Gravity V2 with stochastic optimization
// ---------------------------------------------------------------------------

struct GravityV2Params : GravityParams {
    int    random_restarts   = 5;    // number of random order attempts
    int    jiggle_rounds     = 3;    // jiggle-and-resettle rounds per restart
    double jiggle_mm         = 3.0;  // perturbation magnitude
    int    swap_attempts     = 30;   // pair swap attempts per round
    bool   multi_angle       = true; // try fan of directions
    unsigned seed            = 42;   // RNG seed (0 = random)
};

struct GravityV2Result : GravityResult {
    int    restarts_tried    = 0;
    int    best_restart      = 0;
    int    swaps_accepted    = 0;
    double area_after_basic  = 0;  // area after plain gravity (no tricks)
    double area_after_v2     = 0;  // area after full v2 pipeline
};

static GravityV2Result gravity_settle_v2(
    ArrangePolygons& items,
    const BoundingBox& bed,
    const GravityV2Params& gp)
{
    auto t0 = std::chrono::steady_clock::now();
    GravityV2Result result;

    std::mt19937 rng(gp.seed ? gp.seed : (unsigned)std::chrono::steady_clock::now().time_since_epoch().count());

    int max_plate = -1;
    for (const auto& ap : items)
        if (ap.bed_idx > max_plate) max_plate = ap.bed_idx;
    if (max_plate < 0) return result;

    // Measure starting state
    result.cluster_area_before = total_cluster_area(items);

    // Run basic gravity first to establish baseline
    auto baseline_snap = save_positions(items);

    // --- Phase 1: Basic gravity (v1 style but with binary search) ---
    for (int pi = 0; pi <= max_plate; pi++) {
        Point target;
        switch (gp.target) {
            case GravityTarget::CENTER:   target = bed.center(); break;
            case GravityTarget::CORNER_BL: target = bed.min; break;
            case GravityTarget::CORNER_BR: target = {bed.max.x(), bed.min.y()}; break;
            case GravityTarget::CORNER_TL: target = {bed.min.x(), bed.max.y()}; break;
            case GravityTarget::CORNER_TR: target = bed.max; break;
            case GravityTarget::CENTROID:  target = plate_centroid(items, pi); break;
        }

        // Sort farthest first
        std::vector<int> order;
        for (size_t i = 0; i < items.size(); i++)
            if (items[i].bed_idx == pi) order.push_back((int)i);

        std::sort(order.begin(), order.end(), [&](int a, int b) {
            BoundingBox ba = get_extents(item_footprint(items[a]));
            BoundingBox bb = get_extents(item_footprint(items[b]));
            Point ca = ba.center(), cb = bb.center();
            double da = std::pow(ca.x() - target.x(), 2) + std::pow(ca.y() - target.y(), 2);
            double db = std::pow(cb.x() - target.x(), 2) + std::pow(cb.y() - target.y(), 2);
            return da > db;
        });

        // Multi-pass until convergence
        for (int pass = 0; pass < gp.max_passes; pass++) {
            int moved = 0;
            for (int idx : order) {
                auto sr = gp.multi_angle ?
                    slide_multi_angle(items, idx, target, bed, gp) :
                    slide_binary(items, idx,
                        (double)(target.x() - get_extents(item_footprint(items[idx])).center().x()) /
                            std::max(1.0, std::sqrt(std::pow(target.x() - get_extents(item_footprint(items[idx])).center().x(), 2) +
                                                     std::pow(target.y() - get_extents(item_footprint(items[idx])).center().y(), 2))),
                        (double)(target.y() - get_extents(item_footprint(items[idx])).center().y()) /
                            std::max(1.0, std::sqrt(std::pow(target.x() - get_extents(item_footprint(items[idx])).center().x(), 2) +
                                                     std::pow(target.y() - get_extents(item_footprint(items[idx])).center().y(), 2))),
                        bed, gp);

                if (sr.moved) {
                    items[idx].translation = sr.best_translation;
                    moved++;
                }
            }
            if (moved == 0) break;
        }
    }

    result.area_after_basic = total_cluster_area(items);
    auto best_snap = save_positions(items);
    double best_area = result.area_after_basic;

    // --- Phase 2: Random restarts with jiggle + resettle ---
    for (int restart = 0; restart < gp.random_restarts; restart++) {
        // Start from baseline (pre-gravity) positions
        restore_positions(items, baseline_snap);
        result.restarts_tried++;

        for (int pi = 0; pi <= max_plate; pi++) {
            Point target;
            switch (gp.target) {
                case GravityTarget::CENTER:   target = bed.center(); break;
                case GravityTarget::CORNER_BL: target = bed.min; break;
                case GravityTarget::CORNER_BR: target = {bed.max.x(), bed.min.y()}; break;
                case GravityTarget::CORNER_TL: target = {bed.min.x(), bed.max.y()}; break;
                case GravityTarget::CORNER_TR: target = bed.max; break;
                case GravityTarget::CENTROID:  target = plate_centroid(items, pi); break;
            }

            // Random item order
            std::vector<int> order;
            for (size_t i = 0; i < items.size(); i++)
                if (items[i].bed_idx == pi) order.push_back((int)i);
            std::shuffle(order.begin(), order.end(), rng);

            // Gravity settle with this order
            for (int pass = 0; pass < gp.max_passes; pass++) {
                double area = gravity_pass_v2(items, order, target, bed, gp, pi);
                // Check if anything moved — just re-sort and try
                std::sort(order.begin(), order.end(), [&](int a, int b) {
                    BoundingBox ba = get_extents(item_footprint(items[a]));
                    BoundingBox bb = get_extents(item_footprint(items[b]));
                    Point ca = ba.center(), cb = bb.center();
                    double da = std::pow(ca.x() - target.x(), 2) + std::pow(ca.y() - target.y(), 2);
                    double db = std::pow(cb.x() - target.x(), 2) + std::pow(cb.y() - target.y(), 2);
                    return da > db;
                });
            }

            // Jiggle rounds
            for (int jr = 0; jr < gp.jiggle_rounds; jr++) {
                jiggle_items(items, pi, bed, gp, rng, gp.jiggle_mm);
                // Re-settle after jiggle
                for (int pass = 0; pass < 3; pass++) {
                    std::vector<int> jig_order = order;
                    std::shuffle(jig_order.begin(), jig_order.end(), rng);
                    gravity_pass_v2(items, jig_order, target, bed, gp, pi);
                }
            }

            // Pair swaps
            int swaps = try_pair_swaps(items, pi, bed, gp, rng, gp.swap_attempts);
            result.swaps_accepted += swaps;

            // Final settle after swaps
            if (swaps > 0) {
                for (int pass = 0; pass < 3; pass++)
                    gravity_pass_v2(items, order, target, bed, gp, pi);
            }
        }

        double area = total_cluster_area(items);
        if (area < best_area) {
            best_area = area;
            best_snap = save_positions(items);
            result.best_restart = restart;
            if (gp.verbose)
                std::cout << "  v2 restart " << restart << ": new best " << (int)area << "mm2\n";
        }
    }

    // Restore best result
    restore_positions(items, best_snap);
    result.area_after_v2 = best_area;
    result.cluster_area_after = best_area;
    result.items_moved = -1; // not tracked per-item in v2

    auto t1 = std::chrono::steady_clock::now();
    result.elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return result;
}
