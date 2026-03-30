// Anti-gravity compaction experiment
//
// After normal arrangement, push items outward from bed center to
// maximize contiguous free space in the middle. Then try to repack
// items from the last plate into the freed holes on earlier plates.
//
// This is a post-processing pass that operates on ArrangePolygon
// translations — it doesn't touch the bitmap arranger internals.

#pragma once

#include "libslic3r/Arrange/BitmapArranger.hpp"
#include "libslic3r/Arrange.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ClipperUtils.hpp"

#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace Slic3r;
using namespace Slic3r::arrangement;

struct AntiGravityResult {
    int plates_before;
    int plates_after;
    int items_moved;       // items moved from last plate to earlier plates
    int items_pushed;      // items pushed outward
    double elapsed_ms;
};

// Push items on a plate radially outward from the bed center.
// Uses Clipper intersection to check for collisions.
// Returns number of items that moved.
static int anti_gravity_push(
    ArrangePolygons &items,
    const BoundingBox &bed,
    int plate_idx,
    coord_t step = scaled(1.0))  // 1mm steps
{
    // Collect items on this plate
    std::vector<size_t> plate_items;
    for (size_t i = 0; i < items.size(); i++)
        if (items[i].bed_idx == plate_idx) plate_items.push_back(i);

    if (plate_items.empty()) return 0;

    Point bed_center = bed.center();
    int pushed = 0;

    // Sort by distance from center — push closest-to-center items first
    // so they make room for others to follow
    std::sort(plate_items.begin(), plate_items.end(), [&](size_t a, size_t b) {
        ExPolygon pa = items[a].transformed_poly();
        ExPolygon pb = items[b].transformed_poly();
        Point ca = get_extents(pa).center();
        Point cb = get_extents(pb).center();
        double da = std::sqrt(std::pow((double)(ca.x() - bed_center.x()), 2) +
                              std::pow((double)(ca.y() - bed_center.y()), 2));
        double db = std::sqrt(std::pow((double)(cb.x() - bed_center.x()), 2) +
                              std::pow((double)(cb.y() - bed_center.y()), 2));
        return da < db;
    });

    for (size_t idx : plate_items) {
        ExPolygon placed = items[idx].transformed_poly();
        Point item_center = get_extents(placed).center();

        // Direction: from bed center outward through item center
        double dx = (double)(item_center.x() - bed_center.x());
        double dy = (double)(item_center.y() - bed_center.y());
        double len = std::sqrt(dx * dx + dy * dy);
        if (len < 1.0) continue; // already at center, skip

        // Normalize direction
        double nx = dx / len;
        double ny = dy / len;

        // Try sliding outward in steps until we collide
        Vec2crd best_translation = items[idx].translation;
        bool moved = false;

        for (int s = 1; s < 200; s++) { // max 200mm of push
            Vec2crd candidate = {
                items[idx].translation.x() + (coord_t)(s * step * nx),
                items[idx].translation.y() + (coord_t)(s * step * ny)
            };

            // Check: does the item at this position stay on the bed?
            ExPolygon test_poly = items[idx].poly;
            test_poly.rotate(items[idx].rotation);
            test_poly.translate(candidate.x(), candidate.y());
            BoundingBox test_bb = get_extents(test_poly);
            if (test_bb.min.x() < bed.min.x() || test_bb.max.x() > bed.max.x() ||
                test_bb.min.y() < bed.min.y() || test_bb.max.y() > bed.max.y())
                break; // hit bed edge

            // Check: does it collide with any other item on this plate?
            bool collision = false;
            for (size_t other_idx : plate_items) {
                if (other_idx == idx) continue;
                ExPolygon other = items[other_idx].transformed_poly();
                ExPolygons overlap = intersection_ex(test_poly, other);
                double overlap_area = 0;
                for (const auto &ep : overlap)
                    overlap_area += std::abs(ep.area());
                if (overlap_area > scaled(0.5) * scaled(0.5)) {
                    collision = true;
                    break;
                }
            }

            if (collision) break;

            best_translation = candidate;
            moved = true;
        }

        if (moved) {
            items[idx].translation = best_translation;
            pushed++;
        }
    }

    return pushed;
}

// Try to drop items from the last plate into the center holes on earlier plates.
// Center-out scan: start at bed center, spiral outward looking for a collision-free spot.
// This finds holes that skyline placement misses because skyline fills bottom-up.
static int try_fill_holes(
    ArrangePolygons &items,
    const BoundingBox &bed,
    const ArrangeParams &params)
{
    int max_plate = -1;
    for (auto &it : items)
        if (it.bed_idx > max_plate) max_plate = it.bed_idx;
    if (max_plate <= 0) return 0;

    // Collect last plate items (smallest first — easier to fit in holes)
    std::vector<size_t> last_plate_indices;
    for (size_t i = 0; i < items.size(); i++)
        if (items[i].bed_idx == max_plate) last_plate_indices.push_back(i);
    if (last_plate_indices.empty()) return 0;

    std::sort(last_plate_indices.begin(), last_plate_indices.end(), [&](size_t a, size_t b) {
        return std::abs(items[a].poly.area()) < std::abs(items[b].poly.area());
    });

    Point bed_center = bed.center();
    coord_t step = scaled(2.0); // 2mm scan steps
    int moved = 0;

    for (int pi = 0; pi < max_plate; pi++) {
        // Collect placed polys on this plate for collision checking
        auto get_plate_polys = [&]() {
            std::vector<ExPolygon> plate_polys;
            for (size_t i = 0; i < items.size(); i++)
                if (items[i].bed_idx == pi) plate_polys.push_back(items[i].transformed_poly());
            return plate_polys;
        };

        for (int li = (int)last_plate_indices.size() - 1; li >= 0; li--) {
            size_t idx = last_plate_indices[li];
            ExPolygon candidate_poly = items[idx].poly;
            candidate_poly.rotate(items[idx].rotation);
            BoundingBox item_bb = get_extents(candidate_poly);
            coord_t item_w = item_bb.max.x() - item_bb.min.x();
            coord_t item_h = item_bb.max.y() - item_bb.min.y();

            // Center-out spiral scan
            int max_radius = (int)(std::max(bed.max.x() - bed.min.x(), bed.max.y() - bed.min.y()) / step);
            bool placed = false;

            for (int r = 0; r <= max_radius && !placed; r++) {
                // Scan ring at radius r
                int y_lo = (int)((bed_center.y() - r * step - item_h / 2));
                int y_hi = (int)((bed_center.y() + r * step - item_h / 2));
                int x_lo = (int)((bed_center.x() - r * step - item_w / 2));
                int x_hi = (int)((bed_center.x() + r * step - item_w / 2));

                for (coord_t ty = y_lo; ty <= y_hi && !placed; ty += step) {
                    for (coord_t tx = x_lo; tx <= x_hi && !placed; tx += step) {
                        // Only check ring perimeter (skip interior, already checked)
                        if (r > 0 && tx > x_lo + step && tx < x_hi - step &&
                            ty > y_lo + step && ty < y_hi - step)
                            continue;

                        // Build test polygon at this position
                        ExPolygon test = candidate_poly;
                        coord_t off_x = tx - item_bb.min.x();
                        coord_t off_y = ty - item_bb.min.y();
                        test.translate(off_x, off_y);

                        // Check bed bounds
                        BoundingBox tbb = get_extents(test);
                        if (tbb.min.x() < bed.min.x() || tbb.max.x() > bed.max.x() ||
                            tbb.min.y() < bed.min.y() || tbb.max.y() > bed.max.y())
                            continue;

                        // Check collision with all items on this plate
                        bool collision = false;
                        auto plate_polys = get_plate_polys();
                        for (const auto &pp : plate_polys) {
                            ExPolygons overlap = intersection_ex(test, pp);
                            double area = 0;
                            for (const auto &ep : overlap) area += std::abs(ep.area());
                            if (area > scaled(0.5) * scaled(0.5)) { collision = true; break; }
                        }

                        if (!collision) {
                            // Found a spot! Place it here.
                            items[idx].translation = {off_x, off_y};
                            items[idx].bed_idx = pi;
                            last_plate_indices.erase(last_plate_indices.begin() + li);
                            moved++;
                            placed = true;
                        }
                    }
                }
            }
        }
    }

    return moved;
}

// Full anti-gravity pipeline:
// 1. Normal arrangement (caller does this)
// 2. Push items outward from center on each plate except last
// 3. Try to move last-plate items into the freed holes
static AntiGravityResult anti_gravity_compact(
    ArrangePolygons &items,
    const BoundingBox &bed,
    const ArrangeParams &params)
{
    auto start = std::chrono::high_resolution_clock::now();
    AntiGravityResult result = {};

    // Count plates before
    for (auto &it : items)
        result.plates_before = std::max(result.plates_before, it.bed_idx + 1);

    if (result.plates_before <= 1) {
        // Only one plate — push outward for aesthetics but nothing to compact
        result.items_pushed = anti_gravity_push(items, bed, 0);
        result.plates_after = result.plates_before;
        auto end = std::chrono::high_resolution_clock::now();
        result.elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
        return result;
    }

    // Push items outward on all plates except the last
    int last_plate = result.plates_before - 1;
    for (int pi = 0; pi < last_plate; pi++)
        result.items_pushed += anti_gravity_push(items, bed, pi);

    // Try to move last-plate items into holes
    result.items_moved = try_fill_holes(items, bed, params);

    // Count plates after
    result.plates_after = 0;
    for (auto &it : items)
        result.plates_after = std::max(result.plates_after, it.bed_idx + 1);

    auto end = std::chrono::high_resolution_clock::now();
    result.elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
    return result;
}
