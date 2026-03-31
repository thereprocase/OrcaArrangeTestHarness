// Gravity compaction post-processor for BitmapArranger output.
//
// After placement, slides items toward the bed center (or centroid of placed
// items) until they collide, with configurable padding. Optional rotation
// at collision points: wiggle item to find a rotation where it contacts
// a second neighbor, yielding the tightest possible fit.
//
// This is a DENSIFICATION step — it doesn't change which items are on which
// plate, just makes them pack tighter within each plate.

#pragma once

#include "libslic3r/Arrange.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ClipperUtils.hpp"

#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>

using namespace Slic3r;
using namespace Slic3r::arrangement;

enum class GravityTarget {
    CENTER,     // attract toward bed center
    CORNER_BL,  // attract toward bottom-left corner
    CORNER_BR,  // attract toward bottom-right corner
    CORNER_TL,  // attract toward top-left corner
    CORNER_TR,  // attract toward top-right corner
    CENTROID,   // attract toward centroid of items already placed
};

struct GravityParams {
    double padding_mm       = 2.0;   // min gap between items (applied as dilation)
    bool   allow_rotation   = false; // try rotation at collision points
    double rotation_step_deg = 5.0;  // rotation granularity in degrees
    double step_mm          = 0.25;  // movement step size (fine — seconds per plate is OK)
    int    max_steps        = 1200;  // max iterations per item per pass (~300mm travel)
    int    max_passes       = 15;    // repeat until convergence or max_passes
    GravityTarget target    = GravityTarget::CENTER;
    bool   verbose          = false;
};

struct GravityResult {
    int    items_moved       = 0;
    double total_distance_mm = 0;
    double avg_dist_to_center_before = 0;
    double avg_dist_to_center_after  = 0;
    double cluster_area_before = 0;  // bounding box area of all items
    double cluster_area_after  = 0;
    int    passes_used       = 0;
    double elapsed_ms        = 0;
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Get the transformed polygon for an item (rotated + translated).
static ExPolygon item_footprint(const ArrangePolygon& ap) {
    ExPolygon p = ap.poly;
    p.rotate(ap.rotation);
    p.translate(ap.translation.x(), ap.translation.y());
    return p;
}

// Dilate a polygon by half the padding on each side.
// Returns ExPolygons because offset can split concave shapes.
static ExPolygons dilated_footprint(const ArrangePolygon& ap, double padding_mm) {
    ExPolygon p = item_footprint(ap);
    if (padding_mm <= 0) return {p};
    return offset_ex(p, scaled(padding_mm / 2.0));
}

// Check if a candidate position collides with any other item on the same plate,
// respecting padding. Returns true if collision detected.
static bool collides_with_others(
    const ExPolygon& candidate_poly,
    double padding_mm,
    const ArrangePolygons& items,
    int skip_idx,
    int plate_idx)
{
    // Dilate candidate by half padding
    ExPolygons candidate_dilated = offset_ex(candidate_poly, scaled(padding_mm / 2.0));
    if (candidate_dilated.empty()) return false;

    for (size_t i = 0; i < items.size(); i++) {
        if ((int)i == skip_idx) continue;
        if (items[i].bed_idx != plate_idx) continue;

        // Dilate other item by half padding
        ExPolygon other = item_footprint(items[i]);
        ExPolygons other_dilated = offset_ex(other, scaled(padding_mm / 2.0));

        for (const auto& cd : candidate_dilated) {
            for (const auto& od : other_dilated) {
                ExPolygons overlap = intersection_ex(cd, od);
                double area = 0;
                for (const auto& o : overlap) area += std::abs(o.area());
                if (area > scaled(0.1) * scaled(0.1)) return true;
            }
        }
    }
    return false;
}

// Check if a polygon is fully within bed bounds (with padding margin).
static bool within_bed(const ExPolygon& poly, const BoundingBox& bed, double padding_mm) {
    BoundingBox pb = get_extents(poly);
    coord_t margin = scaled(padding_mm / 2.0);
    return pb.min.x() >= bed.min.x() + margin &&
           pb.min.y() >= bed.min.y() + margin &&
           pb.max.x() <= bed.max.x() - margin &&
           pb.max.y() <= bed.max.y() - margin;
}

// Count how many other items (+ bed edges) are "nearly touching" this item.
// Used to score rotation candidates — more contacts = tighter fit.
static int count_contacts(
    const ExPolygon& poly,
    double padding_mm,
    const ArrangePolygons& items,
    int skip_idx,
    int plate_idx,
    const BoundingBox& bed)
{
    int contacts = 0;
    double touch_threshold_mm = padding_mm * 1.5; // within 150% of padding = "contact"
    coord_t threshold = scaled(touch_threshold_mm / 2.0);

    // Check bed edges
    BoundingBox pb = get_extents(poly);
    if (pb.min.x() - bed.min.x() < threshold) contacts++;
    if (pb.min.y() - bed.min.y() < threshold) contacts++;
    if (bed.max.x() - pb.max.x() < threshold) contacts++;
    if (bed.max.y() - pb.max.y() < threshold) contacts++;

    // Check other items
    ExPolygons poly_expanded = offset_ex(poly, threshold);
    for (size_t i = 0; i < items.size(); i++) {
        if ((int)i == skip_idx) continue;
        if (items[i].bed_idx != plate_idx) continue;

        ExPolygon other = item_footprint(items[i]);
        for (const auto& pe : poly_expanded) {
            ExPolygons overlap = intersection_ex(pe, other);
            double area = 0;
            for (const auto& o : overlap) area += std::abs(o.area());
            if (area > scaled(0.1) * scaled(0.1)) {
                contacts++;
                break;
            }
        }
    }
    return contacts;
}

// Compute centroid of all items on a plate.
static Point plate_centroid(const ArrangePolygons& items, int plate_idx) {
    double sx = 0, sy = 0;
    int count = 0;
    for (const auto& ap : items) {
        if (ap.bed_idx != plate_idx) continue;
        BoundingBox bb = get_extents(item_footprint(ap));
        Point c = bb.center();
        sx += c.x();
        sy += c.y();
        count++;
    }
    if (count == 0) return {0, 0};
    return {(coord_t)(sx / count), (coord_t)(sy / count)};
}

// Compute bounding box of all items on a plate.
static BoundingBox plate_cluster_bb(const ArrangePolygons& items, int plate_idx) {
    BoundingBox bb;
    for (const auto& ap : items) {
        if (ap.bed_idx != plate_idx) continue;
        bb.merge(get_extents(item_footprint(ap)));
    }
    return bb;
}

// ---------------------------------------------------------------------------
// Core: Gravity settle one item toward a target point
// ---------------------------------------------------------------------------

struct SlideResult {
    Vec2crd best_translation;
    double  best_rotation;
    double  distance_moved_mm;
    int     contacts;
    bool    moved;
};

static SlideResult slide_item_toward(
    ArrangePolygons& items,
    int item_idx,
    Point target,
    const BoundingBox& bed,
    const GravityParams& gp)
{
    auto& ap = items[item_idx];
    int plate = ap.bed_idx;

    // Current center
    BoundingBox cur_bb = get_extents(item_footprint(ap));
    Point cur_center = cur_bb.center();

    // Direction vector toward target
    double dx = (double)(target.x() - cur_center.x());
    double dy = (double)(target.y() - cur_center.y());
    double dist = std::sqrt(dx * dx + dy * dy);

    SlideResult result;
    result.best_translation = ap.translation;
    result.best_rotation = ap.rotation;
    result.distance_moved_mm = 0;
    result.contacts = 0;
    result.moved = false;

    if (dist < scaled(0.5)) return result; // already at target

    // Normalize direction
    double nx = dx / dist;
    double ny = dy / dist;
    coord_t step = scaled(gp.step_mm);

    // Step toward target, find last valid position
    Vec2crd last_good = ap.translation;
    int steps_moved = 0;

    for (int s = 1; s <= gp.max_steps; s++) {
        Vec2crd candidate = {
            ap.translation.x() + (coord_t)(nx * step * s),
            ap.translation.y() + (coord_t)(ny * step * s)
        };

        // Build test polygon at candidate position
        ExPolygon test = ap.poly;
        test.rotate(ap.rotation);
        test.translate(candidate.x(), candidate.y());

        // Bounds check
        if (!within_bed(test, bed, gp.padding_mm)) break;

        // Collision check
        Vec2crd saved = ap.translation;
        ap.translation = candidate;
        bool hit = collides_with_others(test, gp.padding_mm, items, item_idx, plate);
        ap.translation = saved;

        if (hit) break;

        last_good = candidate;
        steps_moved = s;
    }

    if (steps_moved == 0) return result;

    // If rotation enabled, try to squeeze further at the collision boundary
    if (gp.allow_rotation && steps_moved > 0) {
        double orig_rotation = ap.rotation;
        double step_rad = gp.rotation_step_deg * M_PI / 180.0;
        int n_tries = (int)(360.0 / gp.rotation_step_deg);

        Vec2crd best_pos = last_good;
        double best_rot = ap.rotation;
        int best_contacts = 0;
        double best_dist = steps_moved * gp.step_mm;

        for (int ri = 1; ri < n_tries; ri++) {
            double test_rot = orig_rotation + ri * step_rad;

            // From last_good, try to slide further with this rotation
            Vec2crd rot_last_good = last_good;
            int rot_steps = steps_moved;

            for (int s = steps_moved + 1; s <= gp.max_steps; s++) {
                Vec2crd candidate = {
                    ap.translation.x() + (coord_t)(nx * step * s),
                    ap.translation.y() + (coord_t)(ny * step * s)
                };

                ExPolygon test = ap.poly;
                test.rotate(test_rot);
                test.translate(candidate.x(), candidate.y());

                if (!within_bed(test, bed, gp.padding_mm)) break;

                Vec2crd saved_t = ap.translation;
                double saved_r = ap.rotation;
                ap.translation = candidate;
                ap.rotation = test_rot;
                bool hit = collides_with_others(test, gp.padding_mm, items, item_idx, plate);
                ap.translation = saved_t;
                ap.rotation = saved_r;

                if (hit) break;

                rot_last_good = candidate;
                rot_steps = s;
            }

            if (rot_steps > steps_moved) {
                // This rotation lets us slide further
                ExPolygon at_pos = ap.poly;
                at_pos.rotate(test_rot);
                at_pos.translate(rot_last_good.x(), rot_last_good.y());

                Vec2crd saved_t = ap.translation;
                double saved_r = ap.rotation;
                ap.translation = rot_last_good;
                ap.rotation = test_rot;
                int contacts = count_contacts(at_pos, gp.padding_mm, items,
                                              item_idx, plate, bed);
                ap.translation = saved_t;
                ap.rotation = saved_r;

                double moved_mm = rot_steps * gp.step_mm;
                // Score: prefer more distance moved, then more contacts
                if (moved_mm > best_dist ||
                    (moved_mm == best_dist && contacts > best_contacts)) {
                    best_pos = rot_last_good;
                    best_rot = test_rot;
                    best_contacts = contacts;
                    best_dist = moved_mm;
                }
            }
        }

        last_good = best_pos;
        result.best_rotation = best_rot;
        result.contacts = best_contacts;
        steps_moved = (int)(best_dist / gp.step_mm);
    }

    result.best_translation = last_good;
    result.distance_moved_mm = steps_moved * gp.step_mm;
    result.moved = true;
    return result;
}

// ---------------------------------------------------------------------------
// Main: Gravity settle all items on all plates
// ---------------------------------------------------------------------------

static GravityResult gravity_settle(
    ArrangePolygons& items,
    const BoundingBox& bed,
    const GravityParams& gp)
{
    auto t0 = std::chrono::steady_clock::now();
    GravityResult result;

    // Find all plates
    int max_plate = -1;
    for (const auto& ap : items)
        if (ap.bed_idx > max_plate) max_plate = ap.bed_idx;
    if (max_plate < 0) return result;

    // Measure before
    double total_dist_before = 0;
    int counted = 0;
    for (int pi = 0; pi <= max_plate; pi++) {
        BoundingBox cbb = plate_cluster_bb(items, pi);
        if (!cbb.defined) continue;
        result.cluster_area_before += unscaled(cbb.size().x()) * unscaled(cbb.size().y());

        Point center = bed.center();
        for (const auto& ap : items) {
            if (ap.bed_idx != pi) continue;
            BoundingBox ibb = get_extents(item_footprint(ap));
            Point ic = ibb.center();
            double d = std::sqrt(std::pow(unscaled(ic.x() - center.x()), 2) +
                                 std::pow(unscaled(ic.y() - center.y()), 2));
            total_dist_before += d;
            counted++;
        }
    }
    if (counted > 0) result.avg_dist_to_center_before = total_dist_before / counted;

    // Multi-pass gravity settle
    for (int pass = 0; pass < gp.max_passes; pass++) {
        int moved_this_pass = 0;

        for (int pi = 0; pi <= max_plate; pi++) {
            // Compute attraction target based on GravityTarget
            Point target;
            switch (gp.target) {
                case GravityTarget::CENTER:   target = bed.center(); break;
                case GravityTarget::CORNER_BL: target = bed.min; break;
                case GravityTarget::CORNER_BR: target = {bed.max.x(), bed.min.y()}; break;
                case GravityTarget::CORNER_TL: target = {bed.min.x(), bed.max.y()}; break;
                case GravityTarget::CORNER_TR: target = bed.max; break;
                case GravityTarget::CENTROID:  target = plate_centroid(items, pi); break;
            }

            // Collect items on this plate, sorted farthest-from-target first
            std::vector<int> plate_items;
            for (size_t i = 0; i < items.size(); i++) {
                if (items[i].bed_idx == pi) plate_items.push_back((int)i);
            }

            std::sort(plate_items.begin(), plate_items.end(),
                [&](int a, int b) {
                    BoundingBox ba = get_extents(item_footprint(items[a]));
                    BoundingBox bb_box = get_extents(item_footprint(items[b]));
                    Point ca = ba.center(), cb = bb_box.center();
                    double da = std::pow(ca.x() - target.x(), 2) + std::pow(ca.y() - target.y(), 2);
                    double db = std::pow(cb.x() - target.x(), 2) + std::pow(cb.y() - target.y(), 2);
                    return da > db; // farthest first
                });

            for (int idx : plate_items) {
                SlideResult sr = slide_item_toward(items, idx, target, bed, gp);
                if (sr.moved) {
                    items[idx].translation = sr.best_translation;
                    items[idx].rotation = sr.best_rotation;
                    result.items_moved++;
                    result.total_distance_mm += sr.distance_moved_mm;
                    moved_this_pass++;
                }
            }
        }

        result.passes_used = pass + 1;
        if (gp.verbose) {
            std::cout << "  gravity pass " << (pass + 1) << ": moved "
                      << moved_this_pass << " items\n";
        }
        if (moved_this_pass == 0) break; // converged
    }

    // Measure after
    double total_dist_after = 0;
    counted = 0;
    for (int pi = 0; pi <= max_plate; pi++) {
        BoundingBox cbb = plate_cluster_bb(items, pi);
        if (!cbb.defined) continue;
        result.cluster_area_after += unscaled(cbb.size().x()) * unscaled(cbb.size().y());

        Point center = bed.center();
        for (const auto& ap : items) {
            if (ap.bed_idx != pi) continue;
            BoundingBox ibb = get_extents(item_footprint(ap));
            Point ic = ibb.center();
            double d = std::sqrt(std::pow(unscaled(ic.x() - center.x()), 2) +
                                 std::pow(unscaled(ic.y() - center.y()), 2));
            total_dist_after += d;
            counted++;
        }
    }
    if (counted > 0) result.avg_dist_to_center_after = total_dist_after / counted;

    auto t1 = std::chrono::steady_clock::now();
    result.elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return result;
}

// ---------------------------------------------------------------------------
// Convenience: check if an arrangement is valid (no overlaps, within bed)
// ---------------------------------------------------------------------------

struct ValidationResult {
    bool   valid = true;
    int    overlaps = 0;
    int    out_of_bounds = 0;
    double max_overlap_area_mm2 = 0;
};

// Accurate validation for concave shapes.
//
// The old validator dilated each polygon by half-padding and checked intersection.
// This created false positives on concave shapes where dilation fills concavities
// (an L-shape dilated becomes nearly rectangular, overlapping a neighbor that fits
// in the concavity). False invalidation is WORSE than false acceptance — it
// disqualifies valid placements.
//
// New approach:
// 1. Check actual polygon intersection (no dilation) for HARD overlap
// 2. Check minimum distance between polygon edges for padding violation
//    using offset_ex with NEGATIVE offset: shrink each poly by padding/2,
//    then check if the shrunken versions still exist (if a poly disappears
//    after shrinking, it's too small for the padding — skip distance check)
// 3. For padding check: offset ONE polygon by -(padding - epsilon) and check
//    intersection with the other. If they intersect, they're closer than padding.
//    This is more accurate than dilating both by half because dilation of concave
//    shapes fills concavities, while erosion preserves them.
//
// Threshold: 1mm² for hard overlap, generous to avoid sub-pixel false positives.

static ValidationResult validate_arrangement(
    const ArrangePolygons& items,
    const BoundingBox& bed,
    double padding_mm,
    int plate_idx = -1) // -1 = all plates
{
    ValidationResult vr;

    // Generous threshold: 1mm² overlap area before flagging
    // (bitmap arranger works at 1mm resolution, so sub-mm overlaps are noise)
    double threshold_mm2 = 1.0;

    for (size_t i = 0; i < items.size(); i++) {
        if (plate_idx >= 0 && items[i].bed_idx != plate_idx) continue;
        if (items[i].bed_idx < 0) continue;

        ExPolygon pi_poly = item_footprint(items[i]);

        // Bounds check (with reduced margin — only flag if actually outside bed)
        BoundingBox pb = get_extents(pi_poly);
        coord_t margin = scaled(0.5); // 0.5mm tolerance for bed edge
        if (pb.min.x() < bed.min.x() - margin ||
            pb.min.y() < bed.min.y() - margin ||
            pb.max.x() > bed.max.x() + margin ||
            pb.max.y() > bed.max.y() + margin) {
            vr.out_of_bounds++;
            vr.valid = false;
        }

        // Pairwise overlap check — actual polygon intersection, no dilation
        for (size_t j = i + 1; j < items.size(); j++) {
            if (items[j].bed_idx != items[i].bed_idx) continue;
            if (items[j].bed_idx < 0) continue;

            ExPolygon pj_poly = item_footprint(items[j]);

            // Step 1: Hard overlap check (actual polygons, no padding)
            ExPolygons hard_overlap = intersection_ex(pi_poly, pj_poly);
            double hard_area_mm2 = 0;
            for (const auto& o : hard_overlap)
                hard_area_mm2 += std::abs(unscaled(unscaled(o.area())));

            if (hard_area_mm2 > threshold_mm2) {
                vr.overlaps++;
                vr.valid = false;
                vr.max_overlap_area_mm2 = std::max(vr.max_overlap_area_mm2, hard_area_mm2);
                continue; // already overlapping, skip padding check
            }

            // Step 2: Padding distance check
            // Expand ONE polygon by (padding - 1mm tolerance) and check intersection.
            // Using 80% of padding as the check distance to be generous —
            // sub-pixel placement can't achieve exact padding at bitmap resolution.
            if (padding_mm > 0.5) {
                double check_dist = padding_mm * 0.8; // 80% of requested padding
                ExPolygons pi_expanded = offset_ex(pi_poly, scaled(check_dist));
                for (const auto& pe : pi_expanded) {
                    ExPolygons pad_overlap = intersection_ex(pe, pj_poly);
                    double pad_area_mm2 = 0;
                    for (const auto& o : pad_overlap)
                        pad_area_mm2 += std::abs(unscaled(unscaled(o.area())));

                    if (pad_area_mm2 > threshold_mm2 * 4) { // 4mm² for padding violations
                        vr.overlaps++;
                        vr.valid = false;
                        vr.max_overlap_area_mm2 = std::max(vr.max_overlap_area_mm2, pad_area_mm2);
                    }
                }
            }
        }
    }
    return vr;
}

// ---------------------------------------------------------------------------
// Keep-plates repack decision: should we re-arrange this plate?
// ---------------------------------------------------------------------------

enum class RepackDecision {
    USE_NEW,          // re-arrange produced a valid tighter pack
    KEEP_ORIGINAL,    // originals were valid, re-arrange wasn't better
    FORCE_REARRANGE,  // originals had overlaps, must use re-arrange
};

struct RepackResult {
    RepackDecision decision;
    double original_cluster_area;
    double new_cluster_area;
    bool   originals_valid;
    bool   new_valid;
};

static RepackResult evaluate_repack(
    const ArrangePolygons& original,
    const ArrangePolygons& rearranged,
    const BoundingBox& bed,
    double padding_mm,
    int plate_idx)
{
    RepackResult rr;

    // Validate originals
    auto orig_v = validate_arrangement(original, bed, padding_mm, plate_idx);
    rr.originals_valid = orig_v.valid;

    // Validate rearranged
    auto new_v = validate_arrangement(rearranged, bed, padding_mm, plate_idx);
    rr.new_valid = new_v.valid;

    // Measure cluster areas
    BoundingBox orig_bb = plate_cluster_bb(original, plate_idx);
    BoundingBox new_bb = plate_cluster_bb(rearranged, plate_idx);
    rr.original_cluster_area = orig_bb.defined ?
        unscaled(orig_bb.size().x()) * unscaled(orig_bb.size().y()) : 0;
    rr.new_cluster_area = new_bb.defined ?
        unscaled(new_bb.size().x()) * unscaled(new_bb.size().y()) : 0;

    // Decision logic
    if (!rr.originals_valid) {
        rr.decision = RepackDecision::FORCE_REARRANGE;
    } else if (rr.new_valid && rr.new_cluster_area < rr.original_cluster_area * 0.95) {
        // New arrangement is ≥5% tighter — use it
        rr.decision = RepackDecision::USE_NEW;
    } else {
        // Originals were valid, new isn't meaningfully tighter — keep originals
        rr.decision = RepackDecision::KEEP_ORIGINAL;
    }
    return rr;
}

// ---------------------------------------------------------------------------
// Coordinate round-trip: global ↔ plate-local
// ---------------------------------------------------------------------------

struct PlateCoordTransform {
    double stride_x;
    double stride_y;
    int cols;

    PlateCoordTransform(double bed_w, double bed_h, int plate_count) {
        const double GAP = 1.0 / 5.0; // LOGICAL_PART_PLATE_GAP
        stride_x = bed_w * (1.0 + GAP);
        stride_y = bed_h * (1.0 + GAP);

        // compute_colum_count
        float value = std::sqrt((float)plate_count);
        float round_value = std::round(value);
        cols = (value > round_value) ? (int)round_value + 1 : (int)round_value;
    }

    // Global → plate-local
    Vec2crd to_local(Vec2crd global_trans, int plate_idx) const {
        int row = plate_idx / cols;
        int col = plate_idx % cols;
        return {
            global_trans.x() - scaled(stride_x * col),
            global_trans.y() + scaled(stride_y * row)
        };
    }

    // Plate-local → global
    Vec2crd to_global(Vec2crd local_trans, int plate_idx) const {
        int row = plate_idx / cols;
        int col = plate_idx % cols;
        return {
            local_trans.x() + scaled(stride_x * col),
            local_trans.y() - scaled(stride_y * row)
        };
    }

    // Check round-trip: global → local → global should be identity
    bool round_trip_ok(Vec2crd global_trans, int plate_idx) const {
        Vec2crd local = to_local(global_trans, plate_idx);
        Vec2crd back = to_global(local, plate_idx);
        return back.x() == global_trans.x() && back.y() == global_trans.y();
    }

    // Simulate grid reflow: what happens when plate_count changes
    int cols_for_count(int count) const {
        float value = std::sqrt((float)count);
        float round_value = std::round(value);
        return (value > round_value) ? (int)round_value + 1 : (int)round_value;
    }
};
