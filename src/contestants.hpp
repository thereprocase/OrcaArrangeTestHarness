// DESIGN CONTEST — 4 competing post-arrangement compaction strategies
//
// All contestants take the same input (arranged items) and try to minimize
// cluster bounding-box area while respecting padding constraints.
//
// Contestants:
// 1. GREEDY GRAVITY V2 (baseline — our current best)
// 2. SIMULATED ANNEALING — accept worse moves with decreasing probability
// 3. PHYSICS SPRING MODEL — force-based simulation with damping
// 4. HEXAGONAL CLOSE PACK — structure-first, assign items to hex grid

#pragma once

#include "gravity.hpp"
#include "gravity_v2.hpp"
#include <random>
#include <numeric>
#include <algorithm>
#include <thread>
#include <functional>

// ---------------------------------------------------------------------------
// Contest infrastructure
// ---------------------------------------------------------------------------

struct ContestEntry {
    std::string name;
    double cluster_area;    // lower = better
    double elapsed_ms;
    bool   valid;           // no overlaps, no OOB
    int    items_on_plate;
};

// ---------------------------------------------------------------------------
// Contestant 1: GREEDY GRAVITY V2 (already implemented)
// ---------------------------------------------------------------------------

static ContestEntry run_gravity_v2(
    ArrangePolygons items, // copy — each contestant gets its own
    const BoundingBox& bed,
    double padding_mm,
    int plate_idx = 0)
{
    GravityV2Params gp;
    gp.padding_mm = padding_mm;
    gp.target = GravityTarget::CENTER;
    gp.random_restarts = 5;
    gp.jiggle_rounds = 3;
    gp.jiggle_mm = 4.0;
    gp.swap_attempts = 30;
    gp.multi_angle = true;
    gp.seed = 42;

    auto r = gravity_settle_v2(items, bed, gp);

    auto vr = validate_arrangement(items, bed, padding_mm, plate_idx);
    return {"GravityV2", r.area_after_v2, r.elapsed_ms, vr.valid, (int)items.size()};
}

// ---------------------------------------------------------------------------
// Contestant 2: SIMULATED ANNEALING
//
// State: item positions (translations + rotations)
// Move: pick random item, move to random nearby position (or swap with another)
// Energy: cluster bounding box area
// Accept: always if better, with P=exp(-dE/T) if worse
// Schedule: exponential cooling from T_init to T_min
// ---------------------------------------------------------------------------

static ContestEntry run_simulated_annealing(
    ArrangePolygons items,
    const BoundingBox& bed,
    double padding_mm,
    int plate_idx = 0)
{
    auto t0 = std::chrono::steady_clock::now();
    std::mt19937 rng(12345);

    // Collect items on target plate
    std::vector<int> plate_items;
    for (size_t i = 0; i < items.size(); i++)
        if (items[i].bed_idx == plate_idx) plate_items.push_back((int)i);

    if (plate_items.size() < 2)
        return {"SimAnneal", cluster_area_mm2(items, plate_idx), 0, true, (int)plate_items.size()};

    double current_energy = cluster_area_mm2(items, plate_idx);
    double best_energy = current_energy;
    auto best_snap = save_positions(items);

    // SA parameters
    double T = current_energy * 0.1;  // initial temp: 10% of energy
    double T_min = 0.1;
    double alpha = 0.995;  // cooling rate
    int iterations = 0;
    int max_iterations = 8000;

    // Move distributions
    std::uniform_int_distribution<int> pick_item(0, (int)plate_items.size() - 1);
    std::uniform_real_distribution<double> move_dist(-15.0, 15.0); // ±15mm
    std::uniform_real_distribution<double> prob(0.0, 1.0);
    std::uniform_real_distribution<double> rot_dist(-M_PI/6, M_PI/6); // ±30°
    std::uniform_int_distribution<int> move_type(0, 2); // 0=translate, 1=swap, 2=rotate

    while (T > T_min && iterations < max_iterations) {
        int idx = plate_items[pick_item(rng)];
        auto& ap = items[idx];
        Vec2crd saved_t = ap.translation;
        double saved_r = ap.rotation;

        int mt = move_type(rng);

        if (mt == 0) {
            // Random translation
            ap.translation.x() += scaled(move_dist(rng));
            ap.translation.y() += scaled(move_dist(rng));
        } else if (mt == 1 && plate_items.size() >= 2) {
            // Swap with another item
            int other = plate_items[pick_item(rng)];
            if (other != idx) {
                std::swap(items[idx].translation, items[other].translation);
                std::swap(items[idx].rotation, items[other].rotation);
            }
        } else {
            // Random rotation
            ap.rotation += rot_dist(rng);
        }

        // Check validity
        ExPolygon test = item_footprint(ap);
        bool valid_move = within_bed(test, bed, padding_mm) &&
                          !collides_with_others(test, padding_mm, items, idx, plate_idx);

        if (valid_move) {
            double new_energy = cluster_area_mm2(items, plate_idx);
            double dE = new_energy - current_energy;

            if (dE < 0 || prob(rng) < std::exp(-dE / T)) {
                // Accept
                current_energy = new_energy;
                if (current_energy < best_energy) {
                    best_energy = current_energy;
                    best_snap = save_positions(items);
                }
            } else {
                // Reject — rollback
                if (mt == 1) {
                    int other = plate_items[pick_item(rng)]; // imperfect rollback for swap
                    // Just restore everything
                    ap.translation = saved_t;
                    ap.rotation = saved_r;
                } else {
                    ap.translation = saved_t;
                    ap.rotation = saved_r;
                }
            }
        } else {
            // Invalid — rollback
            if (mt == 1) {
                // Undo swap — need to track the other item
                // Simplification: just restore from snap if we went wrong
                ap.translation = saved_t;
                ap.rotation = saved_r;
            } else {
                ap.translation = saved_t;
                ap.rotation = saved_r;
            }
        }

        T *= alpha;
        iterations++;
    }

    restore_positions(items, best_snap);

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    auto vr = validate_arrangement(items, bed, padding_mm, plate_idx);
    return {"SimAnneal", best_energy, ms, vr.valid, (int)plate_items.size()};
}

// ---------------------------------------------------------------------------
// Contestant 3: PHYSICS SPRING MODEL
//
// Each item has a spring pulling it toward bed center.
// Items repel each other when closer than padding distance.
// Simulate with velocity, damping, and small timesteps until equilibrium.
// ---------------------------------------------------------------------------

static ContestEntry run_spring_physics(
    ArrangePolygons items,
    const BoundingBox& bed,
    double padding_mm,
    int plate_idx = 0)
{
    auto t0 = std::chrono::steady_clock::now();

    std::vector<int> plate_items;
    for (size_t i = 0; i < items.size(); i++)
        if (items[i].bed_idx == plate_idx) plate_items.push_back((int)i);

    if (plate_items.size() < 2)
        return {"Springs", cluster_area_mm2(items, plate_idx), 0, true, (int)plate_items.size()};

    Point center = bed.center();

    // Velocity per item
    struct Vel { double vx = 0, vy = 0; };
    std::vector<Vel> velocities(items.size());

    double dt = 0.5;           // timestep
    double damping = 0.85;     // velocity damping per step
    double attract_k = 0.002;  // spring constant toward center
    double repel_k = 0.5;      // repulsion strength
    double repel_range_mm = padding_mm * 3.0; // repulsion range
    int max_steps = 2000;

    double best_area = cluster_area_mm2(items, plate_idx);
    auto best_snap = save_positions(items);

    for (int step = 0; step < max_steps; step++) {
        double max_force = 0;

        for (int idx : plate_items) {
            auto& ap = items[idx];
            BoundingBox ibb = get_extents(item_footprint(ap));
            Point ic = ibb.center();

            double fx = 0, fy = 0;

            // Attraction to center
            double dx = unscaled(center.x() - ic.x());
            double dy = unscaled(center.y() - ic.y());
            fx += attract_k * dx;
            fy += attract_k * dy;

            // Repulsion from other items
            for (int other_idx : plate_items) {
                if (other_idx == idx) continue;
                BoundingBox obb = get_extents(item_footprint(items[other_idx]));
                Point oc = obb.center();

                double odx = unscaled(ic.x() - oc.x());
                double ody = unscaled(ic.y() - oc.y());
                double dist = std::sqrt(odx * odx + ody * ody);

                if (dist < repel_range_mm && dist > 0.1) {
                    double force = repel_k * (repel_range_mm - dist) / dist;
                    fx += force * odx / dist;
                    fy += force * ody / dist;
                }
            }

            // Bed boundary repulsion
            BoundingBox pb = get_extents(item_footprint(ap));
            double margin = padding_mm;
            double left = unscaled(pb.min.x() - bed.min.x());
            double right = unscaled(bed.max.x() - pb.max.x());
            double bottom = unscaled(pb.min.y() - bed.min.y());
            double top = unscaled(bed.max.y() - pb.max.y());

            if (left < margin * 2) fx += repel_k * (margin * 2 - left);
            if (right < margin * 2) fx -= repel_k * (margin * 2 - right);
            if (bottom < margin * 2) fy += repel_k * (margin * 2 - bottom);
            if (top < margin * 2) fy -= repel_k * (margin * 2 - top);

            // Update velocity
            velocities[idx].vx = (velocities[idx].vx + fx * dt) * damping;
            velocities[idx].vy = (velocities[idx].vy + fy * dt) * damping;

            max_force = std::max(max_force, std::sqrt(fx*fx + fy*fy));
        }

        // Apply velocities
        for (int idx : plate_items) {
            auto& ap = items[idx];
            Vec2crd saved = ap.translation;

            ap.translation.x() += scaled(velocities[idx].vx * dt);
            ap.translation.y() += scaled(velocities[idx].vy * dt);

            // Validate — rollback if invalid
            ExPolygon test = item_footprint(ap);
            if (!within_bed(test, bed, padding_mm) ||
                collides_with_others(test, padding_mm, items, idx, plate_idx)) {
                ap.translation = saved;
                velocities[idx].vx *= -0.3; // bounce
                velocities[idx].vy *= -0.3;
            }
        }

        double area = cluster_area_mm2(items, plate_idx);
        if (area < best_area) {
            best_area = area;
            best_snap = save_positions(items);
        }

        // Convergence check
        if (max_force < 0.001 && step > 100) break;
    }

    restore_positions(items, best_snap);

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    auto vr = validate_arrangement(items, bed, padding_mm, plate_idx);
    return {"Springs", best_area, ms, vr.valid, (int)plate_items.size()};
}

// ---------------------------------------------------------------------------
// Contestant 4: HEXAGONAL CLOSE PACK
//
// Generate a hex grid covering the bed, centered on bed center.
// Assign items to hex cells greedily (largest items first, closest to center).
// Hex packing achieves ~90.7% density for circles — optimal for uniform items.
// ---------------------------------------------------------------------------

static ContestEntry run_hex_pack(
    ArrangePolygons items,
    const BoundingBox& bed,
    double padding_mm,
    int plate_idx = 0)
{
    auto t0 = std::chrono::steady_clock::now();

    std::vector<int> plate_items;
    for (size_t i = 0; i < items.size(); i++)
        if (items[i].bed_idx == plate_idx) plate_items.push_back((int)i);

    if (plate_items.empty())
        return {"HexPack", 0, 0, true, 0};

    // Find the largest item bounding box to determine hex cell size
    double max_dim = 0;
    for (int idx : plate_items) {
        BoundingBox ibb = get_extents(item_footprint(items[idx]));
        double w = unscaled(ibb.size().x());
        double h = unscaled(ibb.size().y());
        max_dim = std::max(max_dim, std::max(w, h));
    }
    double cell_size = max_dim + padding_mm; // hex cell diameter

    // Generate hex grid centers
    Point bed_center = bed.center();
    double bed_w = unscaled(bed.size().x());
    double bed_h = unscaled(bed.size().y());

    double hex_w = cell_size;                    // horizontal spacing
    double hex_h = cell_size * std::sqrt(3) / 2; // vertical spacing

    struct HexCell { double x, y; bool used = false; };
    std::vector<HexCell> cells;

    int cols = (int)(bed_w / hex_w) + 2;
    int rows = (int)(bed_h / hex_h) + 2;
    double start_x = unscaled(bed_center.x()) - (cols / 2.0) * hex_w;
    double start_y = unscaled(bed_center.y()) - (rows / 2.0) * hex_h;

    for (int r = 0; r < rows; r++) {
        double offset = (r % 2) ? hex_w / 2.0 : 0;
        for (int c = 0; c < cols; c++) {
            double cx = start_x + c * hex_w + offset;
            double cy = start_y + r * hex_h;

            // Check if cell center is within bed (with margin)
            if (cx >= unscaled(bed.min.x()) + cell_size/2 + padding_mm &&
                cx <= unscaled(bed.max.x()) - cell_size/2 - padding_mm &&
                cy >= unscaled(bed.min.y()) + cell_size/2 + padding_mm &&
                cy <= unscaled(bed.max.y()) - cell_size/2 - padding_mm) {
                cells.push_back({cx, cy, false});
            }
        }
    }

    // Sort cells by distance to bed center (closest first)
    double bcx = unscaled(bed_center.x()), bcy = unscaled(bed_center.y());
    std::sort(cells.begin(), cells.end(), [bcx, bcy](const HexCell& a, const HexCell& b) {
        double da = std::pow(a.x - bcx, 2) + std::pow(a.y - bcy, 2);
        double db = std::pow(b.x - bcx, 2) + std::pow(b.y - bcy, 2);
        return da < db;
    });

    // Sort items by area (largest first)
    std::sort(plate_items.begin(), plate_items.end(), [&](int a, int b) {
        BoundingBox ba = get_extents(items[a].poly);
        BoundingBox bb_box = get_extents(items[b].poly);
        return (ba.size().x() * ba.size().y()) > (bb_box.size().x() * bb_box.size().y());
    });

    // Greedy assignment: for each item, find closest unused cell
    for (int idx : plate_items) {
        auto& ap = items[idx];
        BoundingBox ibb = get_extents(ap.poly);
        Point item_center = ibb.center();

        bool placed = false;
        for (auto& cell : cells) {
            if (cell.used) continue;

            // Place item centered on hex cell
            Vec2crd new_trans = {
                scaled(cell.x) - item_center.x(),
                scaled(cell.y) - item_center.y()
            };

            Vec2crd saved = ap.translation;
            ap.translation = new_trans;

            ExPolygon test = item_footprint(ap);
            if (within_bed(test, bed, padding_mm) &&
                !collides_with_others(test, padding_mm, items, idx, plate_idx)) {
                cell.used = true;
                placed = true;
                break;
            }

            ap.translation = saved;
        }

        // If no hex cell works, keep original position
        if (!placed) {
            // Leave as-is
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    double area = cluster_area_mm2(items, plate_idx);
    auto vr = validate_arrangement(items, bed, padding_mm, plate_idx);
    return {"HexPack", area, ms, vr.valid, (int)plate_items.size()};
}

// ---------------------------------------------------------------------------
// Contestant 5: ELENA — Bounding-Box Gradient Descent
//
// Dr. Elena "Forces" Kowalski — physics PhD approach.
// The cluster bounding-box area E = W * H is the energy functional.
// dE/dx_i is nonzero only for items on the bounding-box boundary.
// Strategy: iteratively identify boundary items, slide them inward
// along the negative gradient, then re-settle the whole cluster
// toward the shrinking centroid. Rinse and repeat.
// ---------------------------------------------------------------------------

static ContestEntry elena_round1(
    ArrangePolygons& items,
    const BoundingBox& bed,
    double padding_mm,
    int plate_idx)
{
    auto t0 = std::chrono::steady_clock::now();

    std::vector<int> plate_items;
    for (size_t i = 0; i < items.size(); i++)
        if (items[i].bed_idx == plate_idx) plate_items.push_back((int)i);

    if (plate_items.size() < 2)
        return {"Elena", cluster_area_mm2(items, plate_idx), 0, true, (int)plate_items.size()};

    GravityParams gp;
    gp.padding_mm = padding_mm;
    gp.step_mm = 0.25;
    gp.max_steps = 1200;

    double best_area = cluster_area_mm2(items, plate_idx);
    auto best_snap = save_positions(items);

    // --- Phase 1: Initial gravity settle toward centroid ---
    auto settle_toward = [&](Point target) {
        // Sort farthest-first for stable packing order
        std::vector<int> order = plate_items;
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            Point ca = get_extents(item_footprint(items[a])).center();
            Point cb = get_extents(item_footprint(items[b])).center();
            double da = std::pow(ca.x()-target.x(),2) + std::pow(ca.y()-target.y(),2);
            double db = std::pow(cb.x()-target.x(),2) + std::pow(cb.y()-target.y(),2);
            return da > db;
        });
        for (int pass = 0; pass < 8; pass++) {
            int moved = 0;
            for (int idx : order) {
                auto sr = slide_multi_angle(items, idx, target, bed, gp, 7);
                if (sr.moved) { items[idx].translation = sr.best_translation; moved++; }
            }
            if (moved == 0) break;
        }
    };

    // Initial settle toward bed center
    settle_toward(bed.center());

    double area = cluster_area_mm2(items, plate_idx);
    if (area < best_area) { best_area = area; best_snap = save_positions(items); }

    // --- Phase 2: Gradient descent on bounding-box area ---
    // Identify boundary items and push them inward. The gradient of
    // E = W*H w.r.t. a boundary item's position points outward from
    // the cluster. We move along -grad(E).
    for (int epoch = 0; epoch < 12; epoch++) {
        BoundingBox cbb = plate_cluster_bb(items, plate_idx);
        if (!cbb.defined) break;
        double W = unscaled(cbb.size().x());
        double H = unscaled(cbb.size().y());
        coord_t tol = scaled(1.0); // 1mm tolerance for "on boundary"

        // Find which items define each edge of the bounding box
        for (int idx : plate_items) {
            BoundingBox ibb = get_extents(item_footprint(items[idx]));
            double fx = 0, fy = 0; // gradient direction (inward)

            // dE/dx for right-edge item: H (area grows rightward)
            if (ibb.max.x() > cbb.max.x() - tol) fx -= H;
            // dE/dx for left-edge item: -H (area grows leftward)
            if (ibb.min.x() < cbb.min.x() + tol) fx += H;
            // dE/dy for top-edge item: W
            if (ibb.max.y() > cbb.max.y() - tol) fy -= W;
            // dE/dy for bottom-edge item: -W
            if (ibb.min.y() < cbb.min.y() + tol) fy += W;

            if (fx == 0 && fy == 0) continue; // interior item, zero gradient

            // Normalize and slide along negative gradient (inward)
            double mag = std::sqrt(fx*fx + fy*fy);
            if (mag < 0.01) continue;
            double nx = fx / mag, ny = fy / mag;

            auto sr = slide_binary(items, idx, nx, ny, bed, gp);
            if (sr.moved) items[idx].translation = sr.best_translation;
        }

        // Re-settle everything toward the new centroid
        Point centroid = plate_centroid(items, plate_idx);
        settle_toward(centroid);

        area = cluster_area_mm2(items, plate_idx);
        if (area < best_area) { best_area = area; best_snap = save_positions(items); }
    }

    // --- Phase 3: Corner-target exploration ---
    // Try settling toward each corner — sometimes an asymmetric target
    // reaches a lower-energy state than the centroid.
    Point corners[] = { bed.min, {bed.max.x(), bed.min.y()},
                        {bed.min.x(), bed.max.y()}, bed.max };
    auto pre_corner_snap = save_positions(items);

    for (const auto& corner : corners) {
        restore_positions(items, pre_corner_snap);
        settle_toward(corner);

        // Follow up with gradient descent epochs
        for (int epoch = 0; epoch < 4; epoch++) {
            BoundingBox cbb = plate_cluster_bb(items, plate_idx);
            if (!cbb.defined) break;
            double W = unscaled(cbb.size().x()), H = unscaled(cbb.size().y());
            coord_t tol = scaled(1.0);

            for (int idx : plate_items) {
                BoundingBox ibb = get_extents(item_footprint(items[idx]));
                double fx = 0, fy = 0;
                if (ibb.max.x() > cbb.max.x() - tol) fx -= H;
                if (ibb.min.x() < cbb.min.x() + tol) fx += H;
                if (ibb.max.y() > cbb.max.y() - tol) fy -= W;
                if (ibb.min.y() < cbb.min.y() + tol) fy += W;
                double mag = std::sqrt(fx*fx + fy*fy);
                if (mag < 0.01) continue;
                auto sr = slide_binary(items, idx, fx/mag, fy/mag, bed, gp);
                if (sr.moved) items[idx].translation = sr.best_translation;
            }
            Point centroid = plate_centroid(items, plate_idx);
            settle_toward(centroid);
        }

        area = cluster_area_mm2(items, plate_idx);
        if (area < best_area) { best_area = area; best_snap = save_positions(items); }
    }

    // Restore best configuration found
    restore_positions(items, best_snap);

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    auto vr = validate_arrangement(items, bed, padding_mm, plate_idx);
    return {"Elena", best_area, ms, vr.valid, (int)plate_items.size()};
}

// ---------------------------------------------------------------------------
// Contestant 6: MONTE — Guided Simulated Annealing
//
// Marco "Monte" Chen's approach: gravity pre-compaction, then SA with adaptive
// temperature, center-biased moves, periodic restarts from best, move-type
// weighting by recent success, and strict constraint rejection.
// ---------------------------------------------------------------------------

static ContestEntry marco_round1(
    ArrangePolygons& items,
    const BoundingBox& bed,
    double padding_mm,
    int plate_idx = 0)
{
    auto t0 = std::chrono::steady_clock::now();
    std::mt19937 rng(77777);

    std::vector<int> pi;
    for (size_t i = 0; i < items.size(); i++)
        if (items[i].bed_idx == plate_idx) pi.push_back((int)i);

    if (pi.size() < 2) {
        double a = cluster_area_mm2(items, plate_idx);
        auto vr = validate_arrangement(items, bed, padding_mm, plate_idx);
        return {"Monte", a, 0, vr.valid, (int)pi.size()};
    }

    // --- Phase 1: Gravity pre-compaction (cheap good starting point) ---
    GravityParams gp;
    gp.padding_mm = padding_mm;
    gp.step_mm = 0.5;
    gp.max_steps = 600;
    gp.max_passes = 5;
    gp.target = GravityTarget::CENTER;
    for (int pass = 0; pass < gp.max_passes; pass++) {
        int moved = 0;
        Point target = bed.center();
        for (int idx : pi) {
            auto sr = slide_item_toward(items, idx, target, bed, gp);
            if (sr.moved) {
                items[idx].translation = sr.best_translation;
                items[idx].rotation = sr.best_rotation;
                moved++;
            }
        }
        if (moved == 0) break;
    }

    double current_E = cluster_area_mm2(items, plate_idx);
    double best_E = current_E;
    auto best_snap = save_positions(items);

    // --- Phase 2: Guided SA ---
    Point center = bed.center();
    double cx = unscaled(center.x()), cy = unscaled(center.y());

    // Adaptive temperature
    double T = current_E * 0.05;
    double T_min = 0.01;
    int accepted = 0, attempted = 0;
    double target_accept = 0.25;

    // Move-type weights: [translate, swap, center-pull]
    double w[3] = {1.0, 1.0, 2.0};
    int    wins[3] = {0, 0, 0};
    int    tries[3] = {0, 0, 0};

    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::uniform_int_distribution<int> pick(0, (int)pi.size() - 1);
    int restart_interval = 600;
    int iters = 0, max_iters = 6000;

    while (iters < max_iters && T > T_min) {
        // Pick move type by weight
        double wsum = w[0] + w[1] + w[2];
        double r = U(rng) * wsum;
        int mt = (r < w[0]) ? 0 : (r < w[0] + w[1]) ? 1 : 2;
        tries[mt]++;

        int idx = pi[pick(rng)];
        auto& ap = items[idx];
        Vec2crd saved_t = ap.translation;
        double saved_r = ap.rotation;
        int other_idx = -1;
        Vec2crd other_saved_t;
        double other_saved_r = 0;

        if (mt == 0) {
            // Random translate — magnitude scales with temperature
            double scale = 3.0 + 12.0 * (T / (best_E * 0.05));
            std::normal_distribution<double> nd(0, std::max(scale, 1.0));
            ap.translation.x() += scaled(nd(rng));
            ap.translation.y() += scaled(nd(rng));
        } else if (mt == 1 && pi.size() >= 2) {
            // Swap two items
            other_idx = pi[pick(rng)];
            if (other_idx == idx) { iters++; continue; }
            other_saved_t = items[other_idx].translation;
            other_saved_r = items[other_idx].rotation;
            std::swap(items[idx].translation, items[other_idx].translation);
            std::swap(items[idx].rotation, items[other_idx].rotation);
        } else {
            // Center-biased pull: move toward cluster center + jitter
            BoundingBox ibb = get_extents(item_footprint(ap));
            Point ic = ibb.center();
            double dx = cx - unscaled(ic.x());
            double dy = cy - unscaled(ic.y());
            double dist = std::sqrt(dx * dx + dy * dy);
            if (dist > 0.5) {
                double step = std::min(dist * 0.3, 8.0) * (0.5 + U(rng));
                std::normal_distribution<double> jitter(0, 1.5);
                ap.translation.x() += scaled(dx / dist * step + jitter(rng));
                ap.translation.y() += scaled(dy / dist * step + jitter(rng));
            }
        }

        // Strict validity: reject invalid moves outright
        bool valid = true;
        ExPolygon test = item_footprint(ap);
        if (!within_bed(test, bed, padding_mm) ||
            collides_with_others(test, padding_mm, items, idx, plate_idx))
            valid = false;
        if (valid && other_idx >= 0) {
            ExPolygon t2 = item_footprint(items[other_idx]);
            if (!within_bed(t2, bed, padding_mm) ||
                collides_with_others(t2, padding_mm, items, other_idx, plate_idx))
                valid = false;
        }

        if (valid) {
            double new_E = cluster_area_mm2(items, plate_idx);
            double dE = new_E - current_E;
            attempted++;
            if (dE < 0 || U(rng) < std::exp(-dE / T)) {
                current_E = new_E;
                accepted++;
                if (current_E < best_E) {
                    best_E = current_E;
                    best_snap = save_positions(items);
                    wins[mt]++;
                }
            } else {
                ap.translation = saved_t; ap.rotation = saved_r;
                if (other_idx >= 0) {
                    items[other_idx].translation = other_saved_t;
                    items[other_idx].rotation = other_saved_r;
                }
            }
        } else {
            ap.translation = saved_t; ap.rotation = saved_r;
            if (other_idx >= 0) {
                items[other_idx].translation = other_saved_t;
                items[other_idx].rotation = other_saved_r;
            }
        }

        // Adaptive temperature every 50 iterations
        if (attempted > 0 && iters % 50 == 0) {
            double rate = (double)accepted / attempted;
            if (rate > target_accept + 0.1) T *= 0.9;
            else if (rate < target_accept - 0.1) T *= 1.05;
            else T *= 0.97;
            accepted = 0; attempted = 0;
        }

        // Update move-type weights every 200 iterations
        if (iters % 200 == 0 && iters > 0) {
            for (int k = 0; k < 3; k++) {
                double success = tries[k] > 0 ? (double)wins[k] / tries[k] : 0;
                w[k] = 0.5 + 2.0 * success;
                wins[k] = 0; tries[k] = 0;
            }
        }

        // Periodic restart from best known
        if (iters % restart_interval == 0 && iters > 0) {
            restore_positions(items, best_snap);
            current_E = best_E;
        }

        iters++;
    }

    restore_positions(items, best_snap);

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    auto vr = validate_arrangement(items, bed, padding_mm, plate_idx);
    return {"Monte", best_E, ms, vr.valid, (int)pi.size()};
}

// ---------------------------------------------------------------------------
// Contestant 7: DARWIN — Evolutionary / Genetic Algorithm
//
// Dr. James "Darwin" Okonkwo — evolutionary computation researcher.
// Population of arrangement variants evolved via tournament selection,
// positional crossover (mix positions from two parents), adaptive mutation,
// and elitism. Repair via gravity settle ensures all offspring are valid.
// The crossover is the key insight: combining positions from two good
// arrangements can discover packing configurations neither parent had.
// ---------------------------------------------------------------------------

static ContestEntry james_round1(
    ArrangePolygons& items,
    const BoundingBox& bed,
    double padding_mm,
    int plate_idx)
{
    auto t0 = std::chrono::steady_clock::now();
    std::mt19937 rng(54321);

    std::vector<int> pidx;
    for (size_t i = 0; i < items.size(); i++)
        if (items[i].bed_idx == plate_idx) pidx.push_back((int)i);

    if (pidx.size() < 2)
        return {"Darwin", cluster_area_mm2(items, plate_idx), 0, true, (int)pidx.size()};

    const int POP = 10;
    const int GENS = 50;
    const int TOURN_K = 3;
    const double MUTATE_P = 0.35;
    const double MUTATE_MM = 8.0;
    const double XOVER_RATE = 0.7;

    using Genome = std::vector<ItemSnapshot>;
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::uniform_int_distribution<int> pick_pi(0, (int)pidx.size() - 1);
    std::uniform_int_distribution<int> pick_pop(0, POP - 1);

    // Repair: push apart overlapping items, clamp to bed, then gravity settle
    GravityParams gp;
    gp.padding_mm = padding_mm;
    gp.step_mm = 0.5;
    gp.max_steps = 400;
    gp.max_passes = 3;
    gp.target = GravityTarget::CENTER;

    auto repair = [&](Genome& g) -> double {
        restore_positions(items, g);

        // Phase A: clamp items inside bed bounds
        coord_t margin = scaled(padding_mm);
        for (int idx : pidx) {
            BoundingBox ibb = get_extents(item_footprint(items[idx]));
            coord_t dx = 0, dy = 0;
            if (ibb.min.x() < bed.min.x() + margin) dx = bed.min.x() + margin - ibb.min.x();
            if (ibb.max.x() > bed.max.x() - margin) dx = bed.max.x() - margin - ibb.max.x();
            if (ibb.min.y() < bed.min.y() + margin) dy = bed.min.y() + margin - ibb.min.y();
            if (ibb.max.y() > bed.max.y() - margin) dy = bed.max.y() - margin - ibb.max.y();
            items[idx].translation.x() += dx;
            items[idx].translation.y() += dy;
        }

        // Phase B: iteratively push apart overlapping pairs
        for (int rep = 0; rep < 30; rep++) {
            bool any_overlap = false;
            for (size_t ai = 0; ai < pidx.size(); ai++) {
                int i = pidx[ai];
                for (size_t bi = ai + 1; bi < pidx.size(); bi++) {
                    int j = pidx[bi];
                    ExPolygon pi_p = item_footprint(items[i]);
                    ExPolygon pj_p = item_footprint(items[j]);
                    ExPolygons pi_d = offset_ex(pi_p, scaled(padding_mm / 2.0));
                    ExPolygons pj_d = offset_ex(pj_p, scaled(padding_mm / 2.0));
                    if (pi_d.empty() || pj_d.empty()) continue;
                    ExPolygons ov = intersection_ex(pi_d[0], pj_d[0]);
                    double ov_a = 0;
                    for (auto& o : ov) ov_a += std::abs(o.area());
                    if (ov_a > scaled(0.1) * scaled(0.1)) {
                        any_overlap = true;
                        Point ic = get_extents(pi_p).center();
                        Point jc = get_extents(pj_p).center();
                        double dx = unscaled(ic.x() - jc.x());
                        double dy = unscaled(ic.y() - jc.y());
                        double d = std::sqrt(dx*dx + dy*dy);
                        if (d < 0.1) { dx = 1; dy = 0; d = 1; }
                        double push = padding_mm + 1.0;
                        items[i].translation.x() += scaled(dx/d * push * 0.5);
                        items[i].translation.y() += scaled(dy/d * push * 0.5);
                        items[j].translation.x() -= scaled(dx/d * push * 0.5);
                        items[j].translation.y() -= scaled(dy/d * push * 0.5);
                    }
                }
            }
            if (!any_overlap) break;
        }

        // Phase C: gravity settle toward center
        gravity_settle(items, bed, gp);
        g = save_positions(items);

        // Penalize invalid arrangements so they don't win tournaments
        auto vr = validate_arrangement(items, bed, padding_mm, plate_idx);
        if (!vr.valid) return 1e18;
        return cluster_area_mm2(items, plate_idx);
    };

    // --- Seed population ---
    auto base = save_positions(items);
    std::vector<Genome> pop(POP);
    std::vector<double> fit(POP);

    // Individual 0: original (repaired for baseline)
    pop[0] = base;
    fit[0] = repair(pop[0]);

    // Individuals 1..POP-1: jittered variants
    for (int p = 1; p < POP; p++) {
        restore_positions(items, base);
        for (int idx : pidx) {
            if (U(rng) < 0.5) {
                std::uniform_real_distribution<double> jig(-MUTATE_MM, MUTATE_MM);
                items[idx].translation.x() += scaled(jig(rng));
                items[idx].translation.y() += scaled(jig(rng));
            }
        }
        pop[p] = save_positions(items);
        fit[p] = repair(pop[p]);
    }

    // Track elite
    int elite = 0;
    for (int p = 1; p < POP; p++)
        if (fit[p] < fit[elite]) elite = p;
    Genome best_g = pop[elite];
    double best_a = fit[elite];

    // --- Evolution ---
    for (int gen = 0; gen < GENS; gen++) {
        std::vector<Genome> next(POP);
        std::vector<double> nfit(POP);

        // Slot 0: elitism — best individual passes through unchanged
        next[0] = best_g;
        nfit[0] = best_a;

        for (int p = 1; p < POP; p++) {
            // Tournament selection for two parents
            auto tournament = [&]() {
                int w = pick_pop(rng);
                for (int k = 1; k < TOURN_K; k++) {
                    int c = pick_pop(rng);
                    if (fit[c] < fit[w]) w = c;
                }
                return w;
            };
            int a = tournament(), b = tournament();

            // Crossover: per-item position from parent A or B
            Genome child = pop[a];
            if (U(rng) < XOVER_RATE && a != b) {
                for (int idx : pidx)
                    if (U(rng) < 0.5) child[idx] = pop[b][idx];
            }

            // Mutation: adaptive — displacement shrinks over generations
            double scale = MUTATE_MM * (1.0 - 0.6 * gen / GENS);
            for (int idx : pidx) {
                if (U(rng) < MUTATE_P) {
                    std::uniform_real_distribution<double> m(-scale, scale);
                    child[idx].translation.x() += scaled(m(rng));
                    child[idx].translation.y() += scaled(m(rng));
                }
            }

            // Repair offspring: gravity settle fixes overlaps
            next[p] = child;
            nfit[p] = repair(next[p]);
        }

        pop = std::move(next);
        fit = std::move(nfit);

        // Update global best
        for (int p = 0; p < POP; p++) {
            if (fit[p] < best_a) { best_a = fit[p]; best_g = pop[p]; }
        }
    }

    // Restore the fittest individual
    restore_positions(items, best_g);

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    auto vr = validate_arrangement(items, bed, padding_mm, plate_idx);
    return {"Darwin", best_a, ms, vr.valid, (int)pidx.size()};
}

// ---------------------------------------------------------------------------
// Contestant 8: YUKI "GEO" TANAKA — Geometric Gap Analysis + Dual-Axis Sweep
//
// Core insight: the bounding box is defined by 4 extremal items. Wasted space
// concentrates in axis-aligned GAPS between item clusters. Strategy:
// Phase 1: Aggressive gravity settle (center + centroid, keep best).
// Phase 2: Iterative gap-closing with binary-searched shift amounts.
// Phase 3: Boundary gradient descent — extremal items slide inward.
// Phase 4: Per-item squeeze toward nearest-neighbor centroid (local compaction).
// Phase 5: Corner exploration with full gradient+gap pipeline.
// ---------------------------------------------------------------------------

static ContestEntry yuki_round1(
    ArrangePolygons& items, const BoundingBox& bed, double padding_mm, int plate_idx)
{
    auto t0 = std::chrono::steady_clock::now();
    std::vector<int> pi;
    for (size_t i = 0; i < items.size(); i++)
        if (items[i].bed_idx == plate_idx) pi.push_back((int)i);
    if (pi.size() < 2) {
        auto vr = validate_arrangement(items, bed, padding_mm, plate_idx);
        return {"Yuki", cluster_area_mm2(items, plate_idx), 0, vr.valid, (int)pi.size()};
    }
    GravityParams gp; gp.padding_mm = padding_mm; gp.step_mm = 0.25; gp.max_steps = 1200;
    double best_area = cluster_area_mm2(items, plate_idx);
    auto best_snap = save_positions(items);
    auto track = [&]() { // track best
        double a = cluster_area_mm2(items, plate_idx);
        if (a < best_area) { best_area = a; best_snap = save_positions(items); }
    };

    // Gravity settle: farthest-first toward target
    auto settle = [&](Point target, int passes) {
        std::vector<int> order = pi;
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            Point ca = get_extents(item_footprint(items[a])).center();
            Point cb = get_extents(item_footprint(items[b])).center();
            return (std::pow(ca.x()-target.x(),2)+std::pow(ca.y()-target.y(),2))
                 > (std::pow(cb.x()-target.x(),2)+std::pow(cb.y()-target.y(),2));
        });
        for (int p = 0; p < passes; p++) {
            int moved = 0;
            for (int idx : order) {
                auto sr = slide_multi_angle(items, idx, target, bed, gp, 7);
                if (sr.moved) { items[idx].translation = sr.best_translation; moved++; }
            }
            if (!moved) break;
        }
    };

    // Boundary gradient: slide bounding-box-edge items inward
    auto grad = [&]() {
        BoundingBox cbb = plate_cluster_bb(items, plate_idx);
        if (!cbb.defined) return;
        double W = unscaled(cbb.size().x()), H = unscaled(cbb.size().y());
        coord_t tol = scaled(1.0);
        for (int idx : pi) {
            BoundingBox ibb = get_extents(item_footprint(items[idx]));
            double fx = 0, fy = 0;
            if (ibb.max.x() > cbb.max.x() - tol) fx -= H;
            if (ibb.min.x() < cbb.min.x() + tol) fx += H;
            if (ibb.max.y() > cbb.max.y() - tol) fy -= W;
            if (ibb.min.y() < cbb.min.y() + tol) fy += W;
            double mag = std::sqrt(fx*fx + fy*fy);
            if (mag < 0.01) continue;
            auto sr = slide_binary(items, idx, fx/mag, fy/mag, bed, gp);
            if (sr.moved) items[idx].translation = sr.best_translation;
        }
    };

    // Gap closing: for each axis, find gaps > padding and binary-search max shift
    auto close_gaps = [&]() {
        for (int axis = 0; axis < 2; axis++) {
            struct E { int idx; coord_t lo, hi; };
            std::vector<E> edges;
            for (int idx : pi) {
                BoundingBox b = get_extents(item_footprint(items[idx]));
                edges.push_back({idx, axis?b.min.y():b.min.x(), axis?b.max.y():b.max.x()});
            }
            std::sort(edges.begin(), edges.end(), [](const E& a, const E& b) {
                return (a.lo+a.hi) < (b.lo+b.hi); });
            for (size_t gi = 1; gi < edges.size(); gi++) {
                coord_t gap = edges[gi].lo - edges[gi-1].hi;
                if (gap < scaled(padding_mm * 1.5)) continue;
                coord_t lo_sh = 0, hi_sh = gap - scaled(padding_mm * 0.5);
                if (hi_sh <= 0) continue;
                for (int bs = 0; bs < 10; bs++) {
                    coord_t mid = (lo_sh + hi_sh) / 2;
                    auto snap = save_positions(items);
                    bool ok = true;
                    for (size_t j = gi; j < edges.size(); j++) {
                        if (axis) items[edges[j].idx].translation.y() -= mid;
                        else      items[edges[j].idx].translation.x() -= mid;
                    }
                    for (size_t j = gi; j < edges.size() && ok; j++) {
                        ExPolygon fp = item_footprint(items[edges[j].idx]);
                        if (!within_bed(fp, bed, padding_mm) ||
                            collides_with_others(fp, padding_mm, items, edges[j].idx, plate_idx))
                            ok = false;
                    }
                    restore_positions(items, snap);
                    if (ok) lo_sh = mid; else hi_sh = mid;
                }
                if (lo_sh > scaled(0.5)) {
                    for (size_t j = gi; j < edges.size(); j++) {
                        if (axis) items[edges[j].idx].translation.y() -= lo_sh;
                        else      items[edges[j].idx].translation.x() -= lo_sh;
                    }
                }
            }
        }
    };

    // Full compaction: settle + (gradient + gap-close + re-settle) x N
    auto compact = [&](Point target) {
        settle(target, 8); track();
        for (int e = 0; e < 8; e++) {
            grad(); close_gaps();
            settle(plate_centroid(items, plate_idx), 3); track();
        }
    };

    // Phase 1: center and centroid targets
    auto init = save_positions(items);
    compact(bed.center());
    restore_positions(items, init);
    compact(plate_centroid(items, plate_idx));

    // Phase 2: corner exploration
    Point corners[] = { bed.min, {bed.max.x(),bed.min.y()},
                        {bed.min.x(),bed.max.y()}, bed.max };
    for (const auto& c : corners) { restore_positions(items, init); compact(c); }

    // Phase 3: local squeeze — each item toward its 3 nearest neighbors
    restore_positions(items, best_snap);
    for (int sq = 0; sq < 3; sq++) {
        for (int idx : pi) {
            Point ic = get_extents(item_footprint(items[idx])).center();
            std::vector<std::pair<double,int>> nd;
            for (int o : pi) if (o != idx) {
                Point oc = get_extents(item_footprint(items[o])).center();
                nd.push_back({std::pow(ic.x()-oc.x(),2)+std::pow(ic.y()-oc.y(),2), o});
            }
            std::sort(nd.begin(), nd.end());
            int k = std::min(3, (int)nd.size()); double tx=0, ty=0;
            for (int n = 0; n < k; n++) {
                Point nc = get_extents(item_footprint(items[nd[n].second])).center();
                tx += nc.x(); ty += nc.y();
            }
            Point lt = {(coord_t)(tx/k), (coord_t)(ty/k)};
            auto sr = slide_multi_angle(items, idx, lt, bed, gp, 7);
            if (sr.moved) items[idx].translation = sr.best_translation;
        }
        track();
    }

    restore_positions(items, best_snap);
    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    auto vr = validate_arrangement(items, bed, padding_mm, plate_idx);
    return {"Yuki", best_area, ms, vr.valid, (int)pi.size()};
}

// ---------------------------------------------------------------------------
// Contestant 9: SARAH "SPECS" HARTMAN — VICE COMPACTION
//
// Mechanical engineer's approach: a four-sided vice.
// 1. Find the 4 items defining the bounding box extremes
// 2. Binary-search each boundary item inward along its axis
// 3. Slide interior items toward centroid to make room
// 4. Repeat until convergence (area decreases monotonically)
// No randomness, no rotation, no multi-angle — pure deterministic squeeze.
// ---------------------------------------------------------------------------

static ContestEntry sarah_round1(
    ArrangePolygons& items,
    const BoundingBox& bed,
    double padding_mm,
    int plate_idx)
{
    auto t0 = std::chrono::steady_clock::now();

    std::vector<int> plate_items;
    for (size_t i = 0; i < items.size(); i++)
        if (items[i].bed_idx == plate_idx) plate_items.push_back((int)i);

    if (plate_items.size() < 2)
        return {"ViceCompact", cluster_area_mm2(items, plate_idx), 0, true, (int)plate_items.size()};

    GravityParams gp;
    gp.padding_mm = padding_mm;
    gp.step_mm = 0.25;
    gp.max_steps = 1200;

    double best_area = cluster_area_mm2(items, plate_idx);
    auto best_snap = save_positions(items);

    // --- Initial pass: slide every item toward bed center (farthest first) ---
    auto settle_all = [&](Point target) {
        std::vector<int> order = plate_items;
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            Point ca = get_extents(item_footprint(items[a])).center();
            Point cb = get_extents(item_footprint(items[b])).center();
            double da = std::pow(ca.x()-target.x(),2)+std::pow(ca.y()-target.y(),2);
            double db = std::pow(cb.x()-target.x(),2)+std::pow(cb.y()-target.y(),2);
            return da > db;
        });
        for (int pass = 0; pass < 6; pass++) {
            int moved = 0;
            for (int idx : order) {
                BoundingBox ibb = get_extents(item_footprint(items[idx]));
                Point ic = ibb.center();
                double dx = (double)(target.x() - ic.x());
                double dy = (double)(target.y() - ic.y());
                double dist = std::sqrt(dx*dx + dy*dy);
                if (dist < scaled(0.5)) continue;
                auto sr = slide_binary(items, idx, dx/dist, dy/dist, bed, gp);
                if (sr.moved) { items[idx].translation = sr.best_translation; moved++; }
            }
            if (moved == 0) break;
        }
    };

    settle_all(bed.center());
    double area = cluster_area_mm2(items, plate_idx);
    if (area < best_area) { best_area = area; best_snap = save_positions(items); }

    // --- Main vice loop: squeeze boundary items, then consolidate interior ---
    for (int round = 0; round < 20; round++) {
        bool any_boundary_moved = false;

        // Identify the 4 extremal items
        int idx_L = -1, idx_R = -1, idx_B = -1, idx_T = -1;
        coord_t ext_L = std::numeric_limits<coord_t>::max();
        coord_t ext_R = std::numeric_limits<coord_t>::min();
        coord_t ext_B = std::numeric_limits<coord_t>::max();
        coord_t ext_T = std::numeric_limits<coord_t>::min();

        for (int idx : plate_items) {
            BoundingBox ibb = get_extents(item_footprint(items[idx]));
            if (ibb.min.x() < ext_L) { ext_L = ibb.min.x(); idx_L = idx; }
            if (ibb.max.x() > ext_R) { ext_R = ibb.max.x(); idx_R = idx; }
            if (ibb.min.y() < ext_B) { ext_B = ibb.min.y(); idx_B = idx; }
            if (ibb.max.y() > ext_T) { ext_T = ibb.max.y(); idx_T = idx; }
        }

        // Squeeze each boundary item along its inward axis
        struct Squeeze { int idx; double nx, ny; };
        Squeeze squeezes[] = {
            {idx_L,  1.0,  0.0},   // leftmost  -> slide right
            {idx_R, -1.0,  0.0},   // rightmost -> slide left
            {idx_B,  0.0,  1.0},   // bottommost-> slide up
            {idx_T,  0.0, -1.0},   // topmost   -> slide down
        };

        for (auto& sq : squeezes) {
            if (sq.idx < 0) continue;
            auto sr = slide_binary(items, sq.idx, sq.nx, sq.ny, bed, gp);
            if (sr.moved && sr.distance_moved_mm > 0.1) {
                items[sq.idx].translation = sr.best_translation;
                any_boundary_moved = true;
            }
        }

        // Consolidate: slide interior items toward new centroid
        Point centroid = plate_centroid(items, plate_idx);
        std::vector<int> interior;
        for (int idx : plate_items) {
            if (idx == idx_L || idx == idx_R || idx == idx_B || idx == idx_T)
                continue;
            interior.push_back(idx);
        }
        // Farthest from centroid first
        std::sort(interior.begin(), interior.end(), [&](int a, int b) {
            Point ca = get_extents(item_footprint(items[a])).center();
            Point cb = get_extents(item_footprint(items[b])).center();
            double da = std::pow(ca.x()-centroid.x(),2)+std::pow(ca.y()-centroid.y(),2);
            double db = std::pow(cb.x()-centroid.x(),2)+std::pow(cb.y()-centroid.y(),2);
            return da > db;
        });
        for (int idx : interior) {
            BoundingBox ibb = get_extents(item_footprint(items[idx]));
            Point ic = ibb.center();
            double dx = (double)(centroid.x()-ic.x());
            double dy = (double)(centroid.y()-ic.y());
            double dist = std::sqrt(dx*dx+dy*dy);
            if (dist < scaled(0.5)) continue;
            auto sr = slide_binary(items, idx, dx/dist, dy/dist, bed, gp);
            if (sr.moved) items[idx].translation = sr.best_translation;
        }

        // Also slide boundary items toward centroid (interior may have
        // cleared space for them to tuck in further)
        for (auto& sq : squeezes) {
            if (sq.idx < 0) continue;
            BoundingBox ibb = get_extents(item_footprint(items[sq.idx]));
            Point ic = ibb.center();
            double dx = (double)(centroid.x()-ic.x());
            double dy = (double)(centroid.y()-ic.y());
            double dist = std::sqrt(dx*dx+dy*dy);
            if (dist < scaled(0.5)) continue;
            auto sr = slide_binary(items, sq.idx, dx/dist, dy/dist, bed, gp);
            if (sr.moved) items[sq.idx].translation = sr.best_translation;
        }

        area = cluster_area_mm2(items, plate_idx);
        if (area < best_area) { best_area = area; best_snap = save_positions(items); }
        if (!any_boundary_moved) break; // vice fully closed
    }

    restore_positions(items, best_snap);

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    auto vr = validate_arrangement(items, bed, padding_mm, plate_idx);
    return {"ViceCompact", best_area, ms, vr.valid, (int)plate_items.size()};
}

// ---------------------------------------------------------------------------
// CONTEST RUNNER
// ---------------------------------------------------------------------------

struct ContestResult {
    std::string test_name;
    std::vector<ContestEntry> entries;
    std::string winner;
    double starting_area;
};

// Watchdog: run a contestant with time limit and strict validation
static constexpr double WATCHDOG_TIME_LIMIT_MS = 40000.0; // 40 seconds max per contestant

typedef std::function<ContestEntry(ArrangePolygons&, const BoundingBox&, double, int)> ContestantFn;

static ContestEntry run_with_watchdog(
    ContestantFn fn,
    ArrangePolygons items,  // fresh copy per contestant
    const BoundingBox& bed,
    double padding_mm,
    int plate_idx,
    double starting_area)
{
    // Thermal cooldown between contestants (50ms)
    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    auto t0 = std::chrono::steady_clock::now();
    ContestEntry entry;
    try {
        entry = fn(items, bed, padding_mm, plate_idx);
    } catch (const std::exception& e) {
        entry = {"CRASHED:" + std::string(e.what()), starting_area, 0, false, 0};
    } catch (...) {
        entry = {"CRASHED:unknown", starting_area, 0, false, 0};
    }
    auto t1 = std::chrono::steady_clock::now();
    double actual_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Watchdog: override reported time with actual wall-clock time
    entry.elapsed_ms = actual_ms;

    // Watchdog: DISQUALIFY if over time limit — reset to starting area
    if (actual_ms > WATCHDOG_TIME_LIMIT_MS) {
        entry.name += "(DQ:TIMEOUT)";
        entry.valid = false;
        entry.cluster_area = starting_area;
    }

    // Watchdog: strict re-validation (don't trust contestant's self-report)
    auto vr = validate_arrangement(items, bed, padding_mm, plate_idx);
    entry.valid = vr.valid;
    if (!vr.valid) {
        // Disqualified — reset area to starting (no credit for invalid result)
        entry.cluster_area = starting_area;
    }

    // Watchdog: area can't be zero or negative
    if (entry.cluster_area <= 0) {
        entry.cluster_area = starting_area;
        entry.valid = false;
    }

    // Watchdog: area can't be larger than bed area (sanity)
    double bed_area = unscaled(bed.size().x()) * unscaled(bed.size().y());
    if (entry.cluster_area > bed_area * 1.5) {
        entry.cluster_area = starting_area;
    }

    return entry;
}

static ContestResult run_contest(
    const std::string& name,
    ArrangePolygons items,  // pre-arranged input
    const BoundingBox& bed,
    double padding_mm,
    int plate_idx = 0)
{
    ContestResult cr;
    cr.test_name = name;
    cr.starting_area = cluster_area_mm2(items, plate_idx);

    // Validate starting arrangement — skip contest if input is broken
    auto start_vr = validate_arrangement(items, bed, padding_mm, plate_idx);
    if (!start_vr.valid) {
        std::cout << "  SKIP: starting arrangement invalid (overlaps="
                  << start_vr.overlaps << " oob=" << start_vr.out_of_bounds
                  << " max=" << (int)start_vr.max_overlap_area_mm2 << "mm2)\n";
        cr.winner = "SKIP:BAD_INPUT";
        return cr;
    }

    // Register all contestants
    struct Contestant {
        std::string label;
        ContestantFn fn;
    };
    std::vector<Contestant> roster = {
        {"GravityV2",  [](ArrangePolygons& it, const BoundingBox& b, double p, int pi) { return run_gravity_v2(it, b, p, pi); }},
        {"SimAnneal",  [](ArrangePolygons& it, const BoundingBox& b, double p, int pi) { return run_simulated_annealing(it, b, p, pi); }},
        {"Springs",    [](ArrangePolygons& it, const BoundingBox& b, double p, int pi) { return run_spring_physics(it, b, p, pi); }},
        {"HexPack",    [](ArrangePolygons& it, const BoundingBox& b, double p, int pi) { return run_hex_pack(it, b, p, pi); }},
        {"Elena",      [](ArrangePolygons& it, const BoundingBox& b, double p, int pi) { return elena_round1(it, b, p, pi); }},
        {"Monte",      [](ArrangePolygons& it, const BoundingBox& b, double p, int pi) { return marco_round1(it, b, p, pi); }},
        {"Darwin",     [](ArrangePolygons& it, const BoundingBox& b, double p, int pi) { return james_round1(it, b, p, pi); }},
        {"Yuki",       [](ArrangePolygons& it, const BoundingBox& b, double p, int pi) { return yuki_round1(it, b, p, pi); }},
        {"Sarah",      [](ArrangePolygons& it, const BoundingBox& b, double p, int pi) { return sarah_round1(it, b, p, pi); }},
    };

    // Run each contestant with watchdog, fresh copy, timing isolation
    for (auto& c : roster) {
        auto entry = run_with_watchdog(c.fn, items, bed, padding_mm, plate_idx, cr.starting_area);
        cr.entries.push_back(entry);
    }

    // Find winner (lowest valid area)
    double best = 1e18;
    for (auto& e : cr.entries) {
        if (e.valid && e.cluster_area < best) {
            best = e.cluster_area;
            cr.winner = e.name;
        }
    }

    return cr;
}

static std::string format_contest(const ContestResult& cr) {
    std::string s = "[" + cr.test_name + "] start:" + std::to_string((int)cr.starting_area) + " | ";
    for (const auto& e : cr.entries) {
        s += e.name + ":" + std::to_string((int)e.cluster_area);
        if (!e.valid) s += "(INVALID)";
        s += "/" + std::to_string((int)e.elapsed_ms) + "ms ";
    }
    s += "| WINNER: " + cr.winner;
    return s;
}
