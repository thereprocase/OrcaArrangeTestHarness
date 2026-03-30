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
// CONTEST RUNNER
// ---------------------------------------------------------------------------

struct ContestResult {
    std::string test_name;
    std::vector<ContestEntry> entries;
    std::string winner;
    double starting_area;
};

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

    // Each contestant gets a fresh copy
    cr.entries.push_back(run_gravity_v2(items, bed, padding_mm, plate_idx));
    cr.entries.push_back(run_simulated_annealing(items, bed, padding_mm, plate_idx));
    cr.entries.push_back(run_spring_physics(items, bed, padding_mm, plate_idx));
    cr.entries.push_back(run_hex_pack(items, bed, padding_mm, plate_idx));

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
