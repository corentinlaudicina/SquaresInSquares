#pragma once

// Single-population GA at a fixed box with a ratchet (Phase 3 items 3-5).
// After jcottaar/packing (Santa 2025 1st place, CC BY-SA 4.0); an independent reimplementation for
// squares in C++, with no code copied.
//
// The population shares one box side h, so fitness -- lexicographic (h, residual penalty) -- reduces
// to the residual penalty within a generation. A generation:
//   1. overproduce: `offspring` children, each a parent (uniform) changed by one move (moves.hpp);
//   2. pre-filter: keep the best `prefilter` share by raw penalty, before any relaxation;
//   3. relax rough (big step cap), keep the best `fine_frac` share, relax those fine (small caps);
//   4. select `population` survivors from parents + offspring (diversity.hpp: elite + spread);
//   5. ratchet: if the best survivor's depth and wall protrusion are within `tol`, legalize a copy
//      (legalize.hpp; strict, CPU, double) and keep it if its side beats the best legal side; then
//      shrink h by `ratchet_shrink` for the whole population, rescale everyone and relax them fine.
// Per-move survivor rates are counted every generation. All randomness is drawn on the calling
// thread in a fixed order and every relaxation is independent, so a seed fixes the run exactly,
// whatever the thread count.

#include "pop/diversity.hpp"
#include "pop/moves.hpp"
#include "pop/relax_capped.hpp"
#include "system.hpp"

#include <nlohmann/json.hpp>

#include <functional>

namespace pack {

struct GAConfig {
    int population = 40;
    int offspring = 160;
    double prefilter = 0.6;   // share of offspring relaxed at all
    double fine_frac = 0.5;   // share of those relaxed fine after the rough pass
    int threads = 0;          // relaxation threads; 0 = hardware concurrency (fixed at start)

    int penalty = 1;          // 0 = area2, 1 = depth2
    int history = 3;
    double rough_cap = 0.05;
    long rough_iters = 100;
    double fine_cap = 0.0125;  // first fine stage, each next one 4x smaller
    long fine_iters = 150;
    long fine_stages = 3;
    double tol = 1e-9;         // residual depth / protrusion that counts as feasible for the ratchet

    double ratchet_shrink = 2e-4;  // relative shrink of h at each ratchet
    double init_jiggle = 0.02;     // sigma of the noise that makes the initial population from the start

    MoveParams moves;
    DiversityParams diversity;
};

struct GAOutcome {
    System best;            // best legal packing (the legalized start if nothing beat it)
    double best_side = 0.0;
    nlohmann::json telemetry;
};

// `start` is legalized first. Budget: wall seconds and/or generations (0 = unlimited; one must be set).
// `update`, if set, is called before each generation with the elapsed budget fraction (schedule curves).
GAOutcome run_ga(const System& start, const GAConfig& cfg, unsigned long long seed, double wall_seconds,
                 long max_generations, const std::function<void(double, GAConfig&)>& update = {});

}  // namespace pack
