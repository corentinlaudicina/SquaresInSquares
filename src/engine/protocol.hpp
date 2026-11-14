#pragma once

// The two search protocols of main_shrink.cpp and main_compress.cpp, as library functions driven by
// one config struct. Defaults are those drivers' constants, so a default-constructed config with the
// right mode reproduces them bit for bit for the same seed and a kick/cycle budget. (The old
// drivers are kept unchanged as the reference.)
//
// Differences from the drivers, none of which touch the random stream or the search:
//   * no data/best writes, no trajectory, no fixed log paths -- the caller decides what to write;
//   * a wall-clock budget (std::steady_clock) alongside the kick budget, instead of process CPU;
//   * telemetry is returned as JSON.

#include "moves/lnm.hpp"
#include "system.hpp"

#include <nlohmann/json.hpp>

#include <string>

namespace pack {

struct KickConfig {
    double large_base = 0.1;   // chance of a large move instead of a shake
    double large_max = 0.6;    // ceiling of the escalation
    double large_step = 0.1;   // added per `escalate_every` failures (shrink) / stagnant cycles (cycle)
    long escalate_every = 25;
    long shake_sweeps = 100;
    double shake_T = 0.05;
    double shake_translate = 0.1;
    double shake_rotate = 0.1;
    LnmOptions lnm;            // large-move weights and parameters; log/watch are set by the engine
};

struct EngineConfig {
    std::string mode = "shrink";  // "shrink" (main_shrink) | "cycle" (main_compress)

    // Start: mirror the seed through the main diagonal, relax, compress.
    bool mirror = true;
    double prepare_ds = 0.01;     // main_compress uses 0.005 here and in its cycles
    long prepare_max_relax = 5000;

    KickConfig kick;

    // shrink mode
    double eps = 1e-3, eps_min = 1e-7, eps_max = 0.05, eps_up = 1.2, eps_down = 0.95;

    // cycle mode
    double cycle_ds = 0.005;
    long cycle_max_relax = 5000;
    double side_T = 0.01;         // Metropolis temperature on the side, in units of L / sqrt(n)
    long restart_after = 200;     // stagnant cycles before returning to the best

    bool propagate = false;       // the in-memory per-n pool (symmetry.hpp)
    long propagate_every = 1;

    // Budget: stop at whichever comes first. 0 = unlimited (at least one must be set).
    double wall_seconds = 0.0;
    long max_kicks = 0;
};

struct EngineResult {
    bool found = false;     // the start certified; false means `best` is the unprocessed seed
    System best;            // smallest certified packing found (not yet legalized)
    double best_side = 0.0;
    nlohmann::json telemetry;
};

EngineConfig shrink_defaults();
EngineConfig cycle_defaults();  // main_compress.cpp's constants

EngineResult run_engine(const System& seed, const EngineConfig& cfg, unsigned long long seed_value);

}  // namespace pack
