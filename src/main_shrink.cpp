// Shrinks a packing one step at a time. From the best packing, pick a target box a little smaller
// than it, squeeze the packing into that box, kick it, and relax the overlaps away at that fixed
// box. A relaxation that reaches a valid packing is a new best (certified, then the next target
// is smaller still). A relaxation that does not is thrown away.
//
// Compared with main_compress.cpp:
//   * no compression schedule: one relaxation per kick, at a fixed box, instead of up to 500;
//   * no wall-free relaxation and no box refit before the compression;
//   * no restarts: the walker is always the best packing, so there is nothing to go back to;
//   * the step (eps, the relative shrink below the best) adapts to the success rate: a success
//     raises it, a failure lowers it, so it settles where about one kick in five succeeds.
//
// The search stops when eps falls below kMinEps (no target is worth trying any more), or when the
// kick or CPU budget runs out.

#include "io/io.hpp"
#include "moves/lnm.hpp"
#include "moves/mc.hpp"
#include "moves/relax.hpp"
#include "moves/symmetry.hpp"
#include "system.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <map>
#include <random>
#include <string>

using namespace pack;

namespace {

// Run parameters: edit these and rebuild, or override on the command line where there is an option.
constexpr int kN = 29;  // loads data/init_external/n<kN>.cfg
constexpr unsigned long long kSeed = 1ULL;
constexpr long kKicks = 20;

// The kick (the same values as main_compress.cpp, so the two drivers can be compared).
constexpr double kLargeMoveProb = 0.1;  // base chance of a large move instead of a shake
constexpr double kLargeMoveMax = 0.6;   // ceiling the escalation can reach
constexpr long kEscalateEvery = 25;     // failed kicks in a row per +0.1 on the large-move chance
constexpr long kShakeSweeps = 100;      // Monte Carlo sweeps in one shake
constexpr double kShakeT = 0.05;        // temperature of the shake
constexpr double kShakeStep = 0.1;      // largest translation per move
constexpr double kShakeTurn = 0.1;      // largest rotation per move, radians

// The shrink.
constexpr double kStartEps = 1e-3;  // relative shrink of the target below the best, at the start
constexpr double kMaxEps = 0.05;    // largest relative shrink the adaptation can reach
constexpr double kMinEps = 1e-7;    // below this, stop: no target is worth trying
constexpr double kEpsUp = 1.2;      // eps multiplier after a success
constexpr double kEpsDown = 0.95;   // eps multiplier after a failure

// The compression, used once, for the starting point.
constexpr double kDs = 0.01;     // box shrink per step, as a length
constexpr long kMaxIter = 5000;  // minimisations allowed in that compression

constexpr double kTraceSeconds = 0.2;  // least time between two frames of the same phase
constexpr long kPropagateEvery = 1;    // propagate the pool every this many new bests

static_assert(kN >= 1, "kN must be at least 1");
static_assert(kEpsUp > 1.0 && kEpsDown < 1.0, "a success must grow eps and a failure shrink it");

struct RunOptions {
    int n = kN;
    unsigned long long seed = kSeed;
    long kicks = kKicks;
    double cpu_budget = 0.0;  // process CPU seconds; above zero it replaces `kicks`
    bool save = true;         // write the record and the trajectory
    std::string csv;          // the per-move table; empty picks logs/moves_n<n>.csv
    double symmetrize = 0.0;  // LnmOptions::w_symmetrize; off unless asked for
    bool propagate = false;   // use and feed the per-n pool
    double eps = kStartEps;   // starting relative shrink
};

double option_value(int argc, char** argv, int& k) {
    if (k + 1 >= argc) {
        std::fprintf(stderr, "%s needs a value\n", argv[k]);
        std::exit(2);
    }
    return std::strtod(argv[++k], nullptr);
}

RunOptions parse_options(int argc, char** argv) {
    RunOptions o;
    for (int k = 1; k < argc; ++k) {
        const std::string a = argv[k];
        if (a == "--n") o.n = static_cast<int>(option_value(argc, argv, k));
        else if (a == "--seed") o.seed = static_cast<unsigned long long>(option_value(argc, argv, k));
        else if (a == "--kicks") o.kicks = static_cast<long>(option_value(argc, argv, k));
        else if (a == "--cpu") o.cpu_budget = option_value(argc, argv, k);
        else if (a == "--eps") o.eps = option_value(argc, argv, k);
        else if (a == "--no-save") o.save = false;
        else if (a == "--symmetrize") o.symmetrize = option_value(argc, argv, k);
        else if (a == "--propagate") o.propagate = true;
        else if (a == "--csv" && k + 1 < argc) o.csv = argv[++k];
        else {
            std::fprintf(stderr, "unknown option %s\n", argv[k]);
            std::exit(2);
        }
    }
    return o;
}

double process_cpu_seconds() { return static_cast<double>(std::clock()) / CLOCKS_PER_SEC; }

System load_config(int n) {
    char path[64];
    std::snprintf(path, sizeof path, "data/init_external/n%04d.cfg", n);
    std::string why;
    auto sys = load_record(path, &why);
    if (!sys) {
        std::fprintf(stderr, "%s\n", why.c_str());
        std::exit(1);
    }
    return *sys;
}

std::FILE* open_trace(int n, const System& sys) {
    char path[64];
    std::snprintf(path, sizeof path, "logs/n%04d.traj", n);
    std::FILE* f = open_trajectory(path, sys);
    if (!f) std::fprintf(stderr, "warning: could not open %s\n", path);
    return f;
}

MoveId cycle_row(MoveId kick) {
    switch (kick) {
        case MoveId::Destroy: return MoveId::CycleDestroy;
        case MoveId::Strip: return MoveId::CycleStrip;
        case MoveId::Corner: return MoveId::CycleCorner;
        case MoveId::Loosen: return MoveId::CycleLoosen;
        case MoveId::Symmetrize: return MoveId::CycleSymmetrize;
        default: return MoveId::CycleShake;
    }
}

const char* row_name(MoveId id) {
    switch (id) {
        case MoveId::CycleDestroy: return "destroy";
        case MoveId::CycleStrip: return "strip";
        case MoveId::CycleCorner: return "corner";
        case MoveId::CycleLoosen: return "loosen";
        case MoveId::CycleSymmetrize: return "symmetrize";
        default: return "shake";
    }
}

// Returns the row the whole kick is credited to. `large_prob` is the chance of a large move.
MoveId kick(System& sys, std::mt19937_64& rng, const LnmOptions& lnm, LnmStats& stats,
            MoveLog& log, double large_prob) {
    if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) < large_prob) {
        stats.last = MoveId::Count;
        large_neighbour_move(sys, rng, lnm, &stats);
        return cycle_row(stats.last);
    }
    StepSizes steps;
    steps.translate = kShakeStep;
    steps.rotate = kShakeTurn;
    for (long k = 0; k < kShakeSweeps; ++k) sweep(sys, steps, kShakeT, rng, &log);
    return MoveId::CycleShake;
}

struct KickTally {
    long tries = 0, accepted = 0, improved = 0;
    double cpu = 0.0, gain = 0.0;
};

}  // namespace

int main(int argc, char** argv) {
    const RunOptions run = parse_options(argc, argv);
    std::mt19937_64 rng(run.seed);

    MoveLog log;
    LnmOptions lnm;
    lnm.log = &log;
    lnm.w_symmetrize = run.symmetrize;
    LnmStats lnm_stats;  // outlives the kicks, so the symmetrize bandit can learn across them

    System start = load_config(run.n);
    std::FILE* trace = run.save ? open_trace(run.n, start) : nullptr;
    long kick_no = 0;

    std::string last_label;
    auto last_frame = std::chrono::steady_clock::time_point{};
    FrameWatch watch;
    if (trace) {
        watch = [&](const System& s, const char* label) {
            const auto now = std::chrono::steady_clock::now();
            const std::chrono::duration<double> since = now - last_frame;
            if (label == last_label && since.count() < kTraceSeconds) return;
            last_label = label;
            last_frame = now;
            append_frame(trace, kick_no, s, label);
        };
        watch(start, "initial");
    }
    lnm.watch = watch;

    // --- the starting point: mirror, relax, compress (once; the shrink takes over from here) ---
    mirror_symmetry(start);
    if (watch) watch(start, "mirror");
    RelaxOptions init_relax;
    init_relax.watch = watch;
    relax(start, init_relax);
    const CompressResult initial = anneal(start, kDs, kMaxIter, watch);
    if (!initial.found) {
        std::fprintf(stderr, "n %d: the mirrored configuration never certified\n", run.n);
        return 1;
    }
    std::printf("n %d: mirrored and compressed to side %.9f\n", run.n, initial.side);
    if (watch) watch(start, "result");

    System best = start;
    double best_side = start.L;

    const double cpu_start = process_cpu_seconds();

    // The pool works on its own random stream, so that turning it on does not reshuffle the kicks of the run it is compared against until it actually changes something.
    Pool pool;
    PropagateOptions pull;
    std::mt19937_64 pool_rng(run.seed ^ 0x9e3779b97f4a7c15ULL);
    std::map<int, double> pool_gains;
    const auto add_gains = [&](const PropagateReport& r) {
        for (const auto& [m, gain] : r.improved) pool_gains[m] += gain;
    };
    long new_bests = 0;
    if (run.propagate) pool.offer(best);

    std::map<MoveId, KickTally> tally;
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    double eps = run.eps;
    long fails_in_a_row = 0;
    double best_reached = best_side;  // the smallest side any kick certified

    // --- the shrink loop ---------------------------------------------------------------------
    const auto more = [&] {
        if (eps < kMinEps) return false;  // no target worth trying: the search has converged
        if (run.cpu_budget > 0.0) return process_cpu_seconds() - cpu_start < run.cpu_budget;
        return kick_no < run.kicks;
    };

    for (; more(); ++kick_no) {
        const double t0 = process_cpu_seconds();

        // The target box, and how hard to kick: the longer kicks fail in a row, the bigger the kick.
        const double target = best_side * (1.0 - eps);
        const double large_prob =
            std::min(kLargeMoveMax, kLargeMoveProb + 0.1 * static_cast<double>(fails_in_a_row / kEscalateEvery));

        // 1. squeeze the best into the smaller box: the packing now overlaps, so it is not valid
        System trial = best;
        trial.rescale_box(target);
        if (watch) watch(trial, "squeezed");

        // 2. kick inside that box
        const MoveId row = kick(trial, rng, lnm, lnm_stats, log, large_prob);
        if (watch) watch(trial, "kicked");

        // 3. relax the overlaps at this fixed box: one minimisation, no walls removed
        RelaxOptions ro;
        ro.watch = watch;
        ro.watch_label = "relax";
        const RelaxResult rr = relax(trial, ro);

        // 4. a valid packing in the target box is a new best. Its side is the tight bounding box,
        // which can only be at most the target.
        bool accepted = false;
        double gain = 0.0;
        if (rr.feasible) {
            trial.fit_box();
            if (trial.is_feasible(kCertifyTol) && trial.L < best_side) {
                accepted = true;
                gain = best_side - trial.L;
                best = trial;
                best_side = trial.L;
                best_reached = std::min(best_reached, best_side);
            }
        }

        // 5. bookkeeping
        const double kick_cpu = process_cpu_seconds() - t0;
        log.record(row, accepted, gain, kick_cpu);
        KickTally& k = tally[row];
        ++k.tries;
        k.cpu += kick_cpu;
        if (accepted) {
            ++k.accepted;
            ++k.improved;
            k.gain += gain;
        }

        if (accepted) {
            eps = std::min(kMaxEps, eps * kEpsUp);
            fails_in_a_row = 0;
            std::printf("kick %ld  %-10s  new best side %.9f  (-%.3e)  eps %.2e\n", kick_no, row_name(row),
                        best_side, gain, eps);
            if (watch) watch(best, "best");

            if (run.propagate) {
                pool.offer(best);
                if (++new_bests % kPropagateEvery == 0) {
                    add_gains(propagate(pool, run.n, pull, pool_rng, &log));
                }
            }
        } else {
            eps *= kEpsDown;
            ++fails_in_a_row;
        }

        // 6. with --propagate: take the pool's packing when it holds something better.
        if (run.propagate) {
            const PoolEntry* e = pool.at(run.n);
            if (e && e->side < best_side - 1e-12) {
                best = e->packing;
                best_side = e->side;
                best_reached = std::min(best_reached, best_side);
                fails_in_a_row = 0;
            }
        }
    }

    // --- report ------------------------------------------------------------------------------
    const double cpu_used = process_cpu_seconds() - cpu_start;
    std::printf("\nn %d: best side %.9f after %ld kicks, %.1f CPU s, final eps %.2e\n", run.n, best_side,
                kick_no, cpu_used, eps);
    std::printf("%-10s %7s %7s %7s %9s %12s %14s\n", "kick", "tries", "accept", "improve", "cpu s",
                "gain", "gain per cpu s");
    for (const auto& [id, k] : tally) {
        std::printf("%-10s %7ld %7ld %7ld %9.1f %12.3e %14.3e\n", row_name(id), k.tries, k.accepted,
                    k.improved, k.cpu, k.gain, k.cpu > 0.0 ? k.gain / k.cpu : 0.0);
    }
    int other_improved = 0;
    double other_gain = 0.0;
    for (const auto& [m, g] : pool_gains) {
        std::printf("pool: n=%d improved by %.3e in side\n", m, g);
        if (m != run.n) {
            ++other_improved;
            other_gain += g;
        }
    }

    // The per-move table is written even with --no-save: it is a measurement, not a claim.
    std::string csv = run.csv;
    if (csv.empty()) {
        char path[64];
        std::snprintf(path, sizeof path, "logs/moves_n%04d.csv", run.n);
        csv = path;
    }
    if (!log.write_csv(csv)) std::fprintf(stderr, "warning: could not write %s\n", csv.c_str());

    if (run.save) {
        const LogVerdict v = save_record("data/best", best, best_side, true, kCertifyTol);
        std::printf("record %s: %s (%s)\n", v.kept ? "kept" : "not kept", v.path.c_str(),
                    v.reason.c_str());
    }

    // The one line scripts/benchmark_moves.py reads. Same fields as main_compress.cpp.
    std::printf("result n=%d seed=%llu cycles=%ld cpu=%.3f start=%.12f best_cycle=%.12f final=%.12f "
                "score=%.12f other_n_improved=%d other_n_gain=%.6e\n",
                run.n, run.seed, kick_no, cpu_used, initial.side, best_reached, best_side,
                best_side / std::sqrt(static_cast<double>(run.n)) - 1.0, other_improved, other_gain);

    if (trace) std::fclose(trace);
    return 0;
}
