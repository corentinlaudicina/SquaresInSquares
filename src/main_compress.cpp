// Repeats one cycle from a known configuration: kick it, relax it with no walls, compress it back
// down, and keep the result or not according to its side.
//
// Changes from the earlier driver:
//   * the cycle loop is written out, with a CPU-time budget as well as a cycle count;
//   * the side temperature is an option (--T), and is relative to L / sqrt(n);
//   * the large-move probability escalates with the number of cycles since the last improvement,
//     and a long stagnation sends the walker back to `best`;
//   * the per-n pool is fed from every new best and pulled from when it holds something better;
//   * per-kick statistics (tries, accepts, improvements, CPU seconds, gain per CPU second) are
//     kept here and printed at the end, so the kick mix can be tuned from data.

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
constexpr long kCycles = 20;

// The kick.
constexpr double kLargeMoveProb = 0.1;  // base chance of a large move instead of a shake
constexpr double kLargeMoveMax = 0.6;   // ceiling the escalation can reach
constexpr long kEscalateEvery = 25;     // stagnant cycles per +0.1 on the large-move chance
constexpr long kRestartAfter = 200;     // stagnant cycles before going back to `best`
constexpr long kShakeSweeps = 100;      // Monte Carlo sweeps in one shake
constexpr double kShakeT = 0.05;        // temperature of the shake
constexpr double kShakeStep = 0.1;      // largest translation per move
constexpr double kShakeTurn = 0.1;      // largest rotation per move, radians

// The compression.
constexpr double kDs = 0.005;    // box shrink per step, as a length
constexpr long kMaxIter = 5000;  // minimisations allowed per compression

// Temperature on the side for accepting a cycle's result, in units of L / sqrt(n). 0 accepts only improvements.
constexpr double kSideT = 0.01;

constexpr double kTraceSeconds = 0.2;  // least time between two frames of the same phase
constexpr long kPropagateEvery = 1;    // propagate the pool every this many new bests

static_assert(kN >= 1, "kN must be at least 1");

struct RunOptions {
    int n = kN;
    unsigned long long seed = kSeed;
    long cycles = kCycles;
    double cpu_budget = 0.0;  // process CPU seconds; above zero it replaces `cycles`
    bool save = true;         // write the record and the trajectory
    std::string csv;          // the per-move table; empty picks logs/moves_n<n>.csv
    double symmetrize = 0.0;  // LnmOptions::w_symmetrize; off unless asked for
    bool propagate = false;   // use and feed the per-n pool
    double side_t = kSideT;   // relative temperature on the side
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
        else if (a == "--cycles") o.cycles = static_cast<long>(option_value(argc, argv, k));
        else if (a == "--cpu") o.cpu_budget = option_value(argc, argv, k);
        else if (a == "--no-save") o.save = false;
        else if (a == "--symmetrize") o.symmetrize = option_value(argc, argv, k);
        else if (a == "--propagate") o.propagate = true;
        else if (a == "--T") o.side_t = option_value(argc, argv, k);
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

// Returns the row the whole cycle is credited to. `large_prob` is the chance of a large move.
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
    LnmStats lnm_stats;  // outlives the cycles, so the symmetrize bandit can learn across them

    System start = load_config(run.n);
    std::FILE* trace = run.save ? open_trace(run.n, start) : nullptr;
    long cycle = 0;

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
            append_frame(trace, cycle, s, label);
        };
        watch(start, "initial");
    }
    lnm.watch = watch;

    // --- the starting point: mirror, relax, compress -----------------------------------------
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

    System current = start, best = start;
    double current_side = start.L;
    double best_side = start.L;

    const double cpu_start = process_cpu_seconds();

    // The pool works on its own random stream, so that turning it on does not reshuffle the cycles of the run it is compared against until it actually changes something.
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
    long stagnant = 0;  // cycles since `best` last improved
    double best_cycle = std::numeric_limits<double>::infinity();  // smallest side any cycle reached

    // --- the cycle loop ----------------------------------------------------------------------
    const auto more = [&] {
        if (run.cpu_budget > 0.0) return process_cpu_seconds() - cpu_start < run.cpu_budget;
        return cycle < run.cycles;
    };

    for (; more(); ++cycle) {
        const double t0 = process_cpu_seconds();

        // Escalate the kick the longer nothing has improved; a fixed probability is too timid when stuck and too violent when making progress.
        const double large_prob =
            std::min(kLargeMoveMax, kLargeMoveProb + 0.1 * static_cast<double>(stagnant / kEscalateEvery));

        // 1. kick
        System trial = current;
        const MoveId row = kick(trial, rng, lnm, lnm_stats, log, large_prob);
        if (watch) watch(trial, "kicked");

        // 2. relax with no walls, then refit the box around the result
        RelaxOptions ro;
        ro.watch = watch;
        ro.watch_label = "relax free";
        relax_free(trial, ro);

        // 3. compress back down
        const CompressResult res = anneal(trial, kDs, kMaxIter, watch);

        // 4. accept on the side, Metropolis with a temperature relative to L / sqrt(n)
        bool accepted = false, improved = false;
        double gain = 0.0;
        if (res.found) {
            best_cycle = std::min(best_cycle, res.side);
            const double scale = res.side / std::sqrt(static_cast<double>(run.n));
            const double T = run.side_t * scale;
            const double d = res.side - current_side;
            if (d <= 0.0 || (T > 0.0 && unit(rng) < std::exp(-d / T))) {
                accepted = true;
                current = trial;
                current_side = res.side;
            }
            if (res.side < best_side) {
                improved = true;
                gain = best_side - res.side;
                best = trial;
                best_side = res.side;
                current = trial;  // a new best is always where the walker stands
                current_side = res.side;
                accepted = true;
            }
        }

        // 5. bookkeeping
        const double cycle_cpu = process_cpu_seconds() - t0;
        log.record(row, accepted, gain, cycle_cpu);
        KickTally& k = tally[row];
        ++k.tries;
        k.cpu += cycle_cpu;
        if (accepted) ++k.accepted;
        if (improved) {
            ++k.improved;
            k.gain += gain;
        }

        if (improved) {
            stagnant = 0;
            std::printf("cycle %ld  %-10s  new best side %.9f  (-%.3e)\n", cycle, row_name(row),
                        best_side, gain);
            if (watch) watch(best, "best");

            if (run.propagate) {
                pool.offer(best);
                if (++new_bests % kPropagateEvery == 0) {
                    add_gains(propagate(pool, run.n, pull, pool_rng, &log));
                }
            }
        } else {
            ++stagnant;
        }

        // 6. recovery: pull from the pool if it holds something better; else go back to `best`.
        // Compared with `best`, not `current`: the pool always holds `best`, so comparing with
        // `current` would undo every uphill step the Metropolis test had just accepted.
        if (run.propagate) {
            const PoolEntry* e = pool.at(run.n);
            if (e && e->side < best_side - 1e-12) {
                current = best = e->packing;
                current_side = best_side = e->side;
                stagnant = 0;
            }
        }
        if (stagnant >= kRestartAfter) {
            current = best;
            current_side = best_side;
            stagnant = 0;
        }
    }

    // --- report ------------------------------------------------------------------------------
    const double cpu_used = process_cpu_seconds() - cpu_start;
    std::printf("\nn %d: best side %.9f after %ld cycles, %.1f CPU s\n", run.n, best_side, cycle,
                cpu_used);
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

    // The one line scripts/benchmark_moves.py reads.
    std::printf("result n=%d seed=%llu cycles=%ld cpu=%.3f start=%.12f best_cycle=%.12f final=%.12f "
                "score=%.12f other_n_improved=%d other_n_gain=%.6e\n",
                run.n, run.seed, cycle, cpu_used, initial.side, best_cycle, best_side,
                best_side / std::sqrt(static_cast<double>(run.n)) - 1.0, other_improved, other_gain);

    if (trace) std::fclose(trace);
    return 0;
}