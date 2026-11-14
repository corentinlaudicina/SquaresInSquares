// Repeats one cycle from a known configuration: kick it, relax it with no walls, compress it back down, and keep the result or not according to its side.

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

// Run parameters: edit these and rebuild.
constexpr int kN = 552;                    // loads data/init_external/n<kN>.cfg
constexpr unsigned long long kSeed = 1ULL;
constexpr long kCycles = 20;

// The kick.
constexpr double kLargeMoveProb = 0.1;      // chance of a large move instead of a shake
constexpr long kShakeSweeps = 100;           // Monte Carlo sweeps in one shake
constexpr double kShakeT = 0.05;             // temperature of the shake
constexpr double kShakeStep = 0.1;          // largest translation per move
constexpr double kShakeTurn = 0.1;          // largest rotation per move, radians

// The compression.
constexpr double kDs = 0.01;                // box shrink per step, as a length
constexpr long kMaxIter = 5000;              // minimisations allowed per compression

// Temperature on the side for accepting a cycle's result. 0 accepts only improvements.
constexpr double kSideT = 0.01;

constexpr double kTraceSeconds = 0.2;       // least time between two frames of the same phase

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

// Returns the row the whole cycle is credited to.
MoveId kick(System& sys, std::mt19937_64& rng, const LnmOptions& lnm, LnmStats& stats,
            MoveLog& log) {
    if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) < kLargeMoveProb) {
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
    const FrameWatch no_watch;

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
    double best_side = std::numeric_limits<double>::infinity();  // best result of any cycle
    double smallest = start.L;                                   // what a cycle has to beat

    const double cpu_start = process_cpu_seconds();

    // The pool works on its own random stream, so that turning it on does not reshuffle the cycles of the run it is compared against until it actually changes something.
    Pool pool;
    PropagateOptions pull;
    std::mt19937_64 pool_rng(run.seed ^ 0x9e3779b97f4a7c15ULL);
    std::map<int, double> pool_gains;
    const auto add_gains = [&](const PropagateReport& r) {
        for (const auto& [m, gain] : r.improved) pool_gains[m] += gain;
    };

    return 0;
}
