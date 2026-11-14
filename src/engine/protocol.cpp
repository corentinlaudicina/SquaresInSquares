#include "engine/protocol.hpp"

#include "moves/mc.hpp"
#include "moves/relax.hpp"
#include "moves/symmetry.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <map>
#include <random>

namespace pack {
namespace {

using Clock = std::chrono::steady_clock;
using json = nlohmann::json;

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

// Same draw order as the drivers' kick(): one uniform, then either a large move or the sweeps.
MoveId kick(System& sys, std::mt19937_64& rng, const KickConfig& k, const LnmOptions& lnm,
            LnmStats& stats, MoveLog& log, double large_prob) {
    if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) < large_prob) {
        stats.last = MoveId::Count;
        large_neighbour_move(sys, rng, lnm, &stats);
        return cycle_row(stats.last);
    }
    StepSizes steps;
    steps.translate = k.shake_translate;
    steps.rotate = k.shake_rotate;
    for (long s = 0; s < k.shake_sweeps; ++s) sweep(sys, steps, k.shake_T, rng, &log);
    return MoveId::CycleShake;
}

struct KickTally {
    long tries = 0, accepted = 0, improved = 0;
    double cpu = 0.0, gain = 0.0;
};

double cpu_now() { return static_cast<double>(std::clock()) / CLOCKS_PER_SEC; }

double large_prob(const KickConfig& k, long stuck) {
    return std::min(k.large_max,
                    k.large_base + k.large_step * static_cast<double>(stuck / std::max(1L, k.escalate_every)));
}

}  // namespace

EngineConfig shrink_defaults() { return EngineConfig{}; }

EngineConfig cycle_defaults() {
    EngineConfig c;
    c.mode = "cycle";
    c.prepare_ds = 0.005;
    return c;
}

EngineResult run_engine(const System& seed, const EngineConfig& cfg, unsigned long long seed_value) {
    const auto t_start = Clock::now();
    const double cpu_start0 = cpu_now();
    const auto wall = [&] { return std::chrono::duration<double>(Clock::now() - t_start).count(); };
    const bool shrink = cfg.mode == "shrink";

    EngineResult out;
    json& tel = out.telemetry;
    tel["run"] = {{"n", seed.particles.size()}, {"seed", seed_value}, {"mode", cfg.mode}};

    std::mt19937_64 rng(seed_value);
    MoveLog log;
    LnmOptions lnm = cfg.kick.lnm;
    lnm.log = &log;
    lnm.watch = {};
    LnmStats lnm_stats;

    // --- start: mirror, relax, compress ----------------------------------------------------------
    System start = seed;
    const double seed_side = start.tight_side();
    if (cfg.mirror) mirror_symmetry(start);
    relax(start, RelaxOptions{});
    const CompressResult initial = anneal(start, cfg.prepare_ds, cfg.prepare_max_relax);
    tel["side"]["seed"] = seed_side;
    if (!initial.found) {
        tel["budget"] = {{"stop", "start_never_certified"}, {"wall_s", wall()}};
        out.best = seed;
        out.best_side = seed_side;
        return out;
    }
    out.found = true;
    tel["side"]["start"] = initial.side;

    // Kept as the drivers have it: start.L, the box, not initial.side (audit B3).
    System best = start, current = start;
    double best_side = start.L, current_side = start.L;
    double best_reached = shrink ? best_side : std::numeric_limits<double>::infinity();

    Pool pool;
    PropagateOptions pull;
    std::mt19937_64 pool_rng(seed_value ^ 0x9e3779b97f4a7c15ULL);
    std::map<int, double> pool_gains;
    long new_bests = 0;
    if (cfg.propagate) pool.offer(best);

    std::map<MoveId, KickTally> tally;
    std::map<std::string, long> relax_stops;
    long relax_calls = 0, relax_iters = 0, compress_calls = 0, compress_relax = 0, compress_found = 0;
    long escalations = 0, restarts = 0;
    json trace = json::array();
    double t_best = 0.0;
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    double eps = cfg.eps;
    long stuck = 0;  // shrink: failed kicks in a row; cycle: cycles since the best improved
    long kick_no = 0;
    std::string stop = "budget";

    const double cpu_start = cpu_now();
    const auto more = [&] {
        if (shrink && eps < cfg.eps_min) { stop = "converged"; return false; }
        if (cfg.max_kicks > 0 && kick_no >= cfg.max_kicks) return false;
        if (cfg.wall_seconds > 0.0 && wall() >= cfg.wall_seconds) return false;
        return cfg.max_kicks > 0 || cfg.wall_seconds > 0.0;
    };

    for (; more(); ++kick_no) {
        const double t0 = cpu_now();
        const double p_large = large_prob(cfg.kick, stuck);
        bool accepted = false, improved = false;
        double gain = 0.0;
        MoveId row;

        if (shrink) {
            const double target = best_side * (1.0 - eps);
            System trial = best;
            trial.rescale_box(target);
            row = kick(trial, rng, cfg.kick, lnm, lnm_stats, log, p_large);
            const RelaxResult rr = relax(trial, RelaxOptions{});
            ++relax_calls;
            relax_iters += rr.iterations;
            ++relax_stops[stop_name(rr.stop)];
            if (rr.feasible) {
                trial.fit_box();
                if (trial.is_feasible(kCertifyTol) && trial.L < best_side) {
                    accepted = improved = true;
                    gain = best_side - trial.L;
                    best = trial;
                    best_side = trial.L;
                    best_reached = std::min(best_reached, best_side);
                }
            }
        } else {
            System trial = current;
            row = kick(trial, rng, cfg.kick, lnm, lnm_stats, log, p_large);
            relax_free(trial, RelaxOptions{});
            const CompressResult res = anneal(trial, cfg.cycle_ds, cfg.cycle_max_relax);
            ++compress_calls;
            compress_relax += res.relaxations;
            if (res.found) {
                ++compress_found;
                best_reached = std::min(best_reached, res.side);
                const double T = cfg.side_T * res.side / std::sqrt(static_cast<double>(seed.particles.size()));
                const double d = res.side - current_side;
                if (d <= 0.0 || (T > 0.0 && unit(rng) < std::exp(-d / T))) {
                    accepted = true;
                    current = trial;
                    current_side = res.side;
                }
                if (res.side < best_side) {
                    improved = accepted = true;
                    gain = best_side - res.side;
                    best = trial;
                    best_side = res.side;
                    current = trial;
                    current_side = res.side;
                }
            }
        }

        const double kick_cpu = cpu_now() - t0;
        log.record(row, accepted, gain, kick_cpu);
        KickTally& k = tally[row];
        ++k.tries;
        k.cpu += kick_cpu;
        if (accepted) ++k.accepted;
        if (improved) {
            ++k.improved;
            k.gain += gain;
            t_best = wall();
            if (trace.size() < 200) trace.push_back({t_best, best_side});
        }

        if (shrink) {
            if (improved) {
                eps = std::min(cfg.eps_max, eps * cfg.eps_up);
                stuck = 0;
            } else {
                eps *= cfg.eps_down;
                ++stuck;
            }
        } else {
            stuck = improved ? 0 : stuck + 1;
        }
        if (stuck > 0 && stuck % std::max(1L, cfg.kick.escalate_every) == 0 &&
            large_prob(cfg.kick, stuck) > large_prob(cfg.kick, stuck - 1)) {
            ++escalations;
        }

        if (improved && cfg.propagate) {
            pool.offer(best);
            if (++new_bests % std::max(1L, cfg.propagate_every) == 0) {
                for (const auto& [m, g] : propagate(pool, static_cast<int>(seed.particles.size()), pull,
                                                    pool_rng, &log).improved) {
                    pool_gains[m] += g;
                }
            }
        }
        if (cfg.propagate) {
            const PoolEntry* e = pool.at(static_cast<int>(seed.particles.size()));
            if (e && e->side < best_side - 1e-12) {
                best = e->packing;
                best_side = e->side;
                best_reached = std::min(best_reached, best_side);
                if (!shrink) {
                    current = best;
                    current_side = best_side;
                }
                stuck = 0;
            }
        }
        if (!shrink && stuck >= cfg.restart_after) {
            current = best;
            current_side = best_side;
            stuck = 0;
            ++restarts;
        }
    }

    out.best = best;
    out.best_side = best_side;

    tel["budget"] = {{"wall_s", wall()}, {"cpu_s", cpu_now() - cpu_start0}, {"loop_cpu_s", cpu_now() - cpu_start},
                     {"kicks", kick_no}, {"stop", stop}};
    tel["side"]["best"] = best_side;
    tel["side"]["best_reached"] = best_reached;
    tel["side"]["t_best_s"] = t_best;
    tel["side"]["trace"] = trace;
    tel["final_state"] = {{"energy", best.energy()}, {"max_overlap_area", best.max_overlap()}};
    json kicks = json::object();
    for (const auto& [id, k] : tally) {
        kicks[row_name(id)] = {{"tries", k.tries}, {"acc", k.accepted}, {"impr", k.improved},
                               {"cpu_s", k.cpu}, {"gain", k.gain}};
    }
    tel["kicks"] = kicks;
    json moves = json::object();
    for (int i = 0; i < static_cast<int>(MoveId::Count); ++i) {
        const MoveCounters& c = log.at(static_cast<MoveId>(i));
        if (c.attempts == 0 || i >= static_cast<int>(MoveId::CycleShake)) continue;
        moves[move_name(static_cast<MoveId>(i))] = {{"tries", c.attempts}, {"acc", c.accepted},
                                                    {"cpu_s", c.cpu_seconds}};
    }
    tel["moves"] = moves;
    if (shrink) {
        tel["relax"] = {{"calls", relax_calls}, {"iters", relax_iters}, {"stops", relax_stops}};
    } else {
        tel["compress"] = {{"calls", compress_calls}, {"relaxations", compress_relax}, {"found", compress_found}};
    }
    json pg = json::object();
    for (const auto& [m, g] : pool_gains) pg[std::to_string(m)] = g;
    tel["events"] = {{"escalations", escalations}, {"restarts", restarts}, {"jams", relax_calls - relax_stops["feasible"]},
                     {"eps_final", eps}, {"pool_gains", pg}};
    return out;
}

}  // namespace pack
