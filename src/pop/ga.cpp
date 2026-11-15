#include "pop/ga.hpp"

#include "legalize.hpp"
#include "pop/backend.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>

namespace pack {
namespace {

using Clock = std::chrono::steady_clock;
using json = nlohmann::json;

struct Individual {
    Packing p;
    CostInfo info;
    int move = -1;  // the move that made it; -1 = carried over from an earlier generation
};

CappedParams rough_params(const GAConfig& c) {
    CappedParams p;
    p.penalty = c.penalty == 0 ? Penalty::Area2 : Penalty::Depth2;
    p.history = c.history;
    p.tol = c.tol;
    p.stages = {{c.rough_cap, static_cast<int>(c.rough_iters)}};
    return p;
}

CappedParams fine_params(const GAConfig& c) {
    CappedParams p = rough_params(c);
    p.stages.clear();
    double cap = c.fine_cap;
    for (long k = 0; k < c.fine_stages; ++k, cap *= 0.25) p.stages.push_back({cap, static_cast<int>(c.fine_iters)});
    return p;
}

bool feasible(const CostInfo& c, double tol) { return c.max_depth <= tol && c.max_wall <= tol; }

void rescale(Packing& p, double h) {
    const double f = h / p.h;
    for (size_t i = 0; i < p.xyt.size(); i += 3) {
        p.xyt[i] *= f;
        p.xyt[i + 1] *= f;
    }
    p.h = h;
}

// Relax a set of individuals in place with one backend call.
void relax_all(std::vector<Individual>& v, const std::vector<int>& which, const CappedParams& prm, RelaxBackend& be) {
    std::vector<Packing> batch;
    batch.reserve(which.size());
    for (const int i : which) batch.push_back(std::move(v[static_cast<size_t>(i)].p));
    std::vector<CappedStats> st;
    be.relax_batch(batch, prm, st);
    for (size_t k = 0; k < which.size(); ++k) {
        v[static_cast<size_t>(which[k])].p = std::move(batch[k]);
        v[static_cast<size_t>(which[k])].info = st[k].info;
    }
}

std::vector<int> by_cost(const std::vector<Individual>& v) {
    std::vector<int> o(v.size());
    std::iota(o.begin(), o.end(), 0);
    std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return v[static_cast<size_t>(a)].info.cost < v[static_cast<size_t>(b)].info.cost; });
    return o;
}

struct Island {
    std::vector<Individual> pop;
    double h = 0.0;
    long ratchets = 0, resets = 0, since_progress = 0, longest_stall = 0;
    double mark = std::numeric_limits<double>::infinity();  // best residual at the last progress
};

void sort_pop(std::vector<Individual>& pop) {
    std::stable_sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b) { return a.info.cost < b.info.cost; });
}

// Population in box h: slot k is source k mod S (rescaled to h), jiggled from slot S on, relaxed
// rough then fine. One source: the source plus jiggled copies of it.
void seed_island(Island& is, const std::vector<const System*>& from, double h, const GAConfig& c, std::mt19937_64& rng,
                 RelaxBackend& be) {
    std::vector<Packing> base;
    for (const System* f : from) {
        base.push_back(to_packing(*f));
        rescale(base.back(), h);
    }
    is.h = h;
    is.pop.assign(static_cast<size_t>(std::max(1, c.population)), Individual{});
    for (size_t i = 0; i < is.pop.size(); ++i) {
        is.pop[i].p = base[i % base.size()];
        if (i >= base.size())
            for (double& v : is.pop[i].p.xyt) v += c.init_jiggle * normal(rng);
    }
    std::vector<int> all(is.pop.size());
    std::iota(all.begin(), all.end(), 0);
    relax_all(is.pop, all, rough_params(c), be);
    relax_all(is.pop, all, fine_params(c), be);
    sort_pop(is.pop);
    is.since_progress = 0;
    is.mark = is.pop.front().info.cost;
}

std::vector<int> neighbours(int i, int k, int topology) {
    std::vector<int> nb;
    if (k < 2) return nb;
    if (topology == 1) {
        for (int j = 0; j < k; ++j)
            if (j != i) nb.push_back(j);
    } else {
        nb.push_back((i + k - 1) % k);
        if (k > 2) nb.push_back((i + 1) % k);
    }
    return nb;
}

}  // namespace

GAOutcome run_ga(const System& start, const GAConfig& cfg0, unsigned long long seed, double wall_seconds,
                 long max_generations, const std::function<void(double, GAConfig&)>& update) {
    const auto t0 = Clock::now();
    const auto wall = [&] { return std::chrono::duration<double>(Clock::now() - t0).count(); };
    GAOutcome out;
    json& tel = out.telemetry;
    std::mt19937_64 rng(seed);
    GAConfig c = cfg0;
    CpuBackend be(c.threads);

    // The best legal packing so far, starting with the legalized start.
    out.best = start;
    const LegalizeResult lr0 = legalize(out.best);
    out.best_side = lr0.ok ? out.best.L : std::numeric_limits<double>::infinity();
    tel["side"]["start_legal"] = out.best_side;
    json side_trace = json::array({json::array({wall(), out.best_side})});

    // Sources: the start, then the seeding candidates (item 7), each legal with its own side.
    const std::vector<SeedCandidate> cands = seed_candidates(static_cast<int>(start.particles.size()), c.seeding, rng);
    std::vector<const System*> sources = {&out.best};
    std::vector<std::string> source_names = {"start"};
    json seeds = json::array();
    const System start_legal = out.best;  // sources[0] must not move when out.best improves below
    sources[0] = &start_legal;
    for (const SeedCandidate& sc : cands) {
        sources.push_back(&sc.sys);
        source_names.push_back(sc.source);
        seeds.push_back({{"source", sc.source}, {"side", sc.sys.L}});
        if (sc.sys.L < out.best_side) {
            out.best = sc.sys;
            out.best_side = sc.sys.L;
            side_trace.push_back({wall(), out.best_side});
        }
    }

    const int K = std::max(1, c.islands);
    std::vector<Island> isl(static_cast<size_t>(K));
    std::vector<std::string> island_source(static_cast<size_t>(K), "start");
    for (int i = 0; i < K; ++i) {
        if (K > 1) {
            // Box just under the source's legal side, but at most 5% above the best: ratcheting 0.02%
            // a generation from a much worse lattice would spend the budget getting back here.
            const size_t src = static_cast<size_t>(i) % sources.size();
            island_source[static_cast<size_t>(i)] = source_names[src];
            seed_island(isl[static_cast<size_t>(i)], {sources[src]}, std::min(sources[src]->L, out.best_side * 1.05) * (1.0 - c.ratchet_shrink), c, rng, be);
        } else {
            seed_island(isl[static_cast<size_t>(i)], sources, out.best_side * (1.0 - c.ratchet_shrink), c, rng, be);
        }
    }

    std::array<long, kGaMoves> produced{}, relaxed{}, survived{};
    std::array<double, 6> tsplit{};  // offspring, prefilter, rough, fine, select, ratchet
    long generation = 0, legal_improvements = 0, distance_evals = 0, migrations = 0;
    json gens = json::array();
    std::string stop = "budget";

    const auto more = [&] {
        if (max_generations > 0 && generation >= max_generations) return false;
        if (wall_seconds > 0.0 && wall() >= wall_seconds) return false;
        return max_generations > 0 || wall_seconds > 0.0;
    };
    const auto lap = [&](Clock::time_point& t, int slot) {
        const auto now = Clock::now();
        tsplit[static_cast<size_t>(slot)] += std::chrono::duration<double>(now - t).count();
        t = now;
    };

    for (; more(); ++generation) {
        if (update) {
            double frac = 0.0;
            if (max_generations > 0) frac = static_cast<double>(generation) / static_cast<double>(max_generations);
            if (wall_seconds > 0.0) frac = std::max(frac, wall() / wall_seconds);
            update(std::min(1.0, frac), c);
        }
        // Champions at the start of the generation: what the islands mate with.
        std::vector<Packing> champ;
        if (K > 1)
            for (const Island& is : isl) champ.push_back(is.pop.front().p);

        for (int ii = 0; ii < K; ++ii) {
            Island& is = isl[static_cast<size_t>(ii)];
            std::vector<Individual>& pop = is.pop;
            const std::vector<int> nb = neighbours(ii, K, c.topology);
            auto t = Clock::now();

            // 1. offspring
            const int lambda = std::max(1, c.offspring);
            std::vector<Individual> kids(static_cast<size_t>(lambda));
            Packing mate_buf;
            for (Individual& k : kids) {
                const int a = uniform_int(rng, 0, static_cast<int>(pop.size()) - 1);
                GaMove m = draw_move(c.moves, rng);
                const Packing* mate = nullptr;
                if (m == GaMove::Crossover) {
                    if (!nb.empty()) {  // islands: a connected island's champion, in this island's box
                        mate_buf = champ[static_cast<size_t>(nb[static_cast<size_t>(uniform_int(rng, 0, static_cast<int>(nb.size()) - 1))])];
                        rescale(mate_buf, is.h);
                        mate = &mate_buf;
                        ++migrations;
                    } else if (pop.size() < 2) {
                        m = GaMove::JiggleSmall;
                    } else {
                        int b = uniform_int(rng, 0, static_cast<int>(pop.size()) - 2);
                        if (b >= a) ++b;
                        mate = &pop[static_cast<size_t>(b)].p;
                    }
                }
                k.p = pop[static_cast<size_t>(a)].p;
                k.move = static_cast<int>(m);
                apply_move(k.p, m, c.moves, rng, mate);
                ++produced[static_cast<size_t>(m)];
            }
            lap(t, 0);

            // 2. pre-filter on the raw penalty
            {
                std::vector<Packing> batch;
                batch.reserve(kids.size());
                for (const Individual& k : kids) batch.push_back(k.p);
                std::vector<CostInfo> raw;
                be.cost_batch(batch, rough_params(c).penalty, raw);
                for (size_t i = 0; i < kids.size(); ++i) kids[i].info = raw[i];
            }
            std::vector<int> order = by_cost(kids);
            order.resize(static_cast<size_t>(std::clamp<long>(std::lround(c.prefilter * lambda), 1, lambda)));
            lap(t, 1);

            // 3. relax: rough for the pre-filtered, fine for the best of those
            relax_all(kids, order, rough_params(c), be);
            for (const int i : order) ++relaxed[static_cast<size_t>(kids[static_cast<size_t>(i)].move)];
            lap(t, 2);
            std::stable_sort(order.begin(), order.end(), [&](int x, int y) {
                return kids[static_cast<size_t>(x)].info.cost < kids[static_cast<size_t>(y)].info.cost;
            });
            order.resize(static_cast<size_t>(std::clamp<long>(std::lround(c.fine_frac * static_cast<double>(order.size())), 1,
                                                              static_cast<long>(order.size()))));
            relax_all(kids, order, fine_params(c), be);
            lap(t, 3);

            // 4. selection from parents + fine-relaxed offspring
            std::vector<Individual> cand = std::move(pop);
            for (Individual& p : cand) p.move = -1;
            for (const int i : order) cand.push_back(std::move(kids[static_cast<size_t>(i)]));
            const std::vector<int> corder = by_cost(cand);
            std::vector<Packing> cp;
            cp.reserve(cand.size());
            for (const Individual& x : cand) cp.push_back(x.p);
            ParallelFor parfor;
            if (be.threads() > 1) parfor = [&be](int n, const std::function<void(int, int)>& f) { be.parallel_for(n, f); };
            const std::vector<int> chosen = select_diverse(cp, corder, c.population, c.diversity, &distance_evals, parfor);
            pop.clear();
            std::array<long, kGaMoves> surv_gen{};
            for (const int i : chosen) {
                if (cand[static_cast<size_t>(i)].move >= 0) ++surv_gen[static_cast<size_t>(cand[static_cast<size_t>(i)].move)];
                pop.push_back(std::move(cand[static_cast<size_t>(i)]));
            }
            sort_pop(pop);
            for (int m = 0; m < kGaMoves; ++m) survived[static_cast<size_t>(m)] += surv_gen[static_cast<size_t>(m)];
            lap(t, 4);

            // 5. ratchet, or count the stall; a long stall resets the island
            bool ratcheted = false, reset = false;
            if (feasible(pop.front().info, c.tol)) {
                System legal = to_system(pop.front().p);
                const LegalizeResult lr = legalize(legal);
                if (lr.ok && legal.L < out.best_side) {
                    out.best = legal;
                    out.best_side = legal.L;
                    ++legal_improvements;
                    side_trace.push_back({wall(), out.best_side});
                }
                is.h = std::min(is.h, out.best_side) * (1.0 - c.ratchet_shrink);
                for (Individual& x : pop) rescale(x.p, is.h);
                std::vector<int> all(pop.size());
                std::iota(all.begin(), all.end(), 0);
                relax_all(pop, all, fine_params(c), be);
                sort_pop(pop);
                ++is.ratchets;
                ratcheted = true;
                is.since_progress = 0;
                is.mark = pop.front().info.cost;
            } else if (pop.front().info.cost < is.mark * (1.0 - c.island_threshold)) {
                is.mark = pop.front().info.cost;
                is.since_progress = 0;
            } else {
                is.longest_stall = std::max(is.longest_stall, ++is.since_progress);
                if (c.island_patience > 0 && is.since_progress >= c.island_patience) {
                    seed_island(is, {&out.best}, out.best_side * (1.0 - c.ratchet_shrink), c, rng, be);
                    ++is.resets;
                    reset = true;
                }
            }
            lap(t, 5);

            if (gens.size() < 600 && (generation < 300 || generation % 10 == 0)) {
                json sg = json::object();
                for (int m = 0; m < kGaMoves; ++m)
                    if (surv_gen[static_cast<size_t>(m)]) sg[ga_move_name(static_cast<GaMove>(m))] = surv_gen[static_cast<size_t>(m)];
                json e = {{"g", generation}, {"t", std::round(wall() * 100) / 100}, {"h", is.h},
                          {"best_cost", pop.front().info.cost}, {"best_depth", pop.front().info.max_depth},
                          {"ratchet", ratcheted}, {"survivors", sg}};
                if (K > 1) {
                    e["island"] = ii;
                    e["reset"] = reset;
                }
                gens.push_back(e);
            }
        }
    }

    json moves = json::object();
    for (int m = 0; m < kGaMoves; ++m) {
        const size_t k = static_cast<size_t>(m);
        moves[ga_move_name(static_cast<GaMove>(m))] = {
            {"produced", produced[k]}, {"relaxed", relaxed[k]}, {"survived", survived[k]},
            {"survival_rate", produced[k] ? static_cast<double>(survived[k]) / static_cast<double>(produced[k]) : 0.0}};
    }
    long ratchets = 0, resets = 0, longest_stall = 0;
    double h_final = isl.front().h;
    json islands = json::array();
    for (const Island& is : isl) {
        ratchets += is.ratchets;
        resets += is.resets;
        longest_stall = std::max(longest_stall, is.longest_stall);
        h_final = std::min(h_final, is.h);
        islands.push_back({{"source", island_source[static_cast<size_t>(&is - isl.data())]}, {"h", is.h}, {"ratchets", is.ratchets}, {"resets", is.resets},
                           {"best_cost", is.pop.front().info.cost}, {"longest_stall", is.longest_stall}});
    }
    tel["ga"] = {{"generations", generation}, {"ratchets", ratchets}, {"legal_improvements", legal_improvements},
                 {"h_final", h_final}, {"longest_stall", longest_stall}, {"threads", be.threads()},
                 {"distance_evals", distance_evals}, {"moves", moves},
                 {"seconds", {{"offspring", tsplit[0]}, {"prefilter", tsplit[1]}, {"rough", tsplit[2]},
                              {"fine", tsplit[3]}, {"select", tsplit[4]}, {"ratchet", tsplit[5]}}},
                 {"per_generation", gens}};
    if (K > 1) {
        tel["ga"]["islands"] = islands;
        tel["ga"]["resets"] = resets;
        tel["ga"]["migrations"] = migrations;
        tel["ga"]["topology"] = c.topology == 1 ? "full" : "ring";
    }
    if (!seeds.empty()) tel["ga"]["seeds"] = seeds;
    tel["side"]["trace"] = side_trace;
    tel["side"]["best"] = out.best_side;
    tel["budget"] = {{"wall_s", wall()}, {"kicks", generation}, {"stop", stop}};
    return out;
}

}  // namespace pack
