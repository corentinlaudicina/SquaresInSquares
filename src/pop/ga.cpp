#include "pop/ga.hpp"

#include "legalize.hpp"
#include "pop/backend.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
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

    // Initial population: the start in the first ratchet box, and jiggled copies of it.
    double h = out.best.L * (1.0 - c.ratchet_shrink);
    Packing base = to_packing(out.best);
    rescale(base, h);
    std::vector<Individual> pop(static_cast<size_t>(std::max(1, c.population)));
    for (size_t i = 0; i < pop.size(); ++i) {
        pop[i].p = base;
        if (i > 0)
            for (double& v : pop[i].p.xyt) v += c.init_jiggle * normal(rng);
    }
    {
        std::vector<int> all(pop.size());
        std::iota(all.begin(), all.end(), 0);
        relax_all(pop, all, rough_params(c), be);
        relax_all(pop, all, fine_params(c), be);
    }
    std::stable_sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b) { return a.info.cost < b.info.cost; });

    std::array<long, kGaMoves> produced{}, relaxed{}, survived{};
    std::array<double, 6> tsplit{};  // offspring, prefilter, rough, fine, select, ratchet
    long generation = 0, ratchets = 0, legal_improvements = 0, distance_evals = 0, since_ratchet = 0, longest_stall = 0;
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
        auto t = Clock::now();

        // 1. offspring
        const int lambda = std::max(1, c.offspring);
        std::vector<Individual> kids(static_cast<size_t>(lambda));
        for (Individual& k : kids) {
            const int a = uniform_int(rng, 0, static_cast<int>(pop.size()) - 1);
            GaMove m = draw_move(c.moves, rng);
            const Packing* mate = nullptr;
            if (m == GaMove::Crossover) {
                if (pop.size() < 2) {
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
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            return kids[static_cast<size_t>(a)].info.cost < kids[static_cast<size_t>(b)].info.cost;
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
        std::stable_sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b) { return a.info.cost < b.info.cost; });
        for (int m = 0; m < kGaMoves; ++m) survived[static_cast<size_t>(m)] += surv_gen[static_cast<size_t>(m)];
        lap(t, 4);

        // 5. ratchet
        bool ratcheted = false;
        if (feasible(pop.front().info, c.tol)) {
            System legal = to_system(pop.front().p);
            const LegalizeResult lr = legalize(legal);
            if (lr.ok && legal.L < out.best_side) {
                out.best = legal;
                out.best_side = legal.L;
                ++legal_improvements;
                side_trace.push_back({wall(), out.best_side});
            }
            h = std::min(h, out.best_side) * (1.0 - c.ratchet_shrink);
            for (Individual& x : pop) rescale(x.p, h);
            std::vector<int> all(pop.size());
            std::iota(all.begin(), all.end(), 0);
            relax_all(pop, all, fine_params(c), be);
            std::stable_sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b) { return a.info.cost < b.info.cost; });
            ++ratchets;
            ratcheted = true;
            since_ratchet = 0;
        } else {
            longest_stall = std::max(longest_stall, ++since_ratchet);
        }
        lap(t, 5);

        if (generation < 300 || generation % 10 == 0) {
            json sg = json::object();
            for (int m = 0; m < kGaMoves; ++m)
                if (surv_gen[static_cast<size_t>(m)]) sg[ga_move_name(static_cast<GaMove>(m))] = surv_gen[static_cast<size_t>(m)];
            gens.push_back({{"g", generation}, {"t", std::round(wall() * 100) / 100}, {"h", h},
                            {"best_cost", pop.front().info.cost}, {"best_depth", pop.front().info.max_depth},
                            {"ratchet", ratcheted}, {"survivors", sg}});
        }
    }

    json moves = json::object();
    for (int m = 0; m < kGaMoves; ++m) {
        const size_t k = static_cast<size_t>(m);
        moves[ga_move_name(static_cast<GaMove>(m))] = {
            {"produced", produced[k]}, {"relaxed", relaxed[k]}, {"survived", survived[k]},
            {"survival_rate", produced[k] ? static_cast<double>(survived[k]) / static_cast<double>(produced[k]) : 0.0}};
    }
    tel["ga"] = {{"generations", generation}, {"ratchets", ratchets}, {"legal_improvements", legal_improvements},
                 {"h_final", h}, {"longest_stall", longest_stall}, {"threads", be.threads()},
                 {"distance_evals", distance_evals}, {"moves", moves},
                 {"seconds", {{"offspring", tsplit[0]}, {"prefilter", tsplit[1]}, {"rough", tsplit[2]},
                              {"fine", tsplit[3]}, {"select", tsplit[4]}, {"ratchet", tsplit[5]}}},
                 {"per_generation", gens}};
    tel["side"]["trace"] = side_trace;
    tel["side"]["best"] = out.best_side;
    tel["budget"] = {{"wall_s", wall()}, {"kicks", generation}, {"stop", stop}};
    return out;
}

}  // namespace pack
