// Symmetric genotypes: layouts, gradients through the expansion, round trips, crossover, relaxation,
// and a C4 GA run that is legal and thread independent.

#include "io/io.hpp"
#include "init.hpp"
#include "legalize.hpp"
#include "pop/backend.hpp"
#include "pop/diversity.hpp"
#include "pop/ga.hpp"
#include "pop/symmetric.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>

using namespace pack;

namespace {

// A legal start built here, so the tests need no data files (the seed packings are not in git).
System grid_start(int n) {
    System s = grid_init(make_square(1.0), n);
    legalize(s);
    return s;
}

Packing random_gen(const SymLayout& L, double h, unsigned seed) {
    std::mt19937_64 rng(seed);
    Packing g;
    g.h = h;
    for (int i = 0; i < L.m; ++i) g.xyt.insert(g.xyt.end(), {uniform(rng, 0.3, h - 0.3), uniform(rng, 0.3, h - 0.3), uniform(rng, 0.0, kPi)});
    fix_centre(g, L);
    return g;
}

void test_layouts() {
    assert(make_layout(Sym::C2, 7).m == 4 && make_layout(Sym::C2, 7).centre);
    assert(make_layout(Sym::C2, 8).m == 4 && !make_layout(Sym::C2, 8).centre);
    assert(make_layout(Sym::C4, 17).m == 5 && make_layout(Sym::C4, 17).centre);
    assert(make_layout(Sym::C4, 16).m == 4 && !make_layout(Sym::C4, 16).centre);
    assert(make_layout(Sym::C4, 18).sym == Sym::None && make_layout(Sym::C4, 19).sym == Sym::None);
    const SymLayout L = make_layout(Sym::C4, 17);
    const Packing full = expanded(random_gen(L, 5.0, 1), L);
    assert(full.n() == 17);
}

void test_gradient_through_expansion() {
    Workspace ws;
    int checked = 0;
    for (const auto [sym, n] : {std::pair{Sym::C2, 21}, {Sym::C2, 20}, {Sym::C4, 21}, {Sym::C4, 24}}) {
        const SymLayout L = make_layout(sym, n);
        Packing g = random_gen(L, 4.8, static_cast<unsigned>(n));
        std::vector<double> grad(g.xyt.size());
        cost_symmetric(g, L, Penalty::Depth2, ws, grad.data());
        for (size_t k = 0; k < g.xyt.size(); ++k) {
            const bool centre_pos = L.centre && k >= g.xyt.size() - 3 && k % 3 != 2;
            if (centre_pos) {
                assert(grad[k] == 0.0);  // the centre square's position is not a variable
                continue;
            }
            const double v = g.xyt[k], e = 1e-6;
            g.xyt[k] = v + e;
            const double a = cost_symmetric(g, L, Penalty::Depth2, ws, nullptr).cost;
            g.xyt[k] = v - e;
            const double b = cost_symmetric(g, L, Penalty::Depth2, ws, nullptr).cost;
            g.xyt[k] = v;
            const double fd = (a - b) / (2 * e), scale = std::max({1e-8, std::abs(fd), std::abs(grad[k])});
            if (std::abs(fd - grad[k]) > 1e-5 * scale) {
                std::printf("sym %d n %d var %zu: analytic %.10g fd %.10g\n", static_cast<int>(sym), n, k, grad[k], fd);
                assert(false);
            }
            ++checked;
        }
    }
    std::printf("  gradient through the expansion: %d generator components match finite differences\n", checked);
}

void test_round_trips() {
    for (const auto [sym, n] : {std::pair{Sym::C2, 13}, {Sym::C4, 16}, {Sym::C4, 17}}) {
        const SymLayout L = make_layout(sym, n);
        Packing g = random_gen(L, 6.0, 40u + static_cast<unsigned>(n));
        const Packing full = expanded(g, L);
        // canonicalizing generators does not change the packing they describe
        Packing c = g;
        canonicalize(c, L);
        assert(assign_distance(full, expanded(c, L)) < 1e-9);
        // and the generators can be recovered from the full packing
        const Packing back = symmetrize(full, L);
        assert(back.n() == L.m);
        assert(assign_distance(full, expanded(back, L)) < 1e-9);
        // crossover keeps the count, and the expansion is symmetric by construction
        std::mt19937_64 rng(5);
        Packing child = g;
        crossover_sym(child, random_gen(L, 6.0, 99), MoveParams{}, rng, L);
        assert(child.n() == L.m && expanded(child, L).n() == n);
    }
}

void test_relax_symmetric() {
    // 6x6 grid is C4-symmetric: jiggle its generators, relax, expect a strictly valid full packing
    const SymLayout L = make_layout(Sym::C4, 36);
    Packing full;
    full.h = 6.0001;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) full.xyt.insert(full.xyt.end(), {i + 0.5, j + 0.5, 0.0});
    Packing g = symmetrize(full, L);
    std::mt19937_64 rng(3);
    for (double& v : g.xyt) v += 0.05 * normal(rng);
    Workspace ws;
    const CappedStats st = relax_symmetric(g, L, CappedParams{}, ws);
    System s = to_system(expanded(g, L));
    const StrictReport before = strict_report(s);
    std::printf("  C4 6x6 from jiggled generators: reached=%d after %d iterations, depth %.1e, min gap %.1e\n",
                st.reached, st.iters, st.info.max_depth, before.min_pair_gap);
    // reached = every depth <= tol (1e-12), not strictly apart: legalization finishes it at no real cost
    const LegalizeResult lr = legalize(s);
    assert(st.reached && lr.ok && lr.side_after < 6.0001 + 1e-9);
}

void test_c4_ga() {
    const auto seed = std::optional<System>(grid_start(17));
    GAConfig cfg;
    cfg.population = 10;
    cfg.offspring = 30;
    cfg.symmetry = 2;
    cfg.threads = 1;
    const GAOutcome a = run_ga(*seed, cfg, 7, 0.0, 5);
    cfg.threads = 4;
    const GAOutcome b = run_ga(*seed, cfg, 7, 0.0, 5);
    assert(a.best_side == b.best_side && a.telemetry["ga"]["moves"] == b.telemetry["ga"]["moves"]);
    assert(a.telemetry["ga"]["symmetry"] == "c4" && a.telemetry["ga"]["generators"] == 5);
    assert(strict_report(a.best).ok() && a.best.particles.size() == 17);
    cfg.symmetry = 2;  // impossible for 18: falls back to none and says so
    const auto s18 = std::optional<System>(grid_start(18));
    const GAOutcome c = run_ga(*s18, cfg, 7, 0.0, 2);
    assert(c.telemetry["ga"]["symmetry"].get<std::string>().rfind("none", 0) == 0);
    std::printf("  C4 GA n=17: legal side %.9f, same on 1 and 4 threads; n=18 falls back to none\n", a.best_side);
}

}  // namespace

int main() {
    test_layouts();
    test_gradient_through_expansion();
    test_round_trips();
    test_relax_symmetric();
    test_c4_ga();
    std::printf("symmetric: all tests passed\n");
}
