// GA layer: moves keep the genome well formed, D4 and permutation invariance of the distances,
// diversity selection, and a small GA run that is legal and identical on 1 and 4 threads.

#include "io/io.hpp"
#include "legalize.hpp"
#include "pop/backend.hpp"
#include "pop/ga.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>

using namespace pack;

namespace {

Packing grid(int k, double gap = 0.0) {
    Packing p;
    p.h = k * (1.0 + gap);
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j) p.xyt.insert(p.xyt.end(), {(i + 0.5) * (1 + gap), (j + 0.5) * (1 + gap), 0.0});
    return p;
}

Packing jiggled(Packing p, double sigma, unsigned seed) {
    std::mt19937_64 rng(seed);
    for (double& v : p.xyt) v += sigma * normal(rng);
    return p;
}

void test_moves_keep_shape() {
    std::mt19937_64 rng(1);
    const Packing a = jiggled(grid(5), 0.05, 2), b = jiggled(grid(5), 0.05, 3);
    MoveParams prm;
    for (int m = 0; m < kGaMoves; ++m) {
        for (int rep = 0; rep < 50; ++rep) {
            Packing c = a;
            apply_move(c, static_cast<GaMove>(m), prm, rng, &b);
            assert(c.n() == a.n() && c.h == a.h);
            for (const double v : c.xyt) assert(std::isfinite(v));
        }
    }
    // the move draw follows the weights: a zero weight is never drawn
    prm.weight = {0, 1, 0, 0, 0, 0, 0};
    for (int k = 0; k < 100; ++k) assert(draw_move(prm, rng) == GaMove::JiggleSmall);
}

void test_d4_keeps_packings_valid() {
    // a tilted, valid packing stays valid under all 8 symmetries
    Packing p = jiggled(grid(4, 0.05), 0.01, 7);
    Workspace ws;
    assert(cost_grad(p, Penalty::Depth2, ws, nullptr).cost == 0.0);
    for (int k = 0; k < 8; ++k) {
        Packing q = p;
        apply_d4(q, k);
        assert(cost_grad(q, Penalty::Depth2, ws, nullptr).cost == 0.0);
        assert(strict_report(to_system(q)).ok());
    }
}

void test_distances_are_invariant() {
    const Packing a = jiggled(grid(4, 0.05), 0.03, 11);
    Packing perm = a;  // reversed order: same packing
    for (int i = 0; i < a.n(); ++i)
        for (int k = 0; k < 3; ++k) perm.xyt[3 * static_cast<size_t>(i) + static_cast<size_t>(k)] = a.xyt[3 * static_cast<size_t>(a.n() - 1 - i) + static_cast<size_t>(k)];
    assert(assign_distance(a, perm) < 1e-12);
    for (int g = 0; g < 8; ++g) {
        Packing q = a;
        apply_d4(q, g);
        q.xyt[2] += kPi / 2;  // a quarter turn of one square is the same square
        assert(assign_distance(a, q) < 1e-7);
        assert(hist_distance(a, q) < 1e-9);
    }
    const Packing far = jiggled(grid(4, 0.05), 0.2, 12);
    assert(assign_distance(a, far) > 0.05 && hist_distance(a, far) > 0.0);
    std::printf("  distances: permutation / D4 / quarter-turn invariant; jiggled copy at %.3f (assign)\n",
                assign_distance(a, far));
}

void test_selection_spreads() {
    std::vector<Packing> c;
    const Packing a = jiggled(grid(4, 0.05), 0.03, 21), b = jiggled(grid(4, 0.05), 0.3, 22);
    for (int k = 0; k < 6; ++k) c.push_back(a);  // six copies of a, then b
    c.push_back(b);
    std::vector<int> order = {0, 1, 2, 3, 4, 5, 6};
    DiversityParams prm;
    prm.elite = 0.25;
    const std::vector<int> pick = select_diverse(c, order, 2, prm);
    assert(pick.size() == 2 && pick[0] == 0 && pick[1] == 6);  // elite, then the first far candidate
    const std::vector<int> fill = select_diverse(c, order, 4, prm);  // not enough far ones: fill by fitness
    assert(fill.size() == 4 && fill[1] == 6 && fill[2] == 1);
}

void test_small_ga_is_legal_and_thread_independent() {
    auto seed = load_record("data/init_external/n0010.cfg");
    assert(seed);
    GAConfig cfg;
    cfg.population = 12;
    cfg.offspring = 40;
    cfg.threads = 1;
    const GAOutcome a = run_ga(*seed, cfg, 5, 0.0, 6);
    cfg.threads = 4;
    const GAOutcome b = run_ga(*seed, cfg, 5, 0.0, 6);
    assert(a.best_side == b.best_side);
    assert(a.telemetry["ga"]["ratchets"] == b.telemetry["ga"]["ratchets"]);
    assert(a.telemetry["ga"]["moves"] == b.telemetry["ga"]["moves"]);
    assert(strict_report(a.best).ok());
    assert(a.best_side <= a.telemetry["side"]["start_legal"].get<double>());
    std::printf("  GA n=10, 6 generations: legal side %.9f (start %.9f), %d ratchets, same on 1 and 4 threads\n",
                a.best_side, a.telemetry["side"]["start_legal"].get<double>(), a.telemetry["ga"]["ratchets"].get<int>());
}

}  // namespace

int main() {
    test_moves_keep_shape();
    test_d4_keeps_packings_valid();
    test_distances_are_invariant();
    test_selection_spreads();
    test_small_ga_is_legal_and_thread_independent();
    std::printf("ga: all tests passed\n");
}
