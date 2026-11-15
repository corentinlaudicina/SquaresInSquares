// Population kernel: gradients against finite differences, agreement with the System energy and the
// SAT gap, relaxation to 1e-12, and thread-count independence of the batch backend.

#include "io/io.hpp"
#include "legalize.hpp"
#include "pop/backend.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>

using namespace pack;

namespace {

Packing random_packing(int n, double h, unsigned seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    Packing p;
    p.h = h;
    for (int i = 0; i < n; ++i) {
        p.xyt.push_back(0.3 + u(rng) * (h - 0.6));  // some cross the walls
        p.xyt.push_back(0.3 + u(rng) * (h - 0.6));
        p.xyt.push_back(u(rng) * kPi);
    }
    return p;
}

void test_gradients(Penalty pen, const char* name) {
    Workspace ws;
    int checked = 0, kinks = 0;
    for (unsigned seed = 1; seed <= 5; ++seed) {
        Packing p = random_packing(30, 4.5, seed);
        std::vector<double> g(p.xyt.size());
        cost_grad(p, pen, ws, g.data());
        for (size_t k = 0; k < p.xyt.size(); ++k) {
            const auto fd = [&](double e) {
                const double v = p.xyt[k];
                p.xyt[k] = v + e;
                const double a = cost_grad(p, pen, ws, nullptr).cost;
                p.xyt[k] = v - e;
                const double b = cost_grad(p, pen, ws, nullptr).cost;
                p.xyt[k] = v;
                return (a - b) / (2 * e);
            };
            const double f1 = fd(1e-6), f2 = fd(5e-7);
            const double scale = std::max({1e-8, std::abs(f1), std::abs(g[k])});
            if (std::abs(f1 - f2) > 1e-5 * scale) {  // a kink (contact set or SAT axis switch) in the stencil
                ++kinks;
                continue;
            }
            ++checked;
            if (std::abs(f1 - g[k]) > 1e-5 * scale) {
                std::printf("%s: component %zu (seed %u): analytic %.10g, fd %.10g\n", name, k, seed, g[k], f1);
                assert(false);
            }
        }
    }
    assert(kinks * 50 < checked);  // under 2% skipped
    std::printf("  %s gradient: %d components match finite differences (%d at kinks skipped)\n", name, checked, kinks);
}

void test_area2_matches_system() {
    Packing p = random_packing(40, 5.0, 9);
    System s = to_system(p);
    s.k_pair = s.k_wall = 1.0;
    s.recount();
    Workspace ws;
    const double a = cost_grad(p, Penalty::Area2, ws, nullptr).cost, b = s.total_energy();
    assert(std::abs(a - b) <= 1e-12 * std::max(1.0, b));
}

void test_depth_is_sat_penetration() {
    std::mt19937_64 rng(3);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    Workspace ws;
    int tested = 0;
    for (int t = 0; t < 2000; ++t) {
        Packing p;
        p.h = 10.0;
        p.xyt = {5.0, 5.0, u(rng) * kPi, 5.0 + 0.8 * u(rng), 5.0 + 0.8 * u(rng), u(rng) * kPi};
        const double depth = cost_grad(p, Penalty::Depth2, ws, nullptr).max_depth;
        const double gap = strict_report(to_system(p)).min_pair_gap;  // legalize's independent SAT
        if (gap >= 0.0) {
            assert(depth == 0.0);
            continue;
        }
        assert(std::abs(depth + gap) < 1e-12);
        ++tested;
    }
    assert(tested > 1000);
}

void test_relax_reaches_tolerance() {
    // 6x6 grid jiggled, in its optimal box plus 1e-4: relax must find a packing again. (In the exact
    // box the only solution is the degenerate grid itself; no relaxer converges to that quickly.)
    std::mt19937_64 rng(5);
    std::normal_distribution<double> g(0.0, 0.05);
    Packing p;
    p.h = 6.0001;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) p.xyt.insert(p.xyt.end(), {i + 0.5 + g(rng), j + 0.5 + g(rng), g(rng)});
    Workspace ws;
    for (Penalty pen : {Penalty::Depth2, Penalty::Area2}) {
        Packing q = p;
        CappedParams prm;
        prm.penalty = pen;
        prm.tol = pen == Penalty::Depth2 ? 1e-12 : 1e-20;
        const CappedStats st = relax_capped(q, prm, ws);
        std::printf("  jiggled 6x6, %s: reached=%d after %d iterations, worst %.2e / wall %.2e\n",
                    pen == Penalty::Depth2 ? "depth2" : "area2", st.reached, st.iters, st.info.max_depth,
                    st.info.max_wall);
        if (pen == Penalty::Depth2) {
            assert(st.reached);
            // and then legalization is essentially free
            System s = to_system(q);
            const LegalizeResult lr = legalize(s);
            assert(lr.ok && lr.side_after <= 6.0001 + 1e-12);
        }
    }
}

void test_verlet_list_misses_nothing() {
    std::mt19937_64 rng(17);
    std::normal_distribution<double> step(0.0, 0.03);
    for (Penalty pen : {Penalty::Depth2, Penalty::Area2}) {
        Packing p = random_packing(120, 9.5, 21);
        Workspace listed, fresh;
        listed.verlet = true;
        std::vector<double> ga(p.xyt.size()), gb(p.xyt.size());
        for (int it = 0; it < 300; ++it) {
            for (double& v : p.xyt) v += step(rng);  // ~0.03 per step: a rebuild every few steps
            const CostInfo a = cost_grad(p, pen, listed, ga.data());
            const CostInfo b = cost_grad(p, pen, fresh, gb.data());  // the grid, every time
            assert(std::abs(a.cost - b.cost) <= 1e-12 * std::max(1.0, b.cost));
            assert(a.max_depth == b.max_depth && a.max_wall == b.max_wall);
            for (size_t k = 0; k < ga.size(); ++k) assert(std::abs(ga[k] - gb[k]) <= 1e-10 * std::max(1.0, std::abs(gb[k])));
        }
        std::printf("  verlet list (%s): 300 steps match the grid exactly, %ld rebuilds\n",
                    pen == Penalty::Depth2 ? "depth2" : "area2", listed.verlet_builds);
        assert(listed.verlet_builds > 5 && listed.verlet_builds < 300);
    }
}

void test_backend_is_thread_independent() {
    std::vector<Packing> a;
    for (unsigned s = 0; s < 24; ++s) a.push_back(random_packing(25, 5.4, 100 + s));
    std::vector<Packing> b = a;
    std::vector<CappedStats> sa, sb;
    CappedParams prm;
    CpuBackend(1).relax_batch(a, prm, sa);
    CpuBackend(4).relax_batch(b, prm, sb);
    int reached = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        assert(a[i].xyt == b[i].xyt && sa[i].iters == sb[i].iters);
        reached += sa[i].reached;
    }
    std::printf("  backend: 24 packings identical on 1 and 4 threads (%d reached tol)\n", reached);
}

}  // namespace

int main() {
    test_gradients(Penalty::Depth2, "depth2");
    test_gradients(Penalty::Area2, "area2");
    test_area2_matches_system();
    test_depth_is_sat_penetration();
    test_relax_reaches_tolerance();
    test_verlet_list_misses_nothing();
    test_backend_is_thread_independent();
    std::printf("pop: all tests passed\n");
}
