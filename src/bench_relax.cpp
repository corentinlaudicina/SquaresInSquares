// Micro-benchmark for the population relaxation (Phase 3 item 2): ~1000 perturbed packings at one N,
// relaxed by the CPU backend. Reports kernel cost per evaluation, relax time and success, thread
// scaling, and the classic relax() on a subset for reference.
//
//   build/bench_relax [--N 200] [--count 1000] [--start best/200.json] [--threads 1,2,4,8]

#include "engine/solution_io.hpp"
#include "moves/relax.hpp"
#include "pop/backend.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <sstream>
#include <string>

using namespace pack;
using Clock = std::chrono::steady_clock;

namespace {

double secs(Clock::time_point a) { return std::chrono::duration<double>(Clock::now() - a).count(); }

// GA-like offspring: relocate-ish jiggle of a few squares, plus a small jiggle of everyone.
std::vector<Packing> offspring(const Packing& base, int count, double h, unsigned seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> small(0.0, 0.01), big(0.0, 0.3);
    std::uniform_int_distribution<int> pick(0, base.n() - 1);
    std::vector<Packing> out(static_cast<size_t>(count), base);
    for (Packing& p : out) {
        p.h = h;
        for (double& v : p.xyt) v += small(rng);
        for (int k = 0; k < 5; ++k) {
            const size_t i = 3 * static_cast<size_t>(pick(rng));
            p.xyt[i] += big(rng);
            p.xyt[i + 1] += big(rng);
            p.xyt[i + 2] += big(rng);
        }
    }
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    int n = 200, count = 1000;
    std::string start, threads_arg = "1,2,4,8";
    for (int k = 1; k + 1 < argc; k += 2) {
        const std::string a = argv[k];
        if (a == "--N") n = std::atoi(argv[k + 1]);
        else if (a == "--count") count = std::atoi(argv[k + 1]);
        else if (a == "--start") start = argv[k + 1];
        else if (a == "--threads") threads_arg = argv[k + 1];
    }
    if (start.empty()) {
        char p[64];
        std::snprintf(p, sizeof p, "best/%d.json", n);
        start = p;
    }
    std::string why;
    auto sys = read_solution(start, &why);
    if (!sys) {
        std::fprintf(stderr, "%s\n", why.c_str());
        return 1;
    }
    const Packing base = to_packing(*sys);
    std::printf("bench_relax: N=%d, %d offspring from %s (side %.6f), hardware threads %u\n", n, count, start.c_str(),
                base.h, std::thread::hardware_concurrency());

    // 1. Kernel cost per evaluation, single thread.
    {
        Workspace ws;
        const std::vector<Packing> pop = offspring(base, 50, base.h, 7);
        std::vector<double> g(base.xyt.size());
        for (Penalty pen : {Penalty::Depth2, Penalty::Area2}) {
            const auto t = Clock::now();
            int evals = 0;
            for (int rep = 0; rep < 20; ++rep)
                for (const Packing& p : pop) {
                    cost_grad(p, pen, ws, g.data());
                    ++evals;
                }
            const double s = secs(t);
            std::printf("  kernel %-6s cost+grad: %7.1f us/eval (%.0f evals/s, 1 thread)\n",
                        pen == Penalty::Depth2 ? "depth2" : "area2", 1e6 * s / evals, evals / s);
        }
    }

    // 2. Batch relaxation, two boxes, thread scaling.
    std::vector<int> threads;
    {
        std::stringstream ss(threads_arg);
        for (std::string t; std::getline(ss, t, ',');) threads.push_back(std::atoi(t.c_str()));
    }
    for (const double slack : {1e-3, 0.0}) {
        const double h = base.h * (1.0 + slack);
        const std::vector<Packing> pop = offspring(base, count, h, 11);
        for (Penalty pen : {Penalty::Depth2, Penalty::Area2}) {
            CappedParams prm;
            prm.penalty = pen;
            prm.tol = pen == Penalty::Depth2 ? 1e-9 : 1e-18;  // area 1e-18 ~ depth 1e-9 at a corner
            for (int th : threads) {
                std::vector<Packing> batch = pop;
                std::vector<CappedStats> st;
                CpuBackend be(th);
                const auto t = Clock::now();
                be.relax_batch(batch, prm, st);
                const double s = secs(t);
                long iters = 0;
                int reached = 0;
                for (const CappedStats& c : st) {
                    iters += c.iters;
                    reached += c.reached;
                }
                std::printf("  box +%.0e %-6s threads %d: %6.2f s total, %6.2f ms/relax, %5.0f iters/relax, "
                            "%6.0f relax/s, reached %4d/%d\n",
                            slack, pen == Penalty::Depth2 ? "depth2" : "area2", th, s, 1e3 * s / count,
                            static_cast<double>(iters) / count, count / s, reached, count);
            }
        }
    }

    // 3. Reference: the classic relax() (L-BFGS + line search, k = 100) on 50 of the slack offspring.
    {
        const std::vector<Packing> pop = offspring(base, 50, base.h * 1.001, 11);
        const auto t = Clock::now();
        int ok = 0;
        long iters = 0;
        for (const Packing& p : pop) {
            System s = to_system(p);
            const RelaxResult r = relax(s, RelaxOptions{});
            ok += r.feasible;
            iters += r.iterations;
        }
        const double s = secs(t);
        std::printf("  classic relax() 1 thread: %6.2f ms/relax, %5.0f iters/relax, feasible (area 1e-10) %d/50\n",
                    1e3 * s / 50, static_cast<double>(iters) / 50, ok);
    }
}
