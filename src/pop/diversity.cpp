#include "pop/diversity.hpp"

#include "pop/moves.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace pack {
namespace {

constexpr double kR = 0.7071067811865476;

double wrap90(double t) {  // to [-pi/4, pi/4)
    const double q = 0.5 * kPi;
    return t - q * std::floor(t / q + 0.5);
}

// Hungarian algorithm (shortest augmenting paths, potentials), square cost matrix a[n*n].
double assignment_cost(const std::vector<double>& a, int n) {
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> u(static_cast<size_t>(n) + 1), v(static_cast<size_t>(n) + 1), minv(static_cast<size_t>(n) + 1);
    std::vector<int> p(static_cast<size_t>(n) + 1), way(static_cast<size_t>(n) + 1);
    std::vector<char> used(static_cast<size_t>(n) + 1);
    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::fill(minv.begin(), minv.end(), inf);
        std::fill(used.begin(), used.end(), 0);
        do {
            used[static_cast<size_t>(j0)] = 1;
            const int i0 = p[static_cast<size_t>(j0)];
            double delta = inf;
            int j1 = 0;
            for (int j = 1; j <= n; ++j) {
                if (used[static_cast<size_t>(j)]) continue;
                const double cur = a[static_cast<size_t>(i0 - 1) * static_cast<size_t>(n) + static_cast<size_t>(j - 1)] -
                                   u[static_cast<size_t>(i0)] - v[static_cast<size_t>(j)];
                if (cur < minv[static_cast<size_t>(j)]) {
                    minv[static_cast<size_t>(j)] = cur;
                    way[static_cast<size_t>(j)] = j0;
                }
                if (minv[static_cast<size_t>(j)] < delta) {
                    delta = minv[static_cast<size_t>(j)];
                    j1 = j;
                }
            }
            for (int j = 0; j <= n; ++j) {
                if (used[static_cast<size_t>(j)]) {
                    u[static_cast<size_t>(p[static_cast<size_t>(j)])] += delta;
                    v[static_cast<size_t>(j)] -= delta;
                } else {
                    minv[static_cast<size_t>(j)] -= delta;
                }
            }
            j0 = j1;
        } while (p[static_cast<size_t>(j0)] != 0);
        do {
            const int j1 = way[static_cast<size_t>(j0)];
            p[static_cast<size_t>(j0)] = p[static_cast<size_t>(j1)];
            j0 = j1;
        } while (j0);
    }
    double total = 0.0;
    for (int j = 1; j <= n; ++j)
        total += a[static_cast<size_t>(p[static_cast<size_t>(j)] - 1) * static_cast<size_t>(n) + static_cast<size_t>(j - 1)];
    return total;
}

struct Hist {
    int k = 1;
    std::vector<double> f;  // per cell: count, sum cos 4t, sum sin 4t
};

Hist histogram(const Packing& p) {
    Hist hs;
    hs.k = std::max(2, static_cast<int>(std::lround(std::sqrt(static_cast<double>(p.n())) / 2.0)));
    hs.f.assign(3 * static_cast<size_t>(hs.k * hs.k), 0.0);
    for (int i = 0; i < p.n(); ++i) {
        const size_t q = 3 * static_cast<size_t>(i);
        const auto bin = [&](double v) { return std::clamp(static_cast<int>(v / p.h * hs.k), 0, hs.k - 1); };
        const size_t c = 3 * static_cast<size_t>(bin(p.xyt[q + 1]) * hs.k + bin(p.xyt[q]));
        hs.f[c] += 1.0;
        hs.f[c + 1] += std::cos(4.0 * p.xyt[q + 2]);
        hs.f[c + 2] += std::sin(4.0 * p.xyt[q + 2]);
    }
    return hs;
}

}  // namespace

namespace {

// min over D4 of the assignment distance, stopping early once it is known to be <= stop_below.
double assign_min(const Packing& a, const Packing& b0, double stop_below) {
    const int n = a.n();
    if (n == 0) return 0.0;
    std::vector<double> cost(static_cast<size_t>(n) * static_cast<size_t>(n));
    double best = std::numeric_limits<double>::infinity();
    for (int g = 0; g < 8; ++g) {
        Packing b = b0;
        b.h = a.h;
        apply_d4(b, g);
        for (int i = 0; i < n; ++i) {
            const size_t qi = 3 * static_cast<size_t>(i);
            for (int j = 0; j < n; ++j) {
                const size_t qj = 3 * static_cast<size_t>(j);
                const double dx = a.xyt[qi] - b.xyt[qj], dy = a.xyt[qi + 1] - b.xyt[qj + 1];
                const double dt = kR * wrap90(a.xyt[qi + 2] - b.xyt[qj + 2]);
                cost[static_cast<size_t>(i) * static_cast<size_t>(n) + static_cast<size_t>(j)] = dx * dx + dy * dy + dt * dt;
            }
        }
        best = std::min(best, assignment_cost(cost, n));
        if (std::sqrt(best / n) <= stop_below) break;  // already close: the other symmetries cannot matter
    }
    return std::sqrt(best / n);
}

}  // namespace

double assign_distance(const Packing& a, const Packing& b) { return assign_min(a, b, -1.0); }

double hist_distance(const Packing& a, const Packing& b0) {
    const Hist ha = histogram(a);
    double best = std::numeric_limits<double>::infinity();
    for (int g = 0; g < 8; ++g) {
        Packing b = b0;
        b.h = a.h;
        apply_d4(b, g);
        const Hist hb = histogram(b);
        double s = 0.0;
        for (size_t c = 0; c < ha.f.size(); ++c) s += (ha.f[c] - hb.f[c]) * (ha.f[c] - hb.f[c]);
        best = std::min(best, s);
    }
    return std::sqrt(best / std::max(1, a.n()));
}

std::vector<int> select_diverse(const std::vector<Packing>& cand, const std::vector<int>& order, int count,
                                const DiversityParams& prm, long* distance_evals, const ParallelFor& parfor) {
    count = std::min<int>(count, static_cast<int>(order.size()));
    std::vector<int> chosen;
    std::vector<char> taken(cand.size(), 0);
    const int elite = std::min(count, std::max(1, static_cast<int>(std::lround(prm.elite * count))));
    for (int k = 0; k < elite; ++k) {
        chosen.push_back(order[static_cast<size_t>(k)]);
        taken[static_cast<size_t>(order[static_cast<size_t>(k)])] = 1;
    }
    const int n = cand.empty() ? 0 : cand[0].n();
    const bool assign = prm.metric == Metric::Assign || (prm.metric == Metric::Auto && n <= prm.assign_max_n);
    const double thr = assign ? prm.threshold : prm.hist_threshold;
    for (size_t k = static_cast<size_t>(elite); k < order.size() && static_cast<int>(chosen.size()) < count; ++k) {
        const int c = order[k];
        const auto dist = [&](int s) {
            return assign ? assign_min(cand[static_cast<size_t>(c)], cand[static_cast<size_t>(s)], thr)
                          : hist_distance(cand[static_cast<size_t>(c)], cand[static_cast<size_t>(s)]);
        };
        bool far = true;
        if (parfor && chosen.size() > 1) {
            std::vector<char> near(chosen.size(), 0);
            parfor(static_cast<int>(chosen.size()), [&](int q, int) { near[static_cast<size_t>(q)] = dist(chosen[static_cast<size_t>(q)]) <= thr; });
            if (distance_evals) *distance_evals += static_cast<long>(chosen.size());
            for (const char x : near) far = far && !x;
        } else {
            for (const int s : chosen) {
                if (distance_evals) ++*distance_evals;
                if (dist(s) <= thr) {
                    far = false;
                    break;
                }
            }
        }
        if (far) {
            chosen.push_back(c);
            taken[static_cast<size_t>(c)] = 1;
        }
    }
    for (size_t k = 0; k < order.size() && static_cast<int>(chosen.size()) < count; ++k) {
        if (!taken[static_cast<size_t>(order[k])]) {
            chosen.push_back(order[k]);
            taken[static_cast<size_t>(order[k])] = 1;
        }
    }
    return chosen;
}

}  // namespace pack
