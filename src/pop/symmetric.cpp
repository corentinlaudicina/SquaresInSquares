#include "pop/symmetric.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace pack {
namespace {

// Offset (dx, dy) turned by k steps of the orbit: C2 a half turn, C4 a quarter turn. Exact.
inline void turn(double& dx, double& dy, int k, int orbit) {
    for (int q = 0; q < k; ++q) {
        if (orbit == 2) {
            dx = -dx;
            dy = -dy;
        } else {
            const double t = dx;
            dx = -dy;
            dy = t;
        }
    }
}

// The transpose turn, for gradients.
inline void turn_back(double& gx, double& gy, int k, int orbit) {
    for (int q = 0; q < k; ++q) {
        if (orbit == 2) {
            gx = -gx;
            gy = -gy;
        } else {
            const double t = gx;
            gx = gy;
            gy = -t;
        }
    }
}

bool in_region(double dx, double dy, int orbit) {
    if (orbit == 2) return dy < 0.0 || (dy == 0.0 && dx < 0.0);
    return dx < 0.0 && dy <= 0.0;
}

int orbit_gens(const SymLayout& L) { return L.m - (L.centre ? 1 : 0); }

}  // namespace

SymLayout make_layout(Sym sym, int n) {
    SymLayout L;
    L.n = n;
    if (sym == Sym::None) {
        L.m = n;
        return L;
    }
    const int orbit = sym == Sym::C2 ? 2 : 4;
    const int r = n % orbit;
    if (r > 1) return make_layout(Sym::None, n);
    L.sym = sym;
    L.orbit = orbit;
    L.centre = r == 1;
    L.m = n / orbit + (L.centre ? 1 : 0);
    return L;
}

void expand(const Packing& gen, const SymLayout& L, Packing& full) {
    if (L.sym == Sym::None) {
        full = gen;
        return;
    }
    const double c = 0.5 * gen.h;
    const int g = orbit_gens(L);
    full.h = gen.h;
    full.xyt.resize(3 * static_cast<size_t>(L.n));
    for (int k = 0; k < L.orbit; ++k) {
        for (int i = 0; i < g; ++i) {
            const size_t s = 3 * static_cast<size_t>(i), d = 3 * static_cast<size_t>(k * g + i);
            double dx = gen.xyt[s] - c, dy = gen.xyt[s + 1] - c;
            turn(dx, dy, k, L.orbit);
            full.xyt[d] = c + dx;
            full.xyt[d + 1] = c + dy;
            full.xyt[d + 2] = gen.xyt[s + 2];
        }
    }
    if (L.centre) {
        const size_t s = 3 * static_cast<size_t>(g), d = 3 * static_cast<size_t>(L.n - 1);
        full.xyt[d] = c;
        full.xyt[d + 1] = c;
        full.xyt[d + 2] = gen.xyt[s + 2];
    }
}

Packing expanded(const Packing& gen, const SymLayout& L) {
    Packing f;
    expand(gen, L, f);
    return f;
}

void pull_gradient(const double* gfull, const SymLayout& L, double* ggen) {
    const int g = orbit_gens(L);
    std::fill(ggen, ggen + 3 * static_cast<size_t>(L.m), 0.0);
    for (int k = 0; k < L.orbit; ++k) {
        for (int i = 0; i < g; ++i) {
            const size_t s = 3 * static_cast<size_t>(k * g + i), d = 3 * static_cast<size_t>(i);
            double gx = gfull[s], gy = gfull[s + 1];
            turn_back(gx, gy, k, L.orbit);
            ggen[d] += gx;
            ggen[d + 1] += gy;
            ggen[d + 2] += gfull[s + 2];
        }
    }
    if (L.centre) ggen[3 * static_cast<size_t>(g) + 2] = gfull[3 * static_cast<size_t>(L.n - 1) + 2];
}

void fix_centre(Packing& gen, const SymLayout& L) {
    if (!L.centre) return;
    const size_t s = 3 * static_cast<size_t>(L.m - 1);
    gen.xyt[s] = gen.xyt[s + 1] = 0.5 * gen.h;
}

void canonicalize(Packing& gen, const SymLayout& L) {
    if (L.sym == Sym::None) return;
    const double c = 0.5 * gen.h;
    for (int i = 0; i < orbit_gens(L); ++i) {
        const size_t s = 3 * static_cast<size_t>(i);
        for (int k = 0; k < L.orbit; ++k) {
            double dx = gen.xyt[s] - c, dy = gen.xyt[s + 1] - c;
            turn(dx, dy, k, L.orbit);
            if (in_region(dx, dy, L.orbit)) {
                gen.xyt[s] = c + dx;
                gen.xyt[s + 1] = c + dy;
                break;
            }
        }
    }
    fix_centre(gen, L);
}

Packing symmetrize(const Packing& full, const SymLayout& L) {
    if (L.sym == Sym::None) return full;
    const double c = 0.5 * full.h;
    const int n = full.n(), g = orbit_gens(L);
    Packing canon = full;
    canonicalize(canon, make_layout(L.sym, L.orbit * n));  // every square into the region (no centre)
    std::vector<int> idx(static_cast<size_t>(n));
    std::iota(idx.begin(), idx.end(), 0);
    const auto r2 = [&](int i) {
        const double dx = full.xyt[3 * static_cast<size_t>(i)] - c, dy = full.xyt[3 * static_cast<size_t>(i) + 1] - c;
        return dx * dx + dy * dy;
    };
    int centre_sq = -1;
    if (L.centre) centre_sq = *std::min_element(idx.begin(), idx.end(), [&](int a, int b) { return r2(a) < r2(b); });
    std::vector<int> chosen, spare;
    for (const int i : idx) {
        if (i == centre_sq) continue;
        bool dup = false;
        for (const int j : chosen) {
            const double dx = canon.xyt[3 * static_cast<size_t>(i)] - canon.xyt[3 * static_cast<size_t>(j)];
            const double dy = canon.xyt[3 * static_cast<size_t>(i) + 1] - canon.xyt[3 * static_cast<size_t>(j) + 1];
            if (dx * dx + dy * dy < 0.25) dup = true;
        }
        (dup || static_cast<int>(chosen.size()) >= g ? spare : chosen).push_back(i);
    }
    for (size_t k = 0; static_cast<int>(chosen.size()) < g && k < spare.size(); ++k) chosen.push_back(spare[k]);
    Packing gen;
    gen.h = full.h;
    for (const int i : chosen)
        gen.xyt.insert(gen.xyt.end(), {canon.xyt[3 * static_cast<size_t>(i)], canon.xyt[3 * static_cast<size_t>(i) + 1],
                                       canon.xyt[3 * static_cast<size_t>(i) + 2]});
    if (L.centre) gen.xyt.insert(gen.xyt.end(), {c, c, full.xyt[3 * static_cast<size_t>(centre_sq) + 2]});
    return gen;
}

void crossover_sym(Packing& child, const Packing& mate, const MoveParams& prm, std::mt19937_64& rng, const SymLayout& L) {
    const double h = child.h, c = 0.5 * h;
    canonicalize(child, L);
    Packing b = mate;
    b.h = h;
    if (u01(rng) < prm.crossover_symmetry) {
        // a box symmetry maps a C2 / C4 packing to another one: apply it to the expansion, take generators again
        Packing full = expanded(b, L);
        apply_d4(full, uniform_int(rng, 1, 7));
        b = symmetrize(full, L);
    }
    canonicalize(b, L);
    const double px = uniform(rng, 0.0, L.orbit == 2 ? h : c), py = uniform(rng, 0.0, c);
    const int g = orbit_gens(L);
    std::vector<int> mine, theirs(static_cast<size_t>(g));
    for (int i = 0; i < g; ++i) {
        const double dx = child.xyt[3 * static_cast<size_t>(i)] - px, dy = child.xyt[3 * static_cast<size_t>(i) + 1] - py;
        if (dx * dx + dy * dy <= prm.crossover_radius * prm.crossover_radius) mine.push_back(i);
    }
    if (mine.empty() || g == 0) return;
    std::iota(theirs.begin(), theirs.end(), 0);
    const auto d2 = [&](int j) {
        const double dx = b.xyt[3 * static_cast<size_t>(j)] - px, dy = b.xyt[3 * static_cast<size_t>(j) + 1] - py;
        return dx * dx + dy * dy;
    };
    std::partial_sort(theirs.begin(), theirs.begin() + static_cast<long>(mine.size()), theirs.end(),
                      [&](int x, int y) { return d2(x) < d2(y) || (d2(x) == d2(y) && x < y); });
    for (size_t q = 0; q < mine.size(); ++q)
        for (int k = 0; k < 3; ++k)
            child.xyt[3 * static_cast<size_t>(mine[q]) + static_cast<size_t>(k)] =
                b.xyt[3 * static_cast<size_t>(theirs[q]) + static_cast<size_t>(k)];
    fix_centre(child, L);
}

CostInfo cost_symmetric(const Packing& gen, const SymLayout& L, Penalty pen, Workspace& ws, double* ggen) {
    expand(gen, L, ws.sym_full);
    if (!ggen) return cost_grad(ws.sym_full, pen, ws, nullptr);
    ws.sym_g.resize(ws.sym_full.xyt.size());
    const CostInfo ci = cost_grad(ws.sym_full, pen, ws, ws.sym_g.data());
    pull_gradient(ws.sym_g.data(), L, ggen);
    return ci;
}

}  // namespace pack
