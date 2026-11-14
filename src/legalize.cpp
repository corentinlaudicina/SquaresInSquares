#include "legalize.hpp"

#include "moves/relax.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace pack {
namespace {

// Largest separation of a and b over the edge normals of both (positive = a separating gap).
double sat_gap(const Particle& a, const Particle& b, const Shape& s) {
    double best = -std::numeric_limits<double>::infinity();
    for (const Particle* p : {&a, &b}) {
        for (int k = 0; k < static_cast<int>(s.sat_axes.size()); ++k) {
            const Vec2 axis = world_axis(*p, s, k);
            double alo, ahi, blo, bhi;
            project(a, s, axis, alo, ahi);
            project(b, s, axis, blo, bhi);
            best = std::max(best, std::max(blo - ahi, alo - bhi));
        }
    }
    return best;
}

bool axis_aligned(const Particle& p) { return p.theta() == 0.0; }

// Shift so the lowest corner is a margin above 0 on both axes, and set L to the extent plus that
// margin. The shift is a multiple of q = ulp(P), P a power of two above every coordinate, so every
// x + shift is exact: squares that touched exactly still touch exactly, rather than rounding into
// a 1e-16 overlap. The margin (4 q) absorbs the rounding of tilted corners.
void refit(System& sys) {
    const double inf = std::numeric_limits<double>::infinity();
    double xlo = inf, xhi = -inf, ylo = inf, yhi = -inf;
    for (const Particle& p : sys.particles) {
        double lo, hi;
        project(p, sys.base_shape(), {1.0, 0.0}, lo, hi);
        xlo = std::min(xlo, lo);
        xhi = std::max(xhi, hi);
        project(p, sys.base_shape(), {0.0, 1.0}, lo, hi);
        ylo = std::min(ylo, lo);
        yhi = std::max(yhi, hi);
    }
    // P = the power of two one binade above every coordinate before and after the shift.
    const double bound = std::max({std::abs(xlo), std::abs(xhi), std::abs(ylo), std::abs(yhi), xhi - xlo,
                                   yhi - ylo, 1.0}) * 1.01;
    const double q = std::ldexp(1.0, std::ilogb(bound) + 1 - 52);
    const double m = 4.0 * q;
    const double dx = std::ceil((m - xlo) / q) * q, dy = std::ceil((m - ylo) / q) * q;
    for (Particle& p : sys.particles) p.set_pose(p.x() + dx, p.y() + dy, p.theta());
    sys.L = std::ceil((std::max(xhi + dx, yhi + dy) + m) / q) * q;
    sys.recount();
}

}  // namespace

StrictReport strict_report(const System& sys) {
    StrictReport r;
    r.min_pair_gap = r.min_wall_gap = std::numeric_limits<double>::infinity();
    const Shape& s = sys.base_shape();
    for (const Particle& p : sys.particles) {
        double xlo, xhi, ylo, yhi;
        project(p, s, {1.0, 0.0}, xlo, xhi);
        project(p, s, {0.0, 1.0}, ylo, yhi);
        const double g = std::min({xlo, ylo, sys.L - xhi, sys.L - yhi});
        r.min_wall_gap = std::min(r.min_wall_gap, g);
        if (axis_aligned(p) ? g < 0.0 : g <= 0.0) ++r.bad_walls;
    }
    CellList grid;
    grid.build(sys.particles, sys.L, sys.reach());
    grid.for_each_close_pair([&](int i, int j) {
        const Particle& a = sys.particles[static_cast<size_t>(i)];
        const Particle& b = sys.particles[static_cast<size_t>(j)];
        const double dx = a.x() - b.x(), dy = a.y() - b.y();
        if (dx * dx + dy * dy > 2.0 * s.circumradius * 2.0 * s.circumradius) return;
        const double g = sat_gap(a, b, s);
        r.min_pair_gap = std::min(r.min_pair_gap, g);
        if (axis_aligned(a) && axis_aligned(b) ? g < 0.0 : g <= kSepTol) ++r.bad_pairs;
    });
    return r;
}

LegalizeResult legalize(System& sys, int max_growth) {
    LegalizeResult res;
    res.side_before = sys.tight_side();

    refit(sys);
    if (!strict_report(sys).ok()) {
        System big = sys;
        big.species[0] = make_square(1.0 + kInflate);
        big.recount();
        RelaxOptions ro;
        ro.max_iterations = 3000;
        ro.feasible_tol = 0.0;  // relax until nothing overlaps at all; the strict check decides
        ro.energy_tol = 0.0;
        // Grow geometrically until the inflated copy separates, then bisect between the last box
        // that failed and the first that worked: doubling alone overshoots by up to 2x.
        const auto attempt = [&](System& trial) {
            const RelaxResult rr = relax(trial, ro);
            res.relax_iterations += rr.iterations;
            System real = sys;
            for (size_t i = 0; i < real.particles.size(); ++i) {
                const Particle& p = trial.particles[i];
                real.particles[i].set_pose(p.x(), p.y(), p.theta());
            }
            real.L = trial.L;
            real.recount();
            return real;
        };
        const double L0 = big.L;
        double lo = L0, hi = 0.0;
        System good;
        // Growth is cumulative: each attempt continues from the previous attempt's relaxed state.
        for (int k = 0; k <= max_growth; ++k) {
            if (k > 0) big.rescale_box(big.L * (1.0 + 1e-10 * std::pow(2.0, k - 1)));
            System real = attempt(big);
            if (strict_report(real).ok()) {
                hi = big.L;
                good = real;
                break;
            }
            lo = big.L;
            ++res.growth_steps;
        }
        if (hi > 0.0) {
            for (int b = 0; b < 8 && hi - lo > 1e-12 * L0; ++b) {
                const double mid = 0.5 * (lo + hi);
                System trial = big;
                for (size_t i = 0; i < trial.particles.size(); ++i) {  // start from the good poses
                    const Particle& p = good.particles[i];
                    trial.particles[i].set_pose(p.x(), p.y(), p.theta());
                }
                trial.L = hi;
                trial.rescale_box(mid);
                System real = attempt(trial);
                if (strict_report(real).ok()) {
                    hi = mid;
                    good = real;
                } else {
                    lo = mid;
                }
            }
            sys = good;
        }
    }
    refit(sys);
    res.report = strict_report(sys);
    res.ok = res.report.ok();
    res.side_after = sys.L;
    return res;
}

}  // namespace pack
