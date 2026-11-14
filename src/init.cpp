#include "init.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

namespace pack {

System grid_init(const Shape& shape, int n) {
    assert(n >= 1);
    assert(!shape.vertices.empty());
    const int k = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(n))));

    double xlo = shape.vertices[0].x, xhi = xlo, ylo = shape.vertices[0].y, yhi = ylo;
    for (const Vec2& v : shape.vertices) {
        xlo = std::min(xlo, v.x);
        xhi = std::max(xhi, v.x);
        ylo = std::min(ylo, v.y);
        yhi = std::max(yhi, v.y);
    }
    const double pitch = std::max(xhi - xlo, yhi - ylo);
    const double off_x = 0.5 * (xlo + xhi), off_y = 0.5 * (ylo + yhi);

    System s;
    s.species = {shape};
    s.L = k * pitch;
    s.particles.resize(static_cast<size_t>(n));

    // Interior cells first, then the boundary ring, each in row-major order. Taking the first n of that list is exactly "fill the grid and leave out boundary cells", and the cells it leaves out are the last of the top row.
    int placed = 0;
    for (int pass = 0; pass < 2; ++pass) {
        const bool want_edge = (pass == 1);
        for (int cy = 0; cy < k && placed < n; ++cy) {
            for (int cx = 0; cx < k && placed < n; ++cx) {
                const bool edge = (cx == 0 || cy == 0 || cx == k - 1 || cy == k - 1);
                if (edge != want_edge) continue;
                set_pose(s.particles[static_cast<size_t>(placed++)], (cx + 0.5) * pitch - off_x,
                         (cy + 0.5) * pitch - off_y, 0.0);
            }
        }
    }
    assert(placed == n);

    s.recount();  // exactly 0: every shape sits inside its own cell
    return s;
}

namespace {

// Squares in one staircase block whose cells satisfy a + b <= m, counted from the corner.
long staircase_cells(long m) { return (m < 0) ? 0 : (m + 1) * (m + 2) / 2; }

struct DiagonalPlan {
    double side = 0.0;
    long m = -1;
    long strip = 0;
};

DiagonalPlan diagonal_plan_at(double s) {
    const double r = 1.0 / std::sqrt(2.0);
    DiagonalPlan p;
    p.side = s;
    // A hair of slack in both floors, so a side computed as exactly a breakpoint is not
    // knocked one short by the rounding in computing it.
    p.m = static_cast<long>(std::floor(s - 2.0 - r + 1e-12));
    p.strip = std::max(0L, static_cast<long>(std::floor(s * std::sqrt(2.0) - 1.0 + 1e-12)));
    return p;
}

long diagonal_capacity(const DiagonalPlan& p) { return 2 * staircase_cells(p.m) + p.strip; }

// The smallest side of the family holding n. Capacity only steps up where the staircase or the
// strip gains a square, so the minimum sits on one of those breakpoints: for each staircase
// size, the side where it first fits, or the side where the strip first covers the rest.
DiagonalPlan diagonal_plan(int n) {
    const double r = 1.0 / std::sqrt(2.0);
    DiagonalPlan best;
    best.side = std::numeric_limits<double>::infinity();
    for (long m = -1;; ++m) {
        const double s_m = (m < 0) ? 0.0 : static_cast<double>(m) + 2.0 + r;
        if (s_m >= best.side) break;
        const long rest = std::max(0L, static_cast<long>(n) - 2 * staircase_cells(m));
        const double s_strip = static_cast<double>(rest + 1) * r;  // floor(s sqrt2 - 1) >= rest
        const DiagonalPlan p = diagonal_plan_at(std::max(s_m, s_strip));
        if (diagonal_capacity(p) >= n && p.side < best.side) best = p;
    }
    return best;
}

}  // namespace

System diagonal_init(const Shape& shape, int n) {
    assert(n >= 1);
    // Only a square has the staircase and the strip; anything else keeps the lattice.
    double xlo = shape.vertices[0].x, xhi = xlo, ylo = shape.vertices[0].y, yhi = ylo;
    for (const Vec2& v : shape.vertices) {
        xlo = std::min(xlo, v.x);
        xhi = std::max(xhi, v.x);
        ylo = std::min(ylo, v.y);
        yhi = std::max(yhi, v.y);
    }
    const double a = xhi - xlo;
    // Centred too: the strip turns its squares, and an off-centre origin would carry them
    // off their spots as it turned.
    if (shape.vertices.size() != 4 || std::abs((yhi - ylo) - a) > 1e-12 * a ||
        std::abs(shape.area - a * a) > 1e-12 * a * a ||
        std::abs(xlo + xhi) > 1e-12 * a || std::abs(ylo + yhi) > 1e-12 * a) {
        return grid_init(shape, n);
    }
    const double r = 1.0 / std::sqrt(2.0);

    const DiagonalPlan plan = diagonal_plan(n);
    const double s = plan.side;

    // Centres in units of the side, in the order they are kept: staircase cells nearest the
    // corners first (smallest a + b), then the strip from its middle outwards, so a surplus
    // comes off the strip's ends before it ever thins a staircase.
    std::vector<Pose> stair, strip;
    for (long d = 0; d <= plan.m; ++d) {
        for (long b = 0; b <= d; ++b) {
            const long c = d - b;
            // Lower right: x in [s - c - 1, s - c], y in [b, b + 1]. Upper left: its mirror in
            // the diagonal.
            stair.push_back({s - c - 0.5, b + 0.5, 0.0});
            stair.push_back({b + 0.5, s - c - 0.5, 0.0});
        }
    }
    for (long q = 0; q < plan.strip; ++q) {
        const double u = r + r * static_cast<double>(q);  // centre (u, u), one unit apart
        strip.push_back({u, u, 0.25 * kPi});
    }
    std::stable_sort(strip.begin(), strip.end(), [s](const Pose& p, const Pose& q) {
        return std::abs(p.x - 0.5 * s) < std::abs(q.x - 0.5 * s);
    });

    std::vector<Pose> keep(stair);
    keep.insert(keep.end(), strip.begin(), strip.end());
    assert(static_cast<long>(keep.size()) >= n);
    keep.resize(static_cast<size_t>(n));  // staircases first, then the strip from its middle

    System sys;
    sys.species = {shape};
    sys.L = s * a;
    sys.particles.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const Pose& p = keep[static_cast<size_t>(i)];
        set_pose(sys.particles[static_cast<size_t>(i)], p.x * a, p.y * a, p.theta);
    }
    sys.recount();
    return sys;
}

}  // namespace pack
