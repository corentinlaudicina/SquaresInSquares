#include "moves/regions.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace pack {
namespace {

double uniform(std::mt19937_64& rng, double lo, double hi) {
    if (!(hi > lo)) return lo;
    std::uniform_real_distribution<double> d(lo, hi);
    return d(rng);
}

Vec2 perp(Vec2 v) { return {-v.y, v.x}; }

// Diameter of the shape: the unit the RegionWeights sizes are quoted in, so the same weights describe the same physical region whatever shape or side is in use. Reads the base species (the ordinary unit square) rather than being per-particle aware -- region-picking only ever sizes bands/corners for melting ordinary squares; a merged particle must be dissolved before anything here could touch it.
double diameter(const System& sys) {
    return (sys.base_shape().circumradius > 0.0) ? 2.0 * sys.base_shape().circumradius : 1.0;
}

// Support interval of the box's four corners along `axis`: how far a band offset can travel before the band has left the box entirely.
void box_span(double L, Vec2 axis, double& lo, double& hi) {
    const Vec2 corners[4] = {{0.0, 0.0}, {L, 0.0}, {L, L}, {0.0, L}};
    lo = std::numeric_limits<double>::infinity();
    hi = -lo;
    for (const Vec2 c : corners) {
        const double s = dot(c, axis);
        lo = std::min(lo, s);
        hi = std::max(hi, s);
    }
}

}  // namespace

const char* region_name(Region r) {
    switch (r) {
        case Region::Disc: return "disc";
        case Region::Band: return "band";
        case Region::Corner: return "corner";
    }
    return "?";
}

double contact_pitch(const Shape& s) {
    if (s.vertices.empty() || s.sat_axes.empty()) return 1.0;
    const Vec2 a = s.sat_axes[0];
    double lo = dot(s.vertices[0], a), hi = lo;
    for (size_t k = 1; k < s.vertices.size(); ++k) {
        const double d = dot(s.vertices[k], a);
        lo = std::min(lo, d);
        hi = std::max(hi, d);
    }
    return hi - lo;
}

std::vector<int> select_disc(const System& sys, double cx, double cy, double r) {
    std::vector<int> out;
    if (!(r > 0.0)) return out;
    const double r2 = r * r;
    for (size_t i = 0; i < sys.particles.size(); ++i) {
        const Particle& p = sys.particles[i];
        if (p.frozen) continue;
        const double dx = p.x() - cx, dy = p.y() - cy;
        if (dx * dx + dy * dy <= r2) out.push_back(static_cast<int>(i));
    }
    return out;
}

std::vector<int> select_band(const System& sys, Vec2 dir, double offset, double width) {
    std::vector<int> out;
    if (!(width > 0.0) || !(norm(dir) > 0.0)) return out;
    const Vec2 n = perp(normalized(dir));
    for (size_t i = 0; i < sys.particles.size(); ++i) {
        const Particle& p = sys.particles[i];
        if (p.frozen) continue;
        const double s = p.x() * n.x + p.y() * n.y;
        if (s >= offset && s <= offset + width) out.push_back(static_cast<int>(i));
    }
    return out;
}

std::vector<int> select_corner(const System& sys, int corner, double depth) {
    std::vector<int> out;
    if (!(depth > 0.0) || corner < 0 || corner > 3) return out;
    const double cx = (corner == 1 || corner == 2) ? sys.L : 0.0;
    const double cy = (corner == 2 || corner == 3) ? sys.L : 0.0;
    for (size_t i = 0; i < sys.particles.size(); ++i) {
        const Particle& p = sys.particles[i];
        if (p.frozen) continue;
        if (std::abs(p.x() - cx) <= depth && std::abs(p.y() - cy) <= depth) {
            out.push_back(static_cast<int>(i));
        }
    }
    return out;
}

Pick pick_band(const System& sys, const RegionWeights& w, Vec2 dir, std::mt19937_64& rng) {
    Pick p;
    p.kind = Region::Band;
    p.dir = normalized(dir);
    p.axis = p.dir;  // a band is symmetric in its own midline, and the midline runs along dir

    const double width = uniform(rng, w.band_width_lo, w.band_width_hi) * diameter(sys);
    const Vec2 n = perp(p.dir);
    double lo = 0.0, hi = 0.0;
    box_span(sys.L, n, lo, hi);

    // Keep the whole slab over the box where it can, so a draw is not wasted on a band that mostly hangs outside. A box narrower than the band gives one choice.
    const double offset = (hi - lo > width) ? uniform(rng, lo, hi - width) : lo;
    const double mid = offset + 0.5 * width;

    // Any point on the midline will do; the box centre slid onto it is the one that stays near the squares, which matters because a reflection is taken about this point.
    const Vec2 c0{0.5 * sys.L, 0.5 * sys.L};
    p.centre = c0 + (mid - dot(c0, n)) * n;
    p.extent = width;
    p.members = select_band(sys, p.dir, offset, width);
    return p;
}

Pick pick_region(const System& sys, const RegionWeights& w, std::mt19937_64& rng) {
    const double total = std::max(1e-300, w.disc + w.band + w.corner);
    const double roll = uniform(rng, 0.0, total);

    if (roll < w.disc) {
        Pick p;
        p.kind = Region::Disc;
        p.centre = {uniform(rng, 0.0, sys.L), uniform(rng, 0.0, sys.L)};
        p.extent = uniform(rng, w.disc_radius_lo, w.disc_radius_hi) * diameter(sys);
        // A disc is symmetric in every line through its centre, so the axis is a free draw over a half turn -- a line at a and at a + pi are the same line.
        const double a = uniform(rng, 0.0, kPi);
        p.axis = {std::cos(a), std::sin(a)};
        p.dir = p.axis;
        p.members = select_disc(sys, p.centre.x, p.centre.y, p.extent);
        return p;
    }

    if (roll < w.disc + w.band) {
        // The four lattice directions, diagonals weighted up: see RegionWeights.
        const double aw = std::max(0.0, w.axis_weight), dw = std::max(0.0, w.diagonal_weight);
        const double s = std::max(1e-300, 2.0 * (aw + dw));
        const double r = uniform(rng, 0.0, s);
        const double k = 1.0 / std::sqrt(2.0);
        Vec2 dir{1.0, 0.0};
        if (r < aw) dir = {1.0, 0.0};
        else if (r < aw + dw) dir = {k, k};
        else if (r < 2.0 * aw + dw) dir = {0.0, 1.0};
        else dir = {k, -k};
        return pick_band(sys, w, dir, rng);
    }

    Pick p;
    p.kind = Region::Corner;
    p.corner = static_cast<int>(rng() % 4u);
    p.centre = {(p.corner == 1 || p.corner == 2) ? sys.L : 0.0,
                (p.corner == 2 || p.corner == 3) ? sys.L : 0.0};
    p.extent = uniform(rng, w.corner_depth_lo, w.corner_depth_hi) * diameter(sys);
    // The block's symmetry line is its own diagonal, pointing from the corner into the box.
    const double k = 1.0 / std::sqrt(2.0);
    const double sx = (p.corner == 1 || p.corner == 2) ? -1.0 : 1.0;
    const double sy = (p.corner == 2 || p.corner == 3) ? -1.0 : 1.0;
    p.axis = {sx * k, sy * k};
    p.dir = p.axis;
    p.members = select_corner(sys, p.corner, p.extent);
    return p;
}

Bounds member_bounds(const System& sys, const std::vector<int>& members, double pad) {
    Bounds b;
    if (members.empty()) return b;
    const double inf = std::numeric_limits<double>::infinity();
    b = {inf, -inf, inf, -inf};
    for (const int i : members) {
        const Particle& p = sys.particles[static_cast<size_t>(i)];
        b.xlo = std::min(b.xlo, p.x());
        b.xhi = std::max(b.xhi, p.x());
        b.ylo = std::min(b.ylo, p.y());
        b.yhi = std::max(b.yhi, p.y());
    }
    b.xlo = std::max(0.0, b.xlo - pad);
    b.ylo = std::max(0.0, b.ylo - pad);
    b.xhi = std::min(sys.L, b.xhi + pad);
    b.yhi = std::min(sys.L, b.yhi + pad);
    // A box smaller than the padding can cross the bounds over; an empty interval is a valid answer that the callers draw a single point from, but an inverted one is not.
    if (b.xhi < b.xlo) b.xhi = b.xlo;
    if (b.yhi < b.ylo) b.yhi = b.ylo;
    return b;
}

}  // namespace pack
