// Assert-based tests for the geometry layer. No framework: nothing to install.
// Returns 0 on success, aborts on the first failed assert.

#include "geometry/polygon.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>

using namespace pack;

namespace {

bool close(double a, double b, double tol = 1e-12) { return std::abs(a - b) <= tol; }

double shoelace(const std::vector<Vec2>& v) {
    double s = 0.0;
    const size_t n = v.size();
    for (size_t i = 0; i < n; ++i) s += cross(v[i], v[(i + 1) % n]);
    return 0.5 * s;  // positive iff counter-clockwise
}

// Independent of the clipper: a point is inside a CCW convex polygon iff it is left of every
// edge. Uses cross products only, sharing no logic with the half-plane clipping code.
bool point_inside(const Vec2* v, int n, Vec2 q) {
    for (int i = 0; i < n; ++i) {
        if (cross(v[(i + 1) % n] - v[i], q - v[i]) < 0.0) return false;
    }
    return true;
}

// Monte Carlo estimate of the intersection area, by sampling a's bounding box.
double sampled_overlap_area(const Particle& a, const Particle& b, const Shape& s,
                            std::mt19937_64& rng, int samples) {
    Vec2 va[kMaxVertices], vb[kMaxVertices];
    const int na = world_vertices(a, s, va);
    const int nb = world_vertices(b, s, vb);

    double xlo, xhi, ylo, yhi;
    project(a, s, {1.0, 0.0}, xlo, xhi);
    project(a, s, {0.0, 1.0}, ylo, yhi);

    std::uniform_real_distribution<double> ux(xlo, xhi), uy(ylo, yhi);
    int hits = 0;
    for (int i = 0; i < samples; ++i) {
        const Vec2 q{ux(rng), uy(rng)};
        if (point_inside(va, na, q) && point_inside(vb, nb, q)) ++hits;
    }
    return (xhi - xlo) * (yhi - ylo) * hits / samples;
}

// project() is the most reused primitive here and the wall check depends on it, so test it
// directly rather than only through the measures built on top of it.
void test_project() {
    const Shape s = make_square(1.0);
    const double h = 0.5, r = std::sqrt(2.0) / 2.0;
    Particle p;
    double lo, hi;

    // Axis-aligned and centred: half the side either way.
    set_pose(p, 0.0, 0.0, 0.0);
    project(p, s, {1.0, 0.0}, lo, hi);
    assert(close(lo, -h) && close(hi, h));
    project(p, s, {0.0, 1.0}, lo, hi);
    assert(close(lo, -h) && close(hi, h));

    // Turned 45 degrees: the corners now lead, so the reach is the circumradius.
    set_pose(p, 0.0, 0.0, kPi / 4.0);
    project(p, s, {1.0, 0.0}, lo, hi);
    assert(close(lo, -r, 1e-15) && close(hi, r, 1e-15));

    // A diagonal axis on the unturned square reaches the circumradius too.
    set_pose(p, 0.0, 0.0, 0.0);
    project(p, s, {r, r}, lo, hi);
    assert(close(lo, -r, 1e-15) && close(hi, r, 1e-15));

    // Translation shifts the interval by the centre's projection, width unchanged.
    set_pose(p, 2.0, -3.0, 0.9);
    project(p, s, {1.0, 0.0}, lo, hi);
    assert(close(0.5 * (lo + hi), 2.0, 1e-15));
    const double width = hi - lo;
    set_pose(p, 7.0, -3.0, 0.9);
    project(p, s, {1.0, 0.0}, lo, hi);
    assert(close(0.5 * (lo + hi), 7.0, 1e-15));
    assert(close(hi - lo, width, 1e-15));
}

void test_make_square() {
    const Shape s = make_square(1.0);
    assert(s.vertices.size() == 4);
    assert(close(s.area, 1.0));
    assert(close(s.circumradius, std::sqrt(2.0) / 2.0));
    assert(shoelace(s.vertices) > 0.0);          // counter-clockwise
    assert(close(shoelace(s.vertices), s.area)); // and the stored area agrees
    assert(s.sat_axes.size() == 2);              // deduplicated
    for (const Vec2& a : s.sat_axes) assert(close(norm(a), 1.0));

    const Shape big = make_square(3.0);
    assert(close(big.area, 9.0));
    assert(close(big.circumradius, 1.5 * std::sqrt(2.0)));
}

void test_make_rectangle() {
    const Shape r = make_rectangle(3.0, 1.0);
    assert(r.vertices.size() == 4);
    assert(close(r.area, 3.0));
    assert(close(r.circumradius, std::sqrt(9.0 + 1.0) / 2.0));
    assert(shoelace(r.vertices) > 0.0);
    assert(close(shoelace(r.vertices), r.area));
    assert(r.sat_axes.size() == 2);

    // Degenerates to make_square when both extents agree.
    const Shape sq = make_rectangle(2.0, 2.0);
    const Shape ref = make_square(2.0);
    assert(close(sq.area, ref.area));
    assert(close(sq.circumradius, ref.circumradius));
}

void test_make_regular_polygon() {
    const double R = 1.7;
    const Shape hex = make_regular_polygon(6, R);
    assert(hex.vertices.size() == 6);
    assert(hex.sat_axes.size() == 3);  // opposite edges share an axis
    assert(close(hex.area, 1.5 * std::sqrt(3.0) * R * R, 1e-12));
    assert(close(shoelace(hex.vertices), hex.area, 1e-12));
    for (const Vec2& a : hex.sat_axes) assert(close(norm(a), 1.0));

    const Shape pent = make_regular_polygon(5, R);
    assert(pent.sat_axes.size() == 5);  // odd: no axis is shared
    assert(close(shoelace(pent.vertices), pent.area, 1e-12));

    // Area approaches that of the circumscribed circle as n grows.
    const Shape fine = make_regular_polygon(kMaxVertices, R);
    assert(std::abs(fine.area - kPi * R * R) < 0.03 * kPi * R * R);
}

void test_analytic_overlap() {
    const Shape s = make_square(1.0);
    Particle a, b;
    set_pose(a, 0.0, 0.0, 0.0);

    // Axis-aligned unit squares offset by (dx, dy) overlap on (1-|dx|)(1-|dy|).
    const double offsets[] = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0, 1.5};
    for (double dx : offsets) {
        for (double dy : offsets) {
            set_pose(b, dx, dy, 0.0);
            const double want = std::max(0.0, 1.0 - dx) * std::max(0.0, 1.0 - dy);
            assert(close(overlap_area(a, s, b, s), want, 1e-12));
            set_pose(b, -dx, -dy, 0.0);
            assert(close(overlap_area(a, s, b, s), want, 1e-12));
        }
    }

    // Fully coincident: the whole area.
    set_pose(b, 0.0, 0.0, 0.0);
    assert(close(overlap_area(a, s, b, s), 1.0, 1e-12));

    // Far apart: exactly zero, via the circumradius reject.
    set_pose(b, 5.0, 5.0, 0.7);
    assert(overlap_area(a, s, b, s) == 0.0);
    assert(!overlaps(a, s, b, s));

    // Concentric, one turned 45 degrees: a regular octagon of area 2(sqrt(2)-1).
    set_pose(b, 0.0, 0.0, kPi / 4.0);
    assert(close(overlap_area(a, s, b, s), 2.0 * (std::sqrt(2.0) - 1.0), 1e-12));
}

void test_rigid_invariance() {
    const Shape s = make_square(1.0);
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> upos(-1.2, 1.2), uang(-kPi, kPi);

    for (int trial = 0; trial < 500; ++trial) {
        Particle a, b;
        set_pose(a, upos(rng), upos(rng), uang(rng));
        set_pose(b, upos(rng), upos(rng), uang(rng));
        const double reference = overlap_area(a, s, b, s);

        // Move both by the same rigid transform: rotate by phi about c, then translate by t.
        const double phi = uang(rng), ct = std::cos(phi), st = std::sin(phi);
        const Vec2 c{upos(rng), upos(rng)}, t{upos(rng), upos(rng)};
        Particle a2, b2;
        for (const Particle* src : {&a, &b}) {
            const Vec2 d{src->x() - c.x, src->y() - c.y};
            Particle& dst = (src == &a) ? a2 : b2;
            set_pose(dst, c.x + d.x * ct - d.y * st + t.x,
                     c.y + d.x * st + d.y * ct + t.y, src->theta() + phi);
        }
        assert(close(overlap_area(a2, s, b2, s), reference, 1e-10));
    }
}

// The important one: validate the clipper against point sampling, which shares none of its
// logic. Tolerance is set by the sampler's own statistics, not by the clipper's accuracy.
void test_clipper_against_sampling() {
    const Shape s = make_square(1.0);
    std::mt19937_64 rng(777);
    std::uniform_real_distribution<double> upos(-1.1, 1.1), uang(-kPi, kPi);

    const int samples = 20000;
    const double tol = 0.05;  // ~7 sigma at these sample counts; a broken clipper is off by ~0.1
    int tested = 0;
    for (int trial = 0; trial < 120; ++trial) {
        Particle a, b;
        set_pose(a, 0.0, 0.0, uang(rng));
        set_pose(b, upos(rng), upos(rng), uang(rng));

        const double exact = overlap_area(a, s, b, s);
        const double estimate = sampled_overlap_area(a, b, s, rng, samples);
        assert(std::abs(exact - estimate) < tol);
        if (exact > 1e-6) ++tested;
    }
    // Guard against the test passing because every pair happened to be disjoint.
    assert(tested > 60);

    // Same cross-check for non-square shapes, exercising the many-vertex path.
    const Shape hex = make_regular_polygon(6, 0.6);
    for (int trial = 0; trial < 60; ++trial) {
        Particle a, b;
        set_pose(a, 0.0, 0.0, uang(rng));
        set_pose(b, upos(rng) * 0.6, upos(rng) * 0.6, uang(rng));
        const double exact = overlap_area(a, hex, b, hex);
        const double estimate = sampled_overlap_area(a, b, hex, rng, samples);
        assert(std::abs(exact - estimate) < tol);
    }
}

void test_overlaps_agrees_with_area() {
    const Shape s = make_square(1.0);
    std::mt19937_64 rng(2024);
    std::uniform_real_distribution<double> upos(-1.6, 1.6), uang(-kPi, kPi);

    for (int trial = 0; trial < 20000; ++trial) {
        Particle a, b;
        set_pose(a, 0.0, 0.0, uang(rng));
        set_pose(b, upos(rng), upos(rng), uang(rng));
        const double area = overlap_area(a, s, b, s);
        // Near-contact pairs are skipped on purpose: the boolean test and the area disagree
        // inside a tolerance band there, which is expected rather than a defect.
        if (area > 1e-9) assert(overlaps(a, s, b, s));
    }
}

// Pins the sa.circumradius + sb.circumradius early-reject fix: a pair whose centre distance
// exceeds either shape's OWN circumradius doubled, but not the sum of the two, has to reach the
// clipper rather than being rejected early on the smaller shape's reach alone.
void test_overlap_between_two_different_shapes() {
    const Shape sa = make_square(1.0);           // circumradius sqrt(2)/2 ~ 0.7071
    const Shape sb = make_rectangle(5.0, 1.0);    // circumradius sqrt(26)/2 ~ 2.5495
    assert(sa.circumradius + sb.circumradius > 3.2);
    assert(2.0 * sa.circumradius < 1.5);  // a reach built from sa alone would wrongly reject below

    Particle a, b;
    set_pose(a, 0.0, 0.0, 0.0);
    set_pose(b, 2.0, 0.0, 0.0);  // centre distance 2.0: past 2*sa.circumradius, well within the sum
    assert(overlaps(a, sa, b, sb));
    // a's whole x-range [-0.5, 0.5] sits inside b's [-0.5, 4.5], and both share y in [-0.5, 0.5],
    // so a is entirely covered: the overlap is a's whole area.
    assert(close(overlap_area(a, sa, b, sb), sa.area, 1e-12));

    // Symmetric in argument order.
    assert(close(overlap_area(b, sb, a, sa), sa.area, 1e-12));

    // Far enough apart that even the correct, larger reach rejects them.
    set_pose(b, 20.0, 0.0, 0.0);
    assert(!overlaps(a, sa, b, sb));
    assert(overlap_area(a, sa, b, sb) == 0.0);
}

void test_outside_area() {
    const Shape s = make_square(1.0);
    const double L = 10.0;
    Particle p;

    // Well inside: exactly zero, via the AABB gate.
    set_pose(p, 5.0, 5.0, 0.3);
    assert(outside_area(p, s, L) == 0.0);

    // Centred on the corner: three quarters stick out.
    set_pose(p, 0.0, 0.0, 0.0);
    assert(close(outside_area(p, s, L), 0.75, 1e-12));

    // Same, turned 45 degrees: still three quarters by symmetry.
    set_pose(p, 0.0, 0.0, kPi / 4.0);
    assert(close(outside_area(p, s, L), 0.75, 1e-12));

    // One wall only: a strip of width 0.25 over the full height.
    set_pose(p, 0.25, 5.0, 0.0);
    assert(close(outside_area(p, s, L), 0.25, 1e-12));
    set_pose(p, L - 0.25, 5.0, 0.0);
    assert(close(outside_area(p, s, L), 0.25, 1e-12));
    set_pose(p, 5.0, 0.25, 0.0);
    assert(close(outside_area(p, s, L), 0.25, 1e-12));

    // Entirely outside: the whole area.
    set_pose(p, -3.0, 5.0, 0.4);
    assert(close(outside_area(p, s, L), 1.0, 1e-12));

    // A box smaller than the square: everything but the box sticks out.
    set_pose(p, 0.25, 0.25, 0.0);
    assert(close(outside_area(p, s, 0.5), 0.75, 1e-12));

    // Cross-check the partial cases against point sampling of the shape itself.
    std::mt19937_64 rng(99);
    std::uniform_real_distribution<double> upos(-0.8, 1.8), uang(-kPi, kPi);
    const double small_L = 1.0;
    for (int trial = 0; trial < 200; ++trial) {
        set_pose(p, upos(rng), upos(rng), uang(rng));
        Vec2 v[kMaxVertices];
        const int n = world_vertices(p, s, v);
        double xlo, xhi, ylo, yhi;
        project(p, s, {1.0, 0.0}, xlo, xhi);
        project(p, s, {0.0, 1.0}, ylo, yhi);
        std::uniform_real_distribution<double> ux(xlo, xhi), uy(ylo, yhi);
        const int samples = 20000;
        int hits = 0;
        for (int i = 0; i < samples; ++i) {
            const Vec2 q{ux(rng), uy(rng)};
            const bool in_box = q.x >= 0.0 && q.y >= 0.0 && q.x <= small_L && q.y <= small_L;
            if (point_inside(v, n, q) && !in_box) ++hits;
        }
        const double estimate = (xhi - xlo) * (yhi - ylo) * hits / samples;
        assert(std::abs(outside_area(p, s, small_L) - estimate) < 0.05);
    }
}

void test_shape_from_vertices() {
    // It has to agree with the factories, because a loaded square is used interchangeably with
    // make_square's.
    const Shape want = make_regular_polygon(6, 1.5);
    const std::optional<Shape> got = shape_from_vertices(want.vertices);
    assert(got);
    assert(close(got->area, want.area) && close(got->circumradius, want.circumradius));
    assert(got->sat_axes.size() == 3);  // opposite edges share an axis

    const std::optional<Shape> square = shape_from_vertices(make_square(1.0).vertices);
    assert(square && square->sat_axes.size() == 2 && close(square->area, 1.0));
    Particle a, b;
    b.set_pose(0.7, 0.2, 0.4);
    assert(close(overlap_area(a, *square, b, *square),
                 overlap_area(a, make_square(1.0), b, make_square(1.0))));

    // Everything the clipper cannot be trusted with.
    assert(!shape_from_vertices({}));
    assert(!shape_from_vertices({{0, 0}, {1, 0}}));
    assert(!shape_from_vertices({{0, 0}, {0, 1}, {1, 1}, {1, 0}}));              // clockwise
    assert(!shape_from_vertices({{0, 0}, {2, 0}, {0.5, 0.5}, {0, 2}}));          // dent
    assert(!shape_from_vertices({{0, 0}, {1, 0}, {2, 0}, {1, 1}}));              // three in line
    const double nan = std::nan("");
    assert(!shape_from_vertices({{0, 0}, {1, 0}, {nan, 1}}));

    std::vector<Vec2> pentagram;  // every turn goes the same way, and it is not convex
    for (int k = 0; k < 5; ++k) {
        const double t = 2.0 * kPi * (2 * k) / 5.0;
        pentagram.push_back({std::cos(t), std::sin(t)});
    }
    assert(!shape_from_vertices(pentagram));

    std::vector<Vec2> too_many = make_regular_polygon(kMaxVertices, 1.0).vertices;
    assert(shape_from_vertices(too_many));
    too_many.insert(too_many.begin(), {1.0, -1e-3});
    assert(!shape_from_vertices(too_many));
}

}  // namespace

int main() {
    test_shape_from_vertices();
    test_project();
    test_make_square();
    test_make_rectangle();
    test_make_regular_polygon();
    test_analytic_overlap();
    test_rigid_invariance();
    test_clipper_against_sampling();
    test_overlaps_agrees_with_area();
    test_overlap_between_two_different_shapes();
    test_outside_area();
    std::printf("polygon: all tests passed\n");
    return 0;
}
