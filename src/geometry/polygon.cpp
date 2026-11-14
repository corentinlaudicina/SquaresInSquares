#include "geometry/polygon.hpp"

#include <algorithm>
#include <cassert>

namespace pack {
namespace {  // anonymous namespace (replaces the static declarations of C)

// ---------------------------------------------------------------------------
// Sutherland-Hodgman clipping.
// Idea: the intersection of two convex polygons is what remains of `sub` after
// cutting it with each edge of `clip`, treating each edge as an infinite line
// and keeping the "inside" half-plane. Both polygons must be counter-clockwise.
// Returns the area of what is left.
// ---------------------------------------------------------------------------
double clipped_area(const Vec2* sub, int n, const Vec2* clip, int m) {

    Vec2 buf[2][2 * kMaxVertices + 2];
    int cur = 0, count = n; // cur = buffer holding the current polygon

    for (int i = 0; i < n; ++i) buf[0][i] = sub[i];

    // One pass per clip edge. Stop early if nothing is left.
    for (int e = 0; e < m && count > 0; ++e) {
        const Vec2 p0 = clip[e];
        const Vec2 p1 = clip[(e + 1) % m];           // % m wraps the last edge to vertex 0

        // Outward normal of the edge. For counter-clockwise order the interior is on the left, so rotating the edge direction (dx,dy) to (dy,-dx) points out.

        const Vec2 outward{p1.y - p0.y, p0.x - p1.x};

        const Vec2* in = buf[cur];                   // polygon to clip
        Vec2* out = buf[1 - cur];                    // where the clipped polygon goes
        int k = 0;                                   // vertices written to `out`

        // Walk every edge a -> b of the current polygon.
        for (int i = 0; i < count; ++i) {
            const Vec2 a = in[i];
            const Vec2 b = in[(i + 1) % count];
            // Signed distance (up to scale) from the clip line: <= 0 means inside.
            const double da = dot(a - p0, outward);
            const double db = dot(b - p0, outward);

            if (da <= 0.0) out[k++] = a;             // keep a if it is inside
            if ((da > 0.0) != (db > 0.0)) {          // edge crosses the line
                const double t = da / (da - db);     // fraction along a->b at the crossing
                out[k++] = a + t * (b - a);          // add the crossing point
            }
        }
        assert(k <= 2 * kMaxVertices);
        count = k;
        cur = 1 - cur;                               // swap buffer roles
    }

    if (count < 3) return 0.0;                       // fewer than 3 vertices: no area

    // Shoelace formula: area = 1/2 |sum of cross(v_i, v_{i+1})|.
    const Vec2* v = buf[cur];
    double shoelace = 0.0;
    for (int i = 0; i < count; ++i) shoelace += cross(v[i], v[(i + 1) % count]);
    return 0.5 * std::abs(shoelace);
}

// ---------------------------------------------------------------------------
// Which part of segment p -> q lies inside the convex polygon `clip`?
// Answer as a parameter range: points p + t*(q-p) for t in [t0, t1].
// Returns false if no part survives. (Cyrus-Beck algorithm.)
// ---------------------------------------------------------------------------
bool clip_segment(Vec2 p, Vec2 q, const Vec2* clip, int m, double& t0, double& t1) {
    const Vec2 d = q - p;
    t0 = 0.0; // start with the whole segment
    t1 = 1.0;
    for (int e = 0; e < m; ++e) {
        const Vec2 p0 = clip[e];
        const Vec2 p1 = clip[(e + 1) % m];
        const Vec2 outward{p1.y - p0.y, p0.x - p1.x};

        // Inside this edge means dot(point - p0, outward) <= 0. For the point
        // p + t*d that is   c + t*rate <= 0   with:
        const double c = dot(p - p0, outward);       // value at t = 0
        const double rate = dot(d, outward);         // change per unit t
        if (rate == 0.0) {                           // segment parallel to this edge
            if (c > 0.0) return false;               // ...and entirely outside it
            continue;                                // ...or entirely inside:  no constraint
        }
        const double hit = -c / rate;                // t where the segment crosses the line
        if (rate > 0.0) {
            t1 = std::min(t1, hit);                  // moving outward: it LEAVES here
        } else {
            t0 = std::max(t0, hit);                  // moving inward: it ENTERS here
        }
        if (t0 >= t1) return false;                  // nothing left
    }
    return true;
}

// ---------------------------------------------------------------------------
// Adds one polygon edge's contribution to the pose gradient, counting only the
// part of the edge that lies inside `clip` (the other body or the box).
//   (cx, cy) is the polygon's centre, the pivot for the torque.
// ---------------------------------------------------------------------------
void add_edge_grad(Vec2 p, Vec2 q, const Vec2* clip, int m, double cx, double cy, Pose& g) {
    double t0, t1;
    if (!clip_segment(p, q, clip, m, t0, t1)) return;   // edge not inside: no contribution

    // Outward normal (dy, -dx). Its length equals the edge length, which is
    // exactly the arc-length factor of the integral, so no sqrt is needed.
    const Vec2 d = q - p;
    const Vec2 n{d.y, -d.x};
    const double w = t1 - t0;                        // fraction of the edge that counts

    g.x += n.x * w;                                  // force: integral of the normal
    g.y += n.y * w;

    // Torque: integral of cross(x - centre, normal). That quantity is linear
    // along the edge, so its average over a sub-segment equals its value at the
    // midpoint, and the integral is (midpoint value) * (length factor).
    const Vec2 mid = p + (0.5 * (t0 + t1)) * d;
    g.theta += cross({mid.x - cx, mid.y - cy}, n) * w;
}

}  // namespace (anonymous)

// ---------------------------------------------------------------------------
// Shape factories
// ---------------------------------------------------------------------------
Shape make_square(double side) {
    assert(side > 0.0);
    const double h = 0.5 * side;                     // half side
    Shape s;
    s.vertices = {{-h, -h}, {h, -h}, {h, h}, {-h, h}};  // counter-clockwise
    s.sat_axes = {{1.0, 0.0}, {0.0, 1.0}};           // opposite edges share an axis
    s.circumradius = h * std::sqrt(2.0);             // half the diagonal
    s.area = side * side;
    return s;
}

// An axis-aligned rectangle, centred at origin, theta = 0 -- the shape merge() builds for a
// group of touching squares whose union tiles an exact w-by-h block. Same field layout as
// make_square, generalized to two independent extents.
Shape make_rectangle(double width, double height) {
    assert(width > 0.0 && height > 0.0);
    const double hw = 0.5 * width, hh = 0.5 * height;
    Shape s;
    s.vertices = {{-hw, -hh}, {hw, -hh}, {hw, hh}, {-hw, hh}};  // counter-clockwise
    s.sat_axes = {{1.0, 0.0}, {0.0, 1.0}};
    s.circumradius = std::sqrt(hw * hw + hh * hh);  // half the diagonal
    s.area = width * height;
    return s;
}

std::optional<Shape> shape_from_vertices(const std::vector<Vec2>& vertices) {
    const size_t n = vertices.size();
    if (n < 3 || n > static_cast<size_t>(kMaxVertices)) return std::nullopt;

    Shape s;
    s.vertices = vertices;
    for (size_t i = 0; i < n; ++i) {
        const Vec2 p = vertices[i], edge = vertices[(i + 1) % n] - p;
        // Every other corner strictly left of this edge. Checking only consecutive turns would
        // pass a pentagram, whose turns all go the same way.
        for (size_t k = 2; k < n; ++k) {
            if (!(cross(edge, vertices[(i + k) % n] - p) > 0.0)) return std::nullopt;
        }

        const Vec2 axis = normalized({edge.y, -edge.x});
        bool seen = false;
        for (const Vec2 a : s.sat_axes) seen = seen || std::abs(cross(a, axis)) < 1e-12;
        if (!seen) s.sat_axes.push_back(axis);

        s.area += 0.5 * cross(p, vertices[(i + 1) % n]);
        s.circumradius = std::max(s.circumradius, norm(p));
    }
    return s;
}

Shape make_regular_polygon(int n, double circumradius) {
    assert(n >= 3 && n <= kMaxVertices && circumradius > 0.0);
    Shape s;
    s.vertices.reserve(static_cast<size_t>(n));
    // n points evenly spaced on a circle, counter-clockwise from angle 0.
    for (int k = 0; k < n; ++k) {
        const double a = 2.0 * kPi * k / n;
        s.vertices.push_back({circumradius * std::cos(a), circumradius * std::sin(a)});
    }
    // Even n: opposite edges are parallel (same axis up to sign), so half the
    // normals suffice. Odd n: every edge gives a distinct axis.
    const int n_axes = (n % 2 == 0) ? n / 2 : n;
    s.sat_axes.reserve(static_cast<size_t>(n_axes));
    for (int k = 0; k < n_axes; ++k) {
        const Vec2 d = s.vertices[static_cast<size_t>((k + 1) % n)] -
                       s.vertices[static_cast<size_t>(k)];     // edge direction
        s.sat_axes.push_back(normalized({d.y, -d.x}));         // outward unit normal
    }
    s.circumradius = circumradius;
    s.area = 0.5 * n * circumradius * circumradius * std::sin(2.0 * kPi / n);  // n triangles
    return s;
}

// ---------------------------------------------------------------------------
// Boolean overlap: separating-axis theorem (SAT).
// Two convex shapes are disjoint if and only if some axis exists on which their
// shadows do not overlap. For polygons it is enough to test the edge normals
// of both shapes. One failing axis proves "no overlap" immediately.
// ---------------------------------------------------------------------------
bool overlaps(const Particle& a, const Shape& sa, const Particle& b, const Shape& sb) {
    // Cheap early exit: if the centres are farther apart than the sum of the two circumradii,
    // the bounding circles are disjoint, so the polygons are too. (Compares squared distances
    // to avoid a sqrt.) Reduces to 2 * s.circumradius when sa == sb, today's only case.
    const double dx = b.x() - a.x(), dy = b.y() - a.y();
    const double reach = sa.circumradius + sb.circumradius;
    if (dx * dx + dy * dy > reach * reach) return false;

    // Test the axes of both polygons: first a's own axes, then b's own.
    for (int k = 0; k < static_cast<int>(sa.sat_axes.size()); ++k) {
        const Vec2 axis = world_axis(a, sa, k);
        double alo, ahi, blo, bhi;
        project(a, sa, axis, alo, ahi);               // a's shadow
        project(b, sb, axis, blo, bhi);               // b's shadow
        // Overlap length of the two intervals. <= 0 means a gap (or just
        // touching): a separating axis exists, so no overlap.
        if (std::min(ahi, bhi) - std::max(alo, blo) <= 0.0) return false;
    }
    for (int k = 0; k < static_cast<int>(sb.sat_axes.size()); ++k) {
        const Vec2 axis = world_axis(b, sb, k);
        double alo, ahi, blo, bhi;
        project(a, sa, axis, alo, ahi);
        project(b, sb, axis, blo, bhi);
        if (std::min(ahi, bhi) - std::max(alo, blo) <= 0.0) return false;
    }
    return true;                                     // no separating axis found
}

// ---------------------------------------------------------------------------
// Overlap area: clip one polygon by the other and measure what remains.
// ---------------------------------------------------------------------------
double overlap_area(const Particle& a, const Shape& sa, const Particle& b, const Shape& sb) {
    const double dx = b.x() - a.x(), dy = b.y() - a.y();
    const double reach = sa.circumradius + sb.circumradius;
    if (dx * dx + dy * dy > reach * reach) return 0.0;   // same early exit as above

    Vec2 va[kMaxVertices], vb[kMaxVertices];         // stack buffers for world vertices
    const int na = world_vertices(a, sa, va);
    const int nb = world_vertices(b, sb, vb);
    return clipped_area(va, na, vb, nb);
}

// ---------------------------------------------------------------------------
// Overlap area plus its gradient. Each body's gradient is accumulated from its
// own edges, restricted to the parts inside the other body.
// ---------------------------------------------------------------------------
double overlap_area_grad(const Particle& a, const Shape& sa, const Particle& b, const Shape& sb,
                         Pose* ga, Pose* gb) {
    if (ga) *ga = Pose{};                        // zero the outputs first
    if (gb) *gb = Pose{};

    const double dx = b.x() - a.x(), dy = b.y() - a.y();
    const double reach = sa.circumradius + sb.circumradius;
    if (dx * dx + dy * dy > reach * reach) return 0.0;  // disjoint: area 0, gradient 0

    Vec2 va[kMaxVertices], vb[kMaxVertices];
    const int na = world_vertices(a, sa, va);
    const int nb = world_vertices(b, sb, vb);

    if (ga) {   // a's edges that lie inside b
        for (int e = 0; e < na; ++e) {
            add_edge_grad(va[e], va[(e + 1) % na], vb, nb, a.x(), a.y(), *ga);
        }
    }
    if (gb) {   // b's edges that lie inside a
        for (int e = 0; e < nb; ++e) {
            add_edge_grad(vb[e], vb[(e + 1) % nb], va, na, b.x(), b.y(), *gb);
        }
    }
    return clipped_area(va, na, vb, nb);
}

// ---------------------------------------------------------------------------
// Area outside the box [0,L]^2.
// ---------------------------------------------------------------------------
double outside_area(const Particle& p, const Shape& s, double L) {
    assert(L > 0.0);

    // Quick test: project onto x and y. If the shadows fit in [0, L] the
    // polygon is entirely inside, with no clipping needed.
    double xlo, xhi, ylo, yhi;
    project(p, s, {1.0, 0.0}, xlo, xhi);
    project(p, s, {0.0, 1.0}, ylo, yhi);
    if (xlo >= 0.0 && ylo >= 0.0 && xhi <= L && yhi <= L) return 0.0;

    Vec2 v[kMaxVertices];
    const int n = world_vertices(p, s, v);
    const Vec2 box[4] = {{0.0, 0.0}, {L, 0.0}, {L, L}, {0.0, L}};  // counter-clockwise

    // outside = total area - area inside the box. max(0, ...) guards against
    // tiny negative values from rounding error.
    return std::max(0.0, s.area - clipped_area(v, n, box, 4));
}

// ---------------------------------------------------------------------------
// Area outside the box, plus its gradient.
// ---------------------------------------------------------------------------
double outside_area_grad(const Particle& p, const Shape& s, double L, Pose* g) {
    assert(L > 0.0);
    if (g) *g = Pose{};

    double xlo, xhi, ylo, yhi;
    project(p, s, {1.0, 0.0}, xlo, xhi);
    project(p, s, {0.0, 1.0}, ylo, yhi);
    if (xlo >= 0.0 && ylo >= 0.0 && xhi <= L && yhi <= L) return 0.0;  // inside: area and gradient are 0

    Vec2 v[kMaxVertices];
    const int n = world_vertices(p, s, v);
    const Vec2 box[4] = {{0.0, 0.0}, {L, 0.0}, {L, L}, {0.0, L}};

    if (g) {
        // Gradient of the INSIDE area: the polygon's edges lying inside the box.
        Pose inside;
        for (int e = 0; e < n; ++e) {
            add_edge_grad(v[e], v[(e + 1) % n], box, 4, p.x(), p.y(), inside);
        }
        // outside = constant - inside, so the gradient flips sign.
        g->x = -inside.x;
        g->y = -inside.y;
        g->theta = -inside.theta;
    }
    return std::max(0.0, s.area - clipped_area(v, n, box, 4));
}

}  // namespace pack