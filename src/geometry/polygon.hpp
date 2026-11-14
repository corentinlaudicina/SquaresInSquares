#pragma once

// (Your original file-level comment explains the design; kept as is.)

#include <cassert>
#include <cmath>
#include <optional>
#include <vector>

namespace pack {

// Max vertices per polygon. Clipping an m-gon by an n-gon gives at most m + n
// vertices, so the clipper's fixed-size stack buffers hold 2 * kMaxVertices.
// Fixed arrays on the stack avoid heap allocation in the hot loop.
inline constexpr int kMaxVertices = 16;

inline constexpr double kPi = 3.14159265358979323846;

struct Vec2 {
    double x = 0.0, y = 0.0;
};

// vector algebra. `inline` lets them be defined in a header safely.
inline Vec2 operator+(Vec2 a, Vec2 b) { return {a.x + b.x, a.y + b.y}; }
inline Vec2 operator-(Vec2 a, Vec2 b) { return {a.x - b.x, a.y - b.y}; }
inline Vec2 operator*(double s, Vec2 a) { return {s * a.x, s * a.y}; }  
inline double dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }
inline double cross(Vec2 a, Vec2 b) { return a.x * b.y - a.y * b.x; }
inline double norm(Vec2 a) { return std::sqrt(dot(a, a)); }
inline Vec2 normalized(Vec2 a) { const double n = norm(a); return {a.x / n, a.y / n}; }  

// Geometry of one species.
struct Shape {
    std::vector<Vec2> vertices;  // counter-clockwise, relative to shape center. 
    std::vector<Vec2> sat_axes;  // unit edge normals. Parallel edges give the same axis, so e.g. a square needs 2, not 4.

    double circumradius = 0.0;   // centre-to-farthest-vertex distance
    double area = 0.0;           // polygon area 
};

// Factory functions. Declared here, defined in the .cpp.
Shape make_square(double side = 1.0);                       // axis-aligned at theta = 0
Shape make_regular_polygon(int n, double circumradius);
Shape make_rectangle(double width, double height);           // axis-aligned at theta = 0

// For vertices that came from outside the program. Nothing if they are not 3 to kMaxVertices
// corners of a strictly convex, counter-clockwise polygon: the clipper assumes exactly that, and
// gives wrong areas rather than an error for anything else.
std::optional<Shape> shape_from_vertices(const std::vector<Vec2>& vertices);

// A particle's pose: position (x_, y_) and angle theta_. cos(theta) and sin(theta) are cached because computing them is expensive and  they are needed constantly. Updated with set_pose(double x, double y, double theta), which updates all three together.
class Particle {
public:
    bool frozen = false;  // public: no invariant to protect, so no need to hide it
    int species = 0;       // index into the owning System::species; which Shape this one is

    void set_pose(double x, double y, double theta) {
        x_ = x;
        y_ = y;
        theta_ = theta;
        cos_t_ = std::cos(theta);   // keep the cache consistent
        sin_t_ = std::sin(theta);
    }

    // Read-only accessors
    double x() const { return x_; }
    double y() const { return y_; }
    double theta() const { return theta_; }
    double cos_t() const { return cos_t_; }
    double sin_t() const { return sin_t_; }

private:
    double x_ = 0.0, y_ = 0.0, theta_ = 0.0;
    double cos_t_ = 1.0, sin_t_ = 0.0;   // cos(0), sin(0): consistent defaults
};

// Free-function wrapper so call sites read like the other geometry functions.
inline void set_pose(Particle& p, double x, double y, double theta) {
    p.set_pose(x, y, theta);
}

// World position of vertex k: rotate the body-frame vertex by theta, then
// translate by the centre (x, y).
inline Vec2 world_vertex(const Particle& p, const Shape& s, int k) {
    const Vec2 v = s.vertices[static_cast<size_t>(k)];
    return {p.x() + v.x * p.cos_t() - v.y * p.sin_t(),
            p.y() + v.x * p.sin_t() + v.y * p.cos_t()};
}

// World direction of edge-normal axis k.
inline Vec2 world_axis(const Particle& p, const Shape& s, int k) {
    const Vec2 a = s.sat_axes[static_cast<size_t>(k)];
    return {a.x * p.cos_t() - a.y * p.sin_t(), a.x * p.sin_t() + a.y * p.cos_t()};
}

// Fills `out` with all world-space vertices; returns how many. The caller
// supplies a buffer of at least kMaxVertices entries (usually a stack array).
inline int world_vertices(const Particle& p, const Shape& s, Vec2* out) {
    const int n = static_cast<int>(s.vertices.size());
    assert(n > 0 && n <= kMaxVertices);  // guard against buffer overflow

    for (int k = 0; k < n; ++k) out[k] = world_vertex(p, s, k);

    return n;
}

// Shadow ("projection") of the polygon on the line along `axis`: the interval
// [lo, hi] of positions its points occupy along that line.
//
// Trick: dot(R v, axis) == dot(v, R^T axis). Rotating the single axis into the
// body frame (R^T = rotate by -theta) once is cheaper than rotating every
// vertex into the world frame. The centre's own projection `c` is added at the end.
inline void project(const Particle& p, const Shape& s, Vec2 axis, double& lo, double& hi) {
    assert(!s.vertices.empty());
    const double c = p.x() * axis.x + p.y() * axis.y; // centre's shadow position

    // body-frame
    const Vec2 local{axis.x * p.cos_t() + axis.y * p.sin_t(),
                     -axis.x * p.sin_t() + axis.y * p.cos_t()};

    double a = dot(s.vertices[0], local), b = a; // running min and max

    for (size_t k = 1; k < s.vertices.size(); ++k) {
        const double d = dot(s.vertices[k], local);
        if (d < a) a = d;
        if (d > b) b = d;
    }
    lo = c + a;
    hi = c + b;
}

// A position and angle (x, y, theta). Also the derivative of an area with respect to them: the force (x, y) and torque (theta) that a penalty on that area would produce.
struct Pose {
    double x = 0.0, y = 0.0, theta = 0.0;
};

// True if the interiors overlap. Touching counts as zero overlap. a and b each read their own shape -- sa for a, sb for b -- so two particles of different species compare correctly; the common case sa == sb (same species) costs nothing extra.
bool overlaps(const Particle& a, const Shape& sa, const Particle& b, const Shape& sb);

// Area of the intersection; 0 if disjoint.
double overlap_area(const Particle& a, const Shape& sa, const Particle& b, const Shape& sb);

// Area plus its gradient with respect to a's and b's pose. Pass nullptr for a gradient not needed. Moving a body only changes the overlap area along the part of its boundary that lies inside the other body. So
//     d(area)/d(translation) = integral of the outward normal over that boundary
//     d(area)/d(theta)       = integral of cross(point - centre, normal) over it
// For polygons these are sums over straight pieces. The translation gradients of a and b are equal and opposite.
double overlap_area_grad(const Particle& a, const Shape& sa, const Particle& b, const Shape& sb,
                         Pose* ga, Pose* gb);

// Area of the part of p outside the box [0,L]^2 (0 if fully inside). Requires L > 0.
double outside_area(const Particle& p, const Shape& s, double L);

// Area of the part of p outside the box [0,L]^2 plus the gradient.
double outside_area_grad(const Particle& p, const Shape& s, double L, Pose* g);

}  // namespace pack