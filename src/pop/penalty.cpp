#include "pop/penalty.hpp"

#include <algorithm>
#include <cmath>

namespace pack {
namespace {

constexpr double kReach = 1.4142135623730951;  // centre distance below which two unit squares can touch
constexpr double kSkin = 0.3;                   // Verlet margin; a rebuild once a square moved kSkin / 2

inline double sgn(double v) { return (v > 0.0) - (v < 0.0); }

// Penetration depth of squares i and j and its gradient g[6] = d/d(xi, yi, ti, xj, yj, tj).
// 0 (and g untouched) when some axis separates them.
//
// Along a unit axis n, square P covers [m - r, m + r] with m = n . centre_P and
// r = (|n . a1| + |n . a2|) / 2, a1 = (c, s), a2 = (-s, c) its edge directions. The axes are the
// edge directions of both squares, so n itself turns with its owner O: dn/dtheta_O = perp(n).
double pair_depth(const double* xyt, const Workspace& ws, int i, int j, double* g) {
    const int idx[2] = {i, j};
    // Pass 1: the separation distance along each of the 4 axes, no derivatives. Any separating axis
    // ends it; otherwise the smallest distance (first one on ties) is the depth.
    double best = 1e300;
    int best_o = 0, best_e = 0, best_up = 0;
    for (int o = 0; o < 2; ++o) {              // axis owner
        for (int e = 0; e < 2; ++e) {          // which of its two edge normals
            const double co = ws.c[static_cast<size_t>(idx[o])], so = ws.s[static_cast<size_t>(idx[o])];
            const double nx = e == 0 ? co : -so, ny = e == 0 ? so : co;
            double m[2], r[2];
            for (int q = 0; q < 2; ++q) {
                const size_t k = static_cast<size_t>(idx[q]);
                const double x = xyt[3 * k], y = xyt[3 * k + 1], c = ws.c[k], s = ws.s[k];
                m[q] = nx * x + ny * y;
                const double t1 = nx * c + ny * s, t2 = -nx * s + ny * c;
                r[q] = 0.5 * (std::abs(t1) + std::abs(t2));
            }
            // Distance to separate along n: push one way (hi_0 - lo_1) or the other (hi_1 - lo_0).
            // Equal to the interval overlap unless one interval contains the other.
            const int up = ((m[0] + r[0]) - (m[1] - r[1]) < (m[1] + r[1]) - (m[0] - r[0])) ? 0 : 1;
            const int lo = 1 - up;
            const double ov = (m[up] + r[up]) - (m[lo] - r[lo]);
            if (ov <= 0.0) return 0.0;
            if (ov < best) {
                best = ov;
                best_o = o;
                best_e = e;
                best_up = up;
            }
        }
    }
    if (!g) return best;

    // Pass 2: the gradient of the winning axis only. Along a unit axis n, square P covers
    // [m - r, m + r] with m = n . centre_P and r = (|n . a1| + |n . a2|) / 2, a1 = (c, s), a2 = (-s, c)
    // its edge directions; n turns with its owner O, dn/dtheta_O = perp(n).
    const int o = best_o, e = best_e, up = best_up, lo = 1 - up;
    const double co = ws.c[static_cast<size_t>(idx[o])], so = ws.s[static_cast<size_t>(idx[o])];
    const double nx = e == 0 ? co : -so, ny = e == 0 ? so : co;
    const double px = -ny, py = nx;  // perp(n) = dn/dtheta_owner
    double dm[2][6] = {}, dr[2][6] = {};
    for (int q = 0; q < 2; ++q) {
        const size_t k = static_cast<size_t>(idx[q]);
        const double x = xyt[3 * k], y = xyt[3 * k + 1], c = ws.c[k], s = ws.s[k];
        dm[q][3 * q] = nx;
        dm[q][3 * q + 1] = ny;
        dm[q][3 * o + 2] += px * x + py * y;
        const double t1 = nx * c + ny * s, t2 = -nx * s + ny * c;
        // own rotation: dt1 = t2, dt2 = -t1; owner rotation: perp(n) . a1, perp(n) . a2
        dr[q][3 * q + 2] += 0.5 * (sgn(t1) * t2 - sgn(t2) * t1);
        dr[q][3 * o + 2] += 0.5 * (sgn(t1) * (px * c + py * s) + sgn(t2) * (-px * s + py * c));
    }
    for (int k = 0; k < 6; ++k) g[k] = dm[up][k] + dr[up][k] - dm[lo][k] + dr[lo][k];
    return best;
}

void build_grid(const Packing& p, Workspace& ws, int& nc, double& side) {
    const int n = p.n();
    const int cap = 2 * static_cast<int>(std::ceil(std::sqrt(static_cast<double>(n)))) + 1;
    nc = std::clamp(static_cast<int>(p.h / kReach), 1, cap);
    side = p.h / nc;
    ws.head.assign(static_cast<size_t>(nc) * static_cast<size_t>(nc), -1);  // no realloc once sized
    ws.next.resize(static_cast<size_t>(n));
    ws.cell.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const auto bin = [&](double v) {
            const double b = v / side;
            return b != b ? 0 : std::clamp(static_cast<int>(b), 0, nc - 1);  // NaN-safe (audit B11)
        };
        const int cell = bin(p.xyt[3 * static_cast<size_t>(i) + 1]) * nc + bin(p.xyt[3 * static_cast<size_t>(i)]);
        ws.cell[static_cast<size_t>(i)] = cell;
        ws.next[static_cast<size_t>(i)] = ws.head[static_cast<size_t>(cell)];
        ws.head[static_cast<size_t>(cell)] = i;
    }
}

template <class F>
void for_each_close_pair(const Packing& p, const Workspace& ws, int nc, F&& f) {
    const int n = p.n();
    for (int i = 0; i < n; ++i) {
        const int cx = ws.cell[static_cast<size_t>(i)] % nc, cy = ws.cell[static_cast<size_t>(i)] / nc;
        for (int y = std::max(cy - 1, 0); y <= std::min(cy + 1, nc - 1); ++y)
            for (int x = std::max(cx - 1, 0); x <= std::min(cx + 1, nc - 1); ++x)
                for (int j = ws.head[static_cast<size_t>(y * nc + x)]; j >= 0; j = ws.next[static_cast<size_t>(j)])
                    if (j > i) f(i, j);
    }
}

// The close pairs: from the Verlet list while ws.verlet is on (rebuilt when stale), else from the grid.
template <class F>
void visit_pairs(const Packing& p, Workspace& ws, F&& f) {
    int nc;
    double side;
    const size_t N = static_cast<size_t>(p.n());
    if (!ws.verlet) {
        build_grid(p, ws, nc, side);
        for_each_close_pair(p, ws, nc, f);
        return;
    }
    bool stale = !ws.verlet_valid || ws.verlet_h != p.h || ws.ref.size() != 2 * N;
    for (size_t i = 0; i < N && !stale; ++i) {
        const double dx = p.xyt[3 * i] - ws.ref[2 * i], dy = p.xyt[3 * i + 1] - ws.ref[2 * i + 1];
        stale = dx * dx + dy * dy > 0.25 * kSkin * kSkin;
    }
    if (stale) {
        build_grid(p, ws, nc, side);
        ws.pairs.clear();
        const double cut = (kReach + kSkin) * (kReach + kSkin);
        for_each_close_pair(p, ws, nc, [&](int i, int j) {
            const double dx = p.xyt[3 * static_cast<size_t>(i)] - p.xyt[3 * static_cast<size_t>(j)];
            const double dy = p.xyt[3 * static_cast<size_t>(i) + 1] - p.xyt[3 * static_cast<size_t>(j) + 1];
            if (dx * dx + dy * dy < cut) {
                ws.pairs.push_back(i);
                ws.pairs.push_back(j);
            }
        });
        ws.ref.resize(2 * N);
        for (size_t i = 0; i < N; ++i) {
            ws.ref[2 * i] = p.xyt[3 * i];
            ws.ref[2 * i + 1] = p.xyt[3 * i + 1];
        }
        ws.verlet_h = p.h;
        ws.verlet_valid = true;
        ++ws.verlet_builds;
    }
    for (size_t k = 0; k < ws.pairs.size(); k += 2) f(ws.pairs[k], ws.pairs[k + 1]);
}

}  // namespace

CostInfo cost_grad(const Packing& p, Penalty pen, Workspace& ws, double* grad) {
    const int n = p.n();
    const size_t N = static_cast<size_t>(n);
    const double* xyt = p.xyt.data();
    CostInfo info;
    if (grad) std::fill(grad, grad + 3 * N, 0.0);
    ws.c.resize(N);
    ws.s.resize(N);
    for (size_t i = 0; i < N; ++i) {
        ws.c[i] = std::cos(xyt[3 * i + 2]);
        ws.s[i] = std::sin(xyt[3 * i + 2]);
    }
    if (pen == Penalty::Area2) {
        ws.parts.resize(N);
        for (size_t i = 0; i < N; ++i) ws.parts[i].set_pose(xyt[3 * i], xyt[3 * i + 1], xyt[3 * i + 2]);
        for (size_t i = 0; i < N; ++i) {
            Pose gw;
            const double out = outside_area_grad(ws.parts[i], ws.square, p.h, grad ? &gw : nullptr);
            info.cost += out * out;
            info.max_wall = std::max(info.max_wall, out);
            if (grad && out > 0.0) {
                grad[3 * i] += 2 * out * gw.x;
                grad[3 * i + 1] += 2 * out * gw.y;
                grad[3 * i + 2] += 2 * out * gw.theta;
            }
        }
        visit_pairs(p, ws, [&](int i, int j) {
            Pose gi, gj;
            const size_t a = static_cast<size_t>(i), b = static_cast<size_t>(j);
            const double ar = overlap_area_grad(ws.parts[a], ws.square, ws.parts[b], ws.square,
                                                grad ? &gi : nullptr, grad ? &gj : nullptr);
            if (ar <= 0.0) return;
            info.cost += ar * ar;
            info.max_depth = std::max(info.max_depth, ar);
            if (grad) {
                grad[3 * a] += 2 * ar * gi.x;
                grad[3 * a + 1] += 2 * ar * gi.y;
                grad[3 * a + 2] += 2 * ar * gi.theta;
                grad[3 * b] += 2 * ar * gj.x;
                grad[3 * b + 1] += 2 * ar * gj.y;
                grad[3 * b + 2] += 2 * ar * gj.theta;
            }
        });
        return info;
    }

    // Depth2. Walls: a square covers [x - r, x + r] in x and [y - r, y + r] in y, r = (|c| + |s|) / 2.
    for (size_t i = 0; i < N; ++i) {
        const double c = ws.c[i], s = ws.s[i];
        const double r = 0.5 * (std::abs(c) + std::abs(s));
        const double dr = 0.5 * (-sgn(c) * s + sgn(s) * c);
        for (int axis = 0; axis < 2; ++axis) {
            const double v = xyt[3 * i + static_cast<size_t>(axis)];
            for (int side_hi = 0; side_hi < 2; ++side_hi) {
                const double pr = side_hi ? v + r - p.h : r - v;  // protrusion past h, or past 0
                if (pr <= 0.0) continue;
                info.cost += pr * pr;
                info.max_wall = std::max(info.max_wall, pr);
                if (grad) {
                    grad[3 * i + static_cast<size_t>(axis)] += 2 * pr * (side_hi ? 1.0 : -1.0);
                    grad[3 * i + 2] += 2 * pr * dr;
                }
            }
        }
    }
    visit_pairs(p, ws, [&](int i, int j) {
        const size_t a = static_cast<size_t>(i), b = static_cast<size_t>(j);
        const double dx = xyt[3 * a] - xyt[3 * b], dy = xyt[3 * a + 1] - xyt[3 * b + 1];
        if (dx * dx + dy * dy >= kReach * kReach) return;
        double g[6];
        const double d = pair_depth(xyt, ws, i, j, grad ? g : nullptr);
        if (d <= 0.0) return;
        info.cost += d * d;
        info.max_depth = std::max(info.max_depth, d);
        if (grad) {
            for (int k = 0; k < 3; ++k) {
                grad[3 * a + static_cast<size_t>(k)] += 2 * d * g[k];
                grad[3 * b + static_cast<size_t>(k)] += 2 * d * g[3 + k];
            }
        }
    });
    return info;
}

}  // namespace pack
