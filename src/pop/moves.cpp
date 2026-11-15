#include "pop/moves.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace pack {
namespace {

inline double& X(Packing& p, int i) { return p.xyt[3 * static_cast<size_t>(i)]; }
inline double& Y(Packing& p, int i) { return p.xyt[3 * static_cast<size_t>(i) + 1]; }
inline double& T(Packing& p, int i) { return p.xyt[3 * static_cast<size_t>(i) + 2]; }

std::vector<int> within(const Packing& p, double cx, double cy, double r) {
    std::vector<int> out;
    for (int i = 0; i < p.n(); ++i) {
        const double dx = p.xyt[3 * static_cast<size_t>(i)] - cx, dy = p.xyt[3 * static_cast<size_t>(i) + 1] - cy;
        if (dx * dx + dy * dy <= r * r) out.push_back(i);
    }
    return out;
}

void jiggle(Packing& p, int i, double sigma, std::mt19937_64& rng) {
    X(p, i) += sigma * normal(rng);
    Y(p, i) += sigma * normal(rng);
    T(p, i) += sigma * normal(rng);
}

}  // namespace

double normal(std::mt19937_64& rng) {
    const double u1 = 1.0 - u01(rng);  // (0, 1]
    const double u2 = u01(rng);
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * kPi * u2);
}

const char* ga_move_name(GaMove m) {
    switch (m) {
        case GaMove::Relocate: return "relocate";
        case GaMove::JiggleSmall: return "jiggle_small";
        case GaMove::JiggleLarge: return "jiggle_large";
        case GaMove::Cluster: return "cluster";
        case GaMove::Twist: return "twist";
        case GaMove::Translate: return "translate";
        case GaMove::Crossover: return "crossover";
        case GaMove::Destroy: return "destroy";
        case GaMove::Strip: return "strip";
        case GaMove::Corner: return "corner";
        case GaMove::Count: break;
    }
    return "?";
}

GaMove draw_move(const MoveParams& prm, std::mt19937_64& rng) {
    double total = 0.0;
    for (const double w : prm.weight) total += std::max(0.0, w);
    double r = u01(rng) * total;
    for (int k = 0; k < kGaMoves; ++k) {
        r -= std::max(0.0, prm.weight[static_cast<size_t>(k)]);
        if (r < 0.0) return static_cast<GaMove>(k);
    }
    return GaMove::JiggleSmall;
}

void apply_d4(Packing& p, int k) {
    // k = 0..3: rotation by k quarter turns; 4..7: the same followed by the mirror x -> -x.
    const double c = 0.5 * p.h;
    for (int i = 0; i < p.n(); ++i) {
        double x = X(p, i) - c, y = Y(p, i) - c;
        for (int q = 0; q < (k & 3); ++q) {
            const double t = x;
            x = -y;
            y = t;  // a quarter turn maps a square onto itself: theta unchanged mod 90
        }
        if (k >= 4) {
            x = -x;
            T(p, i) = -T(p, i);
        }
        X(p, i) = x + c;
        Y(p, i) = y + c;
    }
}

void apply_move(Packing& child, GaMove m, const MoveParams& prm, std::mt19937_64& rng, const Packing* mate) {
    const int n = child.n();
    const double h = child.h;
    switch (m) {
        case GaMove::Relocate: {
            const int i = uniform_int(rng, 0, n - 1);
            X(child, i) = uniform(rng, 0.5, h - 0.5);
            Y(child, i) = uniform(rng, 0.5, h - 0.5);
            T(child, i) = u01(rng) < prm.snap_prob ? 0.25 * kPi * uniform_int(rng, 0, 1) : uniform(rng, 0.0, 0.5 * kPi);
            return;
        }
        case GaMove::JiggleSmall:
            jiggle(child, uniform_int(rng, 0, n - 1), prm.jiggle_small, rng);
            return;
        case GaMove::JiggleLarge:
            jiggle(child, uniform_int(rng, 0, n - 1), prm.jiggle_large, rng);
            return;
        case GaMove::Cluster: {
            const int i = uniform_int(rng, 0, n - 1);
            const int k = std::min(n, uniform_int(rng, prm.cluster_min, std::max(prm.cluster_min, prm.cluster_max)));
            std::vector<int> idx(static_cast<size_t>(n));
            std::iota(idx.begin(), idx.end(), 0);
            const double cx = X(child, i), cy = Y(child, i);
            const auto d2 = [&](int j) {
                const double dx = X(child, j) - cx, dy = Y(child, j) - cy;
                return dx * dx + dy * dy;
            };
            std::partial_sort(idx.begin(), idx.begin() + k, idx.end(),
                              [&](int a, int b) { return d2(a) < d2(b) || (d2(a) == d2(b) && a < b); });
            for (int q = 0; q < k; ++q) jiggle(child, idx[static_cast<size_t>(q)], prm.cluster_sigma, rng);
            return;
        }
        case GaMove::Twist: {
            const double cx = uniform(rng, 0.0, h), cy = uniform(rng, 0.0, h);
            const double alpha = uniform(rng, -prm.twist_angle, prm.twist_angle);
            for (const int i : within(child, cx, cy, prm.twist_radius)) {
                const double dx = X(child, i) - cx, dy = Y(child, i) - cy;
                const double phi = alpha * (1.0 - std::sqrt(dx * dx + dy * dy) / prm.twist_radius);
                const double c = std::cos(phi), s = std::sin(phi);
                X(child, i) = cx + c * dx - s * dy;
                Y(child, i) = cy + s * dx + c * dy;
                T(child, i) += phi;
            }
            return;
        }
        case GaMove::Translate: {
            const double cx = uniform(rng, 0.0, h), cy = uniform(rng, 0.0, h);
            const double vx = prm.translate_sigma * normal(rng), vy = prm.translate_sigma * normal(rng);
            for (const int i : within(child, cx, cy, prm.translate_radius)) {
                X(child, i) += vx;
                Y(child, i) += vy;
            }
            return;
        }
        case GaMove::Crossover: {
            // Replace the child's squares near a random point by the mate's squares nearest that point,
            // the mate optionally mapped by a box symmetry first. The count stays exact.
            Packing b = *mate;
            b.h = h;
            if (u01(rng) < prm.crossover_symmetry) apply_d4(b, uniform_int(rng, 1, 7));
            const double px = uniform(rng, 0.0, h), py = uniform(rng, 0.0, h);
            const std::vector<int> mine = within(child, px, py, prm.crossover_radius);
            if (mine.empty()) {
                jiggle(child, uniform_int(rng, 0, n - 1), prm.jiggle_small, rng);
                return;
            }
            std::vector<int> theirs(static_cast<size_t>(n));
            std::iota(theirs.begin(), theirs.end(), 0);
            const auto d2 = [&](int j) {
                const double dx = X(b, j) - px, dy = Y(b, j) - py;
                return dx * dx + dy * dy;
            };
            const size_t k = mine.size();
            std::partial_sort(theirs.begin(), theirs.begin() + static_cast<long>(k), theirs.end(),
                              [&](int a, int c) { return d2(a) < d2(c) || (d2(a) == d2(c) && a < c); });
            for (size_t q = 0; q < k; ++q) {
                const int dst = mine[q], src = theirs[q];
                X(child, dst) = X(b, src);
                Y(child, dst) = Y(b, src);
                T(child, dst) = T(b, src);
            }
            return;
        }
        case GaMove::Destroy: {
            // Scramble a disc: every square in it gets a random pose inside the disc (half of them at a
            // multiple of 45 degrees). The relaxation rebuilds the region from scratch.
            const double r = uniform(rng, prm.destroy_min, std::max(prm.destroy_min, prm.destroy_frac * h));
            const double cx = uniform(rng, 0.0, h), cy = uniform(rng, 0.0, h);
            const std::vector<int> idx = within(child, cx, cy, r);
            if (idx.size() < 2) {
                jiggle(child, uniform_int(rng, 0, n - 1), prm.jiggle_large, rng);
                return;
            }
            for (const int i : idx) {
                double x, y;
                do {
                    x = cx + uniform(rng, -r, r);
                    y = cy + uniform(rng, -r, r);
                } while ((x - cx) * (x - cx) + (y - cy) * (y - cy) > r * r);
                X(child, i) = std::clamp(x, 0.5, h - 0.5);
                Y(child, i) = std::clamp(y, 0.5, h - 0.5);
                T(child, i) = u01(rng) < prm.snap_prob ? 0.25 * kPi * uniform_int(rng, 0, 1) : uniform(rng, 0.0, 0.5 * kPi);
            }
            return;
        }
        case GaMove::Strip: {
            // Re-lay a band as a regular strip: direction 0, 90, 45 or 135 degrees, width 1..strip_width
            // rows; each row's squares at unit pitch along the band, turned to the band's angle.
            const double phi = 0.25 * kPi * uniform_int(rng, 0, 3);
            const double ux = std::cos(phi), uy = std::sin(phi), vx = -uy, vy = ux;
            const int w = uniform_int(rng, 1, std::max(1, prm.strip_width));
            const double px = uniform(rng, 0.0, h), py = uniform(rng, 0.0, h);
            std::vector<std::vector<std::pair<double, int>>> rows(static_cast<size_t>(w));
            for (int i = 0; i < n; ++i) {
                const double dx = X(child, i) - px, dy = Y(child, i) - py;
                const double b = dx * vx + dy * vy + 0.5 * w;
                if (b < 0.0 || b >= w) continue;
                rows[static_cast<size_t>(b)].push_back({dx * ux + dy * uy, i});
            }
            int moved = 0;
            for (int row = 0; row < w; ++row) {
                auto& r = rows[static_cast<size_t>(row)];
                if (r.empty()) continue;
                std::sort(r.begin(), r.end());
                const double mean = [&] { double s = 0; for (auto& e : r) s += e.first; return s / r.size(); }();
                const double a0 = std::round(mean - 0.5 * (r.size() - 1));
                const double off = row - 0.5 * (w - 1);
                for (size_t k = 0; k < r.size(); ++k) {
                    const double a = a0 + static_cast<double>(k);
                    const int i = r[k].second;
                    X(child, i) = std::clamp(px + a * ux + off * vx, 0.5, h - 0.5);
                    Y(child, i) = std::clamp(py + a * uy + off * vy, 0.5, h - 0.5);
                    T(child, i) = phi;
                    ++moved;
                }
            }
            if (moved < 2) jiggle(child, uniform_int(rng, 0, n - 1), prm.jiggle_large, rng);
            return;
        }
        case GaMove::Corner: {
            // Re-lay a corner block as an axis-aligned grid packed into that corner, nearest squares first.
            const int k = uniform_int(rng, 0, 3);
            const double d = uniform(rng, 2.0, std::max(2.0, prm.corner_frac * h));
            const double ox = (k == 1 || k == 2) ? h : 0.0, oy = (k >= 2) ? h : 0.0;
            const double sx = ox == 0.0 ? 1.0 : -1.0, sy = oy == 0.0 ? 1.0 : -1.0;
            std::vector<std::pair<double, int>> m;
            for (int i = 0; i < n; ++i) {
                const double dx = (X(child, i) - ox) * sx, dy = (Y(child, i) - oy) * sy;
                if (dx <= d && dy <= d) m.push_back({dx * dx + dy * dy, i});
            }
            if (m.size() < 2) {
                jiggle(child, uniform_int(rng, 0, n - 1), prm.jiggle_large, rng);
                return;
            }
            std::sort(m.begin(), m.end());
            const int cols = std::max(1, static_cast<int>(std::floor(d)));
            for (size_t q = 0; q < m.size(); ++q) {  // fill the block's diagonals outwards from the corner
                const int a = static_cast<int>(q) % cols, b = static_cast<int>(q) / cols;
                X(child, m[q].second) = ox + sx * (0.5 + a);
                Y(child, m[q].second) = oy + sy * (0.5 + b);
                T(child, m[q].second) = 0.0;
            }
            return;
        }
        case GaMove::Count: return;
    }
}

}  // namespace pack
