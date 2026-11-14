#include "merge/rectangle.hpp"

#include <algorithm>
#include <cmath>
#include <map>

namespace pack {
namespace {

Vec2 rotate(Vec2 v, double cos_t, double sin_t) {
    return {v.x * cos_t - v.y * sin_t, v.x * sin_t + v.y * cos_t};
}

}  // namespace

std::optional<RectangleUnion> detect_rectangle_union(const System& sys,
                                                      const std::vector<int>& members,
                                                      double pos_tol, double theta_tol) {
    const int n = static_cast<int>(sys.particles.size());
    const size_t k = members.size();

    // 1. Count, range, no duplicates.
    if (k < 2) return std::nullopt;
    for (int i : members) {
        if (i < 0 || i >= n) return std::nullopt;
    }
    {
        std::vector<int> sorted = members;
        std::sort(sorted.begin(), sorted.end());
        if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) return std::nullopt;
    }

    // 2. One species, and that species is a square.
    const int sp = sys.particles[static_cast<size_t>(members[0])].species;
    for (int i : members) {
        if (sys.particles[static_cast<size_t>(i)].species != sp) return std::nullopt;
    }
    const Shape& s = sys.species[static_cast<size_t>(sp)];
    if (s.vertices.size() != 4) return std::nullopt;
    double xlo = s.vertices[0].x, xhi = xlo, ylo = s.vertices[0].y, yhi = ylo;
    for (const Vec2& v : s.vertices) {
        xlo = std::min(xlo, v.x);
        xhi = std::max(xhi, v.x);
        ylo = std::min(ylo, v.y);
        yhi = std::max(yhi, v.y);
    }
    const double hx = 0.5 * (xhi - xlo), hy = 0.5 * (yhi - ylo);
    if (std::abs(hx - hy) > pos_tol) return std::nullopt;  // not a square in its own local frame
    const double side = 2.0 * hx;
    if (side <= 0.0) return std::nullopt;
    if (std::abs(side * side - s.area) > pos_tol * std::max(1.0, s.area)) return std::nullopt;

    // 3. One orientation, up to a multiple of 90 degrees.
    const double theta0 = sys.particles[static_cast<size_t>(members[0])].theta();
    constexpr double kHalfPi = kPi / 2.0;
    for (int i : members) {
        const double residual = sys.particles[static_cast<size_t>(i)].theta() - theta0;
        const double quarter = std::round(residual / kHalfPi);
        if (std::abs(residual - quarter * kHalfPi) > theta_tol) return std::nullopt;
    }
    const double c0 = std::cos(theta0), s0 = std::sin(theta0);

    // 4-5. Rotate into the group's local frame (member 0 as origin), snap to the square's own
    // pitch; any residual past pos_tol is a gap or a stray offset.
    const Vec2 origin{sys.particles[static_cast<size_t>(members[0])].x(),
                      sys.particles[static_cast<size_t>(members[0])].y()};
    std::vector<int> cols(k), rows(k);
    for (size_t b = 0; b < k; ++b) {
        const Particle& p = sys.particles[static_cast<size_t>(members[b])];
        const Vec2 d{p.x() - origin.x, p.y() - origin.y};
        // R(-theta0) * d: rotating the world offset into the group's body frame.
        const Vec2 local = rotate(d, c0, -s0);
        const double uc = local.x / side, vc = local.y / side;
        const long col = std::lround(uc), row = std::lround(vc);
        if (std::abs(uc - static_cast<double>(col)) * side > pos_tol) return std::nullopt;
        if (std::abs(vc - static_cast<double>(row)) * side > pos_tol) return std::nullopt;
        cols[b] = static_cast<int>(col);
        rows[b] = static_cast<int>(row);
    }

    // 6. Distinct cells, exactly filling a cols x rows rectangle -- rejects a non-rectangular
    // touching cluster (right count, wrong cell set: an L- or plus-shape).
    std::map<std::pair<int, int>, int> cell_to_member;
    for (size_t b = 0; b < k; ++b) {
        const auto [it, inserted] = cell_to_member.emplace(std::pair{cols[b], rows[b]}, static_cast<int>(b));
        if (!inserted) return std::nullopt;  // two members claim the same cell
    }
    const int col_min = *std::min_element(cols.begin(), cols.end());
    const int col_max = *std::max_element(cols.begin(), cols.end());
    const int row_min = *std::min_element(rows.begin(), rows.end());
    const int row_max = *std::max_element(rows.begin(), rows.end());
    const int ncols = col_max - col_min + 1;
    const int nrows = row_max - row_min + 1;
    if (static_cast<size_t>(ncols) * static_cast<size_t>(nrows) != k) return std::nullopt;

    // 7. Cheap defensive cross-check (redundant given 5-6 for the same-shape case above, but
    // catches a future relaxation of step 2 to "same area, different species").
    const double rect_area = static_cast<double>(ncols) * side * static_cast<double>(nrows) * side;
    const double members_area = static_cast<double>(k) * s.area;
    if (std::abs(rect_area - members_area) > 1e-9 * std::max(1.0, members_area)) {
        return std::nullopt;
    }

    // 8. Build the result: centre from the ideal grid (not a centroid of possibly-jittery actual
    // positions), members reindexed row-major.
    RectangleUnion out;
    out.side = side;
    out.cols = ncols;
    out.rows = nrows;
    out.theta = theta0;
    const Vec2 center_local{(static_cast<double>(col_min + col_max) * 0.5) * side,
                            (static_cast<double>(row_min + row_max) * 0.5) * side};
    const Vec2 center_offset = rotate(center_local, c0, s0);  // R(theta0) * center_local
    out.center = {origin.x + center_offset.x, origin.y + center_offset.y};
    out.ordered_members.assign(static_cast<size_t>(ncols) * static_cast<size_t>(nrows), -1);
    for (size_t b = 0; b < k; ++b) {
        const int c = cols[b] - col_min, r = rows[b] - row_min;
        out.ordered_members[static_cast<size_t>(r * ncols + c)] = members[b];
    }
    return out;
}

}  // namespace pack
