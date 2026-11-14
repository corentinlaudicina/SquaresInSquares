#include "geometry/cells.hpp"

#include <algorithm>
#include <cmath>

namespace pack {

// Sets up the grid and files every particle into its cell.
//   L     : side length of the square box [0, L]^2
//   reach : max distance at which two particles can interact

void CellList::build(const std::vector<Particle>& particles, double L, double reach) {
    const int n = static_cast<int>(particles.size());

    // Upper limit on cells per side
    const int cap = 2 * static_cast<int>(std::ceil(std::sqrt(double(n)))) + 1;

    nc   = std::clamp(static_cast<int>(L / reach), 1, cap);
    side = L / nc;

    // Reset all the arrays
    head.assign(nc * nc, -1);
    next.assign(n, -1);
    prev.assign(n, -1);
    cell_of.assign(n, 0);

    // Assign each particle to its box
    for (int i = 0; i < n; ++i)
        link(i, bin(particles[i].x(), particles[i].y()));
}

int CellList::bin(double x, double y) const {
    auto axis = [&](double v) { return std::clamp(int(v / side), 0, nc - 1); };
    return axis(y) * nc + axis(x);
}

void CellList::move(int i, double x, double y) {
    const int c = bin(x, y);
    if (c == cell_of[static_cast<size_t>(i)]) return;
    unlink(i);
    link(i, c);
}

void CellList::link(int i, int cell) {
    const size_t k = static_cast<size_t>(i);
    int& h = head[static_cast<size_t>(cell)];
    prev[k] = -1;
    next[k] = h;
    if (h >= 0) prev[static_cast<size_t>(h)] = i;
    h = i;
    cell_of[k] = cell;
}

void CellList::unlink(int i) {
    const size_t k = static_cast<size_t>(i);
    if (prev[k] >= 0) {
        next[static_cast<size_t>(prev[k])] = next[k];
    } else {
        head[static_cast<size_t>(cell_of[k])] = next[k];
    }
    if (next[k] >= 0) prev[static_cast<size_t>(next[k])] = prev[k];
}

}  // namespace pack
