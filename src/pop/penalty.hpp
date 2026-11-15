#pragma once

// The population layer's cost: a packing of N unit squares as a flat genome, and a penalty with its
// gradient, with no System, no allocation per call and an O(N) neighbour grid. Squares only.
//
// Penalties (box [0, h]^2 fixed):
//   Area2   sum over pairs of overlap_area^2 + sum over squares of outside_area^2. The same function
//           System::energy() computes at k_pair = k_wall = 1; kept as the reference.
//   Depth2  sum over pairs of depth^2 + sum over squares and walls of protrusion^2, where depth is
//           the pair's penetration depth (smallest interval overlap over the 4 edge normals -- for
//           convex polygons that is the minimum translation distance) and protrusion how far the
//           square crosses a wall. Grows like depth^2 near contact, where Area2 grows like depth^4
//           at a corner, so a relaxation can drive it to ~1e-12 instead of stalling near 1e-5.
//
// Population layer inspired by jcottaar/packing (Santa 2025, 1st place, CC BY-SA 4.0); this is an
// independent reimplementation for squares, no code was copied.

#include "geometry/polygon.hpp"

#include <vector>

namespace pack {

struct Packing {
    double h = 1.0;            // box side
    std::vector<double> xyt;   // x0, y0, theta0, x1, ...
    int n() const { return static_cast<int>(xyt.size() / 3); }
};

enum class Penalty { Area2, Depth2 };

struct CostInfo {
    double cost = 0.0;
    double max_depth = 0.0;  // Depth2: worst pair penetration depth; Area2: worst pair overlap area
    double max_wall = 0.0;   // Depth2: worst protrusion; Area2: worst outside area
};

// Reused buffers: size them once per N and thread, then every call is allocation free.
struct Workspace {
    std::vector<double> c, s;          // cos, sin per square
    std::vector<int> head, next, cell;  // neighbour grid
    std::vector<Particle> parts;       // Area2 only: poses for the polygon clipper
    Shape square = make_square(1.0);
    // relax_capped's buffers (relax_capped.cpp)
    std::vector<double> g, g_old, d, best, hist_s, hist_y, rho, alpha;
    // symmetric genotypes (symmetric.cpp): the expanded packing and its gradient
    Packing sym_full;
    std::vector<double> sym_g;
};

// grad, if not null, has 3N entries: d cost / d (x, y, theta) per square.
CostInfo cost_grad(const Packing& p, Penalty pen, Workspace& ws, double* grad);

}  // namespace pack
