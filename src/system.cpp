#include "system.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

namespace pack {
namespace {

double sq(double v) { return v * v; }

}  // namespace

double System::reach() const { return 2.0 * max_circumradius(); }

EnergyParts System::energy_parts() const {
    EnergyParts parts;
    for (size_t i = 0; i < particles.size(); ++i) {
        parts.wall += sq(outside_area(particles[i], shape_of(static_cast<int>(i)), L));
    }
    CellList grid;
    grid.build(particles, L, reach());
    grid.for_each_close_pair([&](int i, int j) {
        parts.pair += sq(overlap_area(particles[static_cast<size_t>(i)], shape_of(i),
                                      particles[static_cast<size_t>(j)], shape_of(j)));
    });
    return parts;
}

double System::total_energy() const {
    return energy_of(energy_parts(), k_pair, k_wall);
}

EnergyParts System::particle_parts(int i) const {
    const int n = static_cast<int>(particles.size());
    assert(i >= 0 && i < n);
    assert(cells.cell_of.size() == particles.size());  // built at all: recount() does it
    const Particle& pi = particles[static_cast<size_t>(i)];
    const Shape& si = shape_of(i);

    EnergyParts out;
    out.wall = sq(outside_area(pi, si, L));
    cells.for_each_near(cells.bin(pi.x(), pi.y()), [&](int j) {
        if (j == i) return;
        out.pair += sq(overlap_area(pi, si, particles[static_cast<size_t>(j)], shape_of(j)));
    });
    return out;
}

double System::max_overlap() const {
    double worst = 0.0;
    for (size_t i = 0; i < particles.size(); ++i) {
        worst = std::max(worst, outside_area(particles[i], shape_of(static_cast<int>(i)), L));
    }
    CellList grid;
    grid.build(particles, L, reach());
    grid.for_each_close_pair([&](int i, int j) {
        worst = std::max(worst, overlap_area(particles[static_cast<size_t>(i)], shape_of(i),
                                             particles[static_cast<size_t>(j)], shape_of(j)));
    });
    return worst;
}

bool System::is_feasible(double tol) const { return max_overlap() <= tol; }

double System::tight_side() const {
    if (particles.empty()) return 0.0;
    const double inf = std::numeric_limits<double>::infinity();
    double xlo = inf, xhi = -inf, ylo = inf, yhi = -inf;
    for (size_t idx = 0; idx < particles.size(); ++idx) {
        const Particle& p = particles[idx];
        const Shape& s = shape_of(static_cast<int>(idx));
        double lo, hi;
        project(p, s, {1.0, 0.0}, lo, hi);
        xlo = std::min(xlo, lo);
        xhi = std::max(xhi, hi);
        project(p, s, {0.0, 1.0}, lo, hi);
        ylo = std::min(ylo, lo);
        yhi = std::max(yhi, hi);
    }
    return std::max(xhi - xlo, yhi - ylo);
}

int System::mobile_count() const {
    int m = 0;
    for (const Particle& p : particles) {
        if (!p.frozen) ++m;
    }
    return m;
}

void System::set_box(double new_L) {
    assert(new_L > 0.0);
    L = new_L;
    recount();
}

void System::rescale_box(double new_L) {
    assert(new_L > 0.0);
    const double s = new_L / L;
    for (Particle& p : particles) {
        if (p.frozen) continue;
        p.set_pose(p.x() * s, p.y() * s, p.theta());
    }
    L = new_L;

    recount();
}

void System::fit_box() {
    if (particles.empty()) return;
    double xlo = std::numeric_limits<double>::infinity(), ylo = xlo;
    for (size_t i = 0; i < particles.size(); ++i) {
        double lo, hi;
        project(particles[i], shape_of(static_cast<int>(i)), {1.0, 0.0}, lo, hi);
        xlo = std::min(xlo, lo);
        project(particles[i], shape_of(static_cast<int>(i)), {0.0, 1.0}, lo, hi);
        ylo = std::min(ylo, lo);
    }
    for (Particle& p : particles) p.set_pose(p.x() - xlo, p.y() - ylo, p.theta());
    L = tight_side();
    recount();
}

void System::resync_energy() {
    for (size_t i = 0; i < particles.size(); ++i) {
        assert(cells.cell_of[i] == cells.bin(particles[i].x(), particles[i].y()));
    }
    const EnergyParts exact_parts = energy_parts();
    assert(std::abs(parts.pair - exact_parts.pair) <= 1e-6 * std::max(1.0, exact_parts.pair));
    assert(std::abs(parts.wall - exact_parts.wall) <= 1e-6 * std::max(1.0, exact_parts.wall));
    parts = exact_parts;
}

void System::rebuild_cells() { cells.build(particles, L, reach()); }

void System::recount() {
    rebuild_cells();
    parts = energy_parts();
}

}  // namespace pack
