#pragma once

// Holds the configuration and the soft energy; bodies in system.cpp.

#include "geometry/cells.hpp"
#include "geometry/polygon.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <memory>
#include <vector>

namespace pack {

inline constexpr double kFeasibleTol = 1e-5;

// The two configuration-only sums the energy is built from:
//
//     E = k_pair * pair + k_wall * wall
struct EnergyParts {
    double pair = 0.0;  // sum of (pair overlap area)^2
    double wall = 0.0;  // sum of (area outside the box)^2
};

inline double energy_of(const EnergyParts& p, double k_pair, double k_wall) {
    return k_pair * p.pair + k_wall * p.wall;
}

// A merged particle's provenance: what to restore on dissolve(). Immutable once built -- see
// src/merge/merge.hpp -- so System can share it by pointer across copies instead of deep-copying
// it every time a System is copied (which happens constantly: candidates, compress()).
struct MergeRecord;

struct System {
    // One entry per distinct polygon type in use. species[0] is "the base species" -- the one
    // ordinary unit squares use, and the one moves/regions.cpp and moves/lnm.cpp fall back to
    // (see base_shape()) rather than becoming per-particle-aware themselves.
    std::vector<Shape> species;
    std::vector<Particle> particles;
    double L = 1.0;              // box side; a move in a fixed-pressure run
    double k_pair = 100.0;       // stiffness on pair overlap area
    double k_wall = 100.0;       // stiffness on area outside the box

    // Non-null iff particles[i] is a merged rigid group (src/merge/merge.hpp). May be SHORTER than
    // particles: any index at or past merge_of.size() simply isn't merged. That means every caller
    // that builds or resizes particles directly (every test file, init.cpp) needs no changes at all
    // -- merge_of only ever grows when merge()/dissolve() themselves touch it.
    std::vector<std::shared_ptr<const MergeRecord>> merge_of;

    bool is_merged(int i) const {
        return i >= 0 && static_cast<size_t>(i) < merge_of.size() &&
               merge_of[static_cast<size_t>(i)] != nullptr;
    }

    const Shape& shape_of(int i) const {
        return species[static_cast<size_t>(particles[static_cast<size_t>(i)].species)];
    }
    const Shape& base_shape() const { return species.front(); }

    double max_circumradius() const {
        double m = 0.0;
        for (const Shape& s : species) m = std::max(m, s.circumradius);
        return m;
    }

    int add_species(const Shape& s) {
        species.push_back(s);
        return static_cast<int>(species.size()) - 1;
    }

    // Appends one particle of the given species, at the given pose. The caller owes the system
    // a recount()/rebuild_cells() before any energy query -- same contract set_pose(i, ...)
    // already has for an existing particle.
    int add_particle(double x, double y, double theta, int species_index, bool frozen = false) {
        assert(species_index >= 0 && species_index < static_cast<int>(species.size()));
        Particle p;
        p.species = species_index;
        p.frozen = frozen;
        p.set_pose(x, y, theta);
        particles.push_back(p);
        return static_cast<int>(particles.size()) - 1;
    }

    // The energy, unweighed and in two pieces, kept current by every move, so changing a
    // stiffness costs no recount. resync_energy() checks them against the poses, up to the drift
    // of a long sum of increments.
    EnergyParts parts;

    // The one energy there is: the tracked sums at this system's own stiffnesses.
    double energy() const { return energy_of(parts, k_pair, k_wall); }

    // Which particles are near which, so a single-particle energy looks at its neighbours and
    // not at all n. Kept current by set_pose(i, ...) below and rebuilt by recount(), and only
    // particle_energy() relies on it: the whole-system passes build their own, so they stay
    // right even after a caller has written poses behind its back -- which is what lets
    // resync_energy() check this one against them.
    //
    // The contract that buys: every pose write to a particle of a System goes through
    // set_pose(i, ...), or is followed by recount() or rebuild_cells() before the next
    // particle_energy(). resync_energy() asserts it.
    CellList cells;

    // The pose write that keeps `cells` current. O(1). Before the first recount() there is no
    // grid to keep, and particle_energy() refuses to run until there is one.
    void set_pose(int i, double x, double y, double theta) {
        particles[static_cast<size_t>(i)].set_pose(x, y, theta);
        if (cells.cell_of.size() == particles.size()) cells.move(i, x, y);
    }

    void rebuild_cells();

    // One pass over the close pairs, O(n). total_energy() is this weighed by the system's own stiffnesses, so the
    // two cannot disagree.
    EnergyParts energy_parts() const;

    double total_energy() const;

    // Centre distance beyond which two shapes cannot overlap: the cell width the grids need.
    double reach() const;

    // The energy terms that change when particle i alone moves, unweighed.
    EnergyParts particle_parts(int i) const;

    double particle_energy(int i) const { return energy_of(particle_parts(i), k_pair, k_wall); }

    double max_overlap() const;

    bool is_feasible(double tol) const;

    double tight_side() const;

    int mobile_count() const;

    void set_box(double new_L);

    void rescale_box(double new_L);

    // Slides the squares so their bounding box starts at the origin and makes L its tight side.
    void fit_box();

    void resync_energy();

    // Rebuild the grid and recount both from the configuration. For the moves that change everything
    // at once, and for anyone who has just written poses directly and owes the system an energy.
    void recount();
};

inline double enthalpy(const System& s, double pressure) {
    return s.energy() + pressure * s.L * s.L;
}

}  // namespace pack
