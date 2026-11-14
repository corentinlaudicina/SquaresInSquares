#pragma once

#include "system.hpp"

#include <random>
#include <vector>

namespace pack {

enum class Region { Disc, Band, Corner };

const char* region_name(Region r);

// Face-to-face thickness along the shape's first SAT axis: the pitch at which two copies at the same orientation, stacked along that axis, touch exactly. 1.0 for a unit square. Every move that re-lays squares in a row uses it, and it is measured rather than assumed to be 1 so that a different shape or a non-unit side does not silently lay the rows at the wrong spacing.
double contact_pitch(const Shape& s);

// Indices whose centres lie in the region. Frozen squares are skipped, so a region holding one still moves the rest -- unlike the group moves in mc.hpp, which refuse a group outright. A region is a place in the box, not a chosen set of particles, and a pinned square in the middle of it is a reason to work around it rather than to abandon the move.
std::vector<int> select_disc(const System& sys, double cx, double cy, double r);

// `dir` need not be normalized. The band is the slab of normal coordinate [offset, offset + width], infinite along dir: it is clipped by the box, not by this function.
std::vector<int> select_band(const System& sys, Vec2 dir, double offset, double width);

// corner 0..3, counter-clockwise from the origin: (0,0), (L,0), (L,L), (0,L). The block is the depth x depth square at that corner.
std::vector<int> select_corner(const System& sys, int corner, double depth);

// Sizes are in square diameters (2 * circumradius) so they mean the same thing for any shape.
struct RegionWeights {
    // Relative, normalized on use: bands and corners beat discs for the reason at the top.
    double disc = 0.2, band = 0.5, corner = 0.3;

    double disc_radius_lo = 1.0, disc_radius_hi = 3.0;
    double band_width_lo = 1.0, band_width_hi = 3.0;
    double corner_depth_lo = 2.0, corner_depth_hi = 4.0;

    // Band direction, drawn from the four lattice directions a packing of squares actually uses.
    // The diagonals carry twice the weight of the axes: an axis-aligned band is what a squeeze
    // produces by itself, and the diagonal one is the arrangement the sampler has to be handed.
    double axis_weight = 1.0, diagonal_weight = 2.0;
};

// One weighted draw. `members` may come out empty -- a band can miss every square -- and every caller has to treat that as "no move", not as an error.
struct Pick {
    Region kind = Region::Disc;
    std::vector<int> members;

    // The geometry that defined the region, so the move does not have to guess it back:
    //
    //   centre  disc: its centre. band: a point on its midline. corner: the corner itself.
    //   dir     band: along the band, unit. disc, corner: unit, otherwise unused.
    //   axis    a line through `centre` that the region is symmetric about, unit.
    //   extent  disc radius, band width, or corner depth, in length units (not diameters).
    Vec2 centre{0.0, 0.0};
    Vec2 dir{1.0, 0.0};
    Vec2 axis{1.0, 0.0};
    double extent = 0.0;
    int corner = 0;  // only meaningful for Region::Corner
};

Pick pick_region(const System& sys, const RegionWeights& w, std::mt19937_64& rng);

// A band in a named direction, for the moves that need one specific kind rather than a weighted
// draw: the strip move re-lays a diagonal band and has no use for a disc.
Pick pick_band(const System& sys, const RegionWeights& w, Vec2 dir, std::mt19937_64& rng);

// Axis-aligned bounding box of the members' centres, padded by `pad` and clipped to the box.
// Where destroy-and-repair puts the squares it lifted out: inside the hole it just made, not
// anywhere in the box, which is what makes it a local move.
struct Bounds {
    double xlo = 0.0, xhi = 0.0, ylo = 0.0, yhi = 0.0;
};
Bounds member_bounds(const System& sys, const std::vector<int>& members, double pad);

}  // namespace pack
