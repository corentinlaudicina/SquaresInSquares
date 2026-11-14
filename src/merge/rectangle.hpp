#pragma once

// Detects whether a group of touching squares, at their current poses, exactly tiles an axis-aligned (relative to each other) rectangle -- the one shape merge() (merge.hpp) is willing to turn them into. A convex result for a convex input, so the clipping in geometry/polygon.cpp never has to handle anything concave: an arbitrary touching cluster (an L-shape, a staircase) simply isn't eligible, on purpose.

#include "system.hpp"

#include <optional>
#include <vector>

namespace pack {

struct RectangleUnion {
    double side = 0.0;                   // the member squares' own edge length
    int cols = 0, rows = 0;              // grid shape; cols * rows == the member count
    double theta = 0.0;                   // the merged rectangle's orientation
    Vec2 center{0.0, 0.0};                // the rectangle's centre, from the ideal grid -- not
                                           // a centroid of the members' (possibly jittery) poses
    std::vector<int> ordered_members;     // `members`, reordered row-major:
                                           // ordered_members[row * cols + col]
};

// Checks `members` (indices into sys.particles) against every condition a mergeable rectangle has
// to satisfy; returns the first failure's reason via std::nullopt (never mutates sys):
//   - at least 2 members, in range, no duplicates
//   - all one species, and that species is a square
//   - all at one orientation, up to a multiple of 90 degrees (a square at 0 and at 90 is the
//     same shape)
//   - on an exact grid at the square's own pitch, no gaps or stray offsets past pos_tol/theta_tol
//   - the grid cells are all distinct and exactly fill a cols x rows rectangle -- what rejects
//     a non-rectangular touching cluster of the right count but the wrong shape
std::optional<RectangleUnion> detect_rectangle_union(const System& sys,
                                                      const std::vector<int>& members,
                                                      double pos_tol = 1e-6,
                                                      double theta_tol = 1e-6);

}  // namespace pack
