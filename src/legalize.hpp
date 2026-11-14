#pragma once

// Legalization: turn a near-feasible packing of unit squares (overlaps up to kCertifyTol in area,
// i.e. up to ~1e-5 in depth) into one that tools/validate.py accepts -- no positive-area overlap
// and nothing outside [0, side]^2, judged exactly on the float corners.
//
//   1. If every close pair is already separated (strict_report below), only the box is refit.
//   2. Otherwise relax a copy whose squares are enlarged by kInflate (relative), so a relaxed copy
//      leaves the real squares ~kInflate/2 apart -- far above float error. If the copy cannot be
//      separated in its box, the box grows by a tiny geometric step and the copy is relaxed again.
//   3. Refit: shift the squares so the lowest corner sits a few ulps above 0 and set the side to the
//      extent plus the same margin, so float rounding of corners cannot cross the wall.
//
// The validator, not this code, is the judge: this pass only has to make its verdict "valid" with
// overwhelming probability. Phase 3 item 2 is expected to replace step 2 with something smarter.

#include "system.hpp"

namespace pack {

inline constexpr double kInflate = 1e-11;  // relative enlargement of the squares in step 2
inline constexpr double kSepTol = 1e-13;   // gap a tilted pair needs in the engine's float check

struct StrictReport {
    double min_pair_gap = 0.0;  // smallest SAT gap over close pairs (negative = penetration)
    double min_wall_gap = 0.0;  // smallest distance from a corner to the wall (negative = outside)
    int bad_pairs = 0;          // pairs failing the float criterion
    int bad_walls = 0;          // squares failing the float wall criterion
    bool ok() const { return bad_pairs == 0 && bad_walls == 0; }
};

// The engine's own strict check, in double. Axis-aligned pairs (theta exactly 0 for both) pass on
// gap >= 0, since their corners are computed exactly; any other pair needs gap > kSepTol.
StrictReport strict_report(const System& sys);

struct LegalizeResult {
    bool ok = false;            // strict_report passed after the final refit
    double side_before = 0.0;   // tight_side() on entry
    double side_after = 0.0;    // the legal side (sys.L on return)
    int growth_steps = 0;       // box growths step 2 needed
    long relax_iterations = 0;  // L-BFGS iterations spent in step 2
    StrictReport report;
};

// Unit squares only (base species, nothing frozen or merged). On ok, sys is the legal packing with
// sys.L its side; otherwise sys is left at the last attempt and `ok` is false.
LegalizeResult legalize(System& sys, int max_growth = 24);

}  // namespace pack
