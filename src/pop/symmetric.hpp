#pragma once

// Symmetric genotypes (Phase 3 item 8), after jcottaar/packing (CC BY-SA 4.0), reimplemented.
//
// A C2 (180 degree) or C4 (90 degree) symmetric packing about the box centre c = (h/2, h/2) is stored
// as its generators: one square per orbit, plus a centre square when N leaves a remainder of 1 (it
// sits exactly at c with a free angle; a square is unchanged by a quarter or half turn). An image of
// a generator is its centre offset turned k times, with the same angle -- turning a square by 90 or
// 180 degrees gives the same square. So:
//   C2: N = 2m (+1),  any N;     C4: N = 4m (+1),  only N mod 4 in {0, 1}.
// Evaluation expands to the full N squares, takes the usual cost and gradient, and pulls the gradient
// back: an image's position gradient turned back onto its generator, angles summed; the centre
// square's position gets none, so it never moves.
//
// Generators are kept in the fundamental region (C2: the lower half, C4: the lower-left quadrant) for
// crossover, which swaps generators near a point of that region. Order: orbit generators, centre last.

#include "pop/moves.hpp"
#include "pop/relax_capped.hpp"

namespace pack {

enum class Sym { None, C2, C4 };

struct SymLayout {
    Sym sym = Sym::None;
    int n = 0;            // full squares
    int orbit = 1;        // 2 or 4
    int m = 0;            // generators, including the centre one
    bool centre = false;  // the last generator is the fixed centre square
};

// Sym::None in the result when the symmetry cannot hold n squares (C4 with n mod 4 in {2, 3}).
SymLayout make_layout(Sym sym, int n);

void expand(const Packing& gen, const SymLayout& L, Packing& full);
Packing expanded(const Packing& gen, const SymLayout& L);
void pull_gradient(const double* gfull, const SymLayout& L, double* ggen);
void fix_centre(Packing& gen, const SymLayout& L);
void canonicalize(Packing& gen, const SymLayout& L);

// Generators from any packing of L.n squares: the squares' images in the fundamental region, near
// duplicates (orbit mates of an already symmetric packing) merged; the centre square is the one
// nearest the centre. Exact for a packing that already has the symmetry; a starting guess otherwise.
Packing symmetrize(const Packing& full, const SymLayout& L);

// Region-swap crossover on generators (see above); the mate may first be mapped by a box symmetry.
void crossover_sym(Packing& child, const Packing& mate, const MoveParams& prm, std::mt19937_64& rng, const SymLayout& L);

CostInfo cost_symmetric(const Packing& gen, const SymLayout& L, Penalty pen, Workspace& ws, double* ggen);
CappedStats relax_symmetric(Packing& gen, const SymLayout& L, const CappedParams& prm, Workspace& ws);

}  // namespace pack
