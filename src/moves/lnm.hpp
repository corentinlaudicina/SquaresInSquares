#pragma once

// Large-neighbourhood moves. They exist because single-particle moves cannot leave a basin. A
// cold run cannot slide two squares past each other, because every path between the two
// arrangements goes through an overlap it will not pay for. These moves change the arrangement at the scale the good solutions differ at: a diagonal band
// of tilted squares, a band down one side with a corner cut out of it.
//
// All four are built the same way, and the shape is the point:
//
//     copy the system -> perturb it, hard -> relax the perturbation -> accept on enthalpy
//
// The perturbation is allowed to be terrible. It can leave squares deeply overlapped or stacked
// on one another, because relax() (moves/relax.hpp) follows the analytic gradient downhill
// afterwards and the Metropolis test throws away whatever that does not rescue. That division is
// what makes the moves cheap to write: a proposal only has to be good often enough to be worth
// the relaxation.
//
// WHAT IT COSTS, and this is not a detail. relax() is a many-to-one map -- a whole basin goes to
// one minimum -- so the reverse of one of these moves does not have the probability the forward
// one does, and **these moves do not satisfy detailed balance**. A run using them is basin
// hopping: an optimiser over local minima, not a sampler of the fixed-pressure ensemble. For a
// records search that is the right trade -- the object of the exercise is the smallest side
// anyone can certify, not an unbiased average -- but it is a real change in what a run means.
//
// Cost and budget. One attempt relaxes the whole system, which is a few hundred O(n^2) passes:
// on the order of a hundred sweeps, for one proposal. So they are attempted once every few
// hundred sweeps, and the relaxation runs on a proposal budget (LnmOptions::relax_iterations)
// at kFeasibleTol rather than kCertifyTol. Certifying is compress()'s job, at the end, once.
//
// HOW TO READ THE ACCEPTANCE RATES, because the obvious reading is wrong. The proposal has been
// relaxed and the state it is weighed against has not: a state at temperature T holds a
// sampled configuration, overlaps and all, while the proposal has been taken downhill. The
// comparison is loaded in the proposal's favour, so these moves are accepted most of the time --
// 90% and up is ordinary -- and a high rate says nothing about whether the perturbation was any
// good. It says a relaxation pays for most of what a perturbation can break. The informative
// rate is a low one: the corner move sits near 5%, and that is a perturbation even a relaxation
// could not pay for.
//
// It also means each accepted move quenches the state into a local minimum, which the next few
// hundred sweeps at temperature then melt again. The textbook basin-hopping alternative
// is to relax both sides and run the chain on minima, so the Metropolis test compares like with
// like and the rates mean what they look like they mean. That costs a second relaxation per
// attempt and is the obvious next thing to try here.

#include "moves/mc.hpp"
#include "moves/regions.hpp"
#include "moves/relax.hpp"
#include "moves/symmetry.hpp"
#include "system.hpp"

#include <random>

namespace pack {

struct LnmOptions {
    RegionWeights regions;

    // Relative weights for lnm_round's draw. Destroy and repair carries the most because it is
    // the one that works at any count and any density; the strip move is next because the
    // diagonal tilted band is the structure the good solutions have and random moves do not
    // find; loosen is the cheapest to get wrong and the least targeted, so it carries least.
    double w_destroy = 0.4, w_strip = 0.3, w_corner = 0.2, w_loosen = 0.1;

    // Off by default: until a benchmark shows it earns its cost, a symmetrize attempt is a full
    // compress for nothing. At 0 the draw and the random stream are what they were without it.
    double w_symmetrize = 0.0;
    SymmetrizeOptions symmetrize;

    // Relaxation budget for one proposal. Not a certification: see the note above.
    long relax_iterations = 300;

    // Reinsertion: best of this many random poses, scored one square at a time.
    int insert_tries = 12;
    // Chance a reinserted square takes a multiple of 45 degrees rather than a free angle. The
    // two families a packing of squares is built from are 45 degrees apart, so half the draws go
    // straight to one of them instead of landing between the two where nothing fits.
    double snap_prob = 0.5;

    double loosen_grow = 0.03;      // box grows by this fraction before the shake
    long loosen_sweeps = 20;        // sweeps of shaking while the box is loose
    double loosen_soften = 0.1;     // stiffness multiplier for the soft-to-rigid variant
    double loosen_soft_prob = 0.5;  // chance an attempt is the soft-to-rigid kind

    MoveLog* log = nullptr;  // per-move table; never touches the random stream
    FrameWatch watch;        // sees each attempt, its relaxation and the verdict
};

// Per-move-type counters. Kept apart rather than summed because the interesting question about a
// move mix is which move is doing the work, and a move that never fires should be visible as a
// zero rather than hidden inside a total.
struct LnmStats {
    Stats destroy, corner, strip, loosen, symmetrize;
    SymmetrizeBandit bandit;      // lives here so it learns across calls; see large_neighbour_move
    MoveId last = MoveId::Count;  // the move type of the most recent accepted attempt
};
Stats lnm_total(const LnmStats& s);

// Delete every square in a region, put them back into the hole, relax. Usually the most
// productive single large move: it is the only one that does not assume what the answer looks
// like. The region comes from pick_region (regions.hpp), so it is usually a band or a corner
// block rather than a disc -- which is the one piece of knowledge about good packings this move
// does use.
bool try_destroy_repair(System& sys, const LnmOptions& opt, double pressure, double temperature,
                        std::mt19937_64& rng);

// Replace a corner block with a differently shaped one: re-lay its squares as a fresh staircase,
// axis-aligned, at contact pitch, with a randomly drawn row profile. The group that set the
// s(155)-s(264) records chose these corner blocks with an integer program; this samples the
// shape instead, which is the Monte Carlo version of the same idea and costs a draw rather than
// a solver.
bool try_reshape_corner(System& sys, const LnmOptions& opt, double pressure, double temperature,
                        std::mt19937_64& rng);

// Re-lay a band as a tiled strip. Two directions, drawn evenly, because a move that only ever
// made tilted bands would be a one-way street:
//
//   carve       a 45 degree band, every square tilted to match, at contact pitch along the
//               diagonal. This builds the diagonal tilted band by hand rather than hoping the
//               sampler stumbles onto it.
//   straighten  an axis-aligned band, every square square-on. The inverse, for a band that
//               should not have been tilted.
bool try_strip(System& sys, const LnmOptions& opt, double pressure, double temperature,
               std::mt19937_64& rng);

// Grow the box a few percent, shake, squeeze back, relax. Half the attempts also soften both
// stiffnesses for the shake and harden them again before the squeeze, which lets the squares
// pass through each other while they rearrange -- the soft-to-rigid variant, and the only move
// here that can swap two squares across a contact they could never cross at full stiffness.
bool try_loosen(System& sys, const LnmOptions& opt, double pressure, double temperature,
                std::mt19937_64& rng);

// The reinsertion destroy-and-repair uses, for moves that need to put a square back somewhere:
// lift the squares out of the way, then each takes the best of opt.insert_tries poses in `b`.
void park_squares(System& sys, const std::vector<int>& members);
void reinsert_in_bounds(System& sys, int i, const LnmOptions& opt, const Bounds& b,
                        std::mt19937_64& rng);

// One weighted attempt among the moves with positive weight.
bool lnm_round(System& sys, const LnmOptions& opt, LnmStats& stats, double pressure,
               double temperature, std::mt19937_64& rng);

// One large move used as a deliberate kick: always applied, whatever it does to the energy.
// False only if every attempt declined, which leaves sys as it was.
// `stats`, if given, is where the counters go, so a caller can see which move it got.
bool large_neighbour_move(System& sys, std::mt19937_64& rng, const LnmOptions& opt = {},
                          LnmStats* stats = nullptr);

}  // namespace pack
