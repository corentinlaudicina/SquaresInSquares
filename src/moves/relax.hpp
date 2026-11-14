#pragma once

// Local minimisation of the overlap energy with L-BFGS, using the analytic gradient.
//
// Two functions:
//
// relax() : Removes overlaps in a fixed box and answers one question: is this configuration, in this box, a valid packing? It only finds the minimum of the  basin it starts in, so a failure means "not from this start", not "impossible".

// compress() : Calls relax() in a loop, shrinking the box step by step, relaxing from the current 1configuration each time, and halves the step when relax() fails. The last success is a packing with a proven side.

#include "system.hpp"

#include <functional>
#include <vector>

namespace pack {

inline constexpr double kCertifyTol = 1e-10;

// Called with a configuration and what is happening to it, so a run can be watched. Only ever reads, so it cannot change the search.
using FrameWatch = std::function<void(const System&, const char*)>;

double energy_grad(const System& sys, std::vector<Pose>& grad);

struct RelaxOptions {
    long max_iterations = 5000;
    
    int history = 0;
    double feasible_tol = kCertifyTol;  // largest overlap area in a valid packing

    // Stopping tests for a configuration that cannot reach feasible_tol (jammed).
    double grad_tol = 1e-14;    // stop (Converged) when no gradient component exceeds this
    double energy_tol = 1e-12;  // relative energy decrease that counts as no progress
    int stall_patience = 3;     // stop (Stalled) after this many no-progress iterations in a row

    // Every this many iterations, stop (Hopeless) if the recent rate of energy decrease would
    // need more than ten times the remaining budget to reach a valid packing. 0 = off: a large
    // move wants the relaxed state whatever it is, only compress() wants a verdict fast.
    long hopeless_window = 0;

    FrameWatch watch;                  // if set, called every watch_every iterations
    long watch_every = 10;
    const char* watch_label = "relax";
};

struct RelaxResult {
    enum class Stop {
        Feasible,       // overlaps within feasible_tol: a valid packing
        Converged,      // gradient below grad_tol with overlaps left: jammed
        Stalled,        // line search failed, or stall_patience steps without progress
        MaxIterations,  // hit max_iterations
        Hopeless,       // too slow to reach feasible_tol in the budget left; see hopeless_window

        Escaped,
    };

    Stop stop = Stop::MaxIterations;
    bool feasible = false;  // max_overlap <= feasible_tol; the main result
    long iterations = 0;
    long evaluations = 0;   // energy evaluations, including line search trials
    double energy = 0.0;
    double max_overlap = 0.0;
};

const char* stop_name(RelaxResult::Stop stop);

// Minimises the energy in the current box (side sys.L). Leaves sys at the lowest point reached, with sys.parts set to match. Frozen particles never move.
RelaxResult relax(System& sys, const RelaxOptions& opt = {});

struct CompressOptions {
    double shrink = 0.995;        // box multiplier for the first squeeze step (0.5% smaller)
    double grow = 1.05;           // box multiplier in phase one, until a box certifies
    double min_shrink = 1e-7;     // stop once the relative step is at or below this
    long max_relaxations = 500;   // relax() calls for both phases; checked in phase two only
    long max_growth_stages = 40;  // box increases allowed in phase one
    // A shrink that fails is the common case here, so the relaxations give up on a jammed box
    // early instead of spending max_iterations on it.
    RelaxOptions relax{.hopeless_window = 200};
};

struct CompressResult {
    System packing;      // smallest configuration that certified; empty if none did
    bool found = false;  // true if any box certified
    double side = 0.0;   // tight_side() of packing: the claimed side

    double final_step = 0.0;
    bool squeezed_out = false;  // phase two ended on min_shrink, not on max_relaxations
    long stages = 0;            // shrink steps accepted
    long relaxations = 0;       // relax() calls, both phases
};

CompressResult compress(const System& start, const CompressOptions& opt = {});

// Minimises with the wall penalty off, so squares may spread past the box, then refits the box
// around them. Stiffnesses are back as they were on return.
void relax_free(System& sys, const RelaxOptions& opt = {});

// The compression schedule: minimise, and while the result is valid shrink the box by ds (a
// length) and minimise again. A shrink that fails is retried from the last valid packing with
// half the step. Runs at most max_iter minimisations once a valid box is found. On success sys is
// the smallest valid packing; if no box was ever valid, sys is unchanged and found is false.
// Not the Monte Carlo anneal() in moves/mc.hpp, which takes a temperature Schedule.
CompressResult anneal(System& sys, double ds, long max_iter, const FrameWatch& watch = {});

struct BestOfResult {
    CompressResult best;        // the smallest certified side among the candidates
    int winner = -1;            // index of the candidate it came from; -1 if none certified
    std::vector<double> sides;  // per candidate: certified side, or 0 if it certified nothing
};

// Compresses every candidate and keeps the smallest certified side; ties go to the earlier
// candidate. compress() is a local polish -- it cannot change an arrangement, only fit it --
// so which arrangement it is handed decides the answer, and trying several is the cheap way
// to stop a single unlucky pick from deciding it.
//
// The candidates run on up to `threads` threads (<= 0: the hardware's count). compress() is
// deterministic and shares nothing, so the threads cannot change the answer.
BestOfResult compress_best_of(const std::vector<System>& candidates,
                              const CompressOptions& opt = {}, int threads = 0);

}  // namespace pack