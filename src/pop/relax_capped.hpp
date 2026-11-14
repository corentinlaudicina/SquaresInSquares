#pragma once

// Staged L-BFGS for a population: short history, a hard cap on every variable's step, no line
// search. One gradient evaluation per iteration, which is what makes relaxing thousands of
// perturbed packings affordable. Stages run in order with fresh history, a rough one (big cap) to
// untangle a perturbation, then fine ones (smaller caps) to settle the contacts. The best point
// seen is what is returned, so an overshooting step can cost time but never quality.
//
// Variables are x, y and R * theta (R = 1/sqrt(2), the circumradius), so a cap is a length for all
// three. Idea from jcottaar/packing (CC BY-SA 4.0); reimplemented here for squares.

#include "pop/penalty.hpp"

#include <vector>

namespace pack {

struct RelaxStage {
    double cap;  // largest step of any single variable, in length units
    int iters;
};

struct CappedParams {
    Penalty penalty = Penalty::Depth2;
    std::vector<RelaxStage> stages = {{0.05, 150}, {0.0125, 150}, {0.003, 150}, {0.0008, 150}};
    int history = 3;
    // Stop as soon as every pair depth and wall protrusion is below this (Depth2), or every overlap
    // area below it (Area2).
    double tol = 1e-12;
};

struct CappedStats {
    int iters = 0;
    bool reached = false;  // max_depth and max_wall within tol
    CostInfo info;         // at the returned point
};

CappedStats relax_capped(Packing& p, const CappedParams& prm, Workspace& ws);

}  // namespace pack
