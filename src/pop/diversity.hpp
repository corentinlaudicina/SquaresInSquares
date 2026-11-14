#pragma once

// Diversity-aware selection (Phase 3 item 5), after jcottaar/packing (CC BY-SA 4.0), reimplemented.
//
// Distance between two packings in the same box:
//   Assign  squares are interchangeable, so match them: min over the 8 box symmetries of the optimal
//           assignment (Hungarian, O(N^3)) of  |dp|^2 + (R * dtheta)^2,  dtheta wrapped to 90 degrees
//           and R = 1/sqrt(2); reported as the RMS over squares, i.e. "how far a square moves", in
//           length units.
//   Hist    the cheap proxy for large N: per cell of a k x k grid (k ~ sqrt(N)/2), the square count
//           and the mean of (cos 4 theta, sin 4 theta); min over the 8 symmetries of the RMS
//           difference per square. O(N + k^2). Not the same scale as Assign: its own threshold.

#include "pop/penalty.hpp"

#include <functional>
#include <vector>

namespace pack {

enum class Metric { Auto, Assign, Hist };

struct DiversityParams {
    double elite = 0.25;            // share of survivors taken purely by fitness
    double threshold = 0.05;        // Assign: minimum RMS square displacement to count as different
    double hist_threshold = 0.05;   // Hist: minimum histogram distance
    Metric metric = Metric::Auto;   // Auto: Assign up to assign_max_n squares, Hist above
    int assign_max_n = 40;
};

double assign_distance(const Packing& a, const Packing& b);
double hist_distance(const Packing& a, const Packing& b);

// `order` lists candidate indices best first. Returns `count` of them: the elite share first, then
// in fitness order each candidate farther than the threshold from all already chosen; if that runs
// out, the rest in fitness order. `distance_evals`, if given, counts distance computations.
// `parfor`, if set, runs f(0..n-1) in parallel (any order); the distances from one candidate to the
// chosen set then run in parallel, which does not change the result.
using ParallelFor = std::function<void(int, const std::function<void(int, int)>&)>;
std::vector<int> select_diverse(const std::vector<Packing>& cand, const std::vector<int>& order, int count,
                                const DiversityParams& prm, long* distance_evals = nullptr,
                                const ParallelFor& parfor = {});

}  // namespace pack
