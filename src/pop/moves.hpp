#pragma once

// GA offspring moves (Phase 3 item 4). Every move is allowed to leave overlaps: the relaxation that
// follows repairs them, and the pre-filter and selection throw away what it cannot repair. Move set
// after jcottaar/packing (Santa 2025 1st place, CC BY-SA 4.0), reimplemented for squares.
//
// Randomness is drawn with the portable helpers below (bits of mt19937_64, Box-Muller), not with
// std::*_distribution, whose algorithms differ between libc++ and libstdc++ (audit B6): the same seed
// gives the same offspring on the Mac and on Kaggle.

#include "pop/penalty.hpp"

#include <array>
#include <random>

namespace pack {

inline double u01(std::mt19937_64& rng) { return static_cast<double>(rng() >> 11) * 0x1p-53; }
inline double uniform(std::mt19937_64& rng, double lo, double hi) { return lo + (hi - lo) * u01(rng); }
inline int uniform_int(std::mt19937_64& rng, int lo, int hi) {  // inclusive; tiny modulo bias is irrelevant here
    return lo + static_cast<int>(rng() % static_cast<unsigned long long>(hi - lo + 1));
}
double normal(std::mt19937_64& rng);  // standard normal, Box-Muller

enum class GaMove { Relocate, JiggleSmall, JiggleLarge, Cluster, Twist, Translate, Crossover, Count };
inline constexpr int kGaMoves = static_cast<int>(GaMove::Count);
const char* ga_move_name(GaMove m);

struct MoveParams {
    std::array<double, kGaMoves> weight = {0.10, 0.25, 0.15, 0.15, 0.10, 0.10, 0.15};
    double snap_prob = 0.5;        // relocate: chance the new angle is a multiple of 45 degrees
    double jiggle_small = 0.02;    // sigma of the position (length) and angle (radians) noise
    double jiggle_large = 0.15;
    int cluster_min = 2, cluster_max = 8;
    double cluster_sigma = 0.08;
    double twist_radius = 3.0;     // twist: squares within this of a random centre turn about it ...
    double twist_angle = 0.25;     // ... by up to this (radians), fading linearly to 0 at the radius
    double translate_radius = 2.5;
    double translate_sigma = 0.15;
    double crossover_radius = 2.5;
    double crossover_symmetry = 0.5;  // chance the mate is first mapped by a random D4 element
};

GaMove draw_move(const MoveParams& prm, std::mt19937_64& rng);

// Applies m to `child` in place. `mate` is used by Crossover only (and must then be non-null).
void apply_move(Packing& child, GaMove m, const MoveParams& prm, std::mt19937_64& rng, const Packing* mate);

// The 8 symmetries of the box about its centre; k in [0, 8). Angles map exactly mod 90 degrees.
void apply_d4(Packing& p, int k);

}  // namespace pack
