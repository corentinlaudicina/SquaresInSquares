#pragma once

// Seeding the GA from other structures (Phase 3 item 7):
//   up       best/(N+1) with one square removed: the loosest (delete_rattler) and the `up_count`
//            most-contacted ones (delete_for_compress) -- what compress can gain from losing a square
//            is only known by trying;
//   down     best/(N-1) with one square inserted at the best of `down_tries` random poses (depth^2);
//   lattice  grid_init and diagonal_init (init.hpp): the square grid and a stacked tilted strip.
// Every candidate is legalized (legalize.hpp), so each has a valid side to start its ratchet from.
// The main route for large N, where the GA cannot invent a whole structure from one start.

#include "system.hpp"

#include <random>
#include <string>
#include <vector>

namespace pack {

struct SeedParams {
    bool up = false;
    int up_count = 2;
    bool down = false;
    int down_tries = 400;
    bool lattice = false;
    std::string neighbor_dir = "best";  // where N+1 and N-1 are read from (sqpack/1 JSON)
};

struct SeedCandidate {
    std::string source;  // "up:rattler", "up:contact1", "down", "lattice:grid", "lattice:diagonal"
    System sys;          // legal, N squares, sys.L its side
};

// Unit squares only. Draws from rng only when `down` is on (and its neighbour exists).
std::vector<SeedCandidate> seed_candidates(int n, const SeedParams& prm, std::mt19937_64& rng);

// One square added to `sys` (N-1 squares) at the lowest-penalty of `tries` random poses in its box.
System insert_square(const System& sys, int tries, std::mt19937_64& rng);

}  // namespace pack
