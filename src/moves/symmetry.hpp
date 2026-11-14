#pragma once

// The symmetries of the container, the group D4, and what the search does with them.
//
// Every element acts about the box centre (L/2, L/2): the box itself is [0,L]^2, so the origin is a corner and the centre has to be subtracted first. Angles are radians and a square is only defined mod 90 degrees, which is what makes the angle rule so short: a rotation leaves theta alone and any of the four mirrors sends it to -theta. A mirror in the line at angle phi really sends theta to 2 phi - theta, and 2 phi is a multiple of 90 degrees for all four axes of the box. The rule is therefore true for squares and false for anything else, so everything here refuses a configuration holding any other shape (see is_square_config).

#include "moves/move_log.hpp"
#include "moves/relax.hpp"
#include "moves/regions.hpp"
#include "system.hpp"

#include <array>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace pack {

// The 2x2 signed permutation [a b; c d] and the sign theta is multiplied by.
struct GroupElem {
    int a = 1, b = 0, c = 0, d = 1;
    int angle_sign = 1;

    Pose apply(const Pose& p, double L) const;
    GroupElem inverse() const;
    // this after g: (f.compose(g)).apply(p) == f.apply(g.apply(p)).
    GroupElem compose(const GroupElem& g) const;

    bool operator==(const GroupElem&) const = default;
};

enum class D4 { Id, Rot90, Rot180, Rot270, MirrorX, MirrorY, MirrorDiag, MirrorAnti };

// MirrorX negates x (the mirror line is vertical), MirrorY negates y, MirrorDiag swaps x and y, MirrorAnti sends (x, y) to (-y, -x), all about the centre.
GroupElem group_elem(D4 which);
const std::array<GroupElem, 8>& d4_elements();

// Angle wrapped into [-45, 45) degrees, and the periodic distance between two orientations.
double wrap90(double theta);
double angle_distance90(double a, double b);

// True when D4 applies: every particle is the base species and the base species is a square.
bool is_square_config(const System& sys);

// The whole configuration's image, frozen flags travelling with their squares.
System image_of(const System& sys, const GroupElem& g);

// The largest subgroup of D4 all of whose elements map the configuration onto itself: same positions within tol, same angles within tol mod 90 degrees, one-to-one. Always holds at least the identity. Each element is tested on its own and a tolerance is not transitive, so the elements that fit need not be closed under composition -- hence a subgroup search rather than simply returning them.
std::vector<GroupElem> detect_symmetry(const System& sys, double tol);
bool contains(const std::vector<GroupElem>& group, const GroupElem& g);

// The sorted, quantised (x, y, theta) of a configuration, smallest of its eight images. Two configurations related by D4 and a relabelling have equal keys.
using CanonicalKey = std::vector<std::array<long long, 3>>;
CanonicalKey canonical_key(const System& sys);

// The image with that smallest key, squares in key order, angles wrapped.
System canonicalize(const System& sys);

// --- symmetrize_region ----------------------------------------------------------------------
//
// Copy a good part of the box onto a bad part through a group element, repair the seams, and keep the result only if it certifies at a side no worse than the one it started from. It is a basin-hopping move with its own compress, so it is far more expensive than the other large moves and has no detailed balance to preserve.

struct LnmOptions;

// A variant is {copies 1, 2, 3} x {rotation, axis mirror, diagonal mirror}:
//
//   1 copy    one region replaced by the image of another (quadrant or diagonal sector)
//   2 copies  half the box replaced by the image of the other half
//   3 copies  a full orbit from one quadrant, or from one diagonal sector
inline constexpr int kSymmetrizeVariants = 9;
enum class SymFamily { Rotation, AxisMirror, DiagonalMirror };
int symmetrize_variant(int copies, SymFamily family);
std::string symmetrize_variant_name(int variant);

struct SymmetrizeOptions {
    double eps = 0.2;           // chance the regions are drawn at random instead of by waste
    double floor = 0.05;        // bandit: added to every variant's rate, so none starves
    double decay = 0.1;         // bandit: weight of the newest attempt in the recent rate
    double seam_width = 1.5;    // centre distance, in pitches, that counts as across a seam
    long seam_iterations = 300; // window relaxation budget
    double shrink = 0.998;      // the repair's compress: first step
    double grow = 1.01;         // ... box growth while the image does not fit
    long growth_stages = 3;     // ... at most this many growths, then the attempt is abandoned
    long relaxations = 60;      // ... total relax() calls
};

struct SymmetrizeBandit {
    struct Arm {
        long attempts = 0;
        long improvements = 0;
        double cpu_seconds = 0.0;
        double recent = 0.0;  // exponentially weighted improvement rate
    };
    std::array<Arm, kSymmetrizeVariants> arms{};

    int choose(std::mt19937_64& rng, double floor) const;
    void update(int variant, bool improved, double seconds, double decay);
};

// Accepts only a packing that passes is_feasible(kCertifyTol), with exactly the same number of squares, at a side no larger than sys's own. On rejection sys is untouched, bit for bit. Frozen squares never move. `improvement`, if given, gets the decrease of the side.
bool try_symmetrize_region(System& sys, const LnmOptions& opt, SymmetrizeBandit& bandit,
                           std::mt19937_64& rng, double* improvement = nullptr);

// --- symmetry operations called directly ------------------------------------------------------

// A part of the box, by where its centre-test puts a square. Regions are closed, so a square on a boundary belongs to both sides.
enum class RegionKind { Quadrant, Half, DiagonalHalf, Sector };
struct BoxRegion {
    RegionKind kind = RegionKind::Quadrant;
    int index = 0;
};

inline constexpr BoxRegion kBottomLeftQuadrant{RegionKind::Quadrant, 0}, kBottomRightQuadrant{RegionKind::Quadrant, 1}, kTopRightQuadrant{RegionKind::Quadrant, 2}, kTopLeftQuadrant{RegionKind::Quadrant, 3};
inline constexpr BoxRegion kLeftHalf{RegionKind::Half, 0}, kRightHalf{RegionKind::Half, 1}, kBottomHalf{RegionKind::Half, 2}, kTopHalf{RegionKind::Half, 3};
inline constexpr BoxRegion kTopLeftTriangle{RegionKind::DiagonalHalf, 0}, kBottomRightTriangle{RegionKind::DiagonalHalf, 1}, kTopRightTriangle{RegionKind::DiagonalHalf, 2}, kBottomLeftTriangle{RegionKind::DiagonalHalf, 3};
inline constexpr BoxRegion kBottomSector{RegionKind::Sector, 0}, kRightSector{RegionKind::Sector, 1}, kTopSector{RegionKind::Sector, 2}, kLeftSector{RegionKind::Sector, 3};

// Keeps `source` and replaces the squares in each image of it with the copied squares. For example: apply_symmetry(sys, kBottomLeftQuadrant, {D4::Rot90, D4::Rot180, D4::Rot270}) builds C4 from one quadrant, apply_symmetry(sys, kLeftHalf, {D4::MirrorX}) mirrors the left half onto the right. The count stays exact: surplus images are dropped, most overlapped first, and spare squares go back into the destination. It only moves squares and does not minimise, so the result usually overlaps until it is relaxed. Frozen squares never move. False, with sys untouched, if it does not apply.
bool apply_symmetry(System& sys, BoxRegion source, const std::vector<D4>& images);

// The top-left triangle (y >= x) copied onto the bottom-right one through the main diagonal.
bool mirror_symmetry(System& sys);

// --- expand: candidates for other square counts ------------------------------------------------
//
// Each returns a starting point only, for another n; nothing here claims it is a packing. propagate() is what relaxes, compresses and validates them. All of them refuse a configuration with frozen or non-square particles by returning nothing.

// n + 4 when `centre` is false: one square per quadrant in the largest hole of a C4 packing.
// n + 1 when `centre` is true: one square at the box centre. Nothing unless C4 is detected.
std::optional<System> symmetric_grow(const System& sys, bool centre, double tol = 1e-6);

// direction 0 cuts across x, 1 across y. Duplicate shifts one side out by a pitch and fills the gap with a copy of the adjacent strip. Remove is the reverse, though the box only shrinks as far as the other dimension allows, so a tight packing keeps its side and leaves it to compress.
std::optional<System> row_duplicate(const System& sys, int direction, std::mt19937_64& rng);
std::optional<System> row_remove(const System& sys, int direction, std::mt19937_64& rng);

// An L-shaped layer along two walls: n + 2 floor(s) + 1 squares at side s + 1, for unit squares.
std::optional<System> frame_grow(const System& sys);

// n - 1: the square with the fewest contacts, the slackest of those on a tie.
std::optional<System> delete_rattler(const System& sys);

// n - 1, one candidate each for the `count` most-contacted squares: what compress can gain from losing a square is only known by compressing, and a load-bearing square is the likely win.
std::vector<System> delete_for_compress(const System& sys, int count);

struct PoolEntry {
    System packing;
    double side = 0.0;
};

// The best known packing per square count, in memory. Nothing enters without passing the same validation save_record() applies, at kCertifyTol, and then only as a strict improvement.
class Pool {
public:
    // Reads n<NNNN>.cfg from dir for n in [lo, hi]. Returns how many were admitted.
    int load_seeds(const std::string& dir, int lo, int hi);
    // Every *.log record in dir with n in [lo, hi]. Returns how many were admitted.
    int load_records(const std::string& dir, int lo, int hi);

    bool offer(const System& packing);
    const PoolEntry* at(int n) const;

    // Whether this configuration, up to D4 and relabelling, has been seen for its n before.
    bool first_seen(const System& candidate);

private:
    std::map<int, PoolEntry> entries_;
    std::map<int, std::set<CanonicalKey>> seen_;
};

struct PropagateOptions {
    int only_target = 0;      // when positive, only candidates for this n are evaluated
    int cuts = 2;             // row_duplicate and row_remove draws, per direction
    int delete_trials = 3;    // delete_for_compress candidates
    double shrink = 0.998;
    double grow = 1.01;
    long growth_stages = 10;
    long relaxations = 80;
};

struct PropagateReport {
    int tried = 0;
    std::map<int, double> improved;  // n -> decrease of its side
};

// Generates every candidate that applies from the pool's entry for n, then relaxes, compresses and validates each, updating the pool only on a strict improvement of L(m)/sqrt(m). Targets the pool has no entry for are skipped, and so are configurations already tried.
PropagateReport propagate(Pool& pool, int n, const PropagateOptions& opt, std::mt19937_64& rng, MoveLog* log = nullptr);

}  // namespace pack
