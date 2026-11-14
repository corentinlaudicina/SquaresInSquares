// Assert-based tests for src/merge/rectangle.hpp and src/merge/merge.hpp.
//
// Two things matter here: detect_rectangle_union has to say yes to an exact rectangle and no to
// everything else (a gap, a stray tilt, a non-rectangular touching cluster), and merge()/
// dissolve() have to be exact inverses -- a group's energy and every member's pose must survive
// the round trip, because that is the one property the rest of the codebase would notice if it
// broke.

#include "init.hpp"
#include "merge/merge.hpp"
#include "merge/rectangle.hpp"
#include "system.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace pack;

namespace {

bool close(double a, double b, double tol = 1e-9) { return std::abs(a - b) <= tol; }

// A System of unit squares at exactly the poses given, nothing else -- the detector and
// merge()/dissolve() tests build their own layouts by hand rather than relying on grid_init's
// particular placement order.
System squares_at(const std::vector<std::array<double, 3>>& poses) {
    System sys;
    sys.add_species(make_square(1.0));
    for (const auto& p : poses) sys.add_particle(p[0], p[1], p[2], 0);
    sys.recount();
    return sys;
}

std::vector<int> all_indices(const System& sys) {
    std::vector<int> v(sys.particles.size());
    for (size_t i = 0; i < v.size(); ++i) v[i] = static_cast<int>(i);
    return v;
}

void test_detects_a_simple_block() {
    const System sys = squares_at({{0.5, 0.5, 0.0}, {1.5, 0.5, 0.0}, {0.5, 1.5, 0.0}, {1.5, 1.5, 0.0}});
    const auto ru = detect_rectangle_union(sys, all_indices(sys));
    assert(ru.has_value());
    assert(close(ru->side, 1.0));
    assert((ru->cols == 2 && ru->rows == 2));
    assert(close(ru->theta, 0.0));
    assert(close(ru->center.x, 1.0) && close(ru->center.y, 1.0));
    assert(ru->ordered_members.size() == 4);
}

void test_detects_a_rotated_row() {
    // A 1x4 row along local x, then the whole thing rigidly rotated by an arbitrary angle --
    // not a multiple of 90 degrees -- about the first member.
    const double phi = 0.37;
    const double c = std::cos(phi), s = std::sin(phi);
    std::vector<std::array<double, 3>> poses;
    for (int i = 0; i < 4; ++i) {
        poses.push_back({static_cast<double>(i) * c, static_cast<double>(i) * s, phi});
    }
    const System sys = squares_at(poses);
    const auto ru = detect_rectangle_union(sys, all_indices(sys));
    assert(ru.has_value());
    assert(close(ru->side, 1.0));
    assert((ru->cols == 4 && ru->rows == 1) || (ru->cols == 1 && ru->rows == 4));
    assert(close(ru->theta, phi));
}

void test_detects_mod_90_equivalent_orientation() {
    // Same 2x2 block, but one member's own theta is 90 degrees from the others -- a square
    // looks identical at a quarter turn, so this has to pass just as the all-0 block does.
    System sys = squares_at(
        {{0.5, 0.5, 0.0}, {1.5, 0.5, 0.0}, {0.5, 1.5, kPi / 2.0}, {1.5, 1.5, 0.0}});
    const auto ru = detect_rectangle_union(sys, all_indices(sys));
    assert(ru.has_value());
    assert((ru->cols == 2 && ru->rows == 2));
}

void test_rejects_a_gap() {
    // Same block, but one member shifted 1.5 sides out: a gap, not a rectangle.
    const System sys =
        squares_at({{0.5, 0.5, 0.0}, {1.5, 0.5, 0.0}, {0.5, 1.5, 0.0}, {3.0, 1.5, 0.0}});
    assert(!detect_rectangle_union(sys, all_indices(sys)).has_value());
}

void test_rejects_a_tilted_member() {
    const System sys = squares_at(
        {{0.5, 0.5, 0.0}, {1.5, 0.5, 0.17}, {0.5, 1.5, 0.0}, {1.5, 1.5, 0.0}});
    assert(!detect_rectangle_union(sys, all_indices(sys)).has_value());
}

void test_rejects_a_non_rectangular_cluster() {
    // A plus/cross of 5 touching squares: the right count, the wrong cell set.
    const System sys = squares_at({{1.5, 0.5, 0.0},
                                   {0.5, 1.5, 0.0},
                                   {1.5, 1.5, 0.0},
                                   {2.5, 1.5, 0.0},
                                   {1.5, 2.5, 0.0}});
    assert(!detect_rectangle_union(sys, all_indices(sys)).has_value());
}

void test_rejects_a_duplicate_slot() {
    const System sys =
        squares_at({{0.5, 0.5, 0.0}, {0.5, 0.5, 0.0}, {0.5, 1.5, 0.0}, {1.5, 1.5, 0.0}});
    assert(!detect_rectangle_union(sys, all_indices(sys)).has_value());
}

void test_rejects_mixed_species() {
    System sys = squares_at({{0.5, 0.5, 0.0}, {1.5, 0.5, 0.0}, {0.5, 1.5, 0.0}, {1.5, 1.5, 0.0}});
    const int other = sys.add_species(make_square(1.0));
    sys.particles[3].species = other;
    sys.recount();
    assert(!detect_rectangle_union(sys, all_indices(sys)).has_value());
}

void test_rejects_too_few_members() {
    const System sys = squares_at({{0.5, 0.5, 0.0}});
    assert(!detect_rectangle_union(sys, all_indices(sys)).has_value());
    assert(!detect_rectangle_union(sys, {}).has_value());
}

void test_merge_and_dissolve_round_trip() {
    const System original = grid_init(make_square(1.0), 4);  // a 2x2 block, energy 0
    const double e_before = original.total_energy();

    System sys = original;
    const std::optional<int> gi = merge(sys, all_indices(sys));
    assert(gi.has_value());
    assert(sys.particles.size() == 1);
    assert(sys.is_merged(*gi));
    assert(close(sys.species[static_cast<size_t>(sys.particles[static_cast<size_t>(*gi)].species)].area,
                4.0));
    assert(close(sys.total_energy(), e_before, 1e-9));

    const bool ok = dissolve(sys, *gi);
    assert(ok);
    assert(sys.particles.size() == 4);
    assert(!sys.is_merged(0) && !sys.is_merged(1) && !sys.is_merged(2) && !sys.is_merged(3));
    assert(close(sys.total_energy(), e_before, 1e-9));

    // Every original pose comes back, in the group's own member order -- not necessarily the
    // caller's original index order, so match by position rather than by index.
    for (const Particle& orig : original.particles) {
        bool found = false;
        for (const Particle& p : sys.particles) {
            if (close(p.x(), orig.x()) && close(p.y(), orig.y()) && close(p.theta(), orig.theta())) {
                found = true;
                break;
            }
        }
        assert(found);
    }
}

void test_unrelated_particle_keeps_its_pose_through_an_unrelated_merge() {
    // An 8-square system: merge/dissolve a 2x2 block (indices 0-3) and check that particle 4's
    // *pose* survives, checked by value since its index does shift.
    System sys = squares_at({{0.5, 0.5, 0.0},
                             {1.5, 0.5, 0.0},
                             {0.5, 1.5, 0.0},
                             {1.5, 1.5, 0.0},
                             {10.0, 10.0, 0.37},  // unrelated, far away
                             {20.0, 20.0, 0.0},
                             {30.0, 30.0, 0.0},
                             {40.0, 40.0, 0.0}});
    const Vec2 bystander{10.0, 10.0};
    const double bystander_theta = 0.37;

    const std::optional<int> gi = merge(sys, {0, 1, 2, 3});
    assert(gi.has_value());
    bool found = false;
    for (const Particle& p : sys.particles) {
        if (close(p.x(), bystander.x) && close(p.y(), bystander.y) &&
            close(p.theta(), bystander_theta)) {
            found = true;
        }
    }
    assert(found);

    assert(dissolve(sys, *gi));
    found = false;
    for (const Particle& p : sys.particles) {
        if (close(p.x(), bystander.x) && close(p.y(), bystander.y) &&
            close(p.theta(), bystander_theta)) {
            found = true;
        }
    }
    assert(found);
}

void test_merge_rejects_leave_the_system_untouched() {
    // A non-rectangular cluster: merge() must refuse and change nothing.
    System sys = squares_at(
        {{1.5, 0.5, 0.0}, {0.5, 1.5, 0.0}, {1.5, 1.5, 0.0}, {2.5, 1.5, 0.0}, {1.5, 2.5, 0.0}});
    const System before = sys;
    assert(!merge(sys, all_indices(sys)).has_value());
    assert(sys.particles.size() == before.particles.size());
    assert(sys.species.size() == before.species.size());
    for (size_t i = 0; i < sys.particles.size(); ++i) {
        assert(close(sys.particles[i].x(), before.particles[i].x()));
        assert(close(sys.particles[i].y(), before.particles[i].y()));
    }
}

void test_merge_rejects_an_already_merged_member() {
    System sys = squares_at({{0.5, 0.5, 0.0}, {1.5, 0.5, 0.0}, {0.5, 1.5, 0.0}, {1.5, 1.5, 0.0},
                             {3.5, 0.5, 0.0}, {4.5, 0.5, 0.0}, {3.5, 1.5, 0.0}, {4.5, 1.5, 0.0}});
    const std::optional<int> first = merge(sys, {0, 1, 2, 3});
    assert(first.has_value());
    assert(sys.particles.size() == 5);  // 4 removed, 1 appended
    // The second block's squares shifted down into indices 0-3 when the first four were
    // removed. Trying to fold the already-merged group into a second merge, alongside them,
    // has to fail outright rather than merging a merged group.
    const size_t before = sys.particles.size();
    assert(!merge(sys, {*first, 0, 1}).has_value());
    assert(sys.particles.size() == before);  // the rejected attempt changed nothing
}

void test_merge_rejects_mixed_frozen_and_accepts_all_frozen() {
    System mixed = squares_at({{0.5, 0.5, 0.0}, {1.5, 0.5, 0.0}, {0.5, 1.5, 0.0}, {1.5, 1.5, 0.0}});
    mixed.particles[0].frozen = true;
    assert(!merge(mixed, all_indices(mixed)).has_value());

    System frozen = squares_at({{0.5, 0.5, 0.0}, {1.5, 0.5, 0.0}, {0.5, 1.5, 0.0}, {1.5, 1.5, 0.0}});
    for (Particle& p : frozen.particles) p.frozen = true;
    const std::optional<int> gi = merge(frozen, all_indices(frozen));
    assert(gi.has_value());
    assert(frozen.particles[static_cast<size_t>(*gi)].frozen);
}

void test_dissolve_rejects_a_non_merged_index() {
    System sys = squares_at({{0.5, 0.5, 0.0}, {1.5, 0.5, 0.0}});
    assert(!dissolve(sys, 0));
    assert(!dissolve(sys, -1));
    assert(!dissolve(sys, 99));
    assert(sys.particles.size() == 2);  // untouched
}

}  // namespace

int main() {
    test_detects_a_simple_block();
    test_detects_a_rotated_row();
    test_detects_mod_90_equivalent_orientation();
    test_rejects_a_gap();
    test_rejects_a_tilted_member();
    test_rejects_a_non_rectangular_cluster();
    test_rejects_a_duplicate_slot();
    test_rejects_mixed_species();
    test_rejects_too_few_members();
    test_merge_and_dissolve_round_trip();
    test_unrelated_particle_keeps_its_pose_through_an_unrelated_merge();
    test_merge_rejects_leave_the_system_untouched();
    test_merge_rejects_an_already_merged_member();
    test_merge_rejects_mixed_frozen_and_accepts_all_frozen();
    test_dissolve_rejects_a_non_merged_index();
    std::printf("merge: all tests passed\n");
    return 0;
}
