// Assert-based tests for system.cpp and moves/mc.cpp.
// The one that matters is test_incremental_energy_matches: the annealer trusts a running
// energy updated by deltas, and nothing else would notice if those deltas were wrong.

#include "init.hpp"
#include "moves/mc.hpp"
#include "moves/regions.hpp"
#include "system.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace pack;

namespace {

bool close(double a, double b, double tol = 1e-12) { return std::abs(a - b) <= tol; }

// A 3 x 3 lattice of unit squares at unit pitch in a box of side 3, so every selector's answer
// below can be written out by hand instead of computed:
//
//     6 7 8      y = 2.5
//     3 4 5      y = 1.5
//     0 1 2      y = 0.5
//     x = 0.5, 1.5, 2.5
System lattice_nine() {
    System s = grid_init(make_square(1.0), 9);
    s.L = 3.0;
    for (int k = 0; k < 9; ++k) {
        s.set_pose(k, 0.5 + (k % 3), 0.5 + (k / 3), 0.0);
    }
    s.recount();
    return s;
}

bool same(const std::vector<int>& got, const std::vector<int>& want) { return got == want; }

void test_grid_init() {
    for (int n = 1; n <= 12; ++n) {
        const System s = grid_init(make_square(1.0), n);
        assert(static_cast<int>(s.particles.size()) == n);
        // Trivially feasible: nothing overlaps and nothing pokes out of the box.
        assert(s.max_overlap() == 0.0);
        assert(s.is_feasible(0.0));
        assert(close(s.energy(), 0.0));
        assert(close(s.energy(), s.total_energy()));
    }
}

void test_energy_parts_price_any_stiffness() {
    // The energy is linear in each stiffness separately, which is what lets System::energy()
    // weigh two carried sums instead of recounting when a stiffness changes.
    System s = grid_init(make_square(1.0), 6);
    s.rescale_box(2.3);  // squeeze, so there is both pair overlap and wall overlap to measure

    const EnergyParts p = s.energy_parts();
    assert(p.pair > 0.0 && p.wall > 0.0);  // vacuous otherwise
    assert(close(energy_of(p, s.k_pair, s.k_wall), s.total_energy(), 1e-12));

    for (const double f : {0.5, 1.0, 3.0, 17.0}) {
        System t = s;
        t.k_pair = s.k_pair * f;
        t.k_wall = s.k_wall * f;
        assert(close(t.total_energy(), f * s.total_energy(), 1e-9));
        assert(close(energy_of(p, t.k_pair, t.k_wall), t.total_energy(), 1e-9));
    }

    // And the two components are independent: switching one stiffness off leaves exactly the
    // other term, which is what makes them separately recoverable.
    System no_wall = s;
    no_wall.k_wall = 0.0;
    assert(close(no_wall.total_energy(), s.k_pair * p.pair, 1e-12));
    System no_pair = s;
    no_pair.k_pair = 0.0;
    assert(close(no_pair.total_energy(), s.k_wall * p.wall, 1e-12));
}

void test_tight_side() {
    // Two unit squares side by side, touching: 2 wide, 1 tall, so the enclosing square is 2.
    System s = grid_init(make_square(1.0), 2);
    set_pose(s.particles[0], 0.5, 0.5, 0.0);
    set_pose(s.particles[1], 1.5, 0.5, 0.0);
    assert(close(s.tight_side(), 2.0));

    // One square turned 45 degrees reaches sqrt(2) across.
    System one = grid_init(make_square(1.0), 1);
    set_pose(one.particles[0], 0.0, 0.0, kPi / 4.0);
    assert(close(one.tight_side(), std::sqrt(2.0), 1e-15));
}

void test_particle_energy_is_the_whole_delta() {
    // particle_energy(i) must capture every term that changes when i alone moves, otherwise
    // the incremental bookkeeping silently diverges from total_energy().
    System s = grid_init(make_square(1.0), 4);
    s.set_box(2.4);
    std::mt19937_64 rng(5);
    std::uniform_real_distribution<double> d(-0.6, 0.6), ang(-kPi, kPi);

    for (int trial = 0; trial < 400; ++trial) {
        const int i = static_cast<int>(rng() % s.particles.size());
        const double total_before = s.total_energy();
        const double part_before = s.particle_energy(i);

        const Particle& p = s.particles[static_cast<size_t>(i)];
        s.set_pose(i, p.x() + d(rng), p.y() + d(rng), p.theta() + ang(rng));

        const double by_part = s.particle_energy(i) - part_before;
        const double by_total = s.total_energy() - total_before;
        assert(close(by_part, by_total, 1e-9));
    }
}

void test_incremental_energy_matches() {
    // Run real moves and check the running energy against a from-scratch recount.
    System s = grid_init(make_square(1.0), 5);
    s.set_box(2.6);
    StepSizes steps;
    std::mt19937_64 rng(11);

    for (int step = 0; step < 4000; ++step) {
        const int i = static_cast<int>(rng() % s.particles.size());
        if (rng() & 1u) {
            try_rotate(s, i, steps, 0.5, rng);
        } else {
            try_translate(s, i, steps, 0.5, rng);
        }
        if (step % 500 == 0) {
            assert(close(s.energy(), s.total_energy(), 1e-9));
        }
    }
    assert(close(s.energy(), s.total_energy(), 1e-9));
}

void test_rejected_move_restores_the_pose() {
    System s = grid_init(make_square(1.0), 2);
    s.set_box(1.9);  // too small for two unit squares side by side, so moves get rejected
    std::mt19937_64 rng(3);
    StepSizes steps;

    const double x0 = s.particles[0].x(), y0 = s.particles[0].y(), t0 = s.particles[0].theta();
    const double e0 = s.energy();
    int rejected = 0;
    for (int k = 0; k < 200; ++k) {
        const double x = s.particles[0].x(), y = s.particles[0].y(), t = s.particles[0].theta();
        if (!try_translate(s, 0, steps, 0.0, rng)) {  // T = 0: uphill is always rejected
            ++rejected;
            assert(close(s.particles[0].x(), x) && close(s.particles[0].y(), y) &&
                   close(s.particles[0].theta(), t));
        }
    }
    assert(rejected > 0);  // the test is vacuous if everything was accepted
    // At T = 0 the energy can only fall.
    assert(s.energy() <= e0 + 1e-12);
    (void)x0; (void)y0; (void)t0;
}

void test_frozen_particle_never_moves() {
    System s = grid_init(make_square(1.0), 3);
    s.set_box(2.2);
    s.particles[1].frozen = true;
    const double x = s.particles[1].x(), y = s.particles[1].y(), t = s.particles[1].theta();

    Schedule sched;
    sched.sweeps = 200;
    StepSizes steps;
    std::mt19937_64 rng(7);
    anneal(s, sched, steps, rng);

    assert(close(s.particles[1].x(), x) && close(s.particles[1].y(), y) &&
           close(s.particles[1].theta(), t));
}

void test_fit_box() {
    // A grid dragged far from the box: same arrangement afterwards, snug in a box of its own.
    System s = grid_init(make_square(1.0), 4);
    const double side = s.tight_side();
    const double gap_x = s.particles[1].x() - s.particles[0].x();
    for (int i = 0; i < 4; ++i) {
        const Particle& p = s.particles[static_cast<size_t>(i)];
        s.set_pose(i, p.x() - 7.5, p.y() + 3.25, p.theta());
    }
    s.recount();
    assert(s.max_overlap() > 0.5);  // it really is outside the box

    s.fit_box();
    assert(close(s.L, side, 1e-12));
    assert(close(s.tight_side(), side, 1e-12));
    assert(close(s.particles[1].x() - s.particles[0].x(), gap_x, 1e-12));
    assert(s.is_feasible(1e-12));
    assert(close(s.energy(), s.total_energy(), 1e-12));

    // The bounding box counts corners, not centres: a lone square at 45 degrees.
    System one = grid_init(make_square(1.0), 1);
    one.set_pose(0, -4.0, 9.0, kPi / 4.0);
    one.fit_box();
    assert(close(one.L, std::sqrt(2.0), 1e-12));
    assert(close(one.particles[0].x(), std::sqrt(0.5), 1e-12));
    assert(close(one.particles[0].y(), std::sqrt(0.5), 1e-12));
}

void test_anneal_watch_sees_every_sweep_and_changes_nothing() {
    System plain = grid_init(make_square(1.0), 5);
    plain.set_box(2.6);
    System watched = plain;
    Schedule sched;
    sched.sweeps = 300;
    StepSizes a, b;
    std::mt19937_64 ra(3), rb(3);

    long last = -1;
    anneal(plain, sched, a, ra);
    anneal(watched, sched, b, rb, [&](long k, const System& s) {
        assert(k == last + 1);
        last = k;
        assert(close(s.energy(), s.total_energy(), 1e-9));
    });

    assert(last == sched.sweeps - 1);
    for (size_t i = 0; i < plain.particles.size(); ++i) {
        assert(plain.particles[i].x() == watched.particles[i].x());
        assert(plain.particles[i].y() == watched.particles[i].y());
        assert(plain.particles[i].theta() == watched.particles[i].theta());
    }
}

void test_group_moves() {
    // A rigid group move must keep the members' arrangement, and the running energy must
    // still match a from-scratch recount -- that is the claim that lets group_parts double
    // count intra-group pairs.
    System s = grid_init(make_square(1.0), 5);
    s.set_box(2.8);
    const std::vector<int> group{0, 2, 4};
    StepSizes steps;
    std::mt19937_64 rng(42);

    int moved = 0;
    for (int trial = 0; trial < 300; ++trial) {
        // Record every pairwise gap inside the group before the move.
        std::vector<double> gaps;
        for (size_t a = 0; a < group.size(); ++a) {
            for (size_t b = a + 1; b < group.size(); ++b) {
                const Particle& pa = s.particles[static_cast<size_t>(group[a])];
                const Particle& pb = s.particles[static_cast<size_t>(group[b])];
                gaps.push_back(std::hypot(pa.x() - pb.x(), pa.y() - pb.y()));
            }
        }

        const bool ok = (trial % 2) ? try_rotate_group(s, group, steps, 0.5, rng)
                                    : try_translate_group(s, group, steps, 0.5, rng);
        if (ok) ++moved;

        size_t g = 0;
        for (size_t a = 0; a < group.size(); ++a) {
            for (size_t b = a + 1; b < group.size(); ++b, ++g) {
                const Particle& pa = s.particles[static_cast<size_t>(group[a])];
                const Particle& pb = s.particles[static_cast<size_t>(group[b])];
                assert(close(std::hypot(pa.x() - pb.x(), pa.y() - pb.y()), gaps[g], 1e-9));
            }
        }
        assert(close(s.energy(), s.total_energy(), 1e-9));
    }
    assert(moved > 0);  // vacuous if nothing was ever accepted

    // A frozen member blocks the whole group.
    s.particles[2].frozen = true;
    const double x = s.particles[0].x();
    for (int k = 0; k < 50; ++k) {
        assert(!try_translate_group(s, group, steps, 1.0, rng));
        assert(!try_rotate_group(s, group, steps, 1.0, rng));
    }
    assert(close(s.particles[0].x(), x));
    s.particles[2].frozen = false;

    assert(!try_translate_group(s, {}, steps, 1.0, rng));        // empty group
    assert(!try_translate_group(s, {0, 99}, steps, 1.0, rng));   // out-of-range index
}

void test_contact_pitch() {
    // The pitch at which two copies at the same orientation touch: the unit every move that
    // re-lays a row works in, measured off the shape rather than assumed to be 1.
    assert(close(contact_pitch(make_square(1.0)), 1.0));
    assert(close(contact_pitch(make_square(2.5)), 2.5));
}

void test_region_selectors() {
    System s = lattice_nine();

    // A disc tight enough to hold the middle square alone: the nearest other centre is 1 away.
    assert(same(select_disc(s, 1.5, 1.5, 0.6), {4}));
    // Out to just past 1: the four edge neighbours join it, the diagonal ones at 1.414 do not.
    assert(same(select_disc(s, 1.5, 1.5, 1.01), {1, 3, 4, 5, 7}));
    assert(select_disc(s, 1.5, 1.5, 1.5).size() == 9);  // past the diagonals: everything
    assert(select_disc(s, 1.5, 1.5, 0.0).empty());      // a region with no area holds nothing
    assert(select_disc(s, 1.5, 1.5, -1.0).empty());

    // A band is a slab in the coordinate normal to its direction, and the normal is dir turned a
    // quarter turn counter-clockwise -- so for a band along x the normal is +y, and for a band
    // along y it is -x and the offsets come out negative. That sign is the one thing about this
    // interface worth pinning down.
    assert(same(select_band(s, {1.0, 0.0}, 1.0, 1.0), {3, 4, 5}));    // middle row
    assert(same(select_band(s, {1.0, 0.0}, 0.0, 1.0), {0, 1, 2}));    // bottom row
    assert(same(select_band(s, {0.0, 1.0}, -2.0, 1.0), {1, 4, 7}));   // middle column
    assert(select_band(s, {1.0, 0.0}, 1.2, 0.2).empty());             // between two rows

    // The diagonal a tilted band would be carved along: dir = (1,1), normal = (-1,1)/sqrt(2), so
    // the normal coordinate is (y - x)/sqrt(2) and the main diagonal sits at exactly 0.
    assert(same(select_band(s, {1.0, 1.0}, -0.3, 0.6), {0, 4, 8}));
    // The other diagonal is not the mirror image of that call: dir = (1,-1) turns the normal to
    // (1,1)/sqrt(2), so the coordinate is (x + y)/sqrt(2) and the anti-diagonal sits at
    // 3/sqrt(2) rather than at 0. Worth writing down instead of leaving it to be rediscovered.
    const double anti = 3.0 / std::sqrt(2.0);
    assert(same(select_band(s, {1.0, -1.0}, anti - 0.3, 0.6), {2, 4, 6}));
    assert(select_band(s, {0.0, 0.0}, 0.0, 1.0).empty());             // no direction, no band

    // Corners run counter-clockwise from the origin: 0 = (0,0), 1 = (L,0), 2 = (L,L), 3 = (0,L).
    assert(same(select_corner(s, 0, 1.1), {0}));
    assert(same(select_corner(s, 0, 2.1), {0, 1, 3, 4}));
    assert(same(select_corner(s, 1, 1.1), {2}));
    assert(same(select_corner(s, 2, 2.1), {4, 5, 7, 8}));
    assert(same(select_corner(s, 3, 1.1), {6}));
    assert(select_corner(s, 4, 1.1).empty());   // not a corner
    assert(select_corner(s, -1, 1.1).empty());

    // A region is a place in the box, not a chosen set of particles, so a pinned square in the
    // middle of one is worked around rather than being a reason to abandon the move. This is the
    // opposite of what the group moves do with a frozen member, and the difference is deliberate.
    s.particles[4].frozen = true;
    assert(same(select_disc(s, 1.5, 1.5, 1.01), {1, 3, 5, 7}));
    assert(same(select_band(s, {1.0, 0.0}, 1.0, 1.0), {3, 5}));
    assert(same(select_corner(s, 0, 2.1), {0, 1, 3}));
}

void test_pick_region_stays_in_range() {
    System s = lattice_nine();
    s.particles[2].frozen = true;
    RegionWeights w;
    std::mt19937_64 rng(1234);

    int discs = 0, bands = 0, corners = 0, nonempty = 0;
    for (int k = 0; k < 4000; ++k) {
        const Pick p = pick_region(s, w, rng);
        switch (p.kind) {
            case Region::Disc: ++discs; break;
            case Region::Band: ++bands; break;
            case Region::Corner: ++corners; break;
        }
        if (!p.members.empty()) ++nonempty;
        assert(p.extent > 0.0);
        assert(close(norm(p.axis), 1.0, 1e-12));  // the symmetry line is a unit direction
        assert(close(norm(p.dir), 1.0, 1e-12));
        int last = -1;
        for (const int i : p.members) {
            assert(i > last);  // ascending and distinct: a member listed twice would be moved twice
            last = i;
            assert(i >= 0 && i < static_cast<int>(s.particles.size()));
            assert(!s.particles[static_cast<size_t>(i)].frozen);
        }
    }
    // Bands and corners carry most of the weight, for the reason at the top of regions.hpp.
    assert(bands > discs && corners > discs);
    assert(nonempty > 0);
}

void test_a_region_move_keeps_the_region() {
    // The claim the shove move rests on: the rigid motion offered for each region maps that
    // region onto itself, so the same draw afterwards selects the same squares, the reverse move
    // is as likely as the forward one, and the chain keeps detailed balance.
    System s = lattice_nine();
    StepSizes steps;
    steps.translate = 0.4;
    std::mt19937_64 rng(13);

    const Vec2 dir{1.0, 0.0};
    const double offset = 1.0, width = 1.0;
    const std::vector<int> band = select_band(s, dir, offset, width);
    assert(same(band, {3, 4, 5}));

    // A slide runs along the band, and a band is unbounded along its own direction.
    int moved = 0;
    for (int k = 0; k < 200; ++k) {
        if (try_slide_group(s, band, dir, steps, 0.8, rng)) ++moved;
        assert(same(select_band(s, dir, offset, width), band));
    }
    assert(moved > 0);

    // A mirror in the band's own midline mirrors the normal coordinate about the middle of the
    // slab, so it too stays inside it. The midline, not the group's centroid: the centroid is
    // not generally the middle of the slab, and mirroring about it would push members out.
    const Vec2 midline{1.5, offset + 0.5 * width};
    assert(try_reflect_group(s, band, midline, dir, 1e9, rng));
    assert(same(select_band(s, dir, offset, width), band));

    // The negative control, and the reason mc.cpp offers no sideways shove of a band: move the
    // members across the band instead and the same draw no longer finds a single one of them.
    for (const int i : band) {
        const Particle& p = s.particles[static_cast<size_t>(i)];
        s.set_pose(i, p.x(), p.y() + 1.0, p.theta());
    }
    assert(select_band(s, dir, offset, width).empty());
}

void test_shove_region_keeps_the_energy_exact() {
    System s = lattice_nine();
    StepSizes steps;
    steps.translate = 0.4;
    steps.rotate = 0.3;
    std::mt19937_64 rng(77);

    int moved = 0;
    for (int k = 0; k < 600; ++k) {
        if (try_shove_region(s, RegionWeights{}, steps, 0.5, rng)) ++moved;
        // A rigid move of a group leaves every intra-group pair alone, which is what lets the
        // delta be measured on group_parts with those pairs double counted.
        assert(close(s.energy(), s.total_energy(), 1e-9));
        assert(s.particles.size() == 9);
    }
    assert(moved > 0);  // vacuous if nothing was ever accepted
}

void test_reflect_group_is_its_own_inverse() {
    // A reflection is reversible because it is an involution, and that is what the balance
    // argument in mc.hpp rests on. Applied twice at T = 0 with the same line, it has to come
    // back -- and it has to preserve every distance inside the group on the way.
    System s = lattice_nine();
    const std::vector<int> group{0, 1, 3, 4};
    const Vec2 point{1.0, 1.0}, axis{1.0, 1.0};
    std::mt19937_64 rng(5);

    const double d01 = std::hypot(s.particles[0].x() - s.particles[1].x(), s.particles[0].y() - s.particles[1].y());
    assert(try_reflect_group(s, group, point, axis, 1e9, rng));  // hot: always accepted
    assert(close(std::hypot(s.particles[0].x() - s.particles[1].x(), s.particles[0].y() - s.particles[1].y()),
                 d01, 1e-12));
    // Mirrored in the main diagonal through (1,1), square 1 at (1.5,0.5) goes to (0.5,1.5).
    assert(close(s.particles[1].x(), 0.5, 1e-12) && close(s.particles[1].y(), 1.5, 1e-12));

    assert(try_reflect_group(s, group, point, axis, 1e9, rng));
    assert(close(s.particles[1].x(), 1.5, 1e-12) && close(s.particles[1].y(), 0.5, 1e-12));
    assert(close(s.energy(), s.total_energy(), 1e-9));

    // The group rules still apply: a frozen member stops the whole group, unlike a region.
    StepSizes steps;
    s.particles[3].frozen = true;
    assert(!try_reflect_group(s, group, point, axis, 1e9, rng));
    assert(!try_slide_group(s, group, {1.0, 0.0}, steps, 1e9, rng));
    assert(!try_reflect_group(s, group, point, {0.0, 0.0}, 1e9, rng));  // no line, no reflection
}

void test_slide_group_moves_along_one_direction_only() {
    System s = lattice_nine();
    const std::vector<int> group{0, 1, 2};
    StepSizes steps;
    steps.translate = 0.3;
    std::mt19937_64 rng(9);

    int moved = 0;
    for (int k = 0; k < 200; ++k) {
        const double y0 = s.particles[0].y(), y1 = s.particles[1].y();
        if (try_slide_group(s, group, {1.0, 0.0}, steps, 0.6, rng)) ++moved;
        // Along x only: a slide cannot carry a member out of the band it was drawn from.
        assert(close(s.particles[0].y(), y0) && close(s.particles[1].y(), y1));
        assert(close(s.energy(), s.total_energy(), 1e-9));
    }
    assert(moved > 0);
    assert(!try_slide_group(s, group, {0.0, 0.0}, steps, 0.6, rng));  // no direction, no move
}

void test_flip_turns_by_45_degrees() {
    // A square is unchanged by a quarter turn, so the flip has to be 45 degrees: a 90 degree
    // version would be the identity, accepted every time, and would do nothing whatever.
    System s = lattice_nine();
    std::mt19937_64 rng(3);

    int moved = 0;
    for (int k = 0; k < 400; ++k) {
        const double t0 = s.particles[4].theta();
        if (try_flip(s, 4, 1e9, rng)) {  // hot: always accepted
            ++moved;
            assert(close(std::abs(s.particles[4].theta() - t0), 0.25 * kPi, 1e-12));
        }
        assert(close(s.energy(), s.total_energy(), 1e-9));
    }
    assert(moved > 0);

    s.particles[4].frozen = true;
    assert(!try_flip(s, 4, 1e9, rng));
}

void test_member_bounds_stays_in_the_box() {
    System s = lattice_nine();
    const Bounds b = member_bounds(s, {0, 1, 3, 4}, 5.0);  // padded far past the box
    assert(b.xlo == 0.0 && b.ylo == 0.0);
    assert(close(b.xhi, 3.0) && close(b.yhi, 3.0));
    const Bounds tight = member_bounds(s, {4}, 0.0);
    assert(close(tight.xlo, 1.5) && close(tight.xhi, 1.5));
    const Bounds none = member_bounds(s, {}, 1.0);
    assert(none.xlo == 0.0 && none.xhi == 0.0);
}

void test_anneal_resolves_three_squares() {
    // Three unit squares fit in a 2 x 2 box with room to spare (the known optimum is 2), so a
    // short anneal must get the overlaps away entirely.
    System s = grid_init(make_square(1.0), 3);
    // grid_init already IS the 2 x 2 answer, so squeeze the squares together first and give
    // the box back: the anneal then has the overlaps to resolve that this test is about.
    s.rescale_box(1.6);
    s.set_box(2.0);
    const double squeezed = s.energy();
    assert(squeezed > 1.0);  // the squeeze really did create violations

    Schedule sched;
    sched.sweeps = 4000;
    StepSizes steps;
    std::mt19937_64 rng(2024);
    const Stats st = anneal(s, sched, steps, rng);

    assert(st.attempted == 3 * sched.sweeps);
    assert(st.accepted > 0);
    assert(s.energy() < squeezed);
    assert(s.is_feasible(kFeasibleTol));
    assert(s.tight_side() <= 2.0 + kFeasibleTol);
}

void test_anneal_cannot_cheat_an_impossible_box() {
    // Two unit squares cannot fit in a box of side 1.2, so the energy must stay well above
    // zero no matter how long we anneal. This is the test that would catch an energy that
    // forgets a term.
    System s = grid_init(make_square(1.0), 2);
    s.set_box(1.2);

    Schedule sched;
    sched.sweeps = 3000;
    StepSizes steps;
    std::mt19937_64 rng(99);
    anneal(s, sched, steps, rng);

    assert(!s.is_feasible(kFeasibleTol));
    assert(s.max_overlap() > 0.01);
}


void test_incremental_parts_match() {
    // The mirror of test_incremental_energy_matches, and deliberately wider: it drives box
    // moves and group moves as well, which that test never did. A rejected box move that
    // restored the energy but not the parts would leave sums that drift further apart with every
    // rejection.
    System s = grid_init(make_square(1.0), 5);
    s.set_box(2.6);
    StepSizes steps;
    std::mt19937_64 rng(11);
    const std::vector<int> group{0, 1, 2};
    long boxes = 0;

    for (int step = 0; step < 4000; ++step) {
        const int i = static_cast<int>(rng() % s.particles.size());
        switch (rng() % 4) {
            case 0: try_rotate(s, i, steps, 0.5, rng); break;
            case 1: try_translate(s, i, steps, 0.5, rng); break;
            case 2: try_translate_group(s, group, steps, 0.5, rng); break;
            default:
                try_box(s, steps, 3.0, 0.5, rng);
                ++boxes;
                break;
        }
        if (step % 500 == 0) {
            const EnergyParts exact = s.energy_parts();
            assert(close(s.parts.pair, exact.pair, 1e-9));
            assert(close(s.parts.wall, exact.wall, 1e-9));
            assert(close(s.energy(), s.total_energy(), 1e-9));
        }
    }
    assert(boxes > 100);  // vacuous if the box arm never fired
    const EnergyParts exact = s.energy_parts();
    assert(close(s.parts.pair, exact.pair, 1e-9));
    assert(close(s.parts.wall, exact.wall, 1e-9));
}

// The all-pairs sums the cell list replaced, kept here as the reference it is checked against.
EnergyParts brute_parts(const System& s) {
    EnergyParts e;
    for (size_t i = 0; i < s.particles.size(); ++i) {
        const double w = outside_area(s.particles[i], s.shape_of(static_cast<int>(i)), s.L);
        e.wall += w * w;
        for (size_t j = i + 1; j < s.particles.size(); ++j) {
            const double a = overlap_area(s.particles[i], s.shape_of(static_cast<int>(i)), s.particles[j],
                                          s.shape_of(static_cast<int>(j)));
            e.pair += a * a;
        }
    }
    return e;
}

double brute_max_overlap(const System& s) {
    double worst = 0.0;
    for (size_t i = 0; i < s.particles.size(); ++i) {
        worst = std::max(worst, outside_area(s.particles[i], s.shape_of(static_cast<int>(i)), s.L));
        for (size_t j = i + 1; j < s.particles.size(); ++j) {
            worst = std::max(worst, overlap_area(s.particles[i], s.shape_of(static_cast<int>(i)),
                                                 s.particles[j], s.shape_of(static_cast<int>(j))));
        }
    }
    return worst;
}

double brute_particle_energy(const System& s, int i) {
    const Particle& pi = s.particles[static_cast<size_t>(i)];
    const Shape& si = s.shape_of(i);
    const double w = outside_area(pi, si, s.L);
    double e = s.k_wall * w * w;
    for (size_t j = 0; j < s.particles.size(); ++j) {
        if (static_cast<int>(j) == i) continue;
        const double a = overlap_area(pi, si, s.particles[j], s.shape_of(static_cast<int>(j)));
        e += s.k_pair * a * a;
    }
    return e;
}

// Relative, because these sums run from zero to thousands depending on the configuration.
bool close_rel(double a, double b, double tol = 1e-10) {
    return std::abs(a - b) <= tol * std::max(1.0, std::max(std::abs(a), std::abs(b)));
}

void assert_matches_all_pairs(const System& s) {
    const EnergyParts got = s.energy_parts(), want = brute_parts(s);
    assert(close_rel(got.pair, want.pair));
    assert(close_rel(got.wall, want.wall));
    assert(s.max_overlap() == brute_max_overlap(s));  // a max does not care about order
    for (int i = 0; i < static_cast<int>(s.particles.size()); ++i) {
        assert(close_rel(s.particle_energy(i), brute_particle_energy(s, i)));
    }
}

void test_cell_list_matches_all_pairs() {
    // Random configurations at every grid regime: one cell (L below the reach), two per side
    // (where an unclipped stencil would visit cells twice), a real grid, and a box huge
    // compared with its contents (where the cell count is capped). Centres are thrown up to
    // two units outside the box, so the clamped edge cells are exercised too.
    struct Case { int n; double L; };
    const Case cases[] = {{3, 1.0}, {8, 3.0}, {40, 6.0}, {300, 16.0}, {20, 200.0}};
    std::mt19937_64 rng(11);
    std::uniform_real_distribution<double> ang(-kPi, kPi);
    for (const Case c : cases) {
        System s = grid_init(make_square(1.0), c.n);
        s.L = c.L;
        std::uniform_real_distribution<double> pos(-2.0, c.L + 2.0);
        for (Particle& p : s.particles) set_pose(p, pos(rng), pos(rng), ang(rng));
        s.recount();
        assert_matches_all_pairs(s);
    }
}

void test_cell_list_follows_the_moves() {
    // Every move type at a temperature hot enough that most of them land, so particles cross
    // cell boundaries constantly, box moves rebuild the grid and rejected box moves restore the
    // poses into a grid built for the other L. After each sweep the grid must still find every
    // pair -- the per-particle energies against the all-pairs loop, and the running energy
    // against a recount, with the grid's own consistency asserted by resync_energy().
    System s = grid_init(make_square(1.0), 150);
    s.rescale_box(s.L * 0.97);
    StepSizes steps;
    steps.translate = 0.5;
    std::mt19937_64 rng(23);
    for (int k = 0; k < 40; ++k) {
        sweep(s, steps, 5.0, rng);
        try_box(s, steps, 1.0, 5.0, rng);
        const double brute = energy_of(brute_parts(s), s.k_pair, s.k_wall);
        assert(close_rel(s.energy(), brute, 1e-8));
        s.resync_energy();
        for (int i = 0; i < static_cast<int>(s.particles.size()); ++i) {
            assert(close_rel(s.particle_energy(i), brute_particle_energy(s, i)));
        }
    }
}

void test_diagonal_init() {
    // The family's side is k + 1/sqrt(2) wherever the strip is full: 10 and 27 are the known
    // optima, and 973 is the size where it first beats the 32-wide lattice for a thousand-ish.
    const double r = 1.0 / std::sqrt(2.0);
    const struct { int n; double side; } known[] = {{10, 3.0 + r}, {27, 5.0 + r}, {973, 31.0 + r}};
    for (const auto& k : known) {
        const System d = diagonal_init(make_square(1.0), k.n);
        assert(static_cast<int>(d.particles.size()) == k.n);
        assert(close(d.L, k.side, 1e-12));
        assert(d.tight_side() <= d.L + 1e-12);
    }
    // Overlap-free everywhere it is built, beating the lattice or not, and its tracked energy
    // is the configuration's. The tolerance is the certification one: this is a packing.
    for (int n = 1; n <= 160; ++n) {
        const System d = diagonal_init(make_square(1.0), n);
        assert(static_cast<int>(d.particles.size()) == n);
        assert(d.max_overlap() <= 1e-12);
        assert(close(d.energy(), d.total_energy(), 1e-12));
    }
    // A side-2 square scales with it; anything but a square gets the lattice.
    assert(close(diagonal_init(make_square(2.0), 10).L, 2.0 * (3.0 + r), 1e-12));
    const Shape hex = make_regular_polygon(6, 0.6);
    assert(diagonal_init(hex, 10).L == grid_init(hex, 10).L);
}

}  // namespace

int main() {
    test_fit_box();
    test_anneal_watch_sees_every_sweep_and_changes_nothing();
    test_diagonal_init();
    test_grid_init();
    test_energy_parts_price_any_stiffness();
    test_tight_side();
    test_particle_energy_is_the_whole_delta();
    test_incremental_energy_matches();
    test_incremental_parts_match();
    test_rejected_move_restores_the_pose();
    test_frozen_particle_never_moves();
    test_group_moves();
    test_contact_pitch();
    test_region_selectors();
    test_pick_region_stays_in_range();
    test_a_region_move_keeps_the_region();
    test_shove_region_keeps_the_energy_exact();
    test_reflect_group_is_its_own_inverse();
    test_slide_group_moves_along_one_direction_only();
    test_flip_turns_by_45_degrees();
    test_member_bounds_stays_in_the_box();
    test_anneal_resolves_three_squares();
    test_anneal_cannot_cheat_an_impossible_box();
    test_cell_list_matches_all_pairs();
    test_cell_list_follows_the_moves();
    std::printf("system + mc: all tests passed\n");
    return 0;
}
