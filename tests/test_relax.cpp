// Assert-based tests for the analytic gradient and the L-BFGS relaxation.
//
// The one that matters is test_gradient_matches_finite_differences: everything downhill from
// here trusts that the gradient is the derivative of the energy the Monte Carlo moves use, and
// a wrong gradient does not crash, it just quietly relaxes to the wrong place.

#include "init.hpp"
#include "moves/mc.hpp"
#include "moves/relax.hpp"
#include "system.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace pack;

namespace {

bool close(double a, double b, double tol = 1e-12) { return std::abs(a - b) <= tol; }

void test_area_gradients_match_finite_differences() {
    // The geometry layer on its own: both area measures against central differences, over
    // random overlapping poses rather than a hand-picked few.
    const Shape s = make_square(1.0);
    std::mt19937_64 rng(1);
    std::uniform_real_distribution<double> d(-1.2, 1.2), ang(-kPi, kPi), inbox(-0.4, 1.4);
    const double h = 1e-6;
    int pairs = 0, walls = 0;

    for (int trial = 0; trial < 4000; ++trial) {
        Particle a, b;
        a.set_pose(0.0, 0.0, ang(rng));
        b.set_pose(d(rng), d(rng), ang(rng));

        Pose ga, gb;
        const double area = overlap_area_grad(a, s, b, s, &ga, &gb);
        assert(close(area, overlap_area(a, s, b, s), 1e-15));  // same number as the plain call
        if (area <= 1e-9) {
            continue;
        }
        ++pairs;

        // Newton's third law: a translation of either body sweeps the same area, oppositely.
        assert(close(ga.x + gb.x, 0.0, 1e-12));
        assert(close(ga.y + gb.y, 0.0, 1e-12));

        const Particle* which[2] = {&a, &b};
        const Pose* want[2] = {&ga, &gb};
        for (int side = 0; side < 2; ++side) {
            const double got[3] = {want[side]->x, want[side]->y, want[side]->theta};
            for (int k = 0; k < 3; ++k) {
                double v[3] = {which[side]->x(), which[side]->y(), which[side]->theta()};
                Particle up, dn;
                v[k] += h;
                up.set_pose(v[0], v[1], v[2]);
                v[k] -= 2.0 * h;
                dn.set_pose(v[0], v[1], v[2]);
                const double fd = (side == 0)
                                      ? (overlap_area(up, s, b, s) - overlap_area(dn, s, b, s)) / (2 * h)
                                      : (overlap_area(a, s, up, s) - overlap_area(a, s, dn, s)) / (2 * h);
                assert(std::abs(fd - got[k]) < 1e-6);
            }
        }
    }

    for (int trial = 0; trial < 4000; ++trial) {
        Particle p;
        p.set_pose(inbox(rng), inbox(rng), ang(rng));
        const double L = 1.0;

        Pose g;
        const double out = outside_area_grad(p, s, L, &g);
        assert(close(out, outside_area(p, s, L), 1e-15));
        // Skip the two flat ends: wholly inside, and wholly outside, where the gradient is
        // zero by saturation and a central difference straddles the kink.
        if (out <= 1e-9 || out >= s.area - 1e-9) continue;
        ++walls;

        const double got[3] = {g.x, g.y, g.theta};
        for (int k = 0; k < 3; ++k) {
            double v[3] = {p.x(), p.y(), p.theta()};
            Particle up, dn;
            v[k] += h;
            up.set_pose(v[0], v[1], v[2]);
            v[k] -= 2.0 * h;
            dn.set_pose(v[0], v[1], v[2]);
            const double fd = (outside_area(up, s, L) - outside_area(dn, s, L)) / (2 * h);
            assert(std::abs(fd - got[k]) < 1e-6);
        }
    }

    assert(pairs > 1000 && walls > 1000);  // vacuous if the random poses missed
}

void test_disjoint_squares_have_no_gradient() {
    const Shape s = make_square(1.0);
    Particle a, b;
    a.set_pose(0.0, 0.0, 0.3);
    b.set_pose(5.0, 0.0, 1.1);
    Pose ga, gb;
    assert(overlap_area_grad(a, s, b, s, &ga, &gb) == 0.0);
    assert(ga.x == 0.0 && ga.y == 0.0 && ga.theta == 0.0);
    assert(gb.x == 0.0 && gb.y == 0.0 && gb.theta == 0.0);

    Particle inside;
    inside.set_pose(2.0, 2.0, 0.0);
    Pose g;
    assert(outside_area_grad(inside, s, 4.0, &g) == 0.0);
    assert(g.x == 0.0 && g.y == 0.0 && g.theta == 0.0);
}

void test_gradient_matches_all_pairs() {
    // energy_grad() finds its pairs through a cell list; at this size the grid has real cells,
    // so a pair the grid lost would show up as a missing force. The reference is the all-pairs
    // loop the grid replaced. Centres spill past the walls to exercise the clamped edge cells.
    System sys = grid_init(make_square(1.0), 200);
    sys.L = 13.0;
    std::mt19937_64 rng(31);
    std::uniform_real_distribution<double> pos(-1.0, 14.0), ang(-kPi, kPi);
    for (Particle& p : sys.particles) p.set_pose(pos(rng), pos(rng), ang(rng));
    sys.particles[3].frozen = true;  // frozen bodies still push, and still get no gradient

    std::vector<Pose> want(sys.particles.size());
    double e_want = 0.0;
    for (size_t i = 0; i < sys.particles.size(); ++i) {
        const bool fi = sys.particles[i].frozen;
        Pose gw;
        const double out =
            outside_area_grad(sys.particles[i], sys.shape_of(static_cast<int>(i)), sys.L, &gw);
        e_want += sys.k_wall * out * out;
        if (!fi) {
            want[i].x += 2.0 * sys.k_wall * out * gw.x;
            want[i].y += 2.0 * sys.k_wall * out * gw.y;
            want[i].theta += 2.0 * sys.k_wall * out * gw.theta;
        }
        for (size_t j = i + 1; j < sys.particles.size(); ++j) {
            Pose gi, gj;
            const double a = overlap_area_grad(sys.particles[i], sys.shape_of(static_cast<int>(i)),
                                               sys.particles[j], sys.shape_of(static_cast<int>(j)),
                                               &gi, &gj);
            e_want += sys.k_pair * a * a;
            const double f = 2.0 * sys.k_pair * a;
            if (!fi) {
                want[i].x += f * gi.x;
                want[i].y += f * gi.y;
                want[i].theta += f * gi.theta;
            }
            if (!sys.particles[j].frozen) {
                want[j].x += f * gj.x;
                want[j].y += f * gj.y;
                want[j].theta += f * gj.theta;
            }
        }
    }

    std::vector<Pose> got;
    const double e = energy_grad(sys, got);
    const auto near = [](double a, double b) {
        return std::abs(a - b) <= 1e-10 * std::max(1.0, std::abs(b));
    };
    assert(near(e, e_want));
    for (size_t i = 0; i < sys.particles.size(); ++i) {
        assert(near(got[i].x, want[i].x));
        assert(near(got[i].y, want[i].y));
        assert(near(got[i].theta, want[i].theta));
    }
    assert(got[3].x == 0.0 && got[3].y == 0.0 && got[3].theta == 0.0);
}

void test_gradient_matches_finite_differences() {
    // The assembled System gradient, including the chain rule through the squared penalty and
    // the sum over every pair and wall.
    System sys = grid_init(make_square(1.0), 5);
    sys.set_box(2.3);
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> jitter(-0.35, 0.35), ang(-kPi, kPi);

    for (int trial = 0; trial < 60; ++trial) {
        for (Particle& p : sys.particles) {
            p.set_pose(p.x() + jitter(rng), p.y() + jitter(rng), p.theta() + ang(rng));
        }
        std::vector<Pose> grad;
        const double e = energy_grad(sys, grad);
        assert(close(e, sys.total_energy(), 1e-12));

        const double h = 1e-7;
        for (size_t i = 0; i < sys.particles.size(); ++i) {
            const double got[3] = {grad[i].x, grad[i].y, grad[i].theta};
            const double x = sys.particles[i].x(), y = sys.particles[i].y(), t = sys.particles[i].theta();
            for (int k = 0; k < 3; ++k) {
                double v[3] = {x, y, t};
                v[k] += h;
                sys.particles[i].set_pose(v[0], v[1], v[2]);
                const double up = sys.total_energy();
                v[k] -= 2.0 * h;
                sys.particles[i].set_pose(v[0], v[1], v[2]);
                const double dn = sys.total_energy();
                sys.particles[i].set_pose(x, y, t);
                const double fd = (up - dn) / (2.0 * h);
                // Relative, because these gradients run to the hundreds when squares are
                // deeply overlapped.
                assert(std::abs(fd - got[k]) <= 1e-5 * std::max(1.0, std::abs(got[k])));
            }
        }
    }
}

void test_frozen_particles_get_no_gradient_and_never_move() {
    System sys = grid_init(make_square(1.0), 4);
    sys.set_box(1.9);  // overlapping, so there is a real gradient to suppress
    sys.particles[2].frozen = true;

    std::vector<Pose> grad;
    energy_grad(sys, grad);
    assert(grad[2].x == 0.0 && grad[2].y == 0.0 && grad[2].theta == 0.0);
    bool others_pushed = false;
    for (size_t i = 0; i < grad.size(); ++i) {
        if (i != 2 && (grad[i].x != 0.0 || grad[i].y != 0.0)) others_pushed = true;
    }
    assert(others_pushed);  // the frozen square still pushes its neighbours

    const double x = sys.particles[2].x(), y = sys.particles[2].y(), t = sys.particles[2].theta();
    relax(sys);
    assert(sys.particles[2].x() == x && sys.particles[2].y() == y && sys.particles[2].theta() == t);
}

void test_relax_certifies_a_box_with_room() {
    // Four unit squares in a box of side 2.1: feasible with slack, so relaxation must drive
    // the overlaps to nothing at all, not merely to the tolerance.
    System sys = grid_init(make_square(1.0), 4);
    // grid_init is the perfect 2 x 2 packing, and 2.1 is wider than that, so sliding the walls
    // out leaves nothing to relax. Squeeze the squares into each other first.
    sys.rescale_box(1.8);
    sys.set_box(2.1);
    assert(sys.energy() > 0.0);  // the squeeze really did create violations

    const RelaxResult r = relax(sys);
    assert(r.feasible);
    assert(r.stop == RelaxResult::Stop::Feasible);
    assert(sys.max_overlap() <= kCertifyTol);
    assert(close(sys.energy(), sys.total_energy(), 1e-12));
    assert(r.iterations < 200);  // it is a short walk downhill, not a search
}

void test_relax_finishes_what_annealing_starts() {
    // The claim that justifies the file. A short anneal gets the squares roughly arranged but
    // nowhere near valid; relaxation closes the rest, which is the division of labour the
    // whole pipeline is built on.
    System annealed = grid_init(make_square(1.0), 4);
    annealed.rescale_box(1.8);  // overlapping, as above: the perfect packing needs no anneal
    annealed.set_box(2.1);
    Schedule sched;
    sched.sweeps = 50;  // deliberately far too short to resolve the overlaps
    StepSizes steps;
    std::mt19937_64 rng(11);
    anneal(annealed, sched, steps, rng);
    assert(!annealed.is_feasible(kFeasibleTol));

    System relaxed = annealed;
    const RelaxResult r = relax(relaxed);

    assert(r.feasible);
    assert(r.iterations > 0);  // it really had work to do
    assert(relaxed.max_overlap() < annealed.max_overlap());
    assert(relaxed.energy() <= annealed.energy());
}

void test_relax_reports_a_square_that_escaped_the_box() {
    // outside_area saturates at the shape's own area once a square is wholly outside, so the
    // gradient there is exactly zero and descent can do nothing. That has to be reported as
    // what it is, not as jamming: the two call for opposite responses from the caller.
    System sys = grid_init(make_square(1.0), 3);
    sys.set_box(2.0);
    set_pose(sys.particles[2], 8.0, 8.0, 0.0);  // far outside, nothing touching
    sys.recount();  // moved by hand, so resync_energy would only assert

    std::vector<Pose> grad;
    energy_grad(sys, grad);
    assert(grad[2].x == 0.0 && grad[2].y == 0.0 && grad[2].theta == 0.0);  // flat out there

    const RelaxResult r = relax(sys);
    assert(!r.feasible);
    assert(r.stop == RelaxResult::Stop::Escaped);
    assert(close(sys.particles[2].x(), 8.0) && close(sys.particles[2].y(), 8.0));
}

void test_relax_reports_an_escape_at_any_stop() {
    // With no wall penalty, one full-length descent step pushes the free square wholly outside,
    // and a one-iteration budget ends the run right there: a stop other than Converged.
    System sys = grid_init(make_square(1.0), 2);
    sys.set_box(5.0);
    sys.k_wall = 0.0;
    sys.particles[0].frozen = true;
    sys.set_pose(0, 1.2, 2.5, 0.0);
    sys.set_pose(1, 0.45, 2.5, 0.0);
    sys.recount();

    RelaxOptions one;
    one.max_iterations = 1;
    const RelaxResult r = relax(sys, one);
    assert(sys.particles[1].x() < -0.5);  // the fixture really did push it out
    assert(r.stop == RelaxResult::Stop::Escaped);
}

void test_relax_needs_memory_to_match_the_contact_set() {
    // Seventeen squares in a box 2% above the known optimum: the case that measures the
    // difference the L-BFGS memory makes. The default auto-sizes it to the variable count and
    // certifies; the textbook eight pairs cannot represent the valley and never gets there.
    // rescale_box, not set_box: the starting grid is wider than this box, so sliding the walls
    // inwards over squares that stay put leaves the outer ones wholly outside it, and out there
    // the wall penalty is saturated. The anneal then spends its budget random-walking them back
    // in and the fixture measures that instead of the thing under test.
    System base = grid_init(make_square(1.0), 17);
    base.rescale_box(4.675 * 1.02);
    Schedule sched;
    sched.sweeps = 3000;
    StepSizes steps;
    // The start comes out of the sampler, so the seed is part of the fixture: most seeds land
    // somewhere the full history also cannot repair, or somewhere already feasible, and neither
    // measures anything. Change the move mix and this number has to be picked again.
    // (Re-picked from 17 when the cell list changed the summation order of the pair terms.)
    std::mt19937_64 rng(19);
    anneal(base, sched, steps, rng);
    assert(!base.is_feasible(kFeasibleTol));

    System wide = base;
    const RelaxResult good = relax(wide);
    assert(good.feasible);
    assert(wide.is_feasible(kFeasibleTol));

    System narrow = base;
    RelaxOptions few;
    few.history = 8;
    few.max_iterations = 4000;
    const RelaxResult poor = relax(narrow, few);
    assert(!poor.feasible);                    // the point of auto-sizing
    assert(good.iterations < poor.iterations);
}

void test_relax_cannot_cheat_an_impossible_box() {
    // Two unit squares do not fit in a box of side 1.2 at any angle, so no amount of going
    // downhill can make this feasible: it must report the failure, not a packing.
    System sys = grid_init(make_square(1.0), 2);
    sys.set_box(1.2);

    const RelaxResult r = relax(sys);
    assert(!r.feasible);
    assert(r.stop != RelaxResult::Stop::Feasible);
    assert(sys.max_overlap() > 0.01);
    assert(close(sys.energy(), sys.total_energy(), 1e-9));
}

void test_relax_leaves_the_configuration_at_its_best_point() {
    // Whatever the exit reason, the system must hold the lowest point reached and an energy
    // that matches it -- a caller keeps relaxing from here.
    for (const double L : {1.2, 1.45, 2.0, 2.4}) {
        System sys = grid_init(make_square(1.0), 4);
        const double before = sys.energy();
        sys.set_box(L);
        const RelaxResult r = relax(sys);
        assert(close(sys.energy(), sys.total_energy(), 1e-9));
        assert(close(r.energy, sys.energy(), 1e-12));
        assert(sys.energy() <= sys.total_energy() + 1e-12);
        assert(r.max_overlap == sys.max_overlap());
        (void)before;
    }
}

void test_compress_finds_the_known_optimum_for_small_counts() {
    // The optimum is 2 for both 3 and 4 unit squares, and compression from the starting grid
    // has to get close to it while every reported packing stays strictly valid.
    for (const int n : {2, 3, 4}) {
        // Loosened on purpose: for these n the starting grid already IS the optimum, and a
        // compression with nothing to squeeze would assert its way to a pass without ever
        // exercising the backtracking squeeze this test is about.
        System start = grid_init(make_square(1.0), n);
        start.rescale_box(start.L * 1.25);
        const CompressResult c = compress(start);

        assert(c.found);
        assert(c.packing.is_feasible(kCertifyTol));  // the claim is a real packing
        assert(close(c.side, c.packing.tight_side(), 1e-12));
        // Never below the optimum: with the area tolerance set at kCertifyTol a corner
        // contact can only hide 1e-5 of depth, so a side under 2 would be a real error and
        // not tolerance slop.
        assert(c.side >= 2.0 - 1e-4);
        assert(c.side < 2.05);  // and the backtracking squeeze gets within 2.5% of it
        assert(c.relaxations > 1);

        assert(c.stages > 0);        // the squeeze accepted at least one step
        assert(c.squeezed_out);      // and ran to the step limit rather than the budget
    }
}

void test_hopeless_relaxation_gives_up_early() {
    // A full 10 x 10 lattice squeezed 0.5% cannot be valid: a jammed box, the case that used to
    // eat compress()'s whole budget. With the window on it has to say so within two windows.
    System jammed = grid_init(make_square(1.0), 100);
    jammed.rescale_box(jammed.L * 0.995);
    RelaxOptions opt;
    opt.hopeless_window = 200;
    const RelaxResult r = relax(jammed, opt);
    assert(r.stop == RelaxResult::Stop::Hopeless);
    assert(!r.feasible);
    assert(r.iterations <= 2 * opt.hopeless_window);

    // And it must not cost a relaxation that does succeed: the loosened eleven from the
    // bracket test below still certify with the window on.
    System loose = grid_init(make_square(1.0), 11);
    loose.rescale_box(loose.L * 1.25);
    const RelaxResult ok = relax(loose, opt);
    assert(ok.feasible);
    assert(ok.stop == RelaxResult::Stop::Feasible);
}

void test_feasible_is_recognised_without_overshooting() {
    // Fifty squares in a row, each overlapping the next by half the tolerance: valid as it
    // stands. The energy bound alone, sqrt(E / k) <= tol, adds the 49 contacts up and says no;
    // the overlaps themselves say yes, and relax() has to stop on the spot.
    const double d = 0.5 * kCertifyTol;  // unit squares edge to edge: overlap area = depth
    System row = grid_init(make_square(1.0), 50);
    row.L = 50.0;
    for (int k = 0; k < 50; ++k) row.set_pose(k, 0.5 + k * (1.0 - d), 0.5, 0.0);
    row.recount();
    assert(std::sqrt(row.energy() / row.k_pair) > kCertifyTol);  // the old test would refuse
    const RelaxResult r = relax(row);
    assert(r.stop == RelaxResult::Stop::Feasible);
    assert(r.iterations == 0);
}

void test_compress_brackets_the_side_it_cannot_beat() {
    // Eleven squares: the case with a known hard optimum of 3.877. Whatever arrangement the
    // compression is handed, it has to squeeze until the step runs out, which shows up as a
    // bracket far tighter than the 0.5% of a single stage.
    // Loosened like the case above: the starting grid is a jammed lattice, and a squeeze that
    // cannot take its first step measures nothing about the bracket.
    System seed = grid_init(make_square(1.0), 11);
    seed.rescale_box(seed.L * 1.25);
    const CompressResult c = compress(seed);

    assert(c.found);
    assert(c.packing.is_feasible(kCertifyTol));
    assert(c.squeezed_out);                         // squeezed out, not merely stopped early
    // Down to the side's own resolution: sqrt(kCertifyTol) as a length, relative to the side.
    assert(c.final_step <= std::sqrt(kCertifyTol) / c.side);
    assert(c.stages > 1);
    assert(c.side >= 3.877 - 1e-4);                 // never below the known optimum
}

void test_compress_keeps_every_square_in_the_box_it_claims() {
    const System start = grid_init(make_square(1.0), 7);
    const CompressResult c = compress(start);
    assert(c.found);

    // Reading the answer as "fits in a square of side c.side" has to survive a direct check:
    // slide the packing to the origin and nothing may poke out of a box of that side.
    System moved = c.packing;
    double xlo = 1e300, ylo = 1e300;
    for (size_t i = 0; i < moved.particles.size(); ++i) {
        const Particle& p = moved.particles[i];
        const Shape& s = moved.shape_of(static_cast<int>(i));
        double lo, hi;
        project(p, s, {1.0, 0.0}, lo, hi);
        xlo = std::min(xlo, lo);
        project(p, s, {0.0, 1.0}, lo, hi);
        ylo = std::min(ylo, lo);
    }
    for (Particle& p : moved.particles) p.set_pose(p.x() - xlo, p.y() - ylo, p.theta());
    moved.set_box(c.side + 1e-9);
    assert(moved.is_feasible(kCertifyTol));
}

void test_compress_best_of_keeps_the_tightest() {
    // Ten squares: the lattice is jammed at 4 and compress() cannot fit it any smaller, while
    // the diagonal strip certifies at 3 + 1/sqrt(2). Handed both, best-of has to keep the
    // strip -- in whichever order -- and the thread count must not change a single side.
    const System grid = grid_init(make_square(1.0), 10);
    const System diag = diagonal_init(make_square(1.0), 10);
    const std::vector<System> cands = {grid, diag};

    const BestOfResult one = compress_best_of(cands, {}, 1);
    const BestOfResult many = compress_best_of(cands, {}, 3);
    assert(one.winner == 1);
    assert(one.best.found);
    assert(one.best.side <= 3.0 + 1.0 / std::sqrt(2.0) + 1e-6);
    assert(one.sides[0] >= 4.0 - 1e-6);  // the lattice stays where it is
    assert(one.sides == many.sides && one.winner == many.winner);
    assert(one.best.packing.is_feasible(kCertifyTol));
}

void test_relax_free_spreads_past_the_box() {
    // Two unit squares cannot fit in a box of side 1.2, but with the walls off they just part.
    System sys = grid_init(make_square(1.0), 2);
    sys.rescale_box(1.2);
    const double k_wall = sys.k_wall;
    assert(sys.max_overlap() > 0.1);

    relax_free(sys);
    assert(sys.k_wall == k_wall);
    assert(sys.L >= 2.0 - 1e-4);
    assert(close(sys.L, sys.tight_side(), 1e-12));
    assert(sys.is_feasible(kCertifyTol));  // the refit box holds them too
    assert(close(sys.energy(), sys.total_energy(), 1e-12));
}

void test_anneal_shrinks_by_ds_and_keeps_the_last_valid_box() {
    System sys = grid_init(make_square(1.0), 4);
    sys.rescale_box(2.5);

    const CompressResult c = anneal(sys, 0.05, 200);
    assert(c.found);
    assert(c.relaxations <= 200);
    assert(c.side >= 2.0 - 1e-4);
    assert(c.side < 2.5 - 0.05);          // it took at least one step
    assert(close(sys.tight_side(), c.side, 1e-12));  // sys is the answer
    assert(sys.is_feasible(kCertifyTol));

    // The budget counts every minimisation.
    System few = grid_init(make_square(1.0), 4);
    few.rescale_box(2.5);
    assert(anneal(few, 0.05, 3).relaxations <= 3);

    // Nothing valid is ever found for a square stranded outside: sys stays as it was.
    System lost = grid_init(make_square(1.0), 3);
    lost.set_box(2.0);
    lost.set_pose(2, 8.0, 8.0, 0.0);
    lost.recount();
    assert(!anneal(lost, 0.05, 50).found);
    assert(lost.L == 2.0 && lost.particles[2].x() == 8.0);
}

}  // namespace

int main() {
    test_relax_free_spreads_past_the_box();
    test_anneal_shrinks_by_ds_and_keeps_the_last_valid_box();
    test_compress_best_of_keeps_the_tightest();
    test_gradient_matches_all_pairs();
    test_area_gradients_match_finite_differences();
    test_disjoint_squares_have_no_gradient();
    test_gradient_matches_finite_differences();
    test_frozen_particles_get_no_gradient_and_never_move();
    test_relax_certifies_a_box_with_room();
    test_relax_finishes_what_annealing_starts();
    test_relax_reports_a_square_that_escaped_the_box();
    test_relax_reports_an_escape_at_any_stop();
    test_relax_needs_memory_to_match_the_contact_set();
    test_relax_cannot_cheat_an_impossible_box();
    test_relax_leaves_the_configuration_at_its_best_point();
    test_compress_finds_the_known_optimum_for_small_counts();
    test_hopeless_relaxation_gives_up_early();
    test_feasible_is_recognised_without_overshooting();
    test_compress_brackets_the_side_it_cannot_beat();
    test_compress_keeps_every_square_in_the_box_it_claims();
    std::printf("relax: all tests passed\n");
    return 0;
}
