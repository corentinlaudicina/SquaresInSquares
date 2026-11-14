// Assert-based tests for the large-neighbourhood moves in moves/lnm.hpp.
//
// Two matter most. test_every_move_keeps_the_count_and_the_energy, because these moves lift
// squares out of the configuration and re-lay whole bands at once: losing a square or leaving
// the running energy stale are both live possibilities here, and either would still produce
// packings -- of the wrong count, or measured at the wrong energy. And
// test_reject_leaves_no_trace, because all four moves share one copy-and-restore reject path, so
// pinning it once pins it four times.
//
// What is deliberately not tested is that these moves sample anything. They do not: relax() is
// many-to-one, and the note at the top of lnm.hpp says what that costs. The claim under test is
// the one they are actually for -- that they reach packings a gradient cannot -- and that is
// test_large_moves_escape_a_jam.

#include "init.hpp"
#include "moves/lnm.hpp"
#include "moves/regions.hpp"
#include "system.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace pack;

namespace {

bool close(double a, double b, double tol = 1e-12) { return std::abs(a - b) <= tol; }

struct Snap {
    double x, y, theta;
};

std::vector<Snap> poses(const System& s) {
    std::vector<Snap> out;
    out.reserve(s.particles.size());
    for (const Particle& p : s.particles) out.push_back({p.x(), p.y(), p.theta()});
    return out;
}

bool unchanged(const System& s, const std::vector<Snap>& was) {
    if (s.particles.size() != was.size()) return false;
    for (size_t i = 0; i < was.size(); ++i) {
        // Bit-identical, not close: the reject path restores a saved copy, it does not try to
        // undo the arithmetic it did, and a test that allowed a tolerance here would pass for an
        // implementation that drifted a little every time it said no.
        if (s.particles[i].x() != was[i].x || s.particles[i].y() != was[i].y ||
            s.particles[i].theta() != was[i].theta) {
            return false;
        }
    }
    return true;
}

bool any_wholly_outside(const System& s) {
    for (size_t i = 0; i < s.particles.size(); ++i) {
        const Shape& sh = s.shape_of(static_cast<int>(i));
        if (outside_area(s.particles[i], sh, s.L) >= sh.area * (1.0 - 1e-12)) return true;
    }
    return false;
}

// Eleven squares in a box 1% above the known optimum of 3.877, relaxed straight out of the
// squeezed grid. The gradient jams there at a worst overlap of about 2.6e-3 and cannot do
// better from that start, which is exactly the state these moves exist for. No sampler is
// involved, so the fixture does not move when the sweep mix changes.
System jammed_eleven() {
    System s = grid_init(make_square(1.0), 11);
    s.rescale_box(3.877 * 1.01);
    relax(s);
    return s;
}

const double kPressure = 2.0;

using Move = bool (*)(System&, const LnmOptions&, double, double, std::mt19937_64&);
constexpr int kKinds = 4;
const Move kMoves[kKinds] = {try_destroy_repair, try_strip, try_reshape_corner, try_loosen};
const char* const kNames[kKinds] = {"destroy", "strip", "corner", "loosen"};

void test_every_move_keeps_the_count_and_the_energy() {
    for (int k = 0; k < kKinds; ++k) {
        System s = jammed_eleven();
        const size_t n = s.particles.size();
        std::mt19937_64 rng(1000 + k);

        int accepted = 0;
        for (int trial = 0; trial < 60; ++trial) {
            if (kMoves[k](s, LnmOptions{}, kPressure, 0.3, rng)) ++accepted;
            assert(s.particles.size() == n);          // no square lost, none conjured
            assert(close(s.energy(), s.total_energy(), 1e-9));
            assert(s.L > 0.0);
            // propose() refuses an Escaped relaxation, so no accepted configuration can be
            // holding a square out where the wall penalty is flat and nothing can pull it back.
            assert(!any_wholly_outside(s));
        }
        // Vacuous if the move never fired. Printed rather than merely asserted, because a move
        // that is accepted twice in sixty attempts is technically passing and still broken.
        std::printf("  %-8s accepted %2d / 60\n", kNames[k], accepted);
        assert(accepted > 0);
    }
}

void test_reject_leaves_no_trace() {
    // At T = 0 nothing uphill is accepted, and from a relaxed configuration most of these
    // proposals are uphill -- so there are plenty of rejections to inspect.
    for (int k = 0; k < kKinds; ++k) {
        System s = jammed_eleven();
        std::mt19937_64 rng(3000 + k);

        int rejected = 0;
        for (int trial = 0; trial < 40; ++trial) {
            const std::vector<Snap> was = poses(s);
            const double L0 = s.L, e0 = s.energy();
            if (!kMoves[k](s, LnmOptions{}, kPressure, 0.0, rng)) {
                ++rejected;
                assert(unchanged(s, was));
                assert(s.L == L0);
                assert(s.energy() == e0);
                assert(s.k_pair == 100.0 && s.k_wall == 100.0);  // loosen softens both, briefly
            }
        }
        assert(rejected > 0);  // vacuous if everything was accepted
    }
}

void test_loosen_puts_the_box_and_the_stiffness_back() {
    // The move grows the box, and half its attempts also soften both stiffnesses for the shake.
    // Every one of those has to be undone whether the move is accepted or not, and a stiffness
    // left scaled down would quietly make every later energy in the run wrong.
    System s = jammed_eleven();
    const double L0 = s.L, kp = s.k_pair, kw = s.k_wall;
    std::mt19937_64 rng(4242);

    for (int trial = 0; trial < 60; ++trial) {
        try_loosen(s, LnmOptions{}, kPressure, 0.3, rng);
        assert(close(s.L, L0, 1e-12));  // back to the box it started in, accepted or not
        assert(s.k_pair == kp && s.k_wall == kw);
        assert(close(s.energy(), s.total_energy(), 1e-9));
    }
}

void test_frozen_squares_never_move() {
    // A frozen square is skipped by the region selectors and left alone by rescale_box, sweep()
    // and relax(). Four separate reasons it stays put, so it is worth checking that all four
    // hold at once.
    System s = jammed_eleven();
    s.particles[0].frozen = true;
    s.particles[5].frozen = true;
    const double x0 = s.particles[0].x(), y0 = s.particles[0].y(), t0 = s.particles[0].theta();
    const double x5 = s.particles[5].x(), y5 = s.particles[5].y(), t5 = s.particles[5].theta();
    std::mt19937_64 rng(555);

    for (int trial = 0; trial < 40; ++trial) {
        for (int k = 0; k < kKinds; ++k) kMoves[k](s, LnmOptions{}, kPressure, 0.4, rng);
        assert(s.particles[0].x() == x0 && s.particles[0].y() == y0 && s.particles[0].theta() == t0);
        assert(s.particles[5].x() == x5 && s.particles[5].y() == y5 && s.particles[5].theta() == t5);
        assert(close(s.energy(), s.total_energy(), 1e-9));
    }
}

void test_zero_temperature_never_goes_uphill() {
    System s = jammed_eleven();
    std::mt19937_64 rng(606);

    for (int trial = 0; trial < 40; ++trial) {
        for (int k = 0; k < kKinds; ++k) {
            const double h0 = enthalpy(s, kPressure);
            kMoves[k](s, LnmOptions{}, kPressure, 0.0, rng);
            // The enthalpy is the quantity the move decides on, so at T = 0 it is monotone. Not
            // the energy: loosen goes through states at other box sides.
            assert(enthalpy(s, kPressure) <= h0 + 1e-12);
        }
    }
}

void test_moves_are_deterministic() {
    // Same seed, same answer.
    System a = jammed_eleven(), b = jammed_eleven();
    std::mt19937_64 ra(31337), rb(31337);
    LnmStats sa, sb;

    for (int trial = 0; trial < 30; ++trial) {
        lnm_round(a, LnmOptions{}, sa, kPressure, 0.3, ra);
        lnm_round(b, LnmOptions{}, sb, kPressure, 0.3, rb);
    }
    assert(a.L == b.L && a.energy() == b.energy());
    assert(unchanged(a, poses(b)));
    const Stats ta = lnm_total(sa), tb = lnm_total(sb);
    assert(ta.attempted == tb.attempted && ta.accepted == tb.accepted);
}

void test_lnm_round_draws_by_weight() {
    System s = jammed_eleven();
    std::mt19937_64 rng(808);

    // One weight at a time: only that bucket may be touched.
    for (int k = 0; k < kKinds; ++k) {
        LnmOptions opt;
        opt.w_destroy = opt.w_strip = opt.w_corner = opt.w_loosen = 0.0;
        if (k == 0) opt.w_destroy = 1.0;
        if (k == 1) opt.w_strip = 1.0;
        if (k == 2) opt.w_corner = 1.0;
        if (k == 3) opt.w_loosen = 1.0;

        LnmStats st;
        for (int trial = 0; trial < 12; ++trial) lnm_round(s, opt, st, kPressure, 0.3, rng);
        const Stats* each[kKinds] = {&st.destroy, &st.strip, &st.corner, &st.loosen};
        for (int j = 0; j < kKinds; ++j) {
            assert(each[j]->attempted == (j == k ? 12 : 0));
        }
        assert(lnm_total(st).attempted == 12);
    }

    // Every weight zero is "no move", not a crash and not a silent fallback to one of them.
    LnmOptions none;
    none.w_destroy = none.w_strip = none.w_corner = none.w_loosen = 0.0;
    LnmStats st;
    for (int trial = 0; trial < 5; ++trial) {
        assert(!lnm_round(s, none, st, kPressure, 0.3, rng));
    }
    assert(lnm_total(st).attempted == 0);
}

void test_large_moves_escape_a_jam() {
    // The claim the whole file rests on. The gradient jams this configuration: all eleven
    // overlaps would have to shrink together and no downhill direction does that. Relaxing it
    // again buys a few iterations and nothing else -- a stalled relaxation is not quite a true
    // minimum, so a fresh L-BFGS history does shave the last digits, which is why the test
    // measures the overlap against a relative bound and not an exact one. A few dozen large
    // moves, by contrast, find a valid packing in the same box.
    System base = jammed_eleven();
    assert(!base.is_feasible(kFeasibleTol));
    const double jammed_overlap = base.max_overlap();
    System again = base;
    for (int k = 0; k < 5; ++k) {
        const RelaxResult r = relax(again);
        assert(!r.feasible);
    }
    assert(again.max_overlap() >= 0.99 * jammed_overlap);  // downhill has nothing left to give

    System s = base;
    LnmOptions opt;
    LnmStats st;
    std::mt19937_64 rng(1);
    int attempts = 0;
    for (; attempts < 250 && !s.is_feasible(kFeasibleTol); ++attempts) {
        lnm_round(s, opt, st, 0.0, 0.05, rng);
    }
    std::printf("  escaped the jam after %d large moves (worst overlap %.3g -> %.3g)\n", attempts,
                jammed_overlap, s.max_overlap());
    assert(s.is_feasible(kFeasibleTol));
    assert(s.particles.size() == 11);
    // The tight side can sit a little outside the box it was found in, and legitimately so:
    // kFeasibleTol is an area, and at a corner contact an area of 1e-5 is a depth of about
    // sqrt(1e-5), once at each end. That is the sampler's tolerance, not a packing's -- the
    // certified side comes from compress() at kCertifyTol, five orders of magnitude tighter.
    assert(s.tight_side() <= 3.877 * 1.01 + 2.0 * std::sqrt(kFeasibleTol));
}

void test_large_neighbour_move_always_changes_the_configuration() {
    for (unsigned seed = 1; seed <= 10; ++seed) {
        System s = jammed_eleven();
        const std::vector<Snap> was = poses(s);
        std::mt19937_64 rng(seed);
        assert(large_neighbour_move(s, rng));
        assert(!unchanged(s, was));
        assert(s.particles.size() == 11);
        assert(close(s.energy(), s.total_energy(), 1e-9));
    }
}

}  // namespace

int main() {
    test_large_neighbour_move_always_changes_the_configuration();
    test_every_move_keeps_the_count_and_the_energy();
    test_reject_leaves_no_trace();
    test_loosen_puts_the_box_and_the_stiffness_back();
    test_frozen_squares_never_move();
    test_zero_temperature_never_goes_uphill();
    test_moves_are_deterministic();
    test_lnm_round_draws_by_weight();
    test_large_moves_escape_a_jam();
    std::printf("lnm: all tests passed\n");
    return 0;
}
