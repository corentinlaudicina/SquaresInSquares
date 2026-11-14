// legalize(): strict output from near-feasible input, and no cost when the input is already legal.

#include "legalize.hpp"
#include "moves/relax.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>

using namespace pack;

namespace {

System grid(int k) {
    System s;
    s.species.push_back(make_square(1.0));
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j) s.add_particle(i + 0.5, j + 0.5, 0.0, 0);
    s.L = k;
    s.recount();
    return s;
}

void test_touching_grid_is_kept_exact() {
    System s = grid(4);
    assert(strict_report(s).ok());
    const LegalizeResult r = legalize(s);
    assert(r.ok && r.growth_steps == 0 && r.relax_iterations == 0);
    assert(r.side_after >= 4.0 && r.side_after - 4.0 < 1e-13);
    std::printf("  grid 4x4: side %.17g, no relaxation\n", r.side_after);
}

void test_squeezed_grid_is_repaired() {
    // A 5x5 grid squeezed by 1e-6 and relaxed to kCertifyTol, as compress() leaves things: certified,
    // but with real overlaps.
    System s = grid(5);
    s.rescale_box(5.0 * (1.0 - 1e-6));
    for (int i = 0; i < 25; ++i) {
        const Particle& p = s.particles[static_cast<size_t>(i)];
        s.set_pose(i, p.x(), p.y(), 0.02 * std::sin(i));  // tilt them so contacts are corners
    }
    s.recount();
    relax(s, RelaxOptions{});
    assert(!strict_report(s).ok());
    const double before = s.tight_side();
    const LegalizeResult r = legalize(s);
    assert(r.ok);
    assert(strict_report(s).ok());
    assert(s.L - before < 1e-4);  // repair costs about the overlap depth, not a dilation
    std::printf("  squeezed 5x5: %.12f -> %.12f, %d growths, %ld iterations, min gap %.2e\n", before, s.L,
                r.growth_steps, r.relax_iterations, r.report.min_pair_gap);
}

void test_tilted_overlap_is_repaired() {
    // s(5) shape with the tilted centre square 1e-6 too big for its hole.
    System s;
    s.species.push_back(make_square(1.0));
    const double L = 2.0 + 1.0 / std::sqrt(2.0) - 1e-6;
    for (const double x : {0.5, L - 0.5})
        for (const double y : {0.5, L - 0.5}) s.add_particle(x, y, 0.0, 0);
    s.add_particle(L / 2, L / 2, kPi / 4, 0);
    s.L = L;
    s.recount();
    assert(!strict_report(s).ok());
    const LegalizeResult r = legalize(s);
    assert(r.ok);
    assert(s.L > 2.0 + 1.0 / std::sqrt(2.0) - 1e-12);  // cannot beat the optimum
    std::printf("  s(5)-1e-6: legal side %.12f (optimum %.12f)\n", s.L, 2.0 + 1.0 / std::sqrt(2.0));
}

}  // namespace

int main() {
    test_touching_grid_is_kept_exact();
    test_squeezed_grid_is_repaired();
    test_tilted_overlap_is_repaired();
    std::printf("legalize: all tests passed\n");
}
