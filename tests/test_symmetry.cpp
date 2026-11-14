// Assert-based tests for moves/symmetry.hpp.
//
// The group tests come first because everything else trusts them: a wrong element here would
// still produce packings, only never the symmetric ones it claims to.

#include "init.hpp"
#include "io/io.hpp"
#include "moves/lnm.hpp"
#include "moves/mc.hpp"
#include "moves/relax.hpp"
#include "moves/symmetry.hpp"
#include "system.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace pack;

namespace {

const double kDeg = kPi / 180.0;

System config_from(const std::vector<Pose>& poses, double L) {
    System s;
    s.species.push_back(make_square(1.0));
    s.L = L;
    for (const Pose& p : poses) s.add_particle(p.x, p.y, p.theta, 0);
    s.recount();
    return s;
}

System random_config(int n, double L, unsigned seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> pos(0.0, L), ang(-kPi, kPi);
    std::vector<Pose> poses;
    for (int i = 0; i < n; ++i) poses.push_back({pos(rng), pos(rng), ang(rng)});
    return config_from(poses, L);
}

// Same squares as a set, in any order, angles mod 90 degrees.
bool same_set(const System& a, const System& b, double tol) {
    if (a.particles.size() != b.particles.size()) return false;
    std::vector<char> used(b.particles.size(), 0);
    for (const Particle& p : a.particles) {
        bool found = false;
        for (size_t j = 0; j < b.particles.size() && !found; ++j) {
            const Particle& q = b.particles[j];
            if (used[j] || std::hypot(p.x() - q.x(), p.y() - q.y()) > tol) continue;
            if (angle_distance90(p.theta(), q.theta()) > tol) continue;
            used[j] = 1;
            found = true;
        }
        if (!found) return false;
    }
    return true;
}

// `seed` squares and all their images under `group`.
System symmetrised(const std::vector<Pose>& seed, const std::vector<GroupElem>& group, double L) {
    std::vector<Pose> poses;
    for (const Pose& p : seed) {
        for (const GroupElem& g : group) poses.push_back(g.apply(p, L));
    }
    return config_from(poses, L);
}

void test_group_table() {
    const auto& all = d4_elements();
    for (size_t i = 0; i < all.size(); ++i) {
        for (size_t j = i + 1; j < all.size(); ++j) assert(!(all[i] == all[j]));
        assert(all[i].compose(all[i].inverse()) == group_elem(D4::Id));
        assert(all[i].angle_sign == all[i].a * all[i].d - all[i].b * all[i].c);  // sign is det
        for (const GroupElem& g : all) {
            assert(std::find(all.begin(), all.end(), all[i].compose(g)) != all.end());
        }
    }
    const GroupElem r = group_elem(D4::Rot90);
    assert(r.compose(r) == group_elem(D4::Rot180));
    assert(r.compose(r).compose(r) == group_elem(D4::Rot270));
    assert(group_elem(D4::MirrorX).compose(group_elem(D4::MirrorY)) == group_elem(D4::Rot180));
    assert(group_elem(D4::MirrorDiag).compose(group_elem(D4::MirrorAnti)) ==
           group_elem(D4::Rot180));
}

void test_elements_act_about_the_centre() {
    const double L = 10.0;
    const Pose p{1.0, 3.0, 0.2};
    const auto at = [&](D4 w, double x, double y) {
        const Pose q = group_elem(w).apply(p, L);
        assert(std::abs(q.x - x) < 1e-12 && std::abs(q.y - y) < 1e-12);
    };
    at(D4::Rot90, 7.0, 1.0);
    at(D4::Rot180, 9.0, 7.0);
    at(D4::Rot270, 3.0, 9.0);
    at(D4::MirrorX, 9.0, 3.0);
    at(D4::MirrorY, 1.0, 7.0);
    at(D4::MirrorDiag, 3.0, 1.0);
    at(D4::MirrorAnti, 7.0, 9.0);
    assert(group_elem(D4::Rot90).apply(p, L).theta == 0.2);
    assert(group_elem(D4::MirrorDiag).apply(p, L).theta == -0.2);
}

// The geometric claim behind "mirrors send theta to -theta": the image pose has to be the same
// point set as the mirrored vertices, for every element.
void test_angle_rule_matches_the_vertices() {
    const Shape square = make_square(1.0);
    const double L = 8.0;
    const Pose p{2.3, 5.1, 0.37};
    Particle before;
    before.set_pose(p.x, p.y, p.theta);
    for (const GroupElem& g : d4_elements()) {
        const Pose q = g.apply(p, L);
        Particle after;
        after.set_pose(q.x, q.y, q.theta);
        for (int k = 0; k < 4; ++k) {
            const Vec2 v = world_vertex(before, square, k);
            const Pose image = g.apply({v.x, v.y, 0.0}, L);
            bool hit = false;
            for (int m = 0; m < 4; ++m) {
                const Vec2 w = world_vertex(after, square, m);
                hit = hit || std::hypot(w.x - image.x, w.y - image.y) < 1e-12;
            }
            assert(hit);
        }
    }
}

void test_inverse_restores_the_configuration() {
    const System s = random_config(40, 9.0, 11);
    for (const GroupElem& g : d4_elements()) {
        const System there = image_of(s, g);
        const System back = image_of(there, g.inverse());
        assert(same_set(s, back, 1e-12));
    }
}

void test_detect_on_hand_built_symmetry() {
    const double L = 12.0;
    const std::vector<Pose> seed = {{1.3, 2.1, 0.11}, {3.7, 1.2, -0.31}, {4.4, 4.9, 0.27}};
    const auto el = [](D4 w) { return group_elem(w); };

    const std::vector<GroupElem> c4 = {el(D4::Id), el(D4::Rot90), el(D4::Rot180), el(D4::Rot270)};
    const auto got = detect_symmetry(symmetrised(seed, c4, L), 1e-9);
    assert(got.size() == 4);
    for (const GroupElem& g : c4) assert(contains(got, g));

    const auto c2 = detect_symmetry(symmetrised(seed, {el(D4::Id), el(D4::Rot180)}, L), 1e-9);
    assert(c2.size() == 2 && contains(c2, el(D4::Rot180)));

    const auto mx = detect_symmetry(symmetrised(seed, {el(D4::Id), el(D4::MirrorX)}, L), 1e-9);
    assert(mx.size() == 2 && contains(mx, el(D4::MirrorX)));

    const auto md = detect_symmetry(symmetrised(seed, {el(D4::Id), el(D4::MirrorDiag)}, L), 1e-9);
    assert(md.size() == 2 && contains(md, el(D4::MirrorDiag)));

    // The full group, from one generic orbit of eight.
    const auto all = detect_symmetry(
        symmetrised(seed, {d4_elements().begin(), d4_elements().end()}, L), 1e-9);
    assert(all.size() == 8);

    // A nudged copy of the C4 packing fits only at a tolerance above the nudge.
    System nudged = symmetrised(seed, c4, L);
    const Particle& p = nudged.particles[0];
    nudged.set_pose(0, p.x() + 1e-4, p.y(), p.theta());
    assert(detect_symmetry(nudged, 1e-6).size() == 1);
    assert(detect_symmetry(nudged, 1e-3).size() == 4);

    assert(detect_symmetry(random_config(30, 9.0, 5), 1e-6).size() == 1);
}

void test_detect_refuses_non_squares() {
    System s;
    s.species.push_back(make_rectangle(2.0, 1.0));
    s.L = 10.0;
    s.add_particle(5.0, 5.0, 0.0, 0);  // on the centre: every element would look like a fit
    s.recount();
    assert(!is_square_config(s));
    assert(detect_symmetry(s, 1e-6).size() == 1);
}

void test_angle_wraparound_at_45_degrees() {
    assert(std::abs(wrap90(44.0 * kDeg) - 44.0 * kDeg) < 1e-12);
    assert(std::abs(wrap90(46.0 * kDeg) + 44.0 * kDeg) < 1e-12);
    assert(std::abs(wrap90(-46.0 * kDeg) - 44.0 * kDeg) < 1e-12);
    assert(wrap90(45.0 * kDeg) < 0.0);  // half open: +45 is reported as -45
    assert(std::abs(wrap90(45.0 * kDeg) + 45.0 * kDeg) < 1e-12);
    assert(std::abs(wrap90(-45.0 * kDeg) + 45.0 * kDeg) < 1e-12);

    assert(std::abs(angle_distance90(44.5 * kDeg, -44.5 * kDeg) - 1.0 * kDeg) < 1e-12);
    assert(angle_distance90(45.0 * kDeg, -45.0 * kDeg) < 1e-12);
    assert(std::abs(angle_distance90(10.0 * kDeg, 100.0 * kDeg)) < 1e-12);

    // Near +-45 a mirror moves a square across the wrap, and it must still be recognised.
    const double L = 8.0;
    const double t = 45.0 * kDeg - 1e-11;
    const System s = config_from({{2.0, 4.0, t}, {6.0, 4.0, -t}}, L);
    assert(contains(detect_symmetry(s, 1e-9), group_elem(D4::MirrorX)));

    // And +45 and -45 are one orientation to the dedupe key, from either side of the wrap.
    const System up = config_from({{2.0, 3.0, 45.0 * kDeg - 1e-12}}, L);
    const System down = config_from({{2.0, 3.0, -45.0 * kDeg + 1e-12}}, L);
    assert(canonical_key(up) == canonical_key(down));
}

void test_canonical_form_ignores_symmetry_and_labels() {
    const System s = random_config(25, 8.0, 21);
    const CanonicalKey key = canonical_key(s);
    for (const GroupElem& g : d4_elements()) {
        System image = image_of(s, g);
        std::mt19937_64 rng(3);
        std::shuffle(image.particles.begin(), image.particles.end(), rng);
        image.recount();
        assert(canonical_key(image) == key);
    }

    const System c = canonicalize(s);
    assert(canonical_key(c) == key);
    assert(c.particles.size() == s.particles.size());
    assert(canonical_key(canonicalize(c)) == key);  // idempotent
    assert(std::abs(c.total_energy() - s.total_energy()) < 1e-9);

    bool on_some_image = false;
    for (const GroupElem& g : d4_elements()) {
        on_some_image = on_some_image || same_set(c, image_of(s, g), 1e-12);
    }
    assert(on_some_image);

    // Different configurations must not collide.
    assert(canonical_key(random_config(25, 8.0, 22)) != key);
}

void test_canonicalize_keeps_frozen_flags_with_their_squares() {
    System s = random_config(12, 7.0, 31);
    s.particles[4].frozen = true;
    const Particle marked = s.particles[4];
    const System c = canonicalize(s);
    assert(c.mobile_count() == s.mobile_count());
    for (const GroupElem& g : d4_elements()) {
        const Pose q = g.apply({marked.x(), marked.y(), marked.theta()}, s.L);
        for (const Particle& p : c.particles) {
            if (std::hypot(p.x() - q.x, p.y() - q.y) < 1e-12 &&
                angle_distance90(p.theta(), q.theta) < 1e-12) {
                assert(p.frozen);
            }
        }
    }
}

// --- symmetrize_region --------------------------------------------------------------------------

System grid_of(int rows, int cols, double L) {
    std::vector<Pose> poses;
    for (int i = 0; i < cols; ++i) {
        for (int j = 0; j < rows; ++j) poses.push_back({i + 0.5, j + 0.5, 0.0});
    }
    return config_from(poses, L);
}

System compressed_grid(int n) {
    const CompressResult c = compress(grid_init(make_square(1.0), n));
    assert(c.found);
    return c.packing;
}

struct Snap {
    double x, y, theta, L;
};

std::vector<Snap> poses(const System& s) {
    std::vector<Snap> out;
    for (const Particle& p : s.particles) out.push_back({p.x(), p.y(), p.theta(), s.L});
    return out;
}

bool identical(const System& s, const std::vector<Snap>& was) {
    if (s.particles.size() != was.size()) return false;
    for (size_t i = 0; i < was.size(); ++i) {
        const Particle& p = s.particles[i];
        if (p.x() != was[i].x || p.y() != was[i].y || p.theta() != was[i].theta || s.L != was[i].L) {
            return false;
        }
    }
    return true;
}

bool valid_packing(const System& s) {
    return why_not_a_packing(s, s.tight_side(), true, kCertifyTol) == nullptr;
}

void test_variants_are_named_apart() {
    std::vector<std::string> names;
    for (int copies = 1; copies <= 3; ++copies) {
        for (const SymFamily f :
             {SymFamily::Rotation, SymFamily::AxisMirror, SymFamily::DiagonalMirror}) {
            const int v = symmetrize_variant(copies, f);
            assert(v >= 0 && v < kSymmetrizeVariants);
            names.push_back(symmetrize_variant_name(v));
        }
    }
    std::sort(names.begin(), names.end());
    assert(std::unique(names.begin(), names.end()) == names.end());
}

// The contract of the move, for every variant: exactly n squares, and either a valid packing no
// larger than before or the original configuration, bit for bit.
void test_symmetrize_keeps_n_and_certifies_or_restores() {
    for (const int n : {10, 16, 23}) {
        System s = compressed_grid(n);
        SymmetrizeBandit bandit;
        std::mt19937_64 rng(100 + static_cast<unsigned>(n));
        LnmOptions opt;
        int accepted = 0, rejected = 0;
        for (int trial = 0; trial < 20; ++trial) {
            const std::vector<Snap> before = poses(s);
            const double side = s.tight_side();
            double gain = -1.0;
            const bool ok = try_symmetrize_region(s, opt, bandit, rng, &gain);
            assert(static_cast<int>(s.particles.size()) == n);
            if (ok) {
                ++accepted;
                assert(s.is_feasible(kCertifyTol));
                assert(valid_packing(s));
                assert(s.tight_side() <= side);
                assert(gain >= 0.0);
            } else {
                ++rejected;
                assert(identical(s, before));
                assert(gain == 0.0);
            }
        }
        long attempts = 0;
        for (const auto& arm : bandit.arms) attempts += arm.attempts;
        assert(attempts == 20);
        std::printf("  n=%-3d symmetrize accepted %d, rejected %d\n", n, accepted, rejected);
    }
}

void test_symmetrize_runs_every_variant() {
    System s = compressed_grid(16);
    std::mt19937_64 rng(77);
    LnmOptions opt;
    opt.symmetrize.floor = 1.0;  // near-uniform, so that all nine come up quickly
    SymmetrizeBandit bandit;
    for (int trial = 0; trial < 60; ++trial) {
        try_symmetrize_region(s, opt, bandit, rng);
        assert(s.particles.size() == 16);
    }
    for (const auto& arm : bandit.arms) assert(arm.attempts > 0);
}

// Not vacuous: on a packing with room to gain the move has to be able to win.
void test_symmetrize_can_improve_a_loose_packing() {
    System s = compressed_grid(16);
    s.rescale_box(6.0);
    const double loose = s.tight_side();
    assert(valid_packing(s));

    SymmetrizeBandit bandit;
    std::mt19937_64 rng(9);
    bool won = false;
    for (int trial = 0; trial < 30 && !won; ++trial) {
        double gain = 0.0;
        won = try_symmetrize_region(s, LnmOptions{}, bandit, rng, &gain) && gain > 0.0;
    }
    assert(won);
    assert(s.tight_side() < loose);
    assert(valid_packing(s));
}

void test_symmetrize_never_moves_a_frozen_square() {
    System s = compressed_grid(16);
    s.rescale_box(5.0);  // loose, so attempts really are accepted and the check means something
    s.particles[0].frozen = true;
    s.particles[9].frozen = true;
    const Particle a = s.particles[0], b = s.particles[9];

    SymmetrizeBandit bandit;
    std::mt19937_64 rng(13);
    int accepted = 0;
    for (int trial = 0; trial < 20; ++trial) {
        if (try_symmetrize_region(s, LnmOptions{}, bandit, rng)) ++accepted;
        assert(s.particles[0].frozen && s.particles[9].frozen);
        assert(s.mobile_count() == 14);
        assert(s.particles[0].x() == a.x() && s.particles[0].y() == a.y());
        assert(s.particles[9].x() == b.x() && s.particles[9].y() == b.y());
    }
    assert(accepted > 0);
}

void test_symmetrize_refuses_other_shapes() {
    System s;
    s.species.push_back(make_rectangle(2.0, 1.0));
    s.L = 10.0;
    s.add_particle(3.0, 3.0, 0.0, 0);
    s.add_particle(7.0, 7.0, 0.0, 0);
    s.recount();
    const std::vector<Snap> before = poses(s);
    SymmetrizeBandit bandit;
    std::mt19937_64 rng(1);
    assert(!try_symmetrize_region(s, LnmOptions{}, bandit, rng));
    assert(identical(s, before));
}

void test_mirror_symmetry_copies_the_top_left_triangle() {
    System s = grid_of(4, 4, 5.0);
    for (size_t i = 0; i < s.particles.size(); ++i) {
        const Particle& p = s.particles[i];
        if (p.x() == 0.5 && p.y() == 3.5) s.set_pose(static_cast<int>(i), 0.5, 3.9, 0.0);  // breaks the symmetry
    }
    assert(!contains(detect_symmetry(s, 1e-9), group_elem(D4::MirrorDiag)));

    assert(mirror_symmetry(s));
    assert(s.particles.size() == 16);
    assert(contains(detect_symmetry(s, 1e-9), group_elem(D4::MirrorDiag)));
    assert(s.is_feasible(kCertifyTol));
    bool kept = false;
    for (const Particle& p : s.particles) kept = kept || (p.x() == 0.5 && p.y() == 3.9);
    assert(kept);  // the source side is the one that survives

    System unequal = random_config(41, 9.0, 3);  // the two triangles hold different counts
    assert(mirror_symmetry(unequal));
    assert(unequal.particles.size() == 41);
    assert(std::abs(unequal.energy() - unequal.total_energy()) < 1e-9);

    System other;
    other.species.push_back(make_rectangle(2.0, 1.0));
    other.L = 10.0;
    other.add_particle(3.0, 7.0, 0.0, 0);
    other.add_particle(7.0, 3.0, 0.0, 0);
    other.recount();
    assert(!mirror_symmetry(other));
}

// Each destination starts as a nudged copy of its image, so the counts match, the symmetry is broken, and only the operation can restore it.
void test_apply_symmetry_builds_each_operation() {
    const double L = 10.0;
    const std::vector<Pose> seed = {{1.5, 3.2, 0.1}, {3.3, 4.4, 0.2}, {2.0, 4.0, -0.3}};
    struct Op {
        BoxRegion source;
        std::vector<D4> images;
        D4 expected;
    };
    const std::vector<Op> ops = {
        {kBottomLeftQuadrant, {D4::Rot90, D4::Rot180, D4::Rot270}, D4::Rot90},
        {kBottomLeftQuadrant, {D4::MirrorX, D4::MirrorY, D4::Rot180}, D4::MirrorX},
        {kBottomLeftQuadrant, {D4::Rot90}, D4::Id},
        {kLeftHalf, {D4::MirrorX}, D4::MirrorX},
        {kBottomHalf, {D4::Rot180}, D4::Rot180},
        {kTopLeftTriangle, {D4::MirrorDiag}, D4::MirrorDiag},
        {kBottomLeftTriangle, {D4::MirrorAnti}, D4::MirrorAnti},
        {kLeftSector, {D4::MirrorDiag, D4::MirrorAnti, D4::Rot180}, D4::MirrorDiag},
    };
    for (const Op& op : ops) {
        std::vector<Pose> poses = seed;
        for (const Pose& p : seed) {
            for (const D4 w : op.images) {
                const Pose q = group_elem(w).apply(p, L);
                poses.push_back({q.x + 0.05, q.y + 0.05, q.theta});
            }
        }
        System s = config_from(poses, L);
        assert(apply_symmetry(s, op.source, op.images));
        assert(s.particles.size() == poses.size());
        assert(contains(detect_symmetry(s, 1e-9), group_elem(op.expected)));
        for (const Pose& p : seed) {  // the source itself did not move
            bool found = false;
            for (const Particle& q : s.particles) found = found || (q.x() == p.x && q.y() == p.y);
            assert(found);
        }
    }
}

void test_bandit_follows_the_recent_rate() {
    SymmetrizeBandit bandit;
    std::mt19937_64 rng(5);
    const auto frequencies = [&](double floor) {
        std::array<int, kSymmetrizeVariants> count{};
        for (int k = 0; k < 90000; ++k) ++count[static_cast<size_t>(bandit.choose(rng, floor))];
        return count;
    };

    for (const int c : frequencies(0.05)) assert(c > 9000 && c < 11000);  // nothing learned: even

    for (int k = 0; k < 50; ++k) bandit.update(4, true, 0.5, 0.1);
    assert(bandit.arms[4].attempts == 50 && bandit.arms[4].improvements == 50);
    assert(std::abs(bandit.arms[4].cpu_seconds - 25.0) < 1e-9);
    assert(bandit.arms[4].recent > 0.9);
    const auto skewed = frequencies(0.05);
    assert(skewed[4] > 50000);
    for (size_t k = 0; k < skewed.size(); ++k) {
        if (k != 4) assert(skewed[k] > 0);  // the floor: no variant starves
    }

    for (int k = 0; k < 200; ++k) bandit.update(4, false, 0.0, 0.1);
    assert(bandit.arms[4].recent < 0.01);  // recent, not lifetime
    assert(bandit.arms[4].improvements == 50);
}

void test_lnm_round_can_draw_only_symmetrize() {
    System s = compressed_grid(10);
    LnmOptions opt;
    opt.w_destroy = opt.w_strip = opt.w_corner = opt.w_loosen = 0.0;
    opt.w_symmetrize = 1.0;
    LnmStats st;
    std::mt19937_64 rng(3);
    for (int k = 0; k < 8; ++k) lnm_round(s, opt, st, 0.0, 0.1, rng);
    assert(st.symmetrize.attempted == 8);
    assert(st.destroy.attempted == 0 && st.strip.attempted == 0);
    assert(st.corner.attempted == 0 && st.loosen.attempted == 0);
    assert(lnm_total(st).attempted == 8);
}

// With the new move off and the logging on, a fixed seed gives the same run to the last bit.
// The comparison against the code as it was before these moves existed is not repeatable from
// inside the repository; this pins what can be: that neither the option nor the log can change
// the random stream.
void test_disabled_moves_and_logging_change_nothing() {
    const auto run = [](bool instrumented) {
        System s = compressed_grid(23);
        std::mt19937_64 rng(2024);
        StepSizes steps;
        MoveLog log;
        LnmOptions opt;
        LnmStats stats;
        if (instrumented) {
            opt.log = &log;
            opt.w_symmetrize = 0.0;
            opt.symmetrize.eps = 0.9;  // meaningless at weight 0, and must stay so
        }
        for (int cycle = 0; cycle < 20; ++cycle) {
            for (int k = 0; k < 10; ++k) sweep(s, steps, 0.05, rng, instrumented ? &log : nullptr);
            lnm_round(s, opt, stats, 0.0, 0.05, rng);
        }
        if (instrumented) {
            assert(log.at(MoveId::Translate).attempts > 0);
            assert(log.at(MoveId::Symmetrize).attempts == 0);
        }
        return s;
    };
    const System plain = run(false), logged = run(true);
    assert(identical(logged, poses(plain)));
}

// The viewer's hook sees the attempt, then the relaxation, then the verdict, and changes nothing.
void test_watch_sees_each_attempt_and_changes_nothing() {
    const auto run = [](bool watched, std::vector<std::string>* seen) {
        System s = compressed_grid(16);
        s.rescale_box(5.0);
        LnmOptions opt;
        opt.w_destroy = 1.0;
        opt.w_symmetrize = 1.0;
        if (watched) {
            opt.watch = [seen](const System& frame, const char* label) {
                assert(frame.particles.size() == 16);
                seen->push_back(label);
            };
        }
        LnmStats stats;
        std::mt19937_64 rng(6);
        for (int k = 0; k < 12; ++k) lnm_round(s, opt, stats, 0.0, 0.05, rng);
        return s;
    };
    std::vector<std::string> seen;
    const System quiet = run(false, nullptr), watched = run(true, &seen);
    assert(identical(watched, poses(quiet)));

    const auto index_of = [&](const std::string& prefix) {
        for (size_t i = 0; i < seen.size(); ++i) {
            if (seen[i].rfind(prefix, 0) == 0) return static_cast<long>(i);
        }
        return -1L;
    };
    assert(index_of("destroy") >= 0);
    assert(index_of("symmetrize ") >= 0);
    const auto verdict = std::find_if(seen.begin() + index_of("destroy"), seen.end(), [](const std::string& l) {
        return l == "accepted" || l == "rejected";
    });
    assert(verdict != seen.end());
}

void test_frames_carry_their_label() {
    const std::string path = "build/test_label.traj";
    const System s = grid_of(2, 2, 2.0);
    std::FILE* f = open_trajectory(path, s);
    assert(f);
    append_frame(f, 3, s, "mirror");
    append_frame(f, 4, s);
    std::fclose(f);

    f = std::fopen(path.c_str(), "r");
    assert(f);
    std::string text;
    for (int ch = std::fgetc(f); ch != EOF; ch = std::fgetc(f)) text.push_back(static_cast<char>(ch));
    std::fclose(f);
    std::remove(path.c_str());
    assert(text.find("m mirror\nf 3 ") != std::string::npos);
    assert(text.find("f 4 ") != std::string::npos);
    assert(text.find("m mirror", text.find("f 3 ")) == std::string::npos);
}

void test_move_log_counts_and_writes() {
    MoveLog log;
    log.record(MoveId::Flip, true, 0.0, 0.25);
    log.record(MoveId::Flip, false, 0.5, 0.75);
    const MoveCounters& c = log.at(MoveId::Flip);
    assert(c.attempts == 2 && c.accepted == 1 && c.improved == 1);
    assert(c.cpu_seconds == 1.0 && c.improvement == 0.5);

    { MoveTimer off(nullptr, MoveId::Flip); off.finish(true, 1.0); }  // a null log is a no-op

    const std::string path = "build/test_move_log.csv";
    assert(log.write_csv(path));
    std::FILE* f = std::fopen(path.c_str(), "r");
    assert(f);
    int lines = 0;
    for (int ch = std::fgetc(f); ch != EOF; ch = std::fgetc(f)) lines += ch == '\n';
    std::fclose(f);
    std::remove(path.c_str());
    assert(lines == static_cast<int>(MoveId::Count) + 1);  // header, then every row
}

// --- expand -------------------------------------------------------------------------------------

void test_symmetric_grow() {
    assert(!symmetric_grow(random_config(8, 6.0, 1), false));  // no C4, no candidate

    const System corners =
        config_from({{0.5, 0.5, 0.0}, {3.5, 0.5, 0.0}, {3.5, 3.5, 0.0}, {0.5, 3.5, 0.0}}, 4.0);
    const auto four = symmetric_grow(corners, false);
    assert(four && four->particles.size() == 8);
    assert(four->L > corners.L);
    assert(contains(detect_symmetry(*four, 1e-9), group_elem(D4::Rot90)));

    const auto one = symmetric_grow(corners, true);
    assert(one && one->particles.size() == 5);
    const Particle& last = one->particles.back();
    assert(std::abs(last.x() - 0.5 * one->L) < 1e-12 && std::abs(last.y() - 0.5 * one->L) < 1e-12);

    // A square already on the centre leaves no hole there.
    const System full = grid_of(3, 3, 3.0);
    assert(!symmetric_grow(full, true));
}

void test_row_operators() {
    std::mt19937_64 rng(4);
    const System grid = grid_of(3, 3, 3.0);
    for (const int direction : {0, 1}) {
        const auto dup = row_duplicate(grid, direction, rng);
        assert(dup && dup->particles.size() == 12);
        assert(dup->L == 4.0);
        assert(dup->is_feasible(kCertifyTol));

        // Tight in the other dimension, so the row goes but the box cannot shrink.
        const auto rem = row_remove(grid, direction, rng);
        assert(rem && rem->particles.size() == 6);
        assert(rem->L == 3.0);
        assert(rem->is_feasible(kCertifyTol));
    }

    // With the slack a duplicate leaves behind, removal really is its reverse.
    const auto narrower = row_remove(grid_of(3, 4, 4.0), 0, rng);
    assert(narrower && narrower->particles.size() == 9);
    assert(narrower->L == 3.0);
    assert(narrower->is_feasible(kCertifyTol));
}

void test_frame_grow_adds_the_l_shaped_layer() {
    const System grid = grid_of(3, 3, 3.0);
    const auto grown = frame_grow(grid);
    assert(grown && grown->particles.size() == 9 + 2 * 3 + 1);
    assert(grown->L == 4.0);
    assert(grown->is_feasible(kCertifyTol));
    assert(same_set(*grown, grid_of(4, 4, 4.0), 1e-12));
}

void test_delete_operators() {
    const System s = config_from({{0.5, 0.5, 0.0}, {1.5, 0.5, 0.0}, {0.5, 1.5, 0.0}, {1.5, 1.5, 0.0},
                                  {2.0, 3.0, 0.2}},
                                 4.0);
    const auto cut = delete_rattler(s);
    assert(cut && cut->particles.size() == 4);
    for (const Particle& p : cut->particles) assert(std::hypot(p.x() - 2.0, p.y() - 3.0) > 0.5);

    const auto several = delete_for_compress(s, 3);
    assert(several.size() == 3);
    for (const System& c : several) assert(c.particles.size() == 4);
}

void test_expand_refuses_frozen_squares() {
    System s = grid_of(3, 3, 3.0);
    s.particles[0].frozen = true;
    std::mt19937_64 rng(1);
    assert(!frame_grow(s));
    assert(!row_duplicate(s, 0, rng));
    assert(!row_remove(s, 0, rng));
    assert(!delete_rattler(s));
    assert(delete_for_compress(s, 2).empty());
}

void test_pool_admits_only_valid_strict_improvements() {
    Pool pool;
    const System grid = grid_of(2, 2, 2.0);
    assert(pool.offer(grid));
    assert(!pool.offer(grid));  // equal is not an improvement

    const System overlapping = config_from({{0.5, 0.5, 0.0}, {0.9, 0.5, 0.0}}, 2.0);
    assert(!pool.offer(overlapping));
    assert(!pool.at(2));

    System worse = grid;
    worse.rescale_box(3.0);
    assert(!pool.offer(worse));
    assert(pool.at(4)->side == 2.0);

    assert(pool.first_seen(grid));
    assert(!pool.first_seen(image_of(grid, group_elem(D4::Rot90))));  // same packing, turned
}

void test_propagate_updates_on_strict_improvement_only() {
    Pool pool;
    assert(pool.offer(grid_of(2, 2, 2.0)));                                    // n = 4, optimal
    const System poor = config_from({{0.5, 0.5, 0.0}, {2.5, 0.5, 0.0}, {0.5, 2.5, 0.0}}, 3.0);
    assert(pool.offer(poor));                                                   // n = 3, side 3
    assert(pool.offer(grid_of(3, 3, 3.0)));                                    // n = 9, optimal
    const double nine = pool.at(9)->side;

    std::mt19937_64 rng(8);
    MoveLog log;
    const PropagateReport report = propagate(pool, 4, PropagateOptions{}, rng, &log);
    assert(report.tried > 0);
    assert(report.improved.count(3) == 1 && report.improved.at(3) > 0.9);
    assert(std::abs(pool.at(3)->side - 2.0) < 1e-3);  // 3 squares in the 2 x 2 box
    assert(pool.at(9)->side == nine);                // the frame reaches 9 and cannot beat it
    assert(pool.at(4)->side == 2.0);
    assert(log.at(MoveId::DeleteRattler).attempts + log.at(MoveId::DeleteCompress).attempts > 0);

    for (const int n : {3, 4, 9}) assert(valid_packing(pool.at(n)->packing));

    // Filtered to a target nothing can improve, it does nothing to anything else.
    PropagateOptions only;
    only.only_target = 9;
    const double three = pool.at(3)->side;
    assert(propagate(pool, 4, only, rng, nullptr).improved.empty());
    assert(pool.at(3)->side == three);
}

}  // namespace

int main() {
    test_group_table();
    test_elements_act_about_the_centre();
    test_angle_rule_matches_the_vertices();
    test_inverse_restores_the_configuration();
    test_detect_on_hand_built_symmetry();
    test_detect_refuses_non_squares();
    test_angle_wraparound_at_45_degrees();
    test_canonical_form_ignores_symmetry_and_labels();
    test_canonicalize_keeps_frozen_flags_with_their_squares();
    test_variants_are_named_apart();
    test_symmetrize_keeps_n_and_certifies_or_restores();
    test_symmetrize_runs_every_variant();
    test_symmetrize_can_improve_a_loose_packing();
    test_symmetrize_never_moves_a_frozen_square();
    test_symmetrize_refuses_other_shapes();
    test_mirror_symmetry_copies_the_top_left_triangle();
    test_apply_symmetry_builds_each_operation();
    test_bandit_follows_the_recent_rate();
    test_lnm_round_can_draw_only_symmetrize();
    test_disabled_moves_and_logging_change_nothing();
    test_watch_sees_each_attempt_and_changes_nothing();
    test_frames_carry_their_label();
    test_move_log_counts_and_writes();
    test_symmetric_grow();
    test_row_operators();
    test_frame_grow_adds_the_l_shaped_layer();
    test_delete_operators();
    test_expand_refuses_frozen_squares();
    test_pool_admits_only_valid_strict_improvements();
    test_propagate_updates_on_strict_improvement_only();
    std::printf("symmetry: all tests passed\n");
    return 0;
}
