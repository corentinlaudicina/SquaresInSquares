#include "moves/lnm.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace pack {
namespace {

double uniform(std::mt19937_64& rng, double lo, double hi) {
    if (!(hi > lo)) return lo;
    std::uniform_real_distribution<double> d(lo, hi);
    return d(rng);
}

Vec2 perp(Vec2 v) { return {-v.y, v.x}; }

// Support interval of the box's four corners along `axis`: how far along a direction the box reaches, which is what bounds a row of squares and the offset of a cut.
void box_span(double L, Vec2 axis, double& lo, double& hi) {
    const Vec2 corners[4] = {{0.0, 0.0}, {L, 0.0}, {L, L}, {0.0, L}};
    lo = std::numeric_limits<double>::infinity();
    hi = -lo;
    for (const Vec2 c : corners) {
        const double s = dot(c, axis);
        lo = std::min(lo, s);
        hi = std::max(hi, s);
    }
}

// One of the four lattice directions a packing of squares is built on.
Vec2 lattice_direction(std::mt19937_64& rng, bool diagonal) {
    const double k = 1.0 / std::sqrt(2.0);
    const bool second = (rng() & 1u) != 0u;
    if (diagonal) return second ? Vec2{k, -k} : Vec2{k, k};
    return second ? Vec2{0.0, 1.0} : Vec2{1.0, 0.0};
}

// Where the lifted squares wait.
Vec2 park_point(const System& sys) {
    const double far = sys.L + 4.0 * sys.base_shape().circumradius + 1.0;
    return {-far, -far};
}

void park(System& sys, const std::vector<int>& members) {
    const Vec2 p = park_point(sys);
    for (const int i : members) sys.set_pose(i, p.x, p.y, 0.0);
}

// A multiple of 45 degrees with probability snap_prob, a free angle otherwise.
double draw_theta(const LnmOptions& opt, std::mt19937_64& rng) {
    if (uniform(rng, 0.0, 1.0) < opt.snap_prob) {
        return 0.25 * kPi * static_cast<double>(rng() % 8u);
    }
    return uniform(rng, 0.0, 2.0 * kPi);
}

// Best of opt.insert_tries poses from `draw`, scored by particle_energy. Greedy and one square at a time: a proposal does not have to be good, only good often enough to be worth relaxing.
template <class Draw>
void repair_one(System& sys, int i, const LnmOptions& opt, std::mt19937_64& rng, Draw&& draw) {
    const int tries = std::max(1, opt.insert_tries);
    Pose best;
    double best_e = std::numeric_limits<double>::infinity();
    for (int k = 0; k < tries; ++k) {
        const Pose cand = draw(rng);
        sys.set_pose(i, cand.x, cand.y, cand.theta);
        const double e = sys.particle_energy(i);
        if (e < best_e) {
            best_e = e;
            best = cand;
        }
    }
    sys.set_pose(i, best.x, best.y, best.theta);
}

// Copy the system, perturb it, relax, accept on enthalpy. The one place the contract described at the top of lnm.hpp lives, so that no individual move can get the reject path wrong.
template <class Perturb>
bool propose(const char* name, System& sys, const LnmOptions& opt, double pressure,
             double temperature, std::mt19937_64& rng, Perturb&& perturb) {
    // Whole-system copy: the reject path has to be exact, L and the two tracked energy sums and all. A perturbation below is free to leave both stale -- relax() is what makes them true again, and nothing reads them in between.
    const System saved = sys;
    const double h0 = enthalpy(sys, pressure);

    if (!perturb(sys, rng)) {
        sys = saved;  // a perturbation can bail out half-done, so restore rather than trust it
        return false;
    }

    if (opt.watch) opt.watch(sys, name);

    RelaxOptions ro;
    ro.max_iterations = std::max<long>(1, opt.relax_iterations);
    ro.feasible_tol = kFeasibleTol;  // a proposal, not a certification
    ro.watch = opt.watch;
    const RelaxResult r = relax(sys, ro);

    // A square wholly outside the box sits where the wall penalty has saturated: its gradient is exactly zero and no amount of relaxation brings it back (the trap main.cpp documents). That configuration is not a candidate, it is a leak, so it is never accepted.
    if (r.stop != RelaxResult::Stop::Escaped &&
        metropolis(enthalpy(sys, pressure) - h0, temperature, rng)) {
        if (opt.watch) opt.watch(sys, "accepted");
        return true;
    }
    sys = saved;
    if (opt.watch) opt.watch(sys, "rejected");
    return false;
}

// --- the perturbations ------------------------------------------------------------------------

bool destroy_and_repair(System& sys, const LnmOptions& opt, std::mt19937_64& rng) {
    const Pick p = pick_region(sys, opt.regions, rng);
    if (p.members.empty()) return false;

    // Put them back inside the hole they left, not anywhere in the box: that is what makes this a local move, and a local move is one the relaxation can finish.
    const Bounds b = member_bounds(sys, p.members, 0.5 * p.extent);

    std::vector<int> order = p.members;
    // Shuffled, because the first square back has the most room to choose from: a fixed order would hand that advantage to the same member every time.
    std::shuffle(order.begin(), order.end(), rng);

    park(sys, order);
    for (const int i : order) {
        repair_one(sys, i, opt, rng, [&](std::mt19937_64& g) {
            return Pose{uniform(g, b.xlo, b.xhi), uniform(g, b.ylo, b.yhi), draw_theta(opt, g)};
        });
    }
    return true;  // the running energy and parts are stale; relax() recounts both from scratch
}

// Re-lay `members` as a tiled strip: rows along `dir` at `pitch`, each row offset by `pitch` along the normal, from the lowest member in both coordinates, every square at `theta`. A band can hold more squares than fit across the box, so a row that runs out wraps to the next one.
void relay_strip(System& sys, std::vector<int> members, Vec2 dir, double pitch, double theta) {
    if (members.empty() || !(pitch > 0.0)) return;
    const Vec2 n = perp(dir);

    // Order along the strip, so squares keep the neighbours they had instead of being shuffled from one end of the band to the other: the smaller the move, the likelier it is accepted.
    std::sort(members.begin(), members.end(), [&](int a, int c) {
        return dot(Vec2{sys.particles[static_cast<size_t>(a)].x(),
                        sys.particles[static_cast<size_t>(a)].y()},
                   dir) < dot(Vec2{sys.particles[static_cast<size_t>(c)].x(),
                                   sys.particles[static_cast<size_t>(c)].y()},
                              dir);
    });

    double t_min = std::numeric_limits<double>::infinity(), s_min = t_min;
    for (const int i : members) {
        const Vec2 c{sys.particles[static_cast<size_t>(i)].x(), sys.particles[static_cast<size_t>(i)].y()};
        t_min = std::min(t_min, dot(c, dir));
        s_min = std::min(s_min, dot(c, n));
    }

    double t_lo = 0.0, t_hi = 0.0;
    box_span(sys.L, dir, t_lo, t_hi);
    const int per_row = std::max(1, static_cast<int>(std::floor((t_hi - t_min) / pitch)));

    for (size_t k = 0; k < members.size(); ++k) {
        const int row = static_cast<int>(k) / per_row;
        const int col = static_cast<int>(k) % per_row;
        const double t = t_min + col * pitch;
        const double s = s_min + row * pitch;
        sys.set_pose(members[k], t * dir.x + s * n.x,
                                                            t * dir.y + s * n.y, theta);
    }
}

bool carve_or_straighten(System& sys, const LnmOptions& opt, std::mt19937_64& rng) {
    // Half carve, half straighten. A move that only ever tilted bands would be a one-way street, and a chain that can only go one way is not exploring, it is drifting.
    const bool carve = (rng() & 1u) != 0u;
    const Vec2 dir = lattice_direction(rng, carve);
    const Pick p = pick_band(sys, opt.regions, dir, rng);
    if (p.members.size() < 2) return false;

    // Carving tilts every square to match the band, so its faces are perpendicular to the band and consecutive squares meet face to face at contact pitch. Straightening sets them square-on in an axis-aligned band, which is the same construction on an axis.
    const double theta = carve ? std::atan2(dir.y, dir.x) : 0.0;
    relay_strip(sys, p.members, p.dir, contact_pitch(sys.base_shape()), theta);
    return true;
}

bool reshape_corner(System& sys, const LnmOptions& opt, std::mt19937_64& rng) {
    const double diameter =
        (sys.base_shape().circumradius > 0.0) ? 2.0 * sys.base_shape().circumradius : 1.0;
    const int corner = static_cast<int>(rng() % 4u);
    const double depth =
        uniform(rng, opt.regions.corner_depth_lo, opt.regions.corner_depth_hi) * diameter;

    const std::vector<int> members = select_corner(sys, corner, depth);
    if (members.size() < 2) return false;

    const double pitch = contact_pitch(sys.base_shape());
    const int m = static_cast<int>(members.size());

    // Rows and columns the staircase may use, taken from the count rather than from the depth. A well-packed depth x depth block holds close to (depth / pitch)^2 squares, so a staircase allowed only floor(depth / pitch) rows of that length is a square short of being able to hold them, and the move then declines to propose rather than proposing anything. Sizing it from m instead guarantees the room exists.
    const int w = std::max(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(m)))));

    // A staircase is a non-increasing row profile summing to m: draw each row no longer than the one below it. The lower bound on each draw is what keeps the profile finishable -- with `left` rows still to come, none longer than this one, a row shorter than ceil(remaining / left) cannot be completed, and a profile drawn without that bound comes out far too small far too often. Which is exactly what made this move almost never fire.
    std::vector<int> rows;
    int remaining = m, cap = w;
    while (remaining > 0 && static_cast<int>(rows.size()) < w) {
        const int left = w - static_cast<int>(rows.size());
        const int hi = std::min(cap, remaining);
        const int need = (remaining + left - 1) / left;
        const int lo = std::min(hi, std::max(1, need));
        const int len = lo + static_cast<int>(rng() % static_cast<unsigned>(hi - lo + 1));
        rows.push_back(len);
        remaining -= len;
        cap = len;
    }
    if (remaining != 0) return false;  // cannot happen with the bound above; cheap to be sure

    // Lay it out from the corner inwards, axis-aligned to match the two walls that meet there.
    const double cx = (corner == 1 || corner == 2) ? sys.L : 0.0;
    const double cy = (corner == 2 || corner == 3) ? sys.L : 0.0;
    const double sx = (corner == 1 || corner == 2) ? -1.0 : 1.0;
    const double sy = (corner == 2 || corner == 3) ? -1.0 : 1.0;

    size_t k = 0;
    for (size_t r = 0; r < rows.size(); ++r) {
        for (int c = 0; c < rows[r] && k < members.size(); ++c, ++k) {
            sys.set_pose(members[k],
                cx + sx * (c + 0.5) * pitch, cy + sy * (static_cast<double>(r) + 0.5) * pitch,
                0.0);
        }
    }
    return true;
}

bool loosen_and_recompress(System& sys, const LnmOptions& opt, double temperature,
                           std::mt19937_64& rng) {
    if (sys.particles.empty()) return false;
    const double L0 = sys.L;
    const double grow = std::max(0.0, opt.loosen_grow);
    if (!(grow > 0.0)) return false;

    const bool soft = uniform(rng, 0.0, 1.0) < opt.loosen_soft_prob;
    const double k_pair0 = sys.k_pair, k_wall0 = sys.k_wall;

    sys.rescale_box(L0 * (1.0 + grow));  // carries the centres out with the walls, opening gaps
    if (soft) {
        sys.k_pair *= opt.loosen_soften;
        sys.k_wall *= opt.loosen_soften;
    }

    // Shake with the ordinary sweep, at the caller's temperature. The amplitude is set from the
    // slack just created: a shake larger than the gap it has to explore only gets rejected.
    StepSizes shake;
    shake.translate = std::max(1e-9, 0.5 * grow * L0);
    shake.rotate = 0.25 * kPi;
    for (long k = 0; k < std::max<long>(0, opt.loosen_sweeps); ++k) {
        sweep(sys, shake, temperature, rng);
    }

    if (soft) {
        sys.k_pair = k_pair0;
        sys.k_wall = k_wall0;  // hardened again before the squeeze: the rigid half of the move
    }
    sys.rescale_box(L0);  // back to the box it started in, which recounts the energy as well
    return true;
}

}  // namespace

void park_squares(System& sys, const std::vector<int>& members) { park(sys, members); }

void reinsert_in_bounds(System& sys, int i, const LnmOptions& opt, const Bounds& b,
                        std::mt19937_64& rng) {
    repair_one(sys, i, opt, rng, [&](std::mt19937_64& g) {
        return Pose{uniform(g, b.xlo, b.xhi), uniform(g, b.ylo, b.yhi), draw_theta(opt, g)};
    });
}

Stats lnm_total(const LnmStats& s) {
    Stats out;
    for (const Stats* c : {&s.destroy, &s.corner, &s.strip, &s.loosen, &s.symmetrize}) {
        out.attempted += c->attempted;
        out.accepted += c->accepted;
    }
    return out;
}

bool try_destroy_repair(System& sys, const LnmOptions& opt, double pressure, double temperature,
                        std::mt19937_64& rng) {
    return propose("destroy", sys, opt, pressure, temperature, rng,
                   [&](System& s, std::mt19937_64& g) { return destroy_and_repair(s, opt, g); });
}

bool try_reshape_corner(System& sys, const LnmOptions& opt, double pressure, double temperature,
                        std::mt19937_64& rng) {
    return propose("corner", sys, opt, pressure, temperature, rng,
                   [&](System& s, std::mt19937_64& g) { return reshape_corner(s, opt, g); });
}

bool try_strip(System& sys, const LnmOptions& opt, double pressure, double temperature,
               std::mt19937_64& rng) {
    return propose("strip", sys, opt, pressure, temperature, rng,
                   [&](System& s, std::mt19937_64& g) { return carve_or_straighten(s, opt, g); });
}

bool try_loosen(System& sys, const LnmOptions& opt, double pressure, double temperature,
                std::mt19937_64& rng) {
    return propose("loosen", sys, opt, pressure, temperature, rng, [&](System& s, std::mt19937_64& g) {
        return loosen_and_recompress(s, opt, temperature, g);
    });
}

bool lnm_round(System& sys, const LnmOptions& opt, LnmStats& stats, double pressure,
               double temperature, std::mt19937_64& rng) {
    const double base = std::max(0.0, opt.w_destroy) + std::max(0.0, opt.w_strip) +
                        std::max(0.0, opt.w_corner) + std::max(0.0, opt.w_loosen);
    const double total = base + std::max(0.0, opt.w_symmetrize);
    if (!(total > 0.0)) return false;
    const double r = uniform(rng, 0.0, total);

    Stats* bucket = nullptr;
    MoveId id = MoveId::Loosen;
    bool ok = false;
    if (r < opt.w_destroy) {
        bucket = &stats.destroy;
        id = MoveId::Destroy;
        MoveTimer timer(opt.log, id);
        ok = try_destroy_repair(sys, opt, pressure, temperature, rng);
        timer.finish(ok);
    } else if (r < opt.w_destroy + opt.w_strip) {
        bucket = &stats.strip;
        id = MoveId::Strip;
        MoveTimer timer(opt.log, id);
        ok = try_strip(sys, opt, pressure, temperature, rng);
        timer.finish(ok);
    } else if (r < opt.w_destroy + opt.w_strip + opt.w_corner) {
        bucket = &stats.corner;
        id = MoveId::Corner;
        MoveTimer timer(opt.log, id);
        ok = try_reshape_corner(sys, opt, pressure, temperature, rng);
        timer.finish(ok);
    } else if (opt.w_symmetrize > 0.0 && r >= base) {
        bucket = &stats.symmetrize;
        id = MoveId::Symmetrize;
        MoveTimer timer(opt.log, id);
        double gain = 0.0;
        ok = try_symmetrize_region(sys, opt, stats.bandit, rng, &gain);
        timer.finish(ok, gain);
    } else {
        bucket = &stats.loosen;
        MoveTimer timer(opt.log, id);
        ok = try_loosen(sys, opt, pressure, temperature, rng);
        timer.finish(ok);
    }
    ++bucket->attempted;
    if (ok) {
        ++bucket->accepted;
        stats.last = id;
    }
    return ok;
}

bool large_neighbour_move(System& sys, std::mt19937_64& rng, const LnmOptions& opt,
                          LnmStats* stats) {
    const int attempts = 100;           // a move can decline, e.g. a band that misses every square
    const double always_accept = 1e9;   // hot enough that the acceptance test never refuses
    LnmStats local;
    LnmStats& st = stats ? *stats : local;
    for (int k = 0; k < attempts; ++k) {
        if (lnm_round(sys, opt, st, 0.0, always_accept, rng)) return true;
    }
    return false;
}

}  // namespace pack
