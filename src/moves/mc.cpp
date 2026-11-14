#include "moves/mc.hpp"

#include <algorithm>
#include <vector>
#include <cmath>

namespace pack {
namespace {

// A sweep draws its regions from the defaults. There is nothing to tune here -- the sizes are drawn from a range rather than adapted, for the reason in lnm.hpp -- so threading a RegionWeights through sweep() and every one of its callers would buy nothing.
constexpr RegionWeights kSweepRegions{};

double rate_of(const Stats& c) {
    return (c.attempted > 0) ? static_cast<double>(c.accepted) / c.attempted : 0.0;
}

EnergyParts group_parts(const System& sys, const std::vector<int>& group) {
    EnergyParts out;
    for (const int i : group) {
        const EnergyParts one = sys.particle_parts(i);
        out.pair += one.pair;
        out.wall += one.wall;
    }
    return out;
}

double energy_change(const System& sys, const EnergyParts& before, const EnergyParts& after) {
    return energy_of(after, sys.k_pair, sys.k_wall) - energy_of(before, sys.k_pair, sys.k_wall);
}

// Commit energy of move.
void commit(System& sys, const EnergyParts& before, const EnergyParts& after) {
    sys.parts.pair += after.pair - before.pair;
    sys.parts.wall += after.wall - before.wall;
}

std::vector<Pose> snapshot(const System& sys, const std::vector<int>& group) {
    std::vector<Pose> out;
    out.reserve(group.size());
    for (const int i : group) {
        const Particle& p = sys.particles[static_cast<size_t>(i)];
        out.push_back({p.x(), p.y(), p.theta()});
    }
    return out;
}

void restore(System& sys, const std::vector<int>& group, const std::vector<Pose>& saved) {
    for (size_t k = 0; k < group.size(); ++k) {
        sys.set_pose(group[k], saved[k].x, saved[k].y, saved[k].theta);
    }
}

std::vector<Pose> snapshot_all(const System& sys) {
    std::vector<Pose> out;
    out.reserve(sys.particles.size());
    for (const Particle& p : sys.particles) out.push_back({p.x(), p.y(), p.theta()});
    return out;
}

void restore_all(System& sys, const std::vector<Pose>& saved) {
    for (size_t k = 0; k < saved.size(); ++k) {
        sys.set_pose(static_cast<int>(k), saved[k].x, saved[k].y, saved[k].theta);
    }
}

// check that a group does not contain pinned squares.
bool group_movable(const System& sys, const std::vector<int>& group) {
    if (group.empty()) return false;
    for (const int i : group) {
        if (i < 0 || i >= static_cast<int>(sys.particles.size())) return false;
        if (sys.particles[static_cast<size_t>(i)].frozen) return false;
    }
    return true;
}

}  // namespace

bool metropolis(double delta, double temperature, std::mt19937_64& rng) {
    if (delta <= 0.0) return true;
    if (temperature <= 0.0) return false;   // T == 0 is pure descent

    const double exponent = -delta / temperature;
    if (exponent < -700.0) return false;
    std::uniform_real_distribution<double> u(0.0, 1.0);
    return u(rng) < std::exp(exponent);
}

void tune_step(double& step, const Stats& trial, double lo, double hi) {
    if (trial.attempted == 0) return;
    step *= (rate_of(trial) > 0.4) ? 1.05 : 0.95;
    step = std::clamp(step, lo, hi);
}

bool try_translate(System& sys, int i, const StepSizes& steps, double temperature, std::mt19937_64& rng) {
    Particle& p = sys.particles[static_cast<size_t>(i)];
    if (p.frozen) return false;

    const double x0 = p.x(), y0 = p.y(), theta0 = p.theta();
    const EnergyParts before = sys.particle_parts(i);

    std::uniform_real_distribution<double> d(-steps.translate, steps.translate);
    sys.set_pose(i, x0 + d(rng), y0 + d(rng), theta0);

    const EnergyParts after = sys.particle_parts(i);
    if (metropolis(energy_change(sys, before, after), temperature, rng)) {
        commit(sys, before, after);
        return true;
    }
    sys.set_pose(i, x0, y0, theta0);
    return false;
}

bool try_rotate(System& sys, int i, const StepSizes& steps, double temperature,
                std::mt19937_64& rng) {
    Particle& p = sys.particles[static_cast<size_t>(i)];
    if (p.frozen) return false;

    const double x0 = p.x(), y0 = p.y(), theta0 = p.theta();
    const EnergyParts before = sys.particle_parts(i);

    std::uniform_real_distribution<double> d(-steps.rotate, steps.rotate);
    sys.set_pose(i, x0, y0, theta0 + d(rng));

    const EnergyParts after = sys.particle_parts(i);
    if (metropolis(energy_change(sys, before, after), temperature, rng)) {
        commit(sys, before, after);
        return true;
    }
    sys.set_pose(i, x0, y0, theta0);
    return false;
}

bool try_flip(System& sys, int i, double temperature, std::mt19937_64& rng) {
    Particle& p = sys.particles[static_cast<size_t>(i)];
    if (p.frozen) return false;

    const double x0 = p.x(), y0 = p.y(), theta0 = p.theta();
    const EnergyParts before = sys.particle_parts(i);

    // Exactly 45 degrees, either way: the gap between the axis-aligned and the tilted family.
    const double dth = ((rng() & 1u) != 0u) ? 0.25 * kPi : -0.25 * kPi;
    sys.set_pose(i, x0, y0, theta0 + dth);

    const EnergyParts after = sys.particle_parts(i);
    if (metropolis(energy_change(sys, before, after), temperature, rng)) {
        commit(sys, before, after);
        return true;
    }
    sys.set_pose(i, x0, y0, theta0);
    return false;
}

bool try_box(System& sys, const StepSizes& steps, double pressure, double temperature,
             std::mt19937_64& rng) {
    if (sys.particles.empty()) return false;

    const double L0 = sys.L, e0 = sys.energy();
    const EnergyParts p0 = sys.parts;  // the accept side is repriced by rescale_box; this is the other side
    std::uniform_real_distribution<double> d(-steps.box, steps.box);
    const double L1 = L0 + d(rng);
    if (!(L1 > 0.0)) return false;  // outside_area needs a box whose corners wind CCW

    const std::vector<Pose> saved = snapshot_all(sys);
    sys.rescale_box(L1);

    // The fixed-pressure acceptance. Two terms beyond the usual energy difference:
    //
    //   P * (A1 - A0)   the work done against the pressure, A = L^2 in 2D. This is what
    //                   pushes the box inwards.
    //   T * 2m ln(L1/L0)  the Jacobian of the move. The configuration integral runs over
    //                   centres measured in units of L, so rescaling the box rescales the
    //                   volume of configuration space with it, by (L1/L0)^(2m) for m mobile
    //                   squares. Drop it and the box collapses onto the squares: it is the
    //                   entropic push back outwards, and leaving it out is the classic NPT
    //                   bug, silent because the run still produces a packing, just a biased
    //                   one.
    //
    // Folding both into one effective delta keeps a single Metropolis rule for every move
    // type. At T = 0 the Jacobian term vanishes, which is the right limit -- it is entropy,
    // and entropy earns nothing at zero temperature.
    const double delta_h = (sys.energy() - e0) + pressure * (L1 * L1 - L0 * L0);
    const double delta = delta_h - temperature * 2.0 * sys.mobile_count() * std::log(L1 / L0);

    if (metropolis(delta, temperature, rng)) return true;

    restore_all(sys, saved);  // restore, never rescale back: dividing out a float multiply drifts
    sys.L = L0;
    sys.parts = p0;  // exact: restored rather than recounted, or every rejected box move leaks drift
    return false;
}

Stats sweep(System& sys, StepSizes& steps, double temperature, std::mt19937_64& rng,
            MoveLog* log, const SweepMix& mix) {
    const double kWTranslate = mix.translate, kWRotate = mix.rotate, kWFlip = mix.flip, kWShove = mix.shove;
    Stats st;
    const int n = static_cast<int>(sys.particles.size());
    if (n == 0) return st;
    std::uniform_int_distribution<int> pick(0, n - 1);
    std::uniform_real_distribution<double> which(0.0,
                                                 kWTranslate + kWRotate + kWFlip + kWShove);

    Stats ct, cr;
    for (int m = 0; m < n; ++m) {
        // One draw per move rather than cycling through the types, so a frozen or stuck
        // particle cannot lock the run into one kind of move.
        const double r = which(rng);
        bool ok = false;
        if (r < kWTranslate) {
            MoveTimer timer(log, MoveId::Translate);
            ok = try_translate(sys, pick(rng), steps, temperature, rng);
            timer.finish(ok);
            ++ct.attempted;
            if (ok) ++ct.accepted;
        } else if (r < kWTranslate + kWRotate) {
            MoveTimer timer(log, MoveId::Rotate);
            ok = try_rotate(sys, pick(rng), steps, temperature, rng);
            timer.finish(ok);
            ++cr.attempted;
            if (ok) ++cr.accepted;
        } else if (r < kWTranslate + kWRotate + kWFlip) {
            MoveTimer timer(log, MoveId::Flip);
            ok = try_flip(sys, pick(rng), temperature, rng);
            timer.finish(ok);
        } else {
            MoveTimer timer(log, MoveId::Shove);
            ok = try_shove_region(sys, kSweepRegions, steps, temperature, rng);
            timer.finish(ok);
        }
        ++st.attempted;
        if (ok) ++st.accepted;
    }

    // Each step size is tuned from its own bucket only; the flip has none and the shove reports
    // to neither. See the note in mc.hpp.
    tune_step(steps.translate, ct, 1e-9, sys.L);
    tune_step(steps.rotate, cr, 1e-9, kPi);
    return st;
}

bool try_translate_group(System& sys, const std::vector<int>& group,
                         const StepSizes& steps, double temperature, std::mt19937_64& rng) {
    if (!group_movable(sys, group)) return false;

    const std::vector<Pose> saved = snapshot(sys, group);
    const EnergyParts before = group_parts(sys, group);

    std::uniform_real_distribution<double> d(-steps.translate, steps.translate);
    const double dx = d(rng), dy = d(rng);
    for (const int i : group) {
        Particle& p = sys.particles[static_cast<size_t>(i)];
        sys.set_pose(i, p.x() + dx, p.y() + dy, p.theta());
    }

    const EnergyParts after = group_parts(sys, group);
    if (metropolis(energy_change(sys, before, after), temperature, rng)) {
        commit(sys, before, after);
        return true;
    }
    restore(sys, group, saved);  // restore the saved poses, not -dx: undoing a float add drifts
    return false;
}

bool try_slide_group(System& sys, const std::vector<int>& group, Vec2 dir,
                     const StepSizes& steps, double temperature, std::mt19937_64& rng) {
    if (!group_movable(sys, group)) return false;
    const double len = norm(dir);
    if (!(len > 0.0)) return false;
    const Vec2 u{dir.x / len, dir.y / len};

    const std::vector<Pose> saved = snapshot(sys, group);
    const EnergyParts before = group_parts(sys, group);

    std::uniform_real_distribution<double> d(-steps.translate, steps.translate);
    const double t = d(rng);
    for (const int i : group) {
        Particle& p = sys.particles[static_cast<size_t>(i)];
        sys.set_pose(i, p.x() + t * u.x, p.y() + t * u.y, p.theta());
    }

    const EnergyParts after = group_parts(sys, group);
    if (metropolis(energy_change(sys, before, after), temperature, rng)) {
        commit(sys, before, after);
        return true;
    }
    restore(sys, group, saved);
    return false;
}

bool try_rotate_group(System& sys, const std::vector<int>& group, Vec2 pivot,
                      const StepSizes& steps, double temperature, std::mt19937_64& rng) {
    if (!group_movable(sys, group)) return false;

    const std::vector<Pose> saved = snapshot(sys, group);
    const EnergyParts before = group_parts(sys, group);

    std::uniform_real_distribution<double> d(-steps.rotate, steps.rotate);
    const double dth = d(rng), ct = std::cos(dth), stt = std::sin(dth);
    for (size_t k = 0; k < group.size(); ++k) {
        const double ux = saved[k].x - pivot.x, uy = saved[k].y - pivot.y;
        sys.set_pose(group[k],
            pivot.x + ux * ct - uy * stt, pivot.y + ux * stt + uy * ct, saved[k].theta + dth);
    }

    const EnergyParts after = group_parts(sys, group);
    if (metropolis(energy_change(sys, before, after), temperature, rng)) {
        commit(sys, before, after);
        return true;
    }
    restore(sys, group, saved);
    return false;
}

bool try_rotate_group(System& sys, const std::vector<int>& group, const StepSizes& steps,
                      double temperature, std::mt19937_64& rng) {
    if (!group_movable(sys, group)) return false;

    // Turn the whole group about its own centroid, so the members keep their arrangement.
    double cx = 0.0, cy = 0.0;
    for (const int i : group) {
        const Particle& p = sys.particles[static_cast<size_t>(i)];
        cx += p.x();
        cy += p.y();
    }
    const double m = static_cast<double>(group.size());
    return try_rotate_group(sys, group, Vec2{cx / m, cy / m}, steps, temperature, rng);
}

bool try_reflect_group(System& sys, const std::vector<int>& group, Vec2 point, Vec2 axis,
                       double temperature, std::mt19937_64& rng) {
    if (!group_movable(sys, group)) return false;
    const double len = norm(axis);
    if (!(len > 0.0)) return false;
    const Vec2 a{axis.x / len, axis.y / len};

    const std::vector<Pose> saved = snapshot(sys, group);
    const EnergyParts before = group_parts(sys, group);

    // Mirror in the line through `point` along `a`: the offset u goes to 2 (u . a) a - u, and an orientation to 2 * angle(a) - theta. A reflection is an isometry, so every distance and  every angle inside the group survives it -- which is what lets the delta be measured on group_parts, intra-group pairs double counted and all, exactly as for a slide or a turn. The same cancellation carries the parts: these moves are isometries, so every intra-group overlap is the same area before and after, and the double-counted terms subtract out.
    const double twice_phi = 2.0 * std::atan2(a.y, a.x);
    for (size_t k = 0; k < group.size(); ++k) {
        const double ux = saved[k].x - point.x, uy = saved[k].y - point.y;
        const double twice_proj = 2.0 * (ux * a.x + uy * a.y);
        sys.set_pose(group[k], point.x + twice_proj * a.x - ux,
                                                          point.y + twice_proj * a.y - uy,
                                                          twice_phi - saved[k].theta);
    }

    const EnergyParts after = group_parts(sys, group);
    if (metropolis(energy_change(sys, before, after), temperature, rng)) {
        commit(sys, before, after);
        return true;
    }
    restore(sys, group, saved);
    return false;
}

bool try_shove_region(System& sys, const RegionWeights& regions, const StepSizes& steps,
                      double temperature, std::mt19937_64& rng) {
    const Pick p = pick_region(sys, regions, rng);
    if (p.members.size() < 2) return false;  // one square moved rigidly is just a worse translate

    switch (p.kind) {
        case Region::Band:
            // Along the band, or mirrored in its midline: both leave the band's membership
            // alone, and a sideways shove would not. See the note in mc.hpp.
            return ((rng() & 1u) != 0u)
                       ? try_reflect_group(sys, p.members, p.centre, p.axis, temperature, rng)
                       : try_slide_group(sys, p.members, p.dir, steps, temperature, rng);
        case Region::Disc:
            return ((rng() & 1u) != 0u)
                       ? try_reflect_group(sys, p.members, p.centre, p.axis, temperature, rng)
                       : try_rotate_group(sys, p.members, p.centre, steps, temperature, rng);
        case Region::Corner:
            // A corner block has one symmetry line, its diagonal, so there is nothing to choose.
            return try_reflect_group(sys, p.members, p.centre, p.axis, temperature, rng);
    }
    return false;
}

Stats anneal(System& sys, const Schedule& sched, StepSizes& steps, std::mt19937_64& rng,
             const std::function<void(long, const System&)>& watch) {
    Stats st;
    const long n = static_cast<long>(sys.particles.size());
    if (n == 0 || sched.sweeps <= 0) return st;

    // Geometric cooling, so each decade of temperature gets the same number of sweeps.
    const double ratio = (sched.sweeps > 1)
                             ? std::pow(sched.t_cold / sched.t_hot,
                                        1.0 / static_cast<double>(sched.sweeps - 1))
                             : 1.0;
    double T = sched.t_hot;

    for (long k = 0; k < sched.sweeps; ++k) {
        const Stats done = sweep(sys, steps, T, rng);
        st.attempted += done.attempted;
        st.accepted += done.accepted;

        if ((k & 1023) == 1023) sys.resync_energy();
        T *= ratio;
        if (watch) watch(k, sys);
    }

    sys.resync_energy();
    return st;
}

}  // namespace pack
