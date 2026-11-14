#include "moves/relax.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <thread>

namespace pack {
namespace {

// Each mobile particle has three variables: x, y and R * theta, where R is the particle's own
// circumradius. Multiplying theta by R gives all three units of length, so L-BFGS sees them on
// the same scale -- one R per particle, not one for the whole system, so a mix of species each
// get their own correctly-scaled angle variable.
struct Layout {
    std::vector<int> index;     // variable block -> particle index
    std::vector<double> r;      // R per block, parallel to `index`
};

// Lists the mobile (non-frozen) particles and sets each one's R.
Layout plan(const System& sys) {
    Layout lay;
    for (size_t i = 0; i < sys.particles.size(); ++i) {
        if (sys.particles[i].frozen) continue;
        lay.index.push_back(static_cast<int>(i));
        const double rc = sys.shape_of(static_cast<int>(i)).circumradius;
        lay.r.push_back(rc > 0.0 ? rc : 1.0);
    }
    return lay;
}

// Copies the positions and orientations of the mobile particles into u.
void gather(const System& sys, const Layout& lay, std::vector<double>& u) {
    u.resize(3 * lay.index.size());
    for (size_t b = 0; b < lay.index.size(); ++b) {
        const Particle& p = sys.particles[static_cast<size_t>(lay.index[b])];
        u[3 * b + 0] = p.x();
        u[3 * b + 1] = p.y();
        u[3 * b + 2] = lay.r[b] * p.theta();
    }
}

// Writes u back into the poses of the mobile particles.
void scatter(System& sys, const Layout& lay, const std::vector<double>& u) {
    for (size_t b = 0; b < lay.index.size(); ++b) {
        sys.set_pose(lay.index[b], u[3 * b + 0], u[3 * b + 1], u[3 * b + 2] / lay.r[b]);
    }
}

// Computes the energy and its gradient in the scaled variables.
double gradient_at(const System& sys, const Layout& lay, std::vector<Pose>& scratch,
                   std::vector<double>& g) {
    const double e = energy_grad(sys, scratch);
    g.resize(3 * lay.index.size());
    for (size_t b = 0; b < lay.index.size(); ++b) {
        const Pose& pg = scratch[static_cast<size_t>(lay.index[b])];
        g[3 * b + 0] = pg.x;
        g[3 * b + 1] = pg.y;
        // Angle component divided by R because the variable is R * theta.
        g[3 * b + 2] = pg.theta / lay.r[b];
    }
    return e;
}

double dot_of(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

double max_abs(const std::vector<double>& a) {
    double m = 0.0;
    for (const double v : a) m = std::max(m, std::abs(v));
    return m;
}

// Steepest descent with length `trust`. Used when there is no L-BFGS history: on the first iteration and after a reset. Don't scale by min(1, 1/|g|) instead. For |g| < 1 that gives a step of length |g|, and near a valid packing |g| is about 2 * k * overlap, which is tiny. The relaxation then takes thousands of iterations instead of tens. A step that is too long is cheaper: the line search halves it until it is accepted.
void descent_direction(const std::vector<double>& g, double trust, std::vector<double>& d) {
    const double n2 = std::sqrt(dot_of(g, g));
    const double scale = (n2 > 0.0) ? trust / n2 : 0.0;
    d.resize(g.size());
    for (size_t i = 0; i < g.size(); ++i) d[i] = -g[i] * scale;
}

// The largest energy a valid packing can have: every one of the n wall terms and n(n-1)/2 pair
// terms at tol, each weighed by the stiffer of the two constants. Above it, some overlap is
// certainly past tol.
double feasible_energy_ceiling(const System& sys, double tol) {
    const double n = static_cast<double>(sys.particles.size());
    return std::max(sys.k_pair, sys.k_wall) * 0.5 * n * (n + 1.0) * tol * tol;
}

// Whether every overlap is within tol. The energy settles it at either end: below k * tol^2
// no single term can exceed tol, and above the ceiling some term must. In between it asks the
// overlaps themselves. That pass is O(n) with the cell list, and only runs near convergence.
// Before it existed the low bound alone decided, so a relaxation already valid kept iterating
// until every contact sat sqrt(#contacts) below tol -- dozens of wasted iterations per success.
bool proven_feasible(const System& sys, double e, double tol) {
    const double k = std::min(sys.k_pair, sys.k_wall);
    if (k > 0.0 && std::sqrt(e / k) <= tol) return true;
    if (e > feasible_energy_ceiling(sys, tol)) return false;
    return sys.max_overlap() <= tol;
}

// Returns the index of the first square lying entirely outside the box, or -1.
int escaped_index(const System& sys) {
    for (size_t i = 0; i < sys.particles.size(); ++i) {
        const Shape& s = sys.shape_of(static_cast<int>(i));
        const double out = outside_area(sys.particles[i], s, sys.L);
        if (out >= s.area * (1.0 - 1e-12)) return static_cast<int>(i);
    }
    return -1;
}

}  // namespace

// Total energy: k_wall * (area outside the box)^2 for each particle, plus k_pair * (overlap area)^2 for each pair. Writes the gradient into grad. Frozen particles get a zero gradient, but their terms are still counted in the energy.
double energy_grad(const System& sys, std::vector<Pose>& grad) {
    const size_t n = sys.particles.size();
    grad.assign(n, Pose{});
    double e = 0.0;

    for (size_t i = 0; i < n; ++i) {
        const bool fi = sys.particles[i].frozen;
        Pose gw;
        const double out =
            outside_area_grad(sys.particles[i], sys.shape_of(static_cast<int>(i)), sys.L,
                              fi ? nullptr : &gw);
        e += sys.k_wall * out * out;
        if (!fi) {
            // d(k * area^2) = 2 * k * area * d(area). The pair terms below use the same rule.
            const double f = 2.0 * sys.k_wall * out;
            grad[i].x += f * gw.x;
            grad[i].y += f * gw.y;
            grad[i].theta += f * gw.theta;
        }
    }

    // The close pairs only, from a grid built here: L-BFGS has just written every pose, and
    // the pass costs O(n) instead of n^2 / 2 clipper calls.
    CellList grid;
    grid.build(sys.particles, sys.L, sys.reach());
    grid.for_each_close_pair([&](int ii, int jj) {
        const size_t i = static_cast<size_t>(ii), j = static_cast<size_t>(jj);
        const bool fi = sys.particles[i].frozen, fj = sys.particles[j].frozen;
        Pose gi, gj;
        const double a = overlap_area_grad(sys.particles[i], sys.shape_of(ii), sys.particles[j],
                                           sys.shape_of(jj), fi ? nullptr : &gi,
                                           fj ? nullptr : &gj);
        e += sys.k_pair * a * a;
        const double f = 2.0 * sys.k_pair * a;
        if (!fi) {
            grad[i].x += f * gi.x;
            grad[i].y += f * gi.y;
            grad[i].theta += f * gi.theta;
        }
        if (!fj) {
            grad[j].x += f * gj.x;
            grad[j].y += f * gj.y;
            grad[j].theta += f * gj.theta;
        }
    });
    return e;
}

const char* stop_name(RelaxResult::Stop stop) {
    switch (stop) {
        case RelaxResult::Stop::Feasible: return "feasible";
        case RelaxResult::Stop::Converged: return "converged (jammed)";
        case RelaxResult::Stop::Stalled: return "stalled";
        case RelaxResult::Stop::MaxIterations: return "max iterations";
        case RelaxResult::Stop::Hopeless: return "hopeless (jammed, out of budget)";
        case RelaxResult::Stop::Escaped: return "a square escaped the box";
    }
    return "?";
}

// Moves the mobile particles to minimize the energy, using L-BFGS with a backtracking line search. Leaves sys at the final configuration.
RelaxResult relax(System& sys, const RelaxOptions& opt) {
    RelaxResult res;
    const Layout lay = plan(sys);

    // History size. By default, one pair per variable (3 per mobile particle), clamped to [10, 100]. With that many pairs L-BFGS is close to full BFGS. See RelaxOptions::history.
    const size_t m = (opt.history > 0)
                         ? static_cast<size_t>(opt.history)
                         : static_cast<size_t>(std::clamp<long>(
                               3 * static_cast<long>(lay.index.size()), 10, 100));

    std::vector<Pose> scratch;
    std::vector<double> u, g, d, u_try;
    gather(sys, lay, u);

    double e = gradient_at(sys, lay, scratch, g);
    ++res.evaluations;

    // Fills in and returns the result. Whatever ended the run, a square outside the box is reported as Escaped.
    const auto finish = [&](RelaxResult::Stop stop) {
        // L-BFGS drives the poses directly, so the two tracked sums have been stale throughout
        // and this is where they become true again.
        sys.parts = sys.energy_parts();
        if (escaped_index(sys) >= 0) stop = RelaxResult::Stop::Escaped;
        res.stop = stop;
        res.energy = sys.energy();
        res.max_overlap = sys.max_overlap();
        res.feasible = res.max_overlap <= opt.feasible_tol;
        return res;
    };

    // Early exits
    if (proven_feasible(sys, e, opt.feasible_tol)) return finish(RelaxResult::Stop::Feasible);
    if (lay.index.empty()) return finish(RelaxResult::Stop::Converged);
    if (escaped_index(sys) >= 0) return finish(RelaxResult::Stop::Escaped);

    // Length of a steepest-descent step, in scaled variables: one length unit (a square's width), at most half the box.
    const double trust = std::min(1.0, 0.5 * sys.L);

    // L-BFGS history, oldest first: s = step taken, y = change in gradient, rho = 1 / (s . y).
    std::vector<std::vector<double>> hist_s, hist_y;
    std::vector<double> hist_rho;
    std::vector<double> alpha(m, 0.0);
    int stalls = 0;
    double e_window = e;  // energy at the start of the current hopeless_window

    for (res.iterations = 0; res.iterations < opt.max_iterations; ++res.iterations) {
        if (max_abs(g) <= opt.grad_tol) return finish(RelaxResult::Stop::Converged);

        // Two-loop recursion: computes d = -H * g, where H is the inverse Hessian estimate built from the stored pairs.
        d.assign(g.size(), 0.0);
        for (size_t i = 0; i < g.size(); ++i) d[i] = -g[i];
        const size_t k = hist_s.size();
        for (size_t t = k; t-- > 0;) {
            alpha[t] = hist_rho[t] * dot_of(hist_s[t], d);
            for (size_t i = 0; i < d.size(); ++i) d[i] -= alpha[t] * hist_y[t][i];
        }
        if (k > 0) {
            // H0 = gamma * I, with gamma = (s . y) / (y . y) from the newest pair. This is the standard choice and gives the step a reasonable length.
            const double gamma =
                dot_of(hist_s[k - 1], hist_y[k - 1]) / dot_of(hist_y[k - 1], hist_y[k - 1]);
            for (double& v : d) v *= gamma;
        } else {
            // No history (first iteration or after a reset): use steepest descent.
            descent_direction(g, trust, d);
        }
        for (size_t t = 0; t < k; ++t) {
            const double beta = hist_rho[t] * dot_of(hist_y[t], d);
            for (size_t i = 0; i < d.size(); ++i) d[i] += (alpha[t] - beta) * hist_s[t][i];
        }

        double slope = dot_of(g, d);
        if (!(slope < 0.0)) {
            // d points uphill. This happens when the contact set changes and the stored pairs no longer fit the energy. Clear the history and use steepest descent. If that is not downhill either, g is zero or NaN, and we stop.
            hist_s.clear();
            hist_y.clear();
            hist_rho.clear();
            descent_direction(g, trust, d);
            slope = dot_of(g, d);
            if (!(slope < 0.0)) return finish(RelaxResult::Stop::Converged);
        }

        // Backtracking line search: halve the step until the Armijo condition (sufficient decrease) holds, at most 60 times. There is no Wolfe curvature test: it would need a gradient at every trial, and the energy's second derivative jumps where contacts change. Without that test s . y > 0 is not guaranteed, so bad pairs are skipped further down.
        double step = 1.0, e_new = e;
        bool ok = false;
        for (int trial = 0; trial < 60; ++trial) {
            u_try.resize(u.size());
            for (size_t i = 0; i < u.size(); ++i) u_try[i] = u[i] + step * d[i];
            scatter(sys, lay, u_try);
            e_new = sys.total_energy();
            ++res.evaluations;
            if (e_new <= e + 1e-4 * step * slope) {
                ok = true;
                break;
            }
            step *= 0.5;
        }
        if (!ok) {
            scatter(sys, lay, u);  // restore the last accepted configuration
            return finish(RelaxResult::Stop::Stalled);
        }

        // New L-BFGS pair: s = step taken, y = change in gradient.
        std::vector<double> s_new(u.size()), y_new(u.size());
        for (size_t i = 0; i < u.size(); ++i) s_new[i] = u_try[i] - u[i];
        u.swap(u_try);

        const std::vector<double> g_old = g;
        const double e_old = e;
        e = gradient_at(sys, lay, scratch, g);
        ++res.evaluations;
        for (size_t i = 0; i < g.size(); ++i) y_new[i] = g[i] - g_old[i];

        if (opt.watch && (res.iterations + 1) % std::max<long>(1, opt.watch_every) == 0) {
            opt.watch(sys, opt.watch_label);
        }

        if (proven_feasible(sys, e, opt.feasible_tol)) {
            ++res.iterations;
            return finish(RelaxResult::Stop::Feasible);
        }
        // Stop if the relative energy decrease stays below energy_tol for stall_patience
        // iterations in a row. Relative because near a valid packing E is around 1e-8, where a
        // tiny absolute decrease is still real progress. Several in a row so that one poor step
        // does not end the run.
        if (e_old - e <= opt.energy_tol * std::max(e_old, 1e-300)) {
            if (++stalls >= std::max(1, opt.stall_patience)) {
                ++res.iterations;
                return finish(RelaxResult::Stop::Stalled);
            }
        } else {
            stalls = 0;
        }

        // Give up on a relaxation that cannot reach a valid packing within its budget. Every
        // window, the log-energy rate over it says how many more iterations getting down to
        // the feasible ceiling would take at the current pace; past ten times what is left,
        // stop. A jammed box shows the signature from its first few hundred iterations --
        // decay of about 1% per 600 iterations with 40 e-folds still to go -- and compress()
        // used to pay the whole max_iterations for each of them, which was nearly all of its
        // time. A relaxation that is going to succeed converges geometrically and needs a
        // small fraction of its budget, so the factor of ten never comes close to it.
        if (opt.hopeless_window > 0 && (res.iterations + 1) % opt.hopeless_window == 0) {
            const double need = std::log(e / feasible_energy_ceiling(sys, opt.feasible_tol));
            const double rate = std::log(e_window / e) / static_cast<double>(opt.hopeless_window);
            const double left = static_cast<double>(opt.max_iterations - (res.iterations + 1));
            if (need > 0.0 && !(need <= 10.0 * rate * left)) {  // also catches rate <= 0, NaN
                ++res.iterations;
                return finish(RelaxResult::Stop::Hopeless);
            }
            e_window = e;
        }

        // Store the pair only if s . y > 0. Without the Wolfe test a pair can have s . y <= 0,
        // and storing it would make H indefinite. When the history is full, drop the oldest.
        const double sy = dot_of(s_new, y_new);
        if (sy > 1e-30) {
            if (hist_s.size() == m) {
                hist_s.erase(hist_s.begin());
                hist_y.erase(hist_y.begin());
                hist_rho.erase(hist_rho.begin());
            }
            hist_s.push_back(std::move(s_new));
            hist_y.push_back(std::move(y_new));
            hist_rho.push_back(1.0 / sy);
        }
    }
    return finish(RelaxResult::Stop::MaxIterations);
}

// Finds a small box for this configuration. Phase one grows the box until relax() succeeds.
// Phase two shrinks it as far as relax() still succeeds.
CompressResult compress(const System& start, const CompressOptions& opt) {
    CompressResult res;
    System sys = start;
    assert(opt.shrink > 0.0 && opt.shrink < 1.0 && opt.grow > 1.0);

    // Phase one: grow the box by opt.grow until relax() succeeds, at most max_growth_stages
    // times. The input usually comes from a fixed-pressure run and still overlaps, so its own
    // box is too small. Growing too much is fine, since phase two shrinks the box again.
    {
        RelaxResult r = relax(sys, opt.relax);
        ++res.relaxations;
        for (long t = 0; !r.feasible && t < opt.max_growth_stages; ++t) {
            sys.rescale_box(sys.L * opt.grow);
            r = relax(sys, opt.relax);
            ++res.relaxations;
        }
        if (!r.feasible) return res;
        res.found = true;
        res.packing = sys;
        res.side = sys.tight_side();
    }

    // Phase two: shrink the box by (1 - step) and relax. On success, keep the new packing. On
    // failure, halve the step and try again from the last good packing. The step never grows.
    // Stops when step <= min_shrink or after max_relaxations.
    double step = 1.0 - opt.shrink;
    res.final_step = step;
    // The side is only known to sqrt(feasible_tol) as a length -- a corner contact hides the
    // square root of its area as depth -- so a step that moves it by less than that cannot
    // produce a smaller side anyone could claim. Squeezing below it only spends relaxations.
    const double floor = std::max(opt.min_shrink, std::sqrt(opt.relax.feasible_tol) / res.side);
    for (res.stages = 0; step > floor;) {
        if (res.relaxations >= opt.max_relaxations) break;

        System trial = res.packing;
        // Shrink from the space the packing occupies (tight_side), not from the box. The
        // packing often has slack in its box, and starting from tight_side removes it.
        trial.rescale_box(std::min(trial.L, res.side) * (1.0 - step));

        const RelaxResult r = relax(trial, opt.relax);
        ++res.relaxations;
        if (r.feasible) {
            res.packing = trial;
            res.side = trial.tight_side();
            ++res.stages;
        } else {
            step *= 0.5;
            res.final_step = step;
        }
    }
    // True if we stopped because the step got too small, not because of the budget.
    res.squeezed_out = step <= floor;
    return res;
}

BestOfResult compress_best_of(const std::vector<System>& candidates, const CompressOptions& opt,
                              int threads) {
    BestOfResult out;
    const int n = static_cast<int>(candidates.size());
    out.sides.assign(candidates.size(), 0.0);
    if (n == 0) return out;

    std::vector<CompressResult> each(candidates.size());
    if (threads <= 0) threads = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    threads = std::min(threads, n);
    // Static striding rather than a work queue: the answer is fixed by the candidates alone,
    // and the few of them do not need balancing finer than this.
    const auto work = [&](int t) {
        for (int i = t; i < n; i += threads) {
            each[static_cast<size_t>(i)] = compress(candidates[static_cast<size_t>(i)], opt);
        }
    };
    if (threads == 1) {
        work(0);
    } else {
        std::vector<std::thread> pool;
        for (int t = 0; t < threads; ++t) pool.emplace_back(work, t);
        for (std::thread& th : pool) th.join();
    }

    for (int i = 0; i < n; ++i) {
        const CompressResult& c = each[static_cast<size_t>(i)];
        if (!c.found) continue;
        out.sides[static_cast<size_t>(i)] = c.side;
        if (out.winner < 0 || c.side < out.best.side) {
            out.best = c;
            out.winner = i;
        }
    }
    return out;
}

void relax_free(System& sys, const RelaxOptions& opt) {
    const double k_wall = sys.k_wall;
    sys.k_wall = 0.0;
    relax(sys, opt);
    sys.k_wall = k_wall;
    sys.fit_box();
}

CompressResult anneal(System& sys, double ds, long max_iter, const FrameWatch& watch) {
    assert(ds > 0.0 && ds < sys.L);
    CompressOptions opt;
    opt.shrink = 1.0 - ds / sys.L;
    opt.max_relaxations = max_iter;
    opt.relax.watch = watch;
    opt.relax.watch_label = "compress";
    CompressResult c = compress(sys, opt);
    if (c.found) sys = c.packing;
    return c;
}

}  // namespace pack