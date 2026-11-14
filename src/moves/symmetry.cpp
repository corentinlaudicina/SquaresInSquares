
// minimiser that exploits symmetries by constructing reflected / mirrored images etc.

// 4 fold symmetry
// We have squares in the first quadrant at (x, y, θ), and they automatically generate the 3 other ones : (−y, x, θ), (−x, −y, θ), (y, −x, θ). This is nice because we only need to simulate 1/4 of the box.


// mirror symmetries [mid and diagonals]

#include "moves/symmetry.hpp"

#include "io/io.hpp"
#include "moves/lnm.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <limits>

namespace pack {
namespace {

constexpr double kHalfPi = 0.5 * kPi;

// Positions are binned at kPositionStep. The angle step is a whole fraction of the 90 degree period so that the period is an exact number of bins: without that, +45 and -45 degrees -- one orientation, and the commonest tilt in a good packing -- could land one bin apart.
constexpr double kPositionStep = 1e-7;
constexpr long long kAnglePeriod = 1LL << 24;
constexpr double kAngleStep = kHalfPi / static_cast<double>(kAnglePeriod);

long long position_bin(double v) { return std::llround(v / kPositionStep); }

// Wrapped key; `moved` reports that the bin was taken from the +45 end to the -45 end.
long long angle_bin(double theta, bool* moved = nullptr) {
    long long k = std::llround(wrap90(theta) / kAngleStep);
    const bool wrap = 2 * k >= kAnglePeriod;
    if (wrap) k -= kAnglePeriod;
    if (moved) *moved = wrap;
    return k;
}

Pose pose_of(const Particle& p) { return {p.x(), p.y(), p.theta()}; }

bool maps_onto_itself(const System& sys, const GroupElem& g, double tol) {
    const size_t n = sys.particles.size();
    std::vector<char> used(n, 0);
    for (const Particle& p : sys.particles) {
        const Pose q = g.apply(pose_of(p), sys.L);
        bool found = false;
        for (size_t j = 0; j < n && !found; ++j) {
            if (used[j]) continue;
            const Particle& o = sys.particles[j];
            if (std::hypot(q.x - o.x(), q.y - o.y()) > tol) continue;
            if (angle_distance90(q.theta, o.theta()) > tol) continue;
            used[j] = 1;
            found = true;
        }
        if (!found) return false;
    }
    return true;
}

bool closed_under_composition(unsigned mask) {
    const auto& all = d4_elements();
    for (int i = 0; i < 8; ++i) {
        if (!(mask >> i & 1u)) continue;
        for (int j = 0; j < 8; ++j) {
            if (!(mask >> j & 1u)) continue;
            const GroupElem product = all[static_cast<size_t>(i)].compose(all[static_cast<size_t>(j)]);
            const auto it = std::find(all.begin(), all.end(), product);
            if (!(mask >> (it - all.begin()) & 1u)) return false;
        }
    }
    return true;
}

// One image's bins in sorted order; `order` gets the particle index behind each entry.
CanonicalKey image_key(const System& sys, const GroupElem& g, std::vector<int>* order = nullptr) {
    struct Entry {
        std::array<long long, 3> key;
        int index;
    };
    std::vector<Entry> entries;
    entries.reserve(sys.particles.size());
    for (size_t i = 0; i < sys.particles.size(); ++i) {
        const Pose q = g.apply(pose_of(sys.particles[i]), sys.L);
        entries.push_back({{position_bin(q.x), position_bin(q.y), angle_bin(q.theta)},
                           static_cast<int>(i)});
    }
    std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
        return (a.key != b.key) ? a.key < b.key : a.index < b.index;
    });
    CanonicalKey out;
    out.reserve(entries.size());
    if (order) order->clear();
    for (const Entry& e : entries) {
        out.push_back(e.key);
        if (order) order->push_back(e.index);
    }
    return out;
}

}  // namespace

Pose GroupElem::apply(const Pose& p, double L) const {
    const double h = 0.5 * L, u = p.x - h, v = p.y - h;
    return {h + a * u + b * v, h + c * u + d * v, angle_sign * p.theta};
}

// A signed permutation matrix is orthogonal ==> its inverse is its transpose.
GroupElem GroupElem::inverse() const { return {a, c, b, d, angle_sign}; }

GroupElem GroupElem::compose(const GroupElem& g) const {
    return {a * g.a + b * g.c, a * g.b + b * g.d, c * g.a + d * g.c, c * g.b + d * g.d, angle_sign * g.angle_sign};
}

GroupElem group_elem(D4 which) {
    return d4_elements()[static_cast<size_t>(which)];
}

const std::array<GroupElem, 8>& d4_elements() {
    static const std::array<GroupElem, 8> all = {{
        {1, 0, 0, 1, 1},     // Id
        {0, -1, 1, 0, 1},    // Rot90: (x, y) -> (-y, x)
        {-1, 0, 0, -1, 1},   // Rot180
        {0, 1, -1, 0, 1},    // Rot270
        {-1, 0, 0, 1, -1},   // MirrorX
        {1, 0, 0, -1, -1},   // MirrorY
        {0, 1, 1, 0, -1},    // MirrorDiag
        {0, -1, -1, 0, -1},  // MirrorAnti
    }};
    return all;
}

double wrap90(double theta) {
    return theta - kHalfPi * std::floor(theta / kHalfPi + 0.5);
}

double angle_distance90(double a, double b) { return std::abs(wrap90(a - b)); }

bool is_square_config(const System& sys) {
    if (sys.species.empty()) return false;
    const Shape& s = sys.base_shape();
    if (s.vertices.size() != 4 || s.sat_axes.size() != 2) return false;
    if (std::abs(dot(s.sat_axes[0], s.sat_axes[1])) > 1e-12) return false;
    const auto width = [&](Vec2 axis) {
        double lo = dot(s.vertices[0], axis), hi = lo;
        for (const Vec2 v : s.vertices) {
            lo = std::min(lo, dot(v, axis));
            hi = std::max(hi, dot(v, axis));
        }
        return hi - lo;
    };
    if (std::abs(width(s.sat_axes[0]) - width(s.sat_axes[1])) > 1e-12) return false;
    for (size_t i = 0; i < sys.particles.size(); ++i) {
        if (sys.particles[i].species != 0 || sys.is_merged(static_cast<int>(i))) return false;
    }
    return true;
}

System image_of(const System& sys, const GroupElem& g) {
    System out = sys;
    for (Particle& p : out.particles) {
        const Pose q = g.apply(pose_of(p), sys.L);
        p.set_pose(q.x, q.y, q.theta);
    }
    out.recount();
    return out;
}

std::vector<GroupElem> detect_symmetry(const System& sys, double tol) {
    const auto& all = d4_elements();
    unsigned fit = 1u;  // the identity always fits
    if (is_square_config(sys)) {
        for (size_t i = 1; i < all.size(); ++i) {
            if (maps_onto_itself(sys, all[i], tol)) fit |= 1u << i;
        }
    }

    // Eight elements make 256 subsets; the ties that can occur are broken by the smallest mask.
    unsigned best = 1u;
    for (unsigned mask = 1u; mask < 256u; mask += 2u) {  // odd: the identity is always in
        if ((mask & fit) != mask || !closed_under_composition(mask)) continue;
        if (std::popcount(mask) > std::popcount(best)) best = mask;
    }

    std::vector<GroupElem> out;
    for (size_t i = 0; i < all.size(); ++i) {
        if (best >> i & 1u) out.push_back(all[i]);
    }
    return out;
}

bool contains(const std::vector<GroupElem>& group, const GroupElem& g) {
    return std::find(group.begin(), group.end(), g) != group.end();
}

CanonicalKey canonical_key(const System& sys) {
    if (!is_square_config(sys)) return image_key(sys, group_elem(D4::Id));
    CanonicalKey best;
    bool have = false;
    for (const GroupElem& g : d4_elements()) {
        CanonicalKey key = image_key(sys, g);
        if (!have || key < best) {
            best = std::move(key);
            have = true;
        }
    }
    return best;
}

System canonicalize(const System& sys) {
    if (!is_square_config(sys)) return sys;

    const auto& all = d4_elements();
    CanonicalKey best;
    size_t winner = 0;
    for (size_t k = 0; k < all.size(); ++k) {
        CanonicalKey key = image_key(sys, all[k]);
        if (k == 0 || key < best) {
            best = std::move(key);
            winner = k;
        }
    }

    std::vector<int> order;
    image_key(sys, all[winner], &order);
    System out = sys;
    for (size_t k = 0; k < order.size(); ++k) {
        const Particle& from = sys.particles[static_cast<size_t>(order[k])];
        const Pose q = all[winner].apply(pose_of(from), sys.L);
        bool moved = false;
        angle_bin(q.theta, &moved);
        Particle& to = out.particles[k];
        to.frozen = from.frozen;
        to.set_pose(q.x, q.y, wrap90(q.theta) - (moved ? kHalfPi : 0.0));
    }
    out.recount();
    return out;
}

// --- symmetrize_region ----------------------------------------------------------------------

namespace {

// A centre on a boundary belongs to every region that touches it: closed regions, so that a
// square exactly on a mirror line is part of the source rather than of nothing.
constexpr double kBoundaryEps = 1e-6;

bool in_region(const BoxRegion& r, double x, double y, double L) {
    const double dx = x - 0.5 * L, dy = y - 0.5 * L, e = kBoundaryEps;
    switch (r.kind) {
        case RegionKind::Quadrant:  // counter-clockwise from the origin corner
            switch (r.index) {
                case 0: return dx <= e && dy <= e;
                case 1: return dx >= -e && dy <= e;
                case 2: return dx >= -e && dy >= -e;
                default: return dx <= e && dy >= -e;
            }
        case RegionKind::Half:
            switch (r.index) {
                case 0: return dx <= e;
                case 1: return dx >= -e;
                case 2: return dy <= e;
                default: return dy >= -e;
            }
        case RegionKind::DiagonalHalf:
            switch (r.index) {
                case 0: return dy - dx >= -e;
                case 1: return dy - dx <= e;
                case 2: return dy + dx >= -e;
                default: return dy + dx <= e;
            }
        case RegionKind::Sector:  // bottom, right, top, left
            switch (r.index) {
                case 0: return dy <= -std::abs(dx) + e;
                case 1: return dx >= std::abs(dy) - e;
                case 2: return dy >= std::abs(dx) - e;
                default: return dx <= -std::abs(dy) + e;
            }
    }
    return false;
}

double region_area(const BoxRegion& r, double L) {
    const bool quarter = r.kind == RegionKind::Quadrant || r.kind == RegionKind::Sector;
    return (quarter ? 0.25 : 0.5) * L * L;
}

// One way to apply a variant: a source and the elements that carry it onto its destinations.
struct Placement {
    BoxRegion source;
    std::vector<GroupElem> images;
};

std::vector<Placement> placements_of(int variant) {
    const int copies = variant / 3 + 1;
    const SymFamily family = static_cast<SymFamily>(variant % 3);
    const auto el = [](D4 w) { return group_elem(w); };
    const GroupElem r90 = el(D4::Rot90), r180 = el(D4::Rot180), r270 = el(D4::Rot270);
    const GroupElem mx = el(D4::MirrorX), my = el(D4::MirrorY);
    const GroupElem md = el(D4::MirrorDiag), ma = el(D4::MirrorAnti);

    std::vector<Placement> out;
    for (int i = 0; i < 4; ++i) {
        const BoxRegion quadrant{RegionKind::Quadrant, i};
        const BoxRegion sector{RegionKind::Sector, i};
        const BoxRegion half{RegionKind::Half, i};
        const BoxRegion diag{RegionKind::DiagonalHalf, i};
        if (copies == 1) {
            if (family == SymFamily::Rotation) {
                for (const GroupElem& g : {r90, r180, r270}) out.push_back({quadrant, {g}});
            } else if (family == SymFamily::AxisMirror) {
                out.push_back({quadrant, {mx}});
                out.push_back({quadrant, {my}});
            } else {
                out.push_back({sector, {md}});
                out.push_back({sector, {ma}});
            }
        } else if (copies == 2) {
            if (family == SymFamily::Rotation) {
                out.push_back({half, {r180}});
                out.push_back({diag, {r180}});
            } else if (family == SymFamily::AxisMirror) {
                out.push_back({half, {i < 2 ? mx : my}});
            } else {
                out.push_back({diag, {i < 2 ? md : ma}});
            }
        } else if (family == SymFamily::Rotation) {
            out.push_back({quadrant, {r90, r180, r270}});
        } else if (family == SymFamily::AxisMirror) {
            out.push_back({quadrant, {mx, my, r180}});
        } else {
            out.push_back({sector, {md, ma, r180}});
        }
    }
    return out;
}

// Waste per square of a region: free area in units of one square, over the squares inside.
double waste(double area, double unit, int count) {
    if (count == 0) return std::numeric_limits<double>::infinity();
    return (area / unit - count) / count;
}

Placement choose_placement(const System& sys, int variant, double eps, std::mt19937_64& rng) {
    const std::vector<Placement> all = placements_of(variant);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    if (u(rng) < eps) return all[rng() % all.size()];

    const double unit = sys.base_shape().area;
    size_t winner = 0;
    double best = std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < all.size(); ++k) {
        const Placement& pl = all[k];
        const double area = region_area(pl.source, sys.L);
        int inside = 0;
        std::vector<int> dest(pl.images.size(), 0);
        for (const Particle& p : sys.particles) {
            if (in_region(pl.source, p.x(), p.y(), sys.L)) ++inside;
            for (size_t g = 0; g < pl.images.size(); ++g) {
                const Pose back = pl.images[g].inverse().apply({p.x(), p.y(), 0.0}, sys.L);
                if (in_region(pl.source, back.x, back.y, sys.L)) ++dest[g];
            }
        }
        if (inside == 0) continue;  // nothing to copy from here
        double dest_waste = 0.0;
        for (const int d : dest) dest_waste += waste(area, unit, d);
        const double score = waste(area, unit, inside) - dest_waste / dest.size();
        if (score < best) {  // lowest waste at the source, highest at the destination
            best = score;
            winner = k;
        }
    }
    return all[winner];
}

// Area of everything a would-be image overlaps, the walls included.
double image_overlap(const System& sys, const std::vector<char>& leaving, const Pose& pose) {
    const Shape& shape = sys.base_shape();
    Particle candidate;
    candidate.set_pose(pose.x, pose.y, pose.theta);
    const double reach = 2.0 * shape.circumradius;
    double total = outside_area(candidate, shape, sys.L);
    for (size_t j = 0; j < sys.particles.size(); ++j) {
        const Particle& other = sys.particles[j];
        if (leaving[j] || std::hypot(other.x() - pose.x, other.y() - pose.y) > reach) continue;
        total += overlap_area(candidate, shape, other, shape);
    }
    return total;
}

// Writes the images into the destination's particles; `changed` marks every square it wrote.
// Recycling those particles means the array never changes size, so the count cannot drift.
bool apply_placement(System& t, const Placement& pl, const LnmOptions& opt, std::mt19937_64& rng,
                     std::vector<char>& changed) {
    const double L = t.L;
    const size_t n = t.particles.size();

    std::vector<int> source;
    std::vector<int> dest;
    std::vector<char> leaving(n, 0);
    for (size_t i = 0; i < n; ++i) {
        const Particle& p = t.particles[i];
        if (in_region(pl.source, p.x(), p.y(), L)) {
            source.push_back(static_cast<int>(i));
            continue;
        }
        if (p.frozen) continue;  // a frozen square stays, and the images have to work around it
        for (const GroupElem& g : pl.images) {
            const Pose back = g.inverse().apply({p.x(), p.y(), 0.0}, L);
            if (in_region(pl.source, back.x, back.y, L)) {
                dest.push_back(static_cast<int>(i));
                leaving[i] = 1;
                break;
            }
        }
    }
    if (source.empty() || dest.empty()) return false;

    // An image that lands back in the source is on a fixed line of the element, and the original
    // is already there. Two images can coincide the same way.
    std::vector<Pose> images;
    for (const int i : source) {
        const Particle& p = t.particles[static_cast<size_t>(i)];
        for (const GroupElem& g : pl.images) {
            const Pose q = g.apply({p.x(), p.y(), p.theta()}, L);
            if (in_region(pl.source, q.x, q.y, L)) continue;
            const bool repeated = std::any_of(images.begin(), images.end(), [&](const Pose& o) {
                return std::hypot(o.x - q.x, o.y - q.y) < kBoundaryEps;
            });
            if (!repeated) images.push_back(q);
        }
    }

    // Too many images: lose the ones that collide most with what is staying.
    if (images.size() > dest.size()) {
        std::vector<std::pair<double, size_t>> scored;
        scored.reserve(images.size());
        for (size_t k = 0; k < images.size(); ++k) {
            scored.push_back({image_overlap(t, leaving, images[k]), k});
        }
        std::sort(scored.begin(), scored.end());
        std::vector<Pose> kept;
        kept.reserve(dest.size());
        for (size_t k = 0; k < dest.size(); ++k) kept.push_back(images[scored[k].second]);
        images = std::move(kept);
    }

    const Bounds hole = member_bounds(t, dest, 0.5 * contact_pitch(t.base_shape()));
    for (size_t k = 0; k < images.size(); ++k) {
        t.set_pose(dest[k], images[k].x, images[k].y, images[k].theta);
    }
    for (const int i : dest) changed[static_cast<size_t>(i)] = 1;

    // Too few images: the leftovers go back into the destination, one at a time, wherever the
    // reinsertion finds the most room.
    const std::vector<int> extras(dest.begin() + static_cast<long>(images.size()), dest.end());
    park_squares(t, extras);
    for (const int i : extras) reinsert_in_bounds(t, i, opt, hole, rng);

    assert(t.particles.size() == n);
    return true;
}

// The window the repair may touch: squares with a neighbour across the changed/unchanged line.
std::vector<char> seam_window(const System& t, const std::vector<char>& changed, double width) {
    const size_t n = t.particles.size();
    std::vector<char> window(n, 0);
    for (size_t i = 0; i < n; ++i) {
        if (!changed[i]) continue;
        const Particle& a = t.particles[i];
        for (size_t j = 0; j < n; ++j) {
            if (changed[j]) continue;
            const Particle& b = t.particles[j];
            if (std::hypot(a.x() - b.x(), a.y() - b.y()) <= width) {
                window[i] = 1;
                window[j] = 1;
            }
        }
    }
    return window;
}

// relax() honours frozen flags, so a window relaxation is a relaxation with the rest frozen.
void relax_window(System& t, const std::vector<char>& window, long iterations,
                  const FrameWatch& watch) {
    const size_t n = t.particles.size();
    std::vector<char> was(n, 0);
    for (size_t i = 0; i < n; ++i) {
        was[i] = t.particles[i].frozen ? 1 : 0;
        t.particles[i].frozen = t.particles[i].frozen || !window[i];
    }
    RelaxOptions ro;
    ro.max_iterations = std::max<long>(1, iterations);
    ro.feasible_tol = kFeasibleTol;  // a repair, not a certification
    ro.watch = watch;
    ro.watch_label = "seam relax";
    relax(t, ro);
    for (size_t i = 0; i < n; ++i) t.particles[i].frozen = was[i] != 0;
}

}  // namespace

int symmetrize_variant(int copies, SymFamily family) {
    assert(copies >= 1 && copies <= 3);
    return (copies - 1) * 3 + static_cast<int>(family);
}

std::string symmetrize_variant_name(int variant) {
    static const char* const family[3] = {"rotation", "axis", "diagonal"};
    return std::to_string(variant / 3 + 1) + "x" + family[variant % 3];
}

int SymmetrizeBandit::choose(std::mt19937_64& rng, double floor) const {
    double total = 0.0;
    for (const Arm& a : arms) total += a.recent + floor;
    double r = std::uniform_real_distribution<double>(0.0, total)(rng);
    for (int k = 0; k < kSymmetrizeVariants; ++k) {
        r -= arms[static_cast<size_t>(k)].recent + floor;
        if (r < 0.0) return k;
    }
    return kSymmetrizeVariants - 1;  // only by rounding
}

void SymmetrizeBandit::update(int variant, bool improved, double seconds, double decay) {
    Arm& a = arms[static_cast<size_t>(variant)];
    ++a.attempts;
    if (improved) ++a.improvements;
    a.cpu_seconds += seconds;
    a.recent += decay * ((improved ? 1.0 : 0.0) - a.recent);
}

bool try_symmetrize_region(System& sys, const LnmOptions& opt, SymmetrizeBandit& bandit,
                           std::mt19937_64& rng, double* improvement) {
    if (improvement) *improvement = 0.0;
    if (sys.particles.size() < 2 || !is_square_config(sys)) return false;

    const auto t0 = std::chrono::steady_clock::now();
    const SymmetrizeOptions& so = opt.symmetrize;
    const size_t n = sys.particles.size();
    const int variant = bandit.choose(rng, so.floor);
    // A sampled state can overlap, and then its tight side is not a side it has earned.
    const double reference = sys.is_feasible(kCertifyTol) ? sys.tight_side() : sys.L;

    // All the work happens on a copy, which is what makes the reject path exact.
    System trial = sys;
    trial.recount();
    std::vector<char> changed(n, 0);
    const Placement placement = choose_placement(trial, variant, so.eps, rng);

    bool accepted = false;
    double gain = 0.0;
    if (apply_placement(trial, placement, opt, rng, changed)) {
        if (opt.watch) opt.watch(trial, ("symmetrize " + symmetrize_variant_name(variant)).c_str());

        const double pitch = contact_pitch(trial.base_shape());
        relax_window(trial, seam_window(trial, changed, so.seam_width * pitch), so.seam_iterations,
                     opt.watch);

        CompressOptions copt;
        copt.shrink = so.shrink;
        copt.grow = so.grow;
        copt.max_growth_stages = so.growth_stages;
        copt.max_relaxations = so.relaxations;
        copt.relax.watch = opt.watch;
        copt.relax.watch_label = "compress";
        const CompressResult res = compress(trial, copt);

        assert(!res.found || res.packing.particles.size() == n);
        if (res.found && res.packing.is_feasible(kCertifyTol) && res.side <= reference) {
            sys = res.packing;
            accepted = true;
            gain = reference - res.side;
        }
        if (opt.watch) opt.watch(sys, accepted ? "accepted" : "rejected");
    }

    const std::chrono::duration<double> dt = std::chrono::steady_clock::now() - t0;
    bandit.update(variant, gain > 0.0, dt.count(), so.decay);
    if (improvement) *improvement = gain;
    return accepted;
}

bool apply_symmetry(System& sys, BoxRegion source, const std::vector<D4>& images) {
    if (sys.particles.size() < 2 || !is_square_config(sys)) return false;
    Placement placement{source, {}};
    for (const D4 w : images) placement.images.push_back(group_elem(w));
    sys.recount();
    std::vector<char> changed(sys.particles.size(), 0);
    std::mt19937_64 rng(0);  // only places spare squares, and a fixed seed keeps the call repeatable
    const bool done = apply_placement(sys, placement, LnmOptions{}, rng, changed);
    sys.recount();
    return done;
}

bool mirror_symmetry(System& sys) { return apply_symmetry(sys, kTopLeftTriangle, {D4::MirrorDiag}); }

// --- expand -------------------------------------------------------------------------------------

namespace {

bool expandable(const System& sys) {
    return !sys.particles.empty() && is_square_config(sys) &&
           sys.mobile_count() == static_cast<int>(sys.particles.size());
}

// How far a point is from the nearest wall or square: the radius of the hole it sits in.
double clearance_at(const System& sys, double x, double y) {
    double c = std::min({x, sys.L - x, y, sys.L - y});
    const double h = 0.5 * contact_pitch(sys.base_shape());
    for (const Particle& p : sys.particles) {
        const double dx = x - p.x(), dy = y - p.y();
        const double u = dx * p.cos_t() + dy * p.sin_t(), v = -dx * p.sin_t() + dy * p.cos_t();
        c = std::min(c, std::hypot(std::max(std::abs(u) - h, 0.0), std::max(std::abs(v) - h, 0.0)));
    }
    return c;
}

// Number of squares a line at `c` along axis `direction` cuts through.
int straddling(const System& sys, Vec2 axis, double c) {
    int count = 0;
    for (const Particle& p : sys.particles) {
        double lo = 0.0, hi = 0.0;
        project(p, sys.base_shape(), axis, lo, hi);
        if (lo < c - 1e-9 && hi > c + 1e-9) ++count;
    }
    return count;
}

Vec2 axis_of(int direction) { return direction == 0 ? Vec2{1.0, 0.0} : Vec2{0.0, 1.0}; }

// Cut position in [lo, hi] that crosses the fewest squares, summed over `offsets` from it; a
// random one among those tied, so that repeated draws can find different rows.
double cleanest_cut(const System& sys, Vec2 axis, double lo, double hi,
                    const std::vector<double>& offsets, std::mt19937_64& rng) {
    const double step = 0.05 * contact_pitch(sys.base_shape());
    std::vector<double> best_cuts;
    int best = std::numeric_limits<int>::max();
    for (double c = lo; c <= hi; c += step) {
        int count = 0;
        for (const double off : offsets) count += straddling(sys, axis, c + off);
        if (count < best) {
            best = count;
            best_cuts.clear();
        }
        if (count == best) best_cuts.push_back(c);
    }
    return best_cuts.empty() ? lo : best_cuts[rng() % best_cuts.size()];
}

double coordinate(const Particle& p, int direction) { return direction == 0 ? p.x() : p.y(); }

void shift_particle(Particle& p, int direction, double by) {
    if (direction == 0) p.set_pose(p.x() + by, p.y(), p.theta());
    else p.set_pose(p.x(), p.y() + by, p.theta());
}

// Separation of two squares on their four face normals; negative means they overlap.
double sat_gap(const Particle& a, const Particle& b, const Shape& s) {
    double gap = -std::numeric_limits<double>::infinity();
    const auto check = [&](Vec2 axis) {
        double alo = 0.0, ahi = 0.0, blo = 0.0, bhi = 0.0;
        project(a, s, axis, alo, ahi);
        project(b, s, axis, blo, bhi);
        gap = std::max({gap, blo - ahi, alo - bhi});
    };
    for (int k = 0; k < static_cast<int>(s.sat_axes.size()); ++k) {
        check(world_axis(a, s, k));
        check(world_axis(b, s, k));
    }
    return gap;
}

struct Contacts {
    int count = 0;
    double slack = std::numeric_limits<double>::infinity();  // smallest gap to anything
};

std::vector<Contacts> contact_table(const System& sys) {
    const double touching = 1e-6;
    const Shape& s = sys.base_shape();
    const double reach = 2.0 * s.circumradius + touching;
    const size_t n = sys.particles.size();
    std::vector<Contacts> out(n);
    const auto note = [&](Contacts& c, double gap) {
        c.slack = std::min(c.slack, gap);
        if (gap <= touching) ++c.count;
    };
    for (size_t i = 0; i < n; ++i) {
        const Particle& p = sys.particles[i];
        double lo = 0.0, hi = 0.0;
        project(p, s, {1.0, 0.0}, lo, hi);
        note(out[i], lo);
        note(out[i], sys.L - hi);
        project(p, s, {0.0, 1.0}, lo, hi);
        note(out[i], lo);
        note(out[i], sys.L - hi);
        for (size_t j = i + 1; j < n; ++j) {
            const Particle& q = sys.particles[j];
            if (std::hypot(p.x() - q.x(), p.y() - q.y()) > reach) continue;
            const double gap = sat_gap(p, q, s);
            note(out[i], gap);
            note(out[j], gap);
        }
    }
    return out;
}

System without_square(const System& sys, size_t index) {
    System out = sys;
    out.particles.erase(out.particles.begin() + static_cast<long>(index));
    out.recount();
    return out;
}

// compress() from a candidate; the candidate's box may well be too small for it.
void evaluate(Pool& pool, const System& candidate, MoveId op, const PropagateOptions& opt,
              PropagateReport& report, MoveLog* log) {
    const int m = static_cast<int>(candidate.particles.size());
    if (opt.only_target > 0 && m != opt.only_target) return;
    const PoolEntry* entry = pool.at(m);
    if (!entry || !pool.first_seen(candidate)) return;

    MoveTimer timer(log, op);
    ++report.tried;
    const double before = entry->side;

    CompressOptions copt;
    copt.shrink = opt.shrink;
    copt.grow = opt.grow;
    copt.max_growth_stages = opt.growth_stages;
    copt.max_relaxations = opt.relaxations;
    const CompressResult res = compress(candidate, copt);

    bool accepted = false;
    double gain = 0.0;
    if (res.found && pool.offer(res.packing)) {
        accepted = true;
        gain = before - pool.at(m)->side;
        report.improved[m] += gain;
    }
    timer.finish(accepted, gain);
}

}  // namespace

std::optional<System> symmetric_grow(const System& sys, bool centre, double tol) {
    if (!expandable(sys)) return std::nullopt;
    if (!contains(detect_symmetry(sys, tol), group_elem(D4::Rot90))) return std::nullopt;

    const double L = sys.L, h = 0.5 * L;
    const double pitch = contact_pitch(sys.base_shape());
    System out = sys;
    const int n = static_cast<int>(sys.particles.size());
    if (centre) {
        // A C4 packing with a square on the centre has no hole there to fill.
        if (clearance_at(sys, h, h) <= 0.0) return std::nullopt;
        out.add_particle(h, h, 0.0, 0);
    } else {
        // Sampled over one quadrant: the images fill the others. Closer to the centre than
        // pitch / sqrt(2) the four images would land on each other.
        const double step = 0.1 * pitch, too_close = pitch * 0.75;
        double best = 0.0, bx = -1.0, by = -1.0;
        for (double x = 0.5 * step; x <= h; x += step) {
            for (double y = 0.5 * step; y <= h; y += step) {
                if (std::hypot(x - h, y - h) < too_close) continue;
                const double c = clearance_at(sys, x, y);
                if (c > best) {
                    best = c;
                    bx = x;
                    by = y;
                }
            }
        }
        if (!(best > 0.0)) return std::nullopt;
        for (const D4 w : {D4::Id, D4::Rot90, D4::Rot180, D4::Rot270}) {
            const Pose q = group_elem(w).apply({bx, by, 0.0}, L);
            out.add_particle(q.x, q.y, q.theta, 0);
        }
    }

    // The new squares need room from somewhere: open the box by the area they cover.
    const double m = static_cast<double>(out.particles.size());
    out.recount();
    out.rescale_box(L * std::sqrt(m / n));
    return out;
}

std::optional<System> row_duplicate(const System& sys, int direction, std::mt19937_64& rng) {
    if (!expandable(sys)) return std::nullopt;
    const double pitch = contact_pitch(sys.base_shape());
    if (sys.L < pitch) return std::nullopt;

    const double cut = cleanest_cut(sys, axis_of(direction), pitch, sys.L, {0.0}, rng);
    System out = sys;
    std::vector<Pose> strip;
    for (size_t i = 0; i < sys.particles.size(); ++i) {
        const Particle& p = sys.particles[i];
        const double t = coordinate(p, direction);
        if (t >= cut - pitch && t < cut) {
            strip.push_back({direction == 0 ? p.x() + pitch : p.x(),
                             direction == 0 ? p.y() : p.y() + pitch, p.theta()});
        }
        if (t >= cut) shift_particle(out.particles[i], direction, pitch);
    }
    if (strip.empty()) return std::nullopt;
    for (const Pose& q : strip) out.add_particle(q.x, q.y, q.theta, 0);
    out.L = sys.L + pitch;
    out.recount();
    return out;
}

std::optional<System> row_remove(const System& sys, int direction, std::mt19937_64& rng) {
    if (!expandable(sys)) return std::nullopt;
    const double pitch = contact_pitch(sys.base_shape());
    if (sys.L <= 2.0 * pitch) return std::nullopt;

    const double cut = cleanest_cut(sys, axis_of(direction), 0.0, sys.L - pitch, {0.0, pitch}, rng);
    System out = sys;
    out.particles.clear();
    for (const Particle& p : sys.particles) {
        const double t = coordinate(p, direction);
        if (t >= cut && t < cut + pitch) continue;
        out.particles.push_back(p);
        if (t >= cut + pitch) shift_particle(out.particles.back(), direction, -pitch);
    }
    if (out.particles.empty() || out.particles.size() == sys.particles.size()) return std::nullopt;

    // The box is square and only one dimension just lost a row, so it can shrink by a pitch only
    // if the other dimension has that slack: a row_duplicate's leftover, say. Any further and
    // squares would start outside walls that compress could not grow back past.
    double farthest = sys.L - pitch;
    for (const Particle& p : out.particles) {
        double lo = 0.0, hi = 0.0;
        project(p, out.base_shape(), {1.0, 0.0}, lo, hi);
        farthest = std::max(farthest, hi);
        project(p, out.base_shape(), {0.0, 1.0}, lo, hi);
        farthest = std::max(farthest, hi);
    }
    out.L = farthest;
    out.recount();
    return out;
}

std::optional<System> frame_grow(const System& sys) {
    if (!expandable(sys)) return std::nullopt;
    const double p = contact_pitch(sys.base_shape());
    const int run = static_cast<int>(std::floor(sys.L / p)) + 1;

    System out = sys;
    for (Particle& q : out.particles) q.set_pose(q.x() + p, q.y() + p, q.theta());
    for (int k = 0; k < run; ++k) out.add_particle((k + 0.5) * p, 0.5 * p, 0.0, 0);
    for (int k = 1; k < run; ++k) out.add_particle(0.5 * p, (k + 0.5) * p, 0.0, 0);
    out.L = sys.L + p;
    out.recount();
    return out;
}

std::optional<System> delete_rattler(const System& sys) {
    if (!expandable(sys) || sys.particles.size() < 2) return std::nullopt;
    const std::vector<Contacts> table = contact_table(sys);
    size_t pick = 0;
    for (size_t i = 1; i < table.size(); ++i) {
        const bool fewer = table[i].count < table[pick].count;
        const bool tied = table[i].count == table[pick].count;
        if (fewer || (tied && table[i].slack > table[pick].slack)) pick = i;
    }
    return without_square(sys, pick);
}

std::vector<System> delete_for_compress(const System& sys, int count) {
    std::vector<System> out;
    if (!expandable(sys) || sys.particles.size() < 2 || count <= 0) return out;
    const std::vector<Contacts> table = contact_table(sys);
    std::vector<size_t> order(table.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(),
                     [&](size_t a, size_t b) { return table[a].count > table[b].count; });
    const size_t k = std::min(order.size(), static_cast<size_t>(count));
    for (size_t i = 0; i < k; ++i) out.push_back(without_square(sys, order[i]));
    return out;
}

// --- the pool -----------------------------------------------------------------------------------

int Pool::load_seeds(const std::string& dir, int lo, int hi) {
    int admitted = 0;
    for (int n = lo; n <= hi; ++n) {
        char name[32];
        std::snprintf(name, sizeof name, "/n%04d.cfg", n);
        const auto sys = load_record(dir + name);
        if (sys && offer(*sys)) ++admitted;
    }
    return admitted;
}

int Pool::load_records(const std::string& dir, int lo, int hi) {
    int admitted = 0;
    std::error_code ec;
    for (std::filesystem::directory_iterator it(dir, ec), end; !ec && it != end; it.increment(ec)) {
        if (it->path().extension() != ".log") continue;
        const auto sys = load_record(it->path().string());
        if (!sys) continue;
        const int n = static_cast<int>(sys->particles.size());
        if (n >= lo && n <= hi && offer(*sys)) ++admitted;
    }
    return admitted;
}

bool Pool::offer(const System& packing) {
    const double side = packing.tight_side();
    if (why_not_a_packing(packing, side, true, kCertifyTol)) return false;
    const int n = static_cast<int>(packing.particles.size());
    const auto it = entries_.find(n);
    // Same n, so comparing sides is comparing L(n) / sqrt(n).
    if (it != entries_.end() && !(side < it->second.side)) return false;
    entries_[n] = PoolEntry{packing, side};
    return true;
}

const PoolEntry* Pool::at(int n) const {
    const auto it = entries_.find(n);
    return it == entries_.end() ? nullptr : &it->second;
}

bool Pool::first_seen(const System& candidate) {
    const int n = static_cast<int>(candidate.particles.size());
    return seen_[n].insert(canonical_key(candidate)).second;
}

PropagateReport propagate(Pool& pool, int n, const PropagateOptions& opt, std::mt19937_64& rng,
                          MoveLog* log) {
    PropagateReport report;
    const PoolEntry* entry = pool.at(n);
    if (!entry) return report;

    // A copy: a candidate that improves the pool can replace the entry this was read from.
    const System source = entry->packing;
    const auto consider = [&](const std::optional<System>& c, MoveId op) {
        if (c) evaluate(pool, *c, op, opt, report, log);
    };

    consider(delete_rattler(source), MoveId::DeleteRattler);
    for (const System& c : delete_for_compress(source, opt.delete_trials)) {
        evaluate(pool, c, MoveId::DeleteCompress, opt, report, log);
    }
    consider(symmetric_grow(source, false), MoveId::SymmetricGrow);
    consider(symmetric_grow(source, true), MoveId::SymmetricGrow);
    consider(frame_grow(source), MoveId::FrameGrow);
    for (int direction = 0; direction < 2; ++direction) {
        for (int k = 0; k < opt.cuts; ++k) {
            consider(row_duplicate(source, direction, rng), MoveId::RowDuplicate);
            consider(row_remove(source, direction, rng), MoveId::RowRemove);
        }
    }
    return report;
}

}  // namespace pack
