#include "merge/merge.hpp"

#include "merge/rectangle.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>

namespace pack {
namespace {

Vec2 rotate(Vec2 v, double cos_t, double sin_t) {
    return {v.x * cos_t - v.y * sin_t, v.x * sin_t + v.y * cos_t};
}

}  // namespace

std::optional<int> merge(System& sys, const std::vector<int>& members) {
    const int n = static_cast<int>(sys.particles.size());

    // Validate before touching anything: no member already merged, and one shared frozen state.
    if (members.size() < 2) return std::nullopt;
    for (int i : members) {
        if (i < 0 || i >= n) return std::nullopt;
        if (sys.is_merged(i)) return std::nullopt;  // no merging a merged group
    }
    const bool shared_frozen = sys.particles[static_cast<size_t>(members[0])].frozen;
    for (int i : members) {
        if (sys.particles[static_cast<size_t>(i)].frozen != shared_frozen) return std::nullopt;
    }

    const std::optional<RectangleUnion> ru = detect_rectangle_union(sys, members);
    if (!ru) return std::nullopt;

    const Shape rect = make_rectangle(static_cast<double>(ru->cols) * ru->side,
                                      static_cast<double>(ru->rows) * ru->side);
    const double gc = std::cos(ru->theta), gs = std::sin(ru->theta);

    // Record each member's ACTUAL pose relative to the group's centre/orientation -- not the
    // idealized grid slot detect_rectangle_union reports -- so dissolve() reproduces the exact
    // poses even when a member sat fractionally off-grid (within the detector's tolerance).
    MergeRecord rec;
    rec.members.reserve(ru->ordered_members.size());
    for (int idx : ru->ordered_members) {
        const Particle& m = sys.particles[static_cast<size_t>(idx)];
        const Vec2 d{m.x() - ru->center.x, m.y() - ru->center.y};
        const Vec2 local = rotate(d, gc, -gs);  // R(-theta_group) * d
        rec.members.push_back({m.species, m.frozen, local, m.theta() - ru->theta});
    }

    const int new_species = sys.add_species(rect);

    // Remove the members highest-index-first, so an earlier erase never invalidates a later
    // one's index.
    std::vector<int> sorted = members;
    std::sort(sorted.begin(), sorted.end(), std::greater<int>());
    for (int idx : sorted) {
        sys.particles.erase(sys.particles.begin() + idx);
        if (static_cast<size_t>(idx) < sys.merge_of.size()) {
            sys.merge_of.erase(sys.merge_of.begin() + idx);
        }
    }

    const int gi = sys.add_particle(ru->center.x, ru->center.y, ru->theta, new_species,
                                    shared_frozen);
    sys.merge_of.resize(sys.particles.size());  // pads with nullptr up to and including gi
    sys.merge_of[static_cast<size_t>(gi)] = std::make_shared<const MergeRecord>(std::move(rec));

    sys.recount();
    return gi;
}

bool dissolve(System& sys, int group_index) {
    if (!sys.is_merged(group_index)) return false;

    const MergeRecord rec = *sys.merge_of[static_cast<size_t>(group_index)];
    const Particle group = sys.particles[static_cast<size_t>(group_index)];  // copy: erased just below
    const double gc = group.cos_t(), gs = group.sin_t();

    sys.particles.erase(sys.particles.begin() + group_index);
    sys.merge_of.erase(sys.merge_of.begin() + group_index);
    // Every remaining merge_of entry still correctly describes the particle now at its index:
    // erase shifted both vectors by exactly one, together. No padding is needed below -- every
    // member this loop appends is ordinary (merge() never admits an already-merged member), and
    // that is exactly what an index past merge_of.size() already means.

    for (const MergeMember& m : rec.members) {
        const Vec2 world_offset = rotate(m.local_offset, gc, gs);  // R(theta_group) * local
        sys.add_particle(group.x() + world_offset.x, group.y() + world_offset.y,
                         group.theta() + m.theta_offset, m.species, m.frozen);
    }

    sys.recount();
    return true;
}

}  // namespace pack
