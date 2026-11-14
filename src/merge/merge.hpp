#pragma once

// Merges a group of touching squares into one large rectangle, and various methods required to manipulate the new object.

#include "system.hpp"

#include <optional>
#include <vector>

namespace pack {

struct MergeMember {
    int species;
    bool frozen;
    Vec2 local_offset;     // member's position minus the group's, rotated into the group's body frame by -theta (i.e. in the group's own coordinates)

    double theta_offset;   // member's theta minus the group's theta
};

struct MergeRecord {
    std::vector<MergeMember> members;  // fixed order, set at merge time; dissolve() restores in this order (not necessarily the caller's order)
};

std::optional<int> merge(System& sys, const std::vector<int>& members);

// The exact inverse: if sys.is_merged(group_index), removes that one particle and re-inserts its original members at the poses they had the moment they were merged (each with its original species and frozen flag, reconstructed from the group's CURRENT pose plus its stored offsets), calls sys.recount(), and returns true. Returns false, without mutating sys, if group_index is out of range or does not name a merged particle.
bool dissolve(System& sys, int group_index);

}  // namespace pack
