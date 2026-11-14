#include "moves/move_log.hpp"

#include <filesystem>

namespace pack {
namespace {

double per_second(const MoveCounters& c) {
    return (c.cpu_seconds > 0.0) ? c.improvement / c.cpu_seconds : 0.0;
}

}  // namespace

const char* move_name(MoveId id) {
    switch (id) {
        case MoveId::Translate: return "translate";
        case MoveId::Rotate: return "rotate";
        case MoveId::Flip: return "flip";
        case MoveId::Shove: return "shove";
        case MoveId::Destroy: return "destroy";
        case MoveId::Strip: return "strip";
        case MoveId::Corner: return "corner";
        case MoveId::Loosen: return "loosen";
        case MoveId::Symmetrize: return "symmetrize";
        case MoveId::CycleShake: return "cycle_shake";
        case MoveId::CycleDestroy: return "cycle_destroy";
        case MoveId::CycleStrip: return "cycle_strip";
        case MoveId::CycleCorner: return "cycle_corner";
        case MoveId::CycleLoosen: return "cycle_loosen";
        case MoveId::CycleSymmetrize: return "cycle_symmetrize";
        case MoveId::SymmetricGrow: return "symmetric_grow";
        case MoveId::RowDuplicate: return "row_duplicate";
        case MoveId::RowRemove: return "row_remove";
        case MoveId::FrameGrow: return "frame_grow";
        case MoveId::DeleteRattler: return "delete_rattler";
        case MoveId::DeleteCompress: return "delete_compress";
        case MoveId::Count: break;
    }
    return "?";
}

void MoveLog::record(MoveId id, bool accepted, double improvement, double seconds) {
    MoveCounters& c = rows_[static_cast<size_t>(id)];
    ++c.attempts;
    if (accepted) ++c.accepted;
    if (improvement > 0.0) {
        ++c.improved;
        c.improvement += improvement;
    }
    c.cpu_seconds += seconds;
}

void MoveLog::print(std::FILE* out) const {
    std::fprintf(out, "%-18s %11s %11s %9s %11s %12s %13s\n", "move", "attempts", "accepted",
                 "improved", "cpu_s", "improvement", "improve/cpu_s");
    for (size_t k = 0; k < rows_.size(); ++k) {
        const MoveCounters& c = rows_[k];
        std::fprintf(out, "%-18s %11ld %11ld %9ld %11.3f %12.3e %13.3e\n",
                     move_name(static_cast<MoveId>(k)), c.attempts, c.accepted, c.improved,
                     c.cpu_seconds, c.improvement, per_second(c));
    }
}

bool MoveLog::write_csv(const std::string& path) const {
    const std::filesystem::path file(path);
    std::error_code ec;  // a failed mkdir shows up as a failed fopen below
    if (file.has_parent_path()) std::filesystem::create_directories(file.parent_path(), ec);
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    std::fprintf(f, "move,attempts,accepted,improved,cpu_seconds,improvement,improvement_per_cpu\n");
    for (size_t k = 0; k < rows_.size(); ++k) {
        const MoveCounters& c = rows_[k];
        std::fprintf(f, "%s,%ld,%ld,%ld,%.6f,%.12g,%.12g\n", move_name(static_cast<MoveId>(k)),
                     c.attempts, c.accepted, c.improved, c.cpu_seconds, c.improvement,
                     per_second(c));
    }
    std::fclose(f);
    return true;
}

}  // namespace pack
