#pragma once

// What each move type costs and what it earns.
//
// Two levels, because the side only changes in one of them. A sweep move or a large move works
// at fixed L, so it can be accepted but cannot improve the side, and its `improved` stays 0 by
// construction. The cycle_* rows are whole kick -> relax -> compress cycles, credited to the
// kick that started them: that is where a side is won, so those rows and the moves that run
// their own compress (symmetrize, expand) are the ones whose improvement per second means
// anything. Cycle rows contain the CPU of the moves inside them, so the rows do not sum.
//
// Time is a monotonic clock read around each move, because a per-thread CPU clock costs about as
// much as a single-particle move does. On one thread that nothing preempts it is the CPU time.

#include <array>
#include <chrono>
#include <cstdio>
#include <string>

namespace pack {

enum class MoveId {
    Translate,
    Rotate,
    Flip,
    Shove,
    Destroy,
    Strip,
    Corner,
    Loosen,
    Symmetrize,
    CycleShake,
    CycleDestroy,
    CycleStrip,
    CycleCorner,
    CycleLoosen,
    CycleSymmetrize,
    SymmetricGrow,
    RowDuplicate,
    RowRemove,
    FrameGrow,
    DeleteRattler,
    DeleteCompress,
    Count,
};

const char* move_name(MoveId id);

struct MoveCounters {
    long attempts = 0;
    long accepted = 0;
    long improved = 0;          // attempts that made the side strictly smaller
    double cpu_seconds = 0.0;
    double improvement = 0.0;   // total decrease of the side
};

class MoveLog {
public:
    void record(MoveId id, bool accepted, double improvement, double seconds);
    const MoveCounters& at(MoveId id) const { return rows_[static_cast<size_t>(id)]; }

    // Every row is written, zeros included: a move that never fired should be visible as such.
    void print(std::FILE* out) const;
    bool write_csv(const std::string& path) const;

private:
    std::array<MoveCounters, static_cast<size_t>(MoveId::Count)> rows_{};
};

// Times one attempt. With a null log it reads no clock and records nothing.
class MoveTimer {
public:
    MoveTimer(MoveLog* log, MoveId id) : log_(log), id_(id) {
        if (log_) start_ = Clock::now();
    }

    void finish(bool accepted, double improvement = 0.0) {
        if (!log_) return;
        const std::chrono::duration<double> dt = Clock::now() - start_;
        log_->record(id_, accepted, improvement, dt.count());
    }

private:
    using Clock = std::chrono::steady_clock;
    MoveLog* log_;
    MoveId id_;
    Clock::time_point start_{};
};

}  // namespace pack
