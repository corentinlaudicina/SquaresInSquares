#pragma once

// A run's record: not its progress, only its answer.
//
// save_record() writes the one configuration a run is willing to claim -- its final, certified
// packing -- verifies the claim, and keeps it only if it ranks among the kKeepPerN densest already
// held for this square count; otherwise it deletes it.
//
// The written record is a header (n, L, the two stiffnesses, the tolerance), then one block per
// species (`s <index> <nv>` followed by that many `v <x> <y>` lines, the species' local-frame
// vertices), then one `p <x> <y> <theta> <frozen> <species>` line per particle. A merged
// particle's provenance (merge/merge.hpp's MergeRecord) is not written, so a merged particle
// comes back from load_record() as one plain rigid particle that cannot be dissolve()d.

#include "system.hpp"

#include <cstdio>
#include <optional>
#include <string>
#include <vector>

namespace pack {

// How many records survive per square count. Runs are ranked by the side they certified, so
// these are the five densest packings ever found for that n and nothing else.
inline constexpr int kKeepPerN = 5;

// Set this environment variable (to anything but "0" or empty) to keep every record: nothing is
// deleted. A record that verified and ranked is still committed, side and all; one that did not
// stays under its plain name, so it is on disk to be read but no later ranking counts it and
// nothing takes it for a record. For the runs where the failure is the interesting part.
inline constexpr const char* kKeepAllEnv = "SQUARES_KEEP_ALL_LOGS";

// What became of the record. `reason` always says why, kept or not, because "it is gone" is
// the kind of thing a driver has to be able to print.
struct LogVerdict {
    bool kept = false;
    std::string path;                  // the committed name when kept; the plain one when not
    std::string reason;                // why it was kept or dropped
    int rank = 0;                      // 1 = densest record held for this n; 0 when not kept
    int held = 0;                      // records now held for this n, at most kKeepPerN
    std::vector<std::string> removed;  // records this result pushed out
};

// Claims the first free name of the form "<dir>/<YYYY-MM-DD>_<HHMMSS>_n<count>.log", with a _2,
// _3 ... suffix if that one is taken, creating `dir` if needed. Naming and opening are one step
// on purpose -- see the note in io.cpp -- so a batch of runs finishing in the same second cannot
// land on one name.
//
// `certified` is the solver's own verdict -- relax() reaching Feasible, or compress() finding a
// box -- and `tol` the overlap tolerance it reached it at. Both are re-measured, in this order:
//   * the run certified a packing at all
//   * the claimed side is a positive, finite number
//   * every pose in it is finite
//   * the claimed side is the one the poses actually have
//   * no two squares overlap, and none leaves the box, by more than tol
//   * the side is not below sqrt(n * area), which no packing can beat
//
// Kept means the configuration is written in, and the file renamed to "<name>_s<side>.log" --
// the ranking reads the sides straight off the directory, so there is no index file to fall out
// of step with what is on disk -- and any record past kKeepPerN deleted. Anything else means the
// file is deleted. Failing to open a file is reported in the verdict, not fatal.
// The checks listed above, without the file: the reason this is not a packing, or nullptr.
const char* why_not_a_packing(const System& packing, double side, bool certified, double tol);

LogVerdict save_record(const std::string& dir, const System& packing, double side, bool certified,
                       double tol);

// Reads a file in the format above: a saved record, or a starting configuration from
// data/init_external. The returned System is ready to use, with its energy and neighbour grid
// built. Nothing if the file is unreadable or malformed, and then `error`, if given, says where
// and why. It checks the file is well formed, not that the packing is valid: a starting
// configuration is allowed to overlap, so ask is_feasible() if that matters.
std::optional<System> load_record(const std::string& path, std::string* error = nullptr);

// A run's progress, for src/viz/show.py --watch: the record format for the starting configuration,
// then one `f <step> <L> <energy> <x> <y> <theta> ...` line per frame. One frame is one line so a
// reader can tell a finished frame from one still being written. Not a record: load_record()
// refuses the `f` lines.
//
// open_trajectory replaces any file at `path`, creating its folder if needed, and returns nullptr
// if it cannot be opened; append_frame ignores a nullptr, so a run can carry on without one.
// fclose() the result when the run ends.
std::FILE* open_trajectory(const std::string& path, const System& start);
// `label`, if given, says what the frame shows: it goes on an `m <label>` line just before it.
void append_frame(std::FILE* f, long step, const System& sys, const char* label = nullptr);

}  // namespace pack
