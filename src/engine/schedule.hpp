#pragma once

// Schedule files: a search protocol as data (schema sqpack.schedule/1, JSON; the runner converts YAML).
//
//   {"schema": "sqpack.schedule/1", "name": "shrink-baseline", "mode": "shrink",
//    "start":   {"mirror": true},
//    "prepare": {"ds": 0.01, "max_relax": 5000},
//    "kick":    {"large_prob": {"base": 0.1, "max": 0.6, "step": 0.1, "every": 25},
//                "shake": {"sweeps": 100, "T": 0.05, "translate": 0.1, "rotate": 0.1,
//                          "mix": {"translate": 0.45, "rotate": 0.45, "flip": 0.06, "shove": 0.04}},
//                "lnm":   {"weights": {"destroy": 0.4, "strip": 0.3, "corner": 0.2, "loosen": 0.1, "symmetrize": 0},
//                          "relax_iterations": 300, "insert_tries": 12, "snap_prob": 0.5, ...}},
//    "shrink":  {"eps": 1e-3, "eps_min": 1e-7, "eps_max": 0.05, "on_success": 1.2, "on_fail": 0.95},
//    "cycle":   {"ds": 0.005, "max_relax": 5000, "side_T": 0.01, "restart_after": 200},
//    "propagate": {"enabled": false, "every": 1}}
//
// Every key is optional: an omitted key keeps the mode's default (shrink_defaults / cycle_defaults),
// so {} is today's main_shrink. Unknown keys are errors. Any number marked "curve" in field_names()
// may instead be a curve [[t0, v0], [t1, v1], ...]: piecewise linear in t = elapsed fraction of the
// run's budget, constant outside [t0, tlast]. The full key list is `engine --list-schedule-keys`.

#include "engine/protocol.hpp"

#include <nlohmann/json.hpp>

#include <optional>
#include <string>
#include <vector>

namespace pack {

struct ParsedSchedule {
    EngineConfig config;      // curves evaluated at t = 0; config.update applies them per kick
    std::string name;
    std::string hash;         // FNV-1a 64 of the resolved schedule, hex
    nlohmann::json resolved;  // every field with its effective value (curves as given)
};

std::optional<ParsedSchedule> parse_schedule(const nlohmann::json& doc, std::string* error);
std::optional<ParsedSchedule> load_schedule(const std::string& path, std::string* error);

// "kick.shake.T (curve)", ... : the vocabulary, one key per line, for docs and the runner.
std::vector<std::string> field_names();

}  // namespace pack
