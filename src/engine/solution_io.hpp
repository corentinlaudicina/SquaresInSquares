#pragma once

// sqpack/1 solution files (format in tools/validate.py), plus the legacy .cfg/.log records.

#include "system.hpp"

#include <nlohmann/json.hpp>

#include <optional>
#include <string>

namespace pack {

// A .json path is read as sqpack/1, anything else as a legacy record (io.hpp's load_record).
// The System comes back ready to use (unit squares, energy and grid built).
std::optional<System> read_solution(const std::string& path, std::string* error = nullptr);

// Writes sys as sqpack/1 with side = sys.L. Atomic: written to <path>.tmp, then renamed.
bool write_solution(const std::string& path, const System& sys, const nlohmann::json& meta);

}  // namespace pack
