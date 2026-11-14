#include "engine/solution_io.hpp"

#include "io/io.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>

namespace pack {

namespace fs = std::filesystem;
using json = nlohmann::json;

std::optional<System> read_solution(const std::string& path, std::string* error) {
    if (fs::path(path).extension() != ".json") return load_record(path, error);
    const auto fail = [&](const std::string& why) -> std::optional<System> {
        if (error) *error = path + ": " + why;
        return std::nullopt;
    };
    std::ifstream in(path);
    if (!in) return fail("cannot open");
    json doc;
    try {
        doc = json::parse(in);
        if (doc.at("format") != "sqpack/1") return fail("format is not sqpack/1");
        System sys;
        sys.species.push_back(make_square(1.0));
        sys.L = doc.at("side").get<double>();
        const auto& sq = doc.at("squares");
        if (sq.size() != doc.at("n").get<size_t>()) return fail("n does not match the squares");
        for (const auto& s : sq) {
            const double x = s.at(0).get<double>(), y = s.at(1).get<double>(), t = s.at(2).get<double>();
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(t)) return fail("non-finite pose");
            sys.add_particle(x, y, t, 0);
        }
        if (!(sys.L > 0.0) || !std::isfinite(sys.L)) return fail("side is not positive");
        sys.recount();
        return sys;
    } catch (const json::exception& e) {
        return fail(e.what());
    }
}

bool write_solution(const std::string& path, const System& sys, const json& meta) {
    json sq = json::array();
    for (const Particle& p : sys.particles) sq.push_back({p.x(), p.y(), p.theta()});
    const json doc = {{"format", "sqpack/1"}, {"n", sys.particles.size()}, {"side", sys.L},
                      {"squares", sq}, {"meta", meta}};
    const fs::path out(path);
    std::error_code ec;
    if (out.has_parent_path()) fs::create_directories(out.parent_path(), ec);
    const std::string tmp = path + ".tmp";
    {
        std::ofstream f(tmp);
        if (!f) return false;
        f << doc.dump(-1) << '\n';  // nlohmann prints doubles round-trip exact
        if (!f.good()) return false;
    }
    fs::rename(tmp, out, ec);
    return !ec;
}

}  // namespace pack
