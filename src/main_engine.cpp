// engine --N <n> --seed <s> (--time <sec> | --kicks <k>) --out <file.json> [--schedule <file.json>]
//        [--start <file>] [--telemetry <file>] [--mode shrink|cycle]
//
// Runs one search, legalizes its best packing (legalize.hpp) and writes it as sqpack/1, plus a
// telemetry JSON (default: <out minus .json>.telemetry.json) and a one-line summary on stdout.
// Writes nothing else: no data/best, no logs/. The runner validates the output and decides where
// it goes. Exit 0 = wrote a solution, 1 = run failed (no packing), 2 = usage error.

#include "engine/protocol.hpp"
#include "engine/solution_io.hpp"
#include "legalize.hpp"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

#ifndef ENGINE_VERSION
#define ENGINE_VERSION "unknown"
#endif

using namespace pack;
using json = nlohmann::json;

namespace {

[[noreturn]] void usage(const std::string& why) {
    std::fprintf(stderr, "engine: %s\nusage: engine --N n --seed s (--time sec | --kicks k) --out file.json "
                         "[--schedule f.json] [--start file] [--telemetry file] [--mode shrink|cycle]\n",
                 why.c_str());
    std::exit(2);
}

// strtod/strtoull with the whole argument checked, unlike the old drivers (audit B9).
double number(const char* s, const char* flag) {
    char* end = nullptr;
    errno = 0;
    const double v = std::strtod(s, &end);
    if (errno || end == s || *end) usage(std::string("bad value for ") + flag + ": " + s);
    return v;
}
unsigned long long integer(const char* s, const char* flag) {
    char* end = nullptr;
    errno = 0;
    const unsigned long long v = std::strtoull(s, &end, 10);
    if (errno || end == s || *end || *s == '-') usage(std::string("bad value for ") + flag + ": " + s);
    return v;
}

}  // namespace

int main(int argc, char** argv) {
    int n = 0;
    unsigned long long seed = 1;
    double time_s = 0.0;
    long kicks = 0;
    std::string out, schedule, start, telemetry, mode = "shrink";
    for (int k = 1; k < argc; ++k) {
        const std::string a = argv[k];
        if (k + 1 >= argc) usage(a + " needs a value");
        const char* v = argv[++k];
        if (a == "--N") n = static_cast<int>(integer(v, "--N"));
        else if (a == "--seed") seed = integer(v, "--seed");
        else if (a == "--time") time_s = number(v, "--time");
        else if (a == "--kicks") kicks = static_cast<long>(integer(v, "--kicks"));
        else if (a == "--out") out = v;
        else if (a == "--schedule") schedule = v;
        else if (a == "--start") start = v;
        else if (a == "--telemetry") telemetry = v;
        else if (a == "--mode") mode = v;
        else usage("unknown option " + a);
    }
    if (n < 1) usage("--N is required");
    if (out.empty()) usage("--out is required");
    if (time_s <= 0.0 && kicks <= 0) usage("give --time and/or --kicks");
    if (!schedule.empty()) usage("--schedule is not implemented yet");  // step 3
    if (mode != "shrink" && mode != "cycle") usage("--mode must be shrink or cycle");
    if (start.empty()) {
        char p[64];
        std::snprintf(p, sizeof p, "data/init_external/n%04d.cfg", n);
        start = p;
    }
    if (telemetry.empty()) {
        telemetry = out.size() > 5 && out.compare(out.size() - 5, 5, ".json") == 0 ? out.substr(0, out.size() - 5) : out;
        telemetry += ".telemetry.json";
    }

    std::string why;
    auto seed_sys = read_solution(start, &why);
    if (!seed_sys) usage(why);
    if (static_cast<int>(seed_sys->particles.size()) != n) usage(start + " does not hold N squares");

    EngineConfig cfg = mode == "cycle" ? cycle_defaults() : shrink_defaults();
    cfg.wall_seconds = time_s;
    cfg.max_kicks = kicks;

    EngineResult res = run_engine(*seed_sys, cfg, seed);
    json& tel = res.telemetry;
    tel["run"]["engine"] = ENGINE_VERSION;
    tel["run"]["start_file"] = start;

    int rc = 1;
    json leg = {{"ok", false}};
    if (res.found) {
        System legal = res.best;
        const LegalizeResult lr = legalize(legal);
        leg = {{"ok", lr.ok}, {"side_before", lr.side_before}, {"side_after", lr.side_after},
               {"growth_steps", lr.growth_steps}, {"relax_iters", lr.relax_iterations},
               {"min_pair_gap", lr.report.min_pair_gap}, {"min_wall_gap", lr.report.min_wall_gap}};
        const json meta = {{"engine", ENGINE_VERSION}, {"seed", seed}, {"mode", cfg.mode},
                           {"start_file", start}, {"legalized", lr.ok}, {"side_raw", res.best_side},
                           {"wall_s", tel["budget"]["wall_s"]}};
        if (write_solution(out, legal, meta)) rc = 0;
        else std::fprintf(stderr, "engine: could not write %s\n", out.c_str());
        tel["side"]["legal"] = lr.side_after;
    }
    tel["legalize"] = leg;
    {
        std::ofstream f(telemetry);
        f << tel.dump(1) << '\n';
    }

    // The summary: one line a planner can read at a glance.
    const json& s = tel["side"];
    std::string kick_text;
    for (const auto& [name, k] : tel["kicks"].items()) {
        kick_text += " " + name + " " + std::to_string(k["impr"].get<long>()) + "/" + std::to_string(k["tries"].get<long>());
    }
    std::printf("n=%d s%llu %s %.0fs %ld kicks: seed %.9f start %.9f -> best %.9f legal %.9f (%s) |%s | stop %s\n",
                n, seed, cfg.mode.c_str(), tel["budget"].value("wall_s", 0.0), tel["budget"].value("kicks", 0L),
                s.value("seed", 0.0), s.value("start", 0.0), s.value("best", 0.0), s.value("legal", 0.0),
                leg["ok"].get<bool>() ? "LEGAL" : "NOT LEGAL", kick_text.c_str(),
                tel["budget"].value("stop", std::string("?")).c_str());
    return rc;
}
