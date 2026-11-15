// engine --N <n> --seed <s> (--time <sec> | --kicks <k>) --out <file.json> [--schedule <file.json>]
//        [--start <file>] [--telemetry <file>] [--neighbors <dir>] [--live <snapshot.json>] [--mode shrink|cycle|ga]   (--mode only without --schedule)
//        engine --list-schedule-keys
//        engine --check-schedule <file.json>
//        engine --legalize <in> --out <out.json>
//
// Runs one search, legalizes its best packing (legalize.hpp) and writes it as sqpack/1, plus a
// telemetry JSON (default: <out minus .json>.telemetry.json) and a one-line summary on stdout.
// Writes nothing else: no data/best, no logs/. The runner validates the output and decides where
// it goes. Exit 0 = wrote a solution, 1 = run failed (no packing), 2 = usage error.

#include "engine/protocol.hpp"
#include "engine/schedule.hpp"
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
    std::string out, schedule, start, telemetry, mode = "shrink", legalize_only, neighbors = "best", live;
    for (int k = 1; k < argc; ++k) {
        const std::string a = argv[k];
        if (a == "--check-schedule") {  // parse only: exit 0 if the engine accepts the file, 2 if not
            if (k + 1 >= argc) usage("--check-schedule needs a file");
            std::string why;
            const auto ps = load_schedule(argv[k + 1], &why);
            if (!ps) {
                std::fprintf(stderr, "%s\n", why.c_str());
                return 2;
            }
            std::printf("%s %s\n", ps->name.c_str(), ps->hash.c_str());
            return 0;
        }
        if (a == "--list-schedule-keys") {
            for (const std::string& f : field_names()) std::printf("%s\n", f.c_str());
            return 0;
        }
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
        else if (a == "--legalize") legalize_only = v;
        else if (a == "--neighbors") neighbors = v;
        else if (a == "--live") live = v;
        else usage("unknown option " + a);
    }
    // engine --legalize <in> --out <out.json>: no search, only legalization (for imports).
    if (!legalize_only.empty()) {
        if (out.empty()) usage("--out is required");
        std::string why;
        auto sys = read_solution(legalize_only, &why);
        if (!sys) usage(why);
        const LegalizeResult lr = legalize(*sys);
        const json meta = {{"engine", ENGINE_VERSION}, {"legalized", lr.ok}, {"legalized_from", legalize_only},
                           {"side_raw", lr.side_before}};
        if (!write_solution(out, *sys, meta)) return 1;
        std::printf("legalize %s: %.12f -> %.12f (%s, %d growths)\n", legalize_only.c_str(), lr.side_before,
                    lr.side_after, lr.ok ? "LEGAL" : "NOT LEGAL", lr.growth_steps);
        return lr.ok ? 0 : 1;
    }
    if (n < 1) usage("--N is required");
    if (out.empty()) usage("--out is required");
    if (time_s <= 0.0 && kicks <= 0) usage("give --time and/or --kicks");
    if (mode != "shrink" && mode != "cycle" && mode != "ga") usage("--mode must be shrink, cycle or ga");
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

    EngineConfig cfg = mode == "cycle" ? cycle_defaults() : mode == "ga" ? ga_defaults() : shrink_defaults();
    json sched_info = {{"name", "builtin-" + mode}};
    if (!schedule.empty()) {
        auto ps = load_schedule(schedule, &why);
        if (!ps) usage(why);
        cfg = ps->config;
        sched_info = {{"name", ps->name}, {"hash", ps->hash}, {"file", schedule}, {"resolved", ps->resolved}};
    }
    cfg.ga.seeding.neighbor_dir = neighbors;  // GA seeding reads N+1 / N-1 from here (default best/)
    if (!live.empty()) {  // snapshots for tools/live_view.py: written to <live>.tmp, then renamed
        cfg.live = [live](const json& snap) {
            const std::string tmp = live + ".tmp";
            {
                std::ofstream f(tmp);
                f << snap.dump() << '\n';
            }
            std::rename(tmp.c_str(), live.c_str());
        };
    }
    cfg.wall_seconds = time_s;
    cfg.max_kicks = kicks;

    EngineResult res = run_engine(*seed_sys, cfg, seed);
    json& tel = res.telemetry;
    tel["run"]["engine"] = ENGINE_VERSION;
    tel["run"]["start_file"] = start;
    tel["run"]["schedule"] = sched_info;

    int rc = 1;
    json leg = {{"ok", false}};
    if (res.found) {
        System legal = res.best;
        const LegalizeResult lr = legalize(legal);
        leg = {{"ok", lr.ok}, {"side_before", lr.side_before}, {"side_after", lr.side_after},
               {"growth_steps", lr.growth_steps}, {"relax_iters", lr.relax_iterations},
               {"min_pair_gap", lr.report.min_pair_gap}, {"min_wall_gap", lr.report.min_wall_gap}};
        const json meta = {{"engine", ENGINE_VERSION}, {"seed", seed}, {"mode", cfg.mode},
                           {"schedule", sched_info["name"]}, {"schedule_hash", sched_info.value("hash", "")},
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
                n, seed, sched_info["name"].get<std::string>().c_str(), tel["budget"].value("wall_s", 0.0), tel["budget"].value("kicks", 0L),
                s.value("seed", 0.0), s.value("start", 0.0), s.value("best", 0.0), s.value("legal", 0.0),
                leg["ok"].get<bool>() ? "LEGAL" : "NOT LEGAL", kick_text.c_str(),
                tel["budget"].value("stop", std::string("?")).c_str());
    return rc;
}
