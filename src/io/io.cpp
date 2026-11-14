#include "io/io.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <utility>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>

namespace pack {
namespace {

namespace fs = std::filesystem;

// The environment escape hatch: keep every record, verified or not, ranked or not.
bool keep_everything() {
    const char* v = std::getenv(kKeepAllEnv);
    return v && v[0] != '\0' && std::strcmp(v, "0") != 0;
}

// The side as the committed name will carry it. Ranking compares these rather than the raw
// doubles, so a run is never beaten -- or saved -- by digits its name does not record. Six
// decimals is well past the sqrt(tolerance) the side is actually good to.
std::string side_text(double side) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.6f", side);
    return buf;
}

double side_as_named(double side) { return std::strtod(side_text(side).c_str(), nullptr); }

// The side recorded in a committed record's name, for the records of one square count:
//
//     2026-10-28_215204_n11_s3.877084.log   ->  3.877084
//     2026-10-28_215204_n11.log             ->  nothing: never committed
//     2026-10-28_215204_n118_s9.138.log     ->  nothing when asked for n = 11
//
// Anything that does not parse is left alone. A file this does not recognise is never counted
// and never deleted, which is what keeps the pruning off records written before this policy, off
// a run still writing, and off whatever else the directory holds.
bool committed_side(const fs::path& p, std::size_t n, double& side) {
    if (p.extension() != ".log") return false;
    const std::string stem = p.stem().string();

    const std::string tag = "_n" + std::to_string(n);
    const std::size_t at = stem.find(tag);
    if (at == std::string::npos) return false;
    // "_n11" must not match inside "_n118": ours is followed by the clash suffix or by "_s".
    if (at + tag.size() < stem.size() && stem[at + tag.size()] != '_') return false;

    const std::size_t s = stem.rfind("_s");
    if (s == std::string::npos || s < at + tag.size()) return false;
    const std::string digits = stem.substr(s + 2);
    char* end = nullptr;
    const double v = std::strtod(digits.c_str(), &end);
    if (!end || *end != '\0' || digits.empty() || !std::isfinite(v) || v <= 0.0) return false;
    side = v;
    return true;
}

void remove_record(const fs::path& path, std::vector<std::string>* removed = nullptr) {
    std::error_code ec;
    if (fs::remove(path, ec) && removed) removed->push_back(path.string());
}

void write_configuration(std::FILE* f, const System& sys, double tol) {
    std::fprintf(f, "# squares-in-squares record\n");
    std::fprintf(f, "n %zu\n", sys.particles.size());
    std::fprintf(f, "L %.17g\n", sys.L);
    std::fprintf(f, "k_pair %.17g\n", sys.k_pair);
    std::fprintf(f, "k_wall %.17g\n", sys.k_wall);
    std::fprintf(f, "feasible_tol %.17g\n", tol);
    std::fprintf(f, "nspecies %zu\n", sys.species.size());
    for (size_t sp = 0; sp < sys.species.size(); ++sp) {
        const Shape& s = sys.species[sp];
        std::fprintf(f, "s %zu %zu\n", sp, s.vertices.size());
        for (const Vec2& v : s.vertices) std::fprintf(f, "v %.17g %.17g\n", v.x, v.y);
    }
    std::fprintf(f, "# one p <x> <y> <theta> <frozen> <species> per particle\n");
    for (const Particle& p : sys.particles) {
        std::fprintf(f, "p %.17g %.17g %.17g %d %d\n", p.x(), p.y(), p.theta(),
                     p.frozen ? 1 : 0, p.species);
    }
}

}  // namespace

// Re-measures what the run claims. Returns the reason it is not a packing, or nullptr.
const char* why_not_a_packing(const System& packing, double side, bool certified, double tol) {
    if (!certified) return "the run certified no packing";
    if (!std::isfinite(side) || side <= 0.0) {
        return "the claimed side is not a positive number";
    }
    for (const Particle& p : packing.particles) {
        if (!std::isfinite(p.x()) || !std::isfinite(p.y()) || !std::isfinite(p.theta())) {
            return "the result holds a pose that is not a finite number";
        }
    }
    // The side has to be the one these poses have, or the record is about some other
    // configuration than the one it was claimed for.
    const double geometric = packing.tight_side();
    if (std::abs(geometric - side) > 1e-9 * std::max(1.0, geometric)) {
        return "the claimed side is not the side of the poses it was claimed for";
    }
    if (!packing.is_feasible(tol)) return "the result still overlaps";
    // No arrangement of n shapes fits in less area than they cover, so a side below this is
    // arithmetic gone wrong somewhere upstream, not a record. Cheap, and the one check that
    // catches a claim too good to be true.
    double total_area = 0.0;
    for (size_t i = 0; i < packing.particles.size(); ++i) {
        total_area += packing.shape_of(static_cast<int>(i)).area;
    }
    const double area_floor = std::sqrt(total_area);
    if (side < area_floor * (1.0 - 1e-9)) {
        return "the claimed side is below the area floor, so the result cannot be real";
    }
    return nullptr;
}

std::FILE* open_trajectory(const std::string& path, const System& start) {
    const fs::path file(path);
    std::error_code ec;  // a failed mkdir shows up as a failed fopen below
    if (file.has_parent_path()) fs::create_directories(file.parent_path(), ec);
    // Removed, not truncated: a new file is how a viewer following the old one sees a new run.
    fs::remove(file, ec);

    std::FILE* f = std::fopen(path.c_str(), "w");
    if (f) {
        write_configuration(f, start, 0.0);
        std::fflush(f);
    }
    return f;
}

void append_frame(std::FILE* f, long step, const System& sys, const char* label) {
    if (!f) return;
    if (label) std::fprintf(f, "m %s\n", label);
    // Recounted, not sys.energy(): mid-relaxation the tracked sums are stale.
    std::fprintf(f, "f %ld %.17g %.17g", step, sys.L, sys.total_energy());
    for (const Particle& p : sys.particles) {
        std::fprintf(f, " %.10g %.10g %.10g", p.x(), p.y(), p.theta());
    }
    std::fputc('\n', f);
    std::fflush(f);  // so a viewer sees it now, and a crash loses nothing already written
}

LogVerdict save_record(const std::string& dir_name, const System& packing, double side,
                       bool certified, double tol) {
    std::error_code ec;  // a failed mkdir shows up as a failed fopen below
    fs::create_directories(dir_name, ec);

    const std::time_t now = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &now);
#else
    localtime_r(&now, &tm);
#endif
    char stamp[32];
    std::strftime(stamp, sizeof stamp, "%Y-%m-%d_%H%M%S", &tm);

    // The create IS the claim: "wx" fails if the name already exists, so each run walks up
    // the suffixes until one name is its own. Testing for the name first and then opening
    // with a plain "w" is what two runs finishing in the same second lose to -- both see the
    // name free, both open it, and the second one truncates the first one's record.
    const std::size_t n = packing.particles.size();
    const std::string base = dir_name + "/" + stamp + "_n" + std::to_string(n);
    std::FILE* f = nullptr;
    std::string path;
    // The loop ends because every turn tries a different name; the bound is only so that a
    // filesystem reporting EEXIST for a name it does not hold cannot hang the run.
    for (int k = 1; !f && k <= 1000; ++k) {
        path = (k == 1) ? base + ".log" : base + "_" + std::to_string(k) + ".log";
        f = std::fopen(path.c_str(), "wx");
        if (!f && errno != EEXIST) break;  // not a name clash: a real failure to open
    }

    LogVerdict verdict;
    verdict.path = path;
    if (!f) {
        verdict.reason = "no record was opened";
        return verdict;
    }

    // The one place a record is dropped. Under the override it stays on disk, but uncommitted:
    // it keeps its plain name, so no later ranking counts it and nothing reads it as a record.
    const bool keep_all = keep_everything();
    const auto drop = [&](std::string why) {
        verdict.reason = std::move(why);
        if (f) {
            std::fclose(f);
            f = nullptr;
        }
        if (keep_all) {
            verdict.kept = true;
            verdict.reason += "; kept unranked by ";
            verdict.reason += kKeepAllEnv;
        } else {
            remove_record(path);
        }
        return verdict;
    };

    if (const char* bad = why_not_a_packing(packing, side, certified, tol)) return drop(bad);

    // Verified: this is the one configuration the run is willing to claim, so write it now.
    write_configuration(f, packing, tol);
    std::fclose(f);
    f = nullptr;

    // Rank it against the committed records of the same square count, reading their sides out
    // of their names.
    struct Entry {
        double side;
        fs::path path;
        bool mine;
    };
    std::vector<Entry> board;
    const fs::path dir = fs::path(path).parent_path();
    for (fs::directory_iterator it(dir, ec), end; !ec && it != end; it.increment(ec)) {
        double held = 0.0;
        if (!it->is_regular_file(ec)) continue;
        if (committed_side(it->path(), n, held)) board.push_back({held, it->path(), false});
    }
    board.push_back({side_as_named(side), fs::path(path), true});

    // Densest first, and a tie leaves the incumbent in front: equalling the fifth-best side is
    // not a reason to delete a record that is already there.
    std::sort(board.begin(), board.end(), [](const Entry& a, const Entry& b) {
        if (a.side != b.side) return a.side < b.side;
        return !a.mine && b.mine;
    });

    std::size_t rank = 0;
    while (rank < board.size() && !board[rank].mine) ++rank;

    if (rank >= static_cast<std::size_t>(kKeepPerN)) {
        verdict.held = kKeepPerN;
        return drop("not among the " + std::to_string(kKeepPerN) + " densest runs for n " +
                    std::to_string(n) + " (side " + side_text(side) + ", worst kept " +
                    side_text(board[static_cast<std::size_t>(kKeepPerN) - 1].side) + ")");
    }

    // Commit: the side goes into the name, which is what makes the next run's ranking a
    // directory listing rather than a re-read of every record in it.
    fs::path committed = fs::path(path);
    const std::string stem = committed.stem().string();
    // The clash suffix goes BEFORE the side, never after it: the side has to be the last field
    // in the name or committed_side() cannot read it back, and a record it cannot read is one
    // that is never ranked and never pruned -- an orphan that outlives every run.
    //
    // Clashes are not vanishingly rare, either. Committing renames the claimed file out of the
    // way, which frees its plain name for the next run to claim, so two runs of the same n in
    // the same second really can reach this point with the same stem and the same side.
    for (int k = 1; k <= 1000; ++k) {
        const std::string suffix = (k == 1) ? "" : "_" + std::to_string(k);
        const fs::path candidate = dir / (stem + suffix + "_s" + side_text(side) + ".log");
        if (fs::exists(candidate, ec)) continue;
        fs::rename(path, candidate, ec);
        if (!ec) committed = candidate;
        break;
    }
    // A rename that failed leaves the record under its old, uncommitted name rather than losing
    // it; it is then invisible to ranking, which is the safe direction to fail in.
    verdict.kept = true;
    verdict.path = committed.string();
    verdict.rank = static_cast<int>(rank) + 1;

    if (!keep_all) {
        // Everything past the cut is an older record -- ours ranked inside it a moment ago.
        for (std::size_t i = static_cast<std::size_t>(kKeepPerN); i < board.size(); ++i) {
            remove_record(board[i].path, &verdict.removed);
        }
    }
    verdict.held = static_cast<int>(std::min<std::size_t>(board.size(), kKeepPerN));
    verdict.reason = "rank " + std::to_string(verdict.rank) + " of " +
                     std::to_string(verdict.held) + " for n " + std::to_string(n) + " (side " +
                     side_text(side) + ")";
    return verdict;
}

namespace {

// strtod rather than `stream >> double`, which refuses very small values such as 1e-320.
bool to_double(const std::string& word, double& out) {
    char* end = nullptr;
    out = std::strtod(word.c_str(), &end);
    return *end == '\0' && std::isfinite(out);
}

bool to_long(const std::string& word, long& out) {
    char* end = nullptr;
    errno = 0;
    out = std::strtol(word.c_str(), &end, 10);
    return *end == '\0' && errno == 0;
}

}  // namespace

std::optional<System> load_record(const std::string& path, std::string* error) {
    long line_no = 0;
    const auto fail = [&](const std::string& why) -> std::optional<System> {
        if (error) *error = path + ":" + std::to_string(line_no) + ": " + why;
        return std::nullopt;
    };

    std::ifstream in(path);
    if (!in) return fail("cannot open the file");

    System sys;
    double L = -1.0, k_pair = -1.0, k_wall = -1.0;
    long n = -1, nspecies = -1, vertices_left = 0;
    std::vector<Vec2> vertices;

    for (std::string line; std::getline(in, line);) {
        ++line_no;
        std::istringstream words(line);
        std::vector<std::string> f;
        for (std::string w; words >> w;) f.push_back(w);
        if (f.empty() || f[0][0] == '#') continue;
        const std::string& tag = f[0];

        if (vertices_left > 0 && tag != "v") return fail("expected another vertex line");

        if (tag == "n" || tag == "nspecies") {
            long v = 0;
            if (f.size() != 2 || !to_long(f[1], v)) return fail("bad integer after " + tag);
            (tag == "n" ? n : nspecies) = v;
        } else if (tag == "L" || tag == "k_pair" || tag == "k_wall" || tag == "feasible_tol") {
            double v = 0.0;
            if (f.size() != 2 || !to_double(f[1], v)) return fail("bad number after " + tag);
            if (tag == "L") L = v;
            if (tag == "k_pair") k_pair = v;
            if (tag == "k_wall") k_wall = v;
        } else if (tag == "s") {
            long index = 0, count = 0;
            if (f.size() != 3 || !to_long(f[1], index) || !to_long(f[2], count)) {
                return fail("bad species line");
            }
            if (index != static_cast<long>(sys.species.size())) return fail("species out of order");
            if (count < 3 || count > kMaxVertices) return fail("bad vertex count");
            vertices.clear();
            vertices_left = count;
        } else if (tag == "v") {
            Vec2 v;
            if (vertices_left == 0) return fail("vertex outside a species block");
            if (f.size() != 3 || !to_double(f[1], v.x) || !to_double(f[2], v.y)) {
                return fail("bad vertex line");
            }
            vertices.push_back(v);
            if (--vertices_left == 0) {
                const std::optional<Shape> shape = shape_from_vertices(vertices);
                if (!shape) return fail("species is not a convex counter-clockwise polygon");
                sys.add_species(*shape);
            }
        } else if (tag == "p") {
            double x = 0.0, y = 0.0, theta = 0.0;
            long frozen = 0, species = 0;
            if (f.size() != 6 || !to_double(f[1], x) || !to_double(f[2], y) ||
                !to_double(f[3], theta) || !to_long(f[4], frozen) || !to_long(f[5], species)) {
                return fail("bad particle line");
            }
            if (frozen != 0 && frozen != 1) return fail("frozen must be 0 or 1");
            if (species < 0 || species >= static_cast<long>(sys.species.size())) {
                return fail("particle names a species that does not exist");
            }
            sys.add_particle(x, y, theta, static_cast<int>(species), frozen == 1);
        } else {
            return fail("unknown line '" + tag + "'");
        }
    }

    if (vertices_left > 0) return fail("file ends inside a species block");
    if (n < 1) return fail("missing or bad n");
    if (nspecies != static_cast<long>(sys.species.size())) return fail("wrong number of species");
    if (static_cast<long>(sys.particles.size()) != n) return fail("wrong number of particles");
    if (!(L > 0.0)) return fail("missing or non-positive L");
    if (k_pair < 0.0 || k_wall < 0.0) return fail("missing or negative stiffness");

    sys.L = L;
    sys.k_pair = k_pair;
    sys.k_wall = k_wall;
    sys.recount();
    return sys;
}

}  // namespace pack
