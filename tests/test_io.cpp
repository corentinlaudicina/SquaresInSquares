// Assert-based tests for the run record.
//
// Two claims matter here, and nothing else in the suite checks either.
//
// test_each_record_gets_its_own_file: two runs finishing in the same second must not land on
// one name, or one overwrites the other's record and hours of compute go with it.
//
// The retention tests: save_record() deletes records, so every claim it acts on is checked here. A
// bug that keeps one record too many wastes a few kilobytes; a bug that drops a real record, or
// saves an unverified one as a record, loses the run or lies about it.

#include "init.hpp"
#include "io/io.hpp"
#include "system.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace pack;

namespace {

bool close(double a, double b, double tol = 1e-12) { return std::abs(a - b) <= tol; }

// A directory of our own, so a test run can neither read nor clobber the project's records.
std::filesystem::path fresh_dir(const char* name) {
    const std::filesystem::path d =
        std::filesystem::temp_directory_path() / ("sq_test_io_" + std::string(name));
    std::filesystem::remove_all(d);
    return d;
}

std::vector<std::string> read_lines(const std::string& path) {
    std::ifstream in(path);
    assert(in && "the verdict named this file, so it has to be readable");
    std::vector<std::string> out;
    for (std::string line; std::getline(in, line);) out.push_back(line);
    return out;
}

// Every word after the tag, for the first line carrying that tag.
std::vector<std::string> fields(const std::vector<std::string>& lines, const std::string& tag) {
    for (const std::string& line : lines) {
        std::istringstream s(line);
        std::string first;
        s >> first;
        if (first != tag) continue;
        std::vector<std::string> f;
        for (std::string w; s >> w;) f.push_back(w);
        return f;
    }
    assert(false && "tag missing from the record");
    return {};
}

long count_tag(const std::vector<std::string>& lines, const std::string& tag) {
    long k = 0;
    for (const std::string& line : lines) {
        std::istringstream s(line);
        std::string first;
        s >> first;
        if (first == tag) ++k;
    }
    return k;
}

// The tolerance a certified packing is held to: relax.hpp's kCertifyTol, spelled out rather
// than included, because this suite links the core alone and has no business pulling L-BFGS in.
constexpr double kTol = 1e-10;

// `n` unit squares in a row inside a box of side `side`, evenly spaced with the outer two
// against the walls: overlap-free, inside the box, and its tight side is exactly `side`. The
// loose grid grid_init() builds is none of those -- it hangs half a square outside the box --
// and the one knob a ranking test needs is how dense a run gets to claim to be.
System in_a_row(int n, double side) {
    assert(n >= 2 && side >= static_cast<double>(n));
    System sys = grid_init(make_square(1.0), n);
    sys.set_box(side);
    const double pitch = 1.0 + (side - n) / (n - 1);
    for (int k = 0; k < n; ++k) {
        sys.particles[static_cast<size_t>(k)].set_pose(0.5 + k * pitch, 0.5, 0.0);
    }
    sys.recount();  // not resync_energy: the poses moved, so it has drifted
    assert(sys.is_feasible(kTol) && close(sys.tight_side(), side));
    return sys;
}

System two_in_a_row(double side) { return in_a_row(2, side); }

LogVerdict run_claiming(const std::filesystem::path& dir, const System& packing, double side,
                        bool certified = true) {
    return save_record(dir.string(), packing, side, certified, kTol);
}

// Every side the directory is holding a committed record for, densest first. A record without
// a side in its name is not committed -- it was never ranked and is never counted here.
std::vector<double> kept_sides(const std::filesystem::path& dir) {
    std::vector<double> out;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (entry.path().extension() != ".log") continue;
        const std::string stem = entry.path().stem().string();
        const size_t at = stem.rfind("_s");
        if (at != std::string::npos) out.push_back(std::stod(stem.substr(at + 2)));
    }
    std::sort(out.begin(), out.end());
    return out;
}

long count_files(const std::filesystem::path& dir) {
    long k = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        (void)entry;
        ++k;
    }
    return k;
}

void test_a_kept_record_holds_the_certified_configuration() {
    // What lands on disk has to be exactly the packing that was claimed.
    const auto dir = fresh_dir("record");
    const System packing = in_a_row(4, 5.0);

    const LogVerdict v = run_claiming(dir, packing, 5.0);
    assert(v.kept);

    const auto lines = read_lines(v.path);
    assert(fields(lines, "n").at(0) == "4");
    assert(close(std::stod(fields(lines, "L").at(0)), packing.L));
    assert(close(std::stod(fields(lines, "k_pair").at(0)), packing.k_pair));
    assert(close(std::stod(fields(lines, "k_wall").at(0)), packing.k_wall));
    assert(close(std::stod(fields(lines, "feasible_tol").at(0)), kTol));
    // One species block (the packing is homogeneous): "s 0 4" then 4 vertex lines.
    assert(fields(lines, "nspecies").at(0) == "1");
    assert(fields(lines, "s").at(0) == "0");  // index
    assert(fields(lines, "s").at(1) == "4");  // vertex count
    assert(count_tag(lines, "v") == 4);

    // One pose per square, in order, matching the packing that was actually claimed.
    std::vector<std::vector<std::string>> poses;
    for (const std::string& line : lines) {
        std::istringstream s(line);
        std::string tag;
        s >> tag;
        std::vector<std::string> f;
        for (std::string w; s >> w;) f.push_back(w);
        if (tag == "p") poses.push_back(f);
    }
    assert(poses.size() == 4);
    for (size_t i = 0; i < poses.size(); ++i) {
        assert(poses[i].size() == 5);  // x, y, theta, frozen, species
        assert(close(std::stod(poses[i][0]), packing.particles[i].x()));
        assert(close(std::stod(poses[i][1]), packing.particles[i].y()));
        assert(close(std::stod(poses[i][2]), packing.particles[i].theta()));
        assert(poses[i][3] == "0");
        assert(poses[i][4] == "0");  // species 0: the packing's only one
    }

    std::filesystem::remove_all(dir);
}

void test_each_record_gets_its_own_file() {
    // Five identical runs saved at once is the shape of a batch finishing together: same n,
    // same side, same second.
    const auto dir = fresh_dir("unique");
    const System packing = two_in_a_row(2.25);

    std::vector<std::string> paths;
    for (int k = 0; k < 5; ++k) {
        const LogVerdict v = run_claiming(dir, packing, 2.25);
        assert(v.kept);
        paths.push_back(v.path);
    }

    std::sort(paths.begin(), paths.end());
    assert(std::unique(paths.begin(), paths.end()) == paths.end());  // no name used twice

    long files = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        (void)entry;
        ++files;
    }
    assert(files == 5);  // and none of them truncated another

    std::filesystem::remove_all(dir);
}

void test_a_log_it_cannot_open_is_reported_not_fatal() {
    // Logging is optional. A path that cannot be created has to be reported rather than crash,
    // because the driver has already printed its answer.
    System sys = grid_init(make_square(1.0), 2);
    const LogVerdict v = run_claiming("/dev/null/not-a-directory", sys, sys.tight_side());  // /dev/null is not a directory
    assert(!v.kept && v.reason.find("no record") != std::string::npos);
}

void test_a_verified_run_is_committed_with_its_side_in_its_name() {
    // The side in the name is what the next run ranks against, so it has to be there, has to
    // be the side that was claimed, and has to name a file that exists.
    const auto dir = fresh_dir("commit");
    const System packing = two_in_a_row(2.25);

    const LogVerdict v = run_claiming(dir, packing, 2.25);
    assert(v.kept);
    assert(v.rank == 1 && v.held == 1);
    assert(v.removed.empty());
    assert(std::filesystem::exists(v.path));
    assert(count_files(dir) == 1);  // renamed, not copied: no plain-named file left behind
    assert(v.path.find("_n2_s2.250000.log") != std::string::npos);
    // The configuration survives the rename: it is the same file, under the name it earned.
    assert(count_tag(read_lines(v.path), "p") == 2);

    std::filesystem::remove_all(dir);
}

void test_only_the_densest_five_logs_survive() {
    const auto dir = fresh_dir("keep5");

    // Five runs, worst first, so each one in turn is the one on the edge.
    for (double side : {2.50, 2.10, 2.40, 2.00, 2.30}) {
        assert(run_claiming(dir, two_in_a_row(side), side).kept);
    }
    assert(kept_sides(dir) == std::vector<double>({2.00, 2.10, 2.30, 2.40, 2.50}));
    assert(count_files(dir) == 5);  // committed records and nothing else

    // A sixth run, denser than the loosest kept: in, and the loosest out.
    const LogVerdict in = run_claiming(dir, two_in_a_row(2.20), 2.20);
    assert(in.kept && in.rank == 3 && in.held == 5);
    assert(in.removed.size() == 1);
    assert(in.removed[0].find("_s2.500000") != std::string::npos);
    assert(kept_sides(dir) == std::vector<double>({2.00, 2.10, 2.20, 2.30, 2.40}));

    // A looser run than all five: dropped, and it takes nothing with it.
    const LogVerdict out = run_claiming(dir, two_in_a_row(2.60), 2.60);
    assert(!out.kept && out.removed.empty());
    assert(out.reason.find("densest") != std::string::npos);
    assert(!std::filesystem::exists(out.path));
    assert(kept_sides(dir) == std::vector<double>({2.00, 2.10, 2.20, 2.30, 2.40}));

    // Equalling the fifth-best is not beating it: the record already on disk stays.
    const LogVerdict tie = run_claiming(dir, two_in_a_row(2.40), 2.40);
    assert(!tie.kept && tie.removed.empty());
    assert(kept_sides(dir) == std::vector<double>({2.00, 2.10, 2.20, 2.30, 2.40}));

    // Another square count keeps its own five: a dense n = 2 run must not crowd out n = 3.
    const System three = in_a_row(3, 3.5);
    assert(run_claiming(dir, three, 3.5).kept);
    assert(kept_sides(dir).size() == 6);
    assert(count_files(dir) == 6);

    std::filesystem::remove_all(dir);
}

void test_an_unverified_result_is_never_saved() {
    // Every one of these would otherwise go on disk as a record, and a record is read later as
    // a fact about how dense a packing of n squares can be.
    const auto dir = fresh_dir("verify");
    const System packing = two_in_a_row(2.20);

    // The run itself certified nothing: a candidate, not a packing.
    assert(!run_claiming(dir, packing, 2.20, false).kept);

    // A side denser than the poses it is claimed for.
    assert(!run_claiming(dir, packing, 2.00, true).kept);

    // Squares that still overlap at the claimed side.
    System overlapping = packing;
    overlapping.particles[1].set_pose(0.9, 0.5, 0.0);
    overlapping.recount();
    assert(overlapping.max_overlap() > kTol);
    assert(!run_claiming(dir, overlapping, overlapping.tight_side(), true).kept);

    // A pose that is not a number. tight_side() is then NaN too, so this has to be caught as
    // the pose it is rather than slipping through a comparison that is false either way.
    System nan_pose = packing;
    nan_pose.particles[0].set_pose(std::numeric_limits<double>::quiet_NaN(), 0.5, 0.0);
    assert(!run_claiming(dir, nan_pose, 2.20, true).kept);

    // Denser than the squares' own area allows. Here it takes an inconsistent shape to force
    // it; upstream it would be a miscounted n or a side from the wrong configuration. Either
    // way the claim is impossible, and impossible is not a record.
    System impossible = packing;
    impossible.species[0].area = 100.0;
    assert(!run_claiming(dir, impossible, 2.20, true).kept);

    // Nothing above reached the disk: the directory holds no files at all.
    assert(count_files(dir) == 0);

    std::filesystem::remove_all(dir);
}

void test_the_keep_all_override_deletes_nothing() {
    // The escape hatch for the runs where the failure is the interesting part.
    const auto dir = fresh_dir("keepall");
    const System packing = two_in_a_row(2.20);
    setenv(kKeepAllEnv, "1", 1);

    const LogVerdict uncertified = run_claiming(dir, packing, 2.20, false);
    assert(uncertified.kept && std::filesystem::exists(uncertified.path));

    for (double side : {2.00, 2.10, 2.30, 2.40, 2.50, 2.60}) {
        const LogVerdict v = run_claiming(dir, two_in_a_row(side), side);
        assert(v.kept && std::filesystem::exists(v.path));
        assert(v.removed.empty());  // the override deletes nothing, ever
    }
    // Five committed, because the override does not change what ranks -- only what is
    // deleted. The sixth run and the unverified one keep their records under their plain
    // names, where nothing can mistake them for records.
    assert(kept_sides(dir) == std::vector<double>({2.00, 2.10, 2.30, 2.40, 2.50}));
    assert(count_files(dir) == 7);
    assert(uncertified.path.find("_s") == std::string::npos);
    assert(uncertified.reason.find(kKeepAllEnv) != std::string::npos);

    unsetenv(kKeepAllEnv);
    std::filesystem::remove_all(dir);
}

void test_a_saved_record_loads_back_exactly() {
    const auto dir = fresh_dir("roundtrip");
    System packing = in_a_row(4, 5.0);
    packing.k_pair = 7.0;
    packing.k_wall = 11.0;
    packing.particles[0].frozen = true;
    packing.set_pose(1, packing.particles[1].x(), 0.7, 0.1);  // an angle with no short decimal
    packing.recount();

    const LogVerdict v = run_claiming(dir, packing, 5.0);
    assert(v.kept);

    std::string why;
    const std::optional<System> loaded = load_record(v.path, &why);
    assert(loaded && why.empty());
    assert(loaded->L == packing.L);
    assert(loaded->k_pair == 7.0 && loaded->k_wall == 11.0);
    assert(loaded->species.size() == 1);
    assert(close(loaded->base_shape().area, 1.0));
    assert(loaded->particles.size() == 4);
    for (size_t i = 0; i < 4; ++i) {
        // Bit-identical: a start that differs in the last digit is a different trajectory.
        assert(loaded->particles[i].x() == packing.particles[i].x());
        assert(loaded->particles[i].y() == packing.particles[i].y());
        assert(loaded->particles[i].theta() == packing.particles[i].theta());
        assert(loaded->particles[i].frozen == packing.particles[i].frozen);
    }
    // Ready to use straight away: the energy and the neighbour grid are built.
    assert(close(loaded->energy(), loaded->total_energy()));
    assert(close(loaded->particle_energy(1), packing.particle_energy(1)));
    assert(loaded->is_feasible(kTol));

    std::filesystem::remove_all(dir);
}

void test_a_malformed_record_is_refused_with_a_reason() {
    const auto dir = fresh_dir("malformed");
    std::filesystem::create_directories(dir);
    const std::string path = (dir / "x.cfg").string();
    const auto write = [&](const std::string& text) { std::ofstream(path) << text; };

    const std::string good =
        "# comment\nn 2\nL 3\nk_pair 100\nk_wall 100\nfeasible_tol 0\nnspecies 1\ns 0 4\n"
        "v -0.5 -0.5\nv 0.5 -0.5\nv 0.5 0.5\nv -0.5 0.5\np 0.5 0.5 0 0 0\np 1.5 0.5 0 0 0\n";
    write(good);
    assert(load_record(path));

    const struct { const char* from; const char* to; } broken[] = {
        {"n 2\n", "n 3\n"},                                // more squares promised than listed
        {"n 2\n", "n two\n"},
        {"n 2\n", ""},
        {"L 3\n", "L -3\n"},
        {"L 3\n", "L nan\n"},
        {"L 3\n", ""},
        {"k_wall 100", "k_wall -1"},
        {"k_pair 100\n", ""},
        {"nspecies 1", "nspecies 2"},
        {"s 0 4", "s 1 4"},
        {"s 0 4", "s 0 2"},
        {"s 0 4", "s 0 5"},                                // the block runs into the particles
        {"v 0.5 0.5\n", ""},
        {"v 0.5 -0.5\nv 0.5 0.5\nv -0.5 0.5\n", "v -0.5 0.5\nv 0.5 0.5\nv 0.5 -0.5\n"},  // clockwise
        {"p 1.5 0.5 0 0 0", "p 1.5 0.5 0 0 1"},            // no such species
        {"p 1.5 0.5 0 0 0", "p 1.5 0.5 0 2 0"},            // frozen is a flag
        {"p 1.5 0.5 0 0 0", "p inf 0.5 0 0 0"},
        {"p 1.5 0.5 0 0 0", "p 1.5 0.5 0 0"},
        {"p 1.5 0.5 0 0 0", "q 1.5 0.5 0 0 0"},            // unknown line
        {"p 1.5 0.5 0 0 0\n", ""},                         // fewer squares than promised
    };
    for (const auto& b : broken) {
        std::string text = good;
        const size_t at = text.find(b.from);
        assert(at != std::string::npos);
        text.replace(at, std::string(b.from).size(), b.to);
        write(text);
        std::string why;
        assert(!load_record(path, &why));
        assert(why.find(path) == 0);
    }

    write("");
    assert(!load_record(path));
    std::string why;
    assert(!load_record((dir / "absent.cfg").string(), &why));
    assert(why.find("cannot open") != std::string::npos);

    std::filesystem::remove_all(dir);
}

void test_a_trajectory_holds_one_line_per_frame() {
    const auto dir = fresh_dir("trajectory");
    const std::string path = (dir / "sub" / "n0004.traj").string();
    const System s = in_a_row(4, 5.0);

    std::FILE* f = open_trajectory(path, s);  // creates the folder too
    assert(f);
    append_frame(f, 0, s);
    append_frame(f, 7, s);
    std::fclose(f);

    const auto lines = read_lines(path);
    assert(count_tag(lines, "p") == 4);  // the starting configuration, for the viewer's shapes
    assert(count_tag(lines, "f") == 2);
    const auto first = fields(lines, "f");
    assert(first.size() == 3 + 3 * 4);
    assert(first[0] == "0");
    assert(close(std::stod(first[1]), s.L));
    assert(close(std::stod(first[3]), s.particles[0].x(), 1e-8));
    assert(!load_record(path));  // progress, not a record

    // A new run at the same path starts over.
    f = open_trajectory(path, s);
    assert(f);
    std::fclose(f);
    assert(count_tag(read_lines(path), "f") == 0);

    assert(open_trajectory("/dev/null/not-a-directory/x.traj", s) == nullptr);
    std::filesystem::remove_all(dir);
}

}  // namespace

int main() {
    test_a_trajectory_holds_one_line_per_frame();
    test_a_saved_record_loads_back_exactly();
    test_a_malformed_record_is_refused_with_a_reason();
    test_a_kept_record_holds_the_certified_configuration();
    test_each_record_gets_its_own_file();
    test_a_log_it_cannot_open_is_reported_not_fatal();
    test_a_verified_run_is_committed_with_its_side_in_its_name();
    test_only_the_densest_five_logs_survive();
    test_an_unverified_result_is_never_saved();
    test_the_keep_all_override_deletes_nothing();
    std::printf("io: all tests passed\n");
    return 0;
}
