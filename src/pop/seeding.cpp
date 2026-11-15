#include "pop/seeding.hpp"

#include "engine/solution_io.hpp"
#include "init.hpp"
#include "legalize.hpp"
#include "moves/relax.hpp"
#include "moves/symmetry.hpp"
#include "pop/backend.hpp"
#include "pop/moves.hpp"

#include <filesystem>

namespace pack {
namespace {

std::optional<System> neighbour(const SeedParams& prm, int m) {
    if (m < 1) return std::nullopt;
    const std::string path = prm.neighbor_dir + "/" + std::to_string(m) + ".json";
    if (!std::filesystem::exists(path)) return std::nullopt;
    auto s = read_solution(path);
    if (!s || static_cast<int>(s->particles.size()) != m) return std::nullopt;
    return s;
}

void add_legal(std::vector<SeedCandidate>& out, const std::string& source, System sys) {
    if (legalize(sys).ok) out.push_back({source, std::move(sys)});
}

}  // namespace

System insert_square(const System& sys, int tries, std::mt19937_64& rng) {
    Packing p = to_packing(sys);
    p.xyt.insert(p.xyt.end(), {0.0, 0.0, 0.0});
    const size_t k = p.xyt.size() - 3;
    Workspace ws;
    double best = std::numeric_limits<double>::infinity();
    std::array<double, 3> pose{0.5 * p.h, 0.5 * p.h, 0.0};
    for (int t = 0; t < std::max(1, tries); ++t) {
        p.xyt[k] = uniform(rng, 0.5, p.h - 0.5);
        p.xyt[k + 1] = uniform(rng, 0.5, p.h - 0.5);
        p.xyt[k + 2] = u01(rng) < 0.5 ? 0.25 * kPi * uniform_int(rng, 0, 1) : uniform(rng, 0.0, 0.5 * kPi);
        const double c = cost_grad(p, Penalty::Depth2, ws, nullptr).cost;
        if (c < best) {
            best = c;
            pose = {p.xyt[k], p.xyt[k + 1], p.xyt[k + 2]};
        }
    }
    p.xyt[k] = pose[0];
    p.xyt[k + 1] = pose[1];
    p.xyt[k + 2] = pose[2];
    return to_system(p);
}

std::vector<SeedCandidate> seed_candidates(int n, const SeedParams& prm, std::mt19937_64& rng) {
    std::vector<SeedCandidate> out;
    if (prm.up) {
        if (auto up = neighbour(prm, n + 1)) {
            if (auto r = delete_rattler(*up)) add_legal(out, "up:rattler", *r);
            int k = 0;
            for (System& s : delete_for_compress(*up, prm.up_count)) add_legal(out, "up:contact" + std::to_string(++k), std::move(s));
        }
    }
    if (prm.down) {
        if (auto down = neighbour(prm, n - 1)) {
            // The inserted square needs room the legalizer's small growth cannot give: compress()
            // grows the box until it certifies, then squeezes it back down.
            const CompressResult cr = compress(insert_square(*down, prm.down_tries, rng));
            if (cr.found) add_legal(out, "down", cr.packing);
        }
    }
    if (prm.lattice) {
        const Shape sq = make_square(1.0);
        add_legal(out, "lattice:grid", grid_init(sq, n));
        add_legal(out, "lattice:diagonal", diagonal_init(sq, n));
    }
    return out;
}

}  // namespace pack
