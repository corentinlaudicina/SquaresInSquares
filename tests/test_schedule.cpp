// Schedule parsing: defaults, unknown keys, types, curves, the resolved round trip.

#include "engine/schedule.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>

using namespace pack;
using json = nlohmann::json;

namespace {

ParsedSchedule ok(const json& j) {
    std::string why;
    auto ps = parse_schedule(j, &why);
    if (!ps) std::printf("unexpected error: %s\n", why.c_str());
    assert(ps);
    return *ps;
}

std::string err(const json& j) {
    std::string why;
    assert(!parse_schedule(j, &why));
    return why;
}

void test_empty_is_defaults() {
    const EngineConfig d = shrink_defaults(), c = ok(json::object()).config;
    assert(c.mode == "shrink" && c.mirror == d.mirror && c.eps == d.eps && c.kick.shake_T == d.kick.shake_T);
    assert(c.kick.lnm.insert_tries == 12 && !c.update);
    const EngineConfig cy = ok({{"mode", "cycle"}}).config;
    assert(cy.prepare_ds == cycle_defaults().prepare_ds && cy.mode == "cycle");
}

void test_values_and_types() {
    const EngineConfig c = ok({{"start", {{"mirror", false}}}, {"kick", {{"shake", {{"sweeps", 7.0}}}}}}).config;
    assert(!c.mirror && c.kick.shake_sweeps == 7);
    assert(err({{"start", {{"mirror", 1}}}}).find("true or false") != std::string::npos);
    assert(err({{"kick", {{"shake", {{"T", "hot"}}}}}}).find("number or a curve") != std::string::npos);
    assert(err({{"nope", 1}}).find("unknown schedule key nope") != std::string::npos);
    assert(err({{"mode", "anneal"}}).find("mode") != std::string::npos);
    assert(err({{"shrink", {{"on_fail", 1.5}}}}).find("on_fail") != std::string::npos);
}

void test_curves() {
    ParsedSchedule ps = ok({{"kick", {{"shake", {{"T", {{0.0, 0.1}, {0.5, 0.01}, {1.0, 0.01}}}}}}}});
    assert(ps.config.kick.shake_T == 0.1 && ps.config.update);
    EngineConfig c = ps.config;
    ps.config.update(0.25, c);
    assert(std::abs(c.kick.shake_T - 0.055) < 1e-15);
    ps.config.update(0.9, c);
    assert(c.kick.shake_T == 0.01);
    assert(err({{"prepare", {{"ds", {{0.0, 0.1}, {1.0, 0.2}}}}}}).find("cannot be a curve") != std::string::npos);
    assert(err({{"cycle", {{"ds", {{0.5, 0.1}, {0.2, 0.2}}}}}}).find("ascending") != std::string::npos);
}

void test_resolved_round_trip() {
    const ParsedSchedule a = ok({{"name", "x"}, {"kick", {{"shake", {{"T", {{0.0, 0.1}, {1.0, 0.0}}}}}}}});
    const ParsedSchedule b = ok(a.resolved);
    assert(a.hash == b.hash && a.resolved == b.resolved);
    assert(ok({{"name", "x"}}).hash != a.hash);
}

}  // namespace

int main() {
    test_empty_is_defaults();
    test_values_and_types();
    test_curves();
    test_resolved_round_trip();
    std::printf("schedule: all tests passed\n");
}
