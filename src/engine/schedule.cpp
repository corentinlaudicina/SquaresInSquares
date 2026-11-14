#include "engine/schedule.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <memory>

namespace pack {
namespace {

using json = nlohmann::json;

// One entry per schedule key. Numbers go through get/set as doubles (integers are rounded on set);
// `curve` marks the ones read per kick, so a curve on them means something.
struct Field {
    const char* path;
    enum Kind { Number, Integer, Bool } kind;
    bool curve;
    std::function<double(const EngineConfig&)> get;
    std::function<void(EngineConfig&, double)> set;
};

#define NUM(path, member, curve) \
    Field{path, Field::Number, curve, [](const EngineConfig& c) { return static_cast<double>(c.member); }, \
          [](EngineConfig& c, double v) { c.member = v; }}
#define INT(path, member, curve) \
    Field{path, Field::Integer, curve, [](const EngineConfig& c) { return static_cast<double>(c.member); }, \
          [](EngineConfig& c, double v) { c.member = static_cast<decltype(c.member)>(std::llround(v)); }}
#define BOOL(path, member) \
    Field{path, Field::Bool, false, [](const EngineConfig& c) { return c.member ? 1.0 : 0.0; }, \
          [](EngineConfig& c, double v) { c.member = v != 0.0; }}

const std::vector<Field>& fields() {
    static const std::vector<Field> f = {
        BOOL("start.mirror", mirror),
        NUM("prepare.ds", prepare_ds, false),
        INT("prepare.max_relax", prepare_max_relax, false),
        NUM("kick.large_prob.base", kick.large_base, true),
        NUM("kick.large_prob.max", kick.large_max, true),
        NUM("kick.large_prob.step", kick.large_step, true),
        INT("kick.large_prob.every", kick.escalate_every, true),
        INT("kick.shake.sweeps", kick.shake_sweeps, true),
        NUM("kick.shake.T", kick.shake_T, true),
        NUM("kick.shake.translate", kick.shake_translate, true),
        NUM("kick.shake.rotate", kick.shake_rotate, true),
        NUM("kick.shake.mix.translate", kick.mix.translate, true),
        NUM("kick.shake.mix.rotate", kick.mix.rotate, true),
        NUM("kick.shake.mix.flip", kick.mix.flip, true),
        NUM("kick.shake.mix.shove", kick.mix.shove, true),
        NUM("kick.lnm.weights.destroy", kick.lnm.w_destroy, true),
        NUM("kick.lnm.weights.strip", kick.lnm.w_strip, true),
        NUM("kick.lnm.weights.corner", kick.lnm.w_corner, true),
        NUM("kick.lnm.weights.loosen", kick.lnm.w_loosen, true),
        NUM("kick.lnm.weights.symmetrize", kick.lnm.w_symmetrize, true),
        INT("kick.lnm.relax_iterations", kick.lnm.relax_iterations, true),
        INT("kick.lnm.insert_tries", kick.lnm.insert_tries, true),
        NUM("kick.lnm.snap_prob", kick.lnm.snap_prob, true),
        NUM("kick.lnm.loosen.grow", kick.lnm.loosen_grow, true),
        INT("kick.lnm.loosen.sweeps", kick.lnm.loosen_sweeps, true),
        NUM("kick.lnm.loosen.soften", kick.lnm.loosen_soften, true),
        NUM("kick.lnm.loosen.soft_prob", kick.lnm.loosen_soft_prob, true),
        NUM("kick.lnm.regions.disc", kick.lnm.regions.disc, true),
        NUM("kick.lnm.regions.band", kick.lnm.regions.band, true),
        NUM("kick.lnm.regions.corner", kick.lnm.regions.corner, true),
        NUM("kick.lnm.regions.axis_weight", kick.lnm.regions.axis_weight, true),
        NUM("kick.lnm.regions.diagonal_weight", kick.lnm.regions.diagonal_weight, true),
        NUM("shrink.eps", eps, false),
        NUM("shrink.eps_min", eps_min, false),
        NUM("shrink.eps_max", eps_max, true),
        NUM("shrink.on_success", eps_up, true),
        NUM("shrink.on_fail", eps_down, true),
        NUM("cycle.ds", cycle_ds, true),
        INT("cycle.max_relax", cycle_max_relax, true),
        NUM("cycle.side_T", side_T, true),
        INT("cycle.restart_after", restart_after, true),
        BOOL("propagate.enabled", propagate),
        INT("propagate.every", propagate_every, true),
    };
    return f;
}

#undef NUM
#undef INT
#undef BOOL

struct Curve {
    std::vector<std::pair<double, double>> pts;
    double at(double t) const {
        if (t <= pts.front().first) return pts.front().second;
        for (size_t i = 1; i < pts.size(); ++i) {
            if (t <= pts[i].first) {
                const auto [t0, v0] = pts[i - 1];
                const auto [t1, v1] = pts[i];
                return t1 > t0 ? v0 + (v1 - v0) * (t - t0) / (t1 - t0) : v1;
            }
        }
        return pts.back().second;
    }
};

std::optional<Curve> as_curve(const json& v) {
    if (!v.is_array() || v.empty()) return std::nullopt;
    Curve c;
    for (const json& p : v) {
        if (!p.is_array() || p.size() != 2 || !p[0].is_number() || !p[1].is_number()) return std::nullopt;
        const double t = p[0].get<double>();
        if (!(t >= 0.0 && t <= 1.0) || (!c.pts.empty() && t < c.pts.back().first)) return std::nullopt;
        c.pts.emplace_back(t, p[1].get<double>());
    }
    return c;
}

// Flattens nested objects into "a.b.c" -> leaf.
void flatten(const json& j, const std::string& prefix, std::vector<std::pair<std::string, json>>& out) {
    for (const auto& [k, v] : j.items()) {
        const std::string path = prefix.empty() ? k : prefix + "." + k;
        if (v.is_object()) flatten(v, path, out);
        else out.emplace_back(path, v);
    }
}

std::string fnv1a(const std::string& s) {
    unsigned long long h = 1469598103934665603ULL;
    for (const unsigned char ch : s) {
        h ^= ch;
        h *= 1099511628211ULL;
    }
    char buf[17];
    std::snprintf(buf, sizeof buf, "%016llx", h);
    return buf;
}

}  // namespace

std::vector<std::string> field_names() {
    std::vector<std::string> out = {"schema", "name", "description", "mode (shrink | cycle)"};
    for (const Field& f : fields()) out.push_back(std::string(f.path) + (f.curve ? " (curve)" : ""));
    return out;
}

std::optional<ParsedSchedule> parse_schedule(const json& doc, std::string* error) {
    const auto fail = [&](const std::string& why) -> std::optional<ParsedSchedule> {
        if (error) *error = why;
        return std::nullopt;
    };
    if (!doc.is_object()) return fail("a schedule is a JSON object");
    if (doc.contains("schema") && doc["schema"] != "sqpack.schedule/1") return fail("schema is not sqpack.schedule/1");
    const std::string mode = doc.value("mode", std::string("shrink"));
    if (mode != "shrink" && mode != "cycle") return fail("mode must be shrink or cycle");

    ParsedSchedule ps;
    ps.config = mode == "cycle" ? cycle_defaults() : shrink_defaults();
    ps.name = doc.value("name", std::string("unnamed"));

    std::vector<std::pair<std::string, json>> leaves;
    json body = doc;
    for (const char* k : {"schema", "name", "description", "mode"}) body.erase(k);
    flatten(body, "", leaves);

    auto curves = std::make_shared<std::vector<std::pair<const Field*, Curve>>>();
    for (const auto& [path, v] : leaves) {
        const Field* f = nullptr;
        for (const Field& g : fields()) {
            if (path == g.path) f = &g;
        }
        if (!f) return fail("unknown schedule key " + path);
        if (f->kind == Field::Bool) {
            if (!v.is_boolean()) return fail(path + " must be true or false");
            f->set(ps.config, v.get<bool>() ? 1.0 : 0.0);
        } else if (v.is_number()) {
            const double x = v.get<double>();
            if (!std::isfinite(x)) return fail(path + " is not finite");
            f->set(ps.config, x);
        } else if (auto c = as_curve(v)) {
            if (!f->curve) return fail(path + " cannot be a curve: it is read once, before the loop");
            f->set(ps.config, c->at(0.0));
            curves->emplace_back(f, *c);
        } else {
            return fail(path + " must be a number or a curve [[t, v], ...] with t ascending in [0, 1]");
        }
    }
    if (ps.config.prepare_ds <= 0.0 || ps.config.cycle_ds <= 0.0) return fail("ds must be positive");
    if (ps.config.eps_up <= 1.0 || ps.config.eps_down >= 1.0) return fail("shrink.on_success must be > 1 and on_fail < 1");

    if (!curves->empty()) {
        ps.config.update = [curves](double t, EngineConfig& c) {
            for (const auto& [f, curve] : *curves) f->set(c, curve.at(t));
        };
    }

    // The resolved schedule: every field's effective value, curves as given.
    ps.resolved = {{"schema", "sqpack.schedule/1"}, {"name", ps.name}, {"mode", mode}};
    for (const Field& f : fields()) {
        json* node = &ps.resolved;
        std::string p = f.path;
        for (size_t dot; (dot = p.find('.')) != std::string::npos; p = p.substr(dot + 1)) node = &(*node)[p.substr(0, dot)];
        json value;
        if (f.kind == Field::Bool) value = f.get(ps.config) != 0.0;
        else if (f.kind == Field::Integer) value = std::llround(f.get(ps.config));
        else value = f.get(ps.config);
        for (const auto& [g, c] : *curves) {
            if (g == &f) {
                value = json::array();
                for (const auto& [t, v] : c.pts) value.push_back({t, v});
            }
        }
        (*node)[p] = value;
    }
    ps.hash = fnv1a(ps.resolved.dump());
    return ps;
}

std::optional<ParsedSchedule> load_schedule(const std::string& path, std::string* error) {
    std::ifstream in(path);
    if (!in) {
        if (error) *error = path + ": cannot open";
        return std::nullopt;
    }
    try {
        auto ps = parse_schedule(json::parse(in), error);
        if (!ps && error) *error = path + ": " + *error;
        return ps;
    } catch (const json::exception& e) {
        if (error) *error = path + ": " + e.what();
        return std::nullopt;
    }
}

}  // namespace pack
