#include "pop/relax_capped.hpp"
#include "pop/symmetric.hpp"

#include <algorithm>
#include <cmath>

namespace pack {
namespace {

constexpr double kR = 0.7071067811865476;  // theta is scaled by the circumradius

double dot(const double* a, const double* b, size_t n) {
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

bool within(const CostInfo& c, double tol) { return c.max_depth <= tol && c.max_wall <= tol; }

}  // namespace

// The relaxation over a variable vector x (x, y, theta per square), whatever it encodes. `raw(g)`
// returns the cost at x and writes its gradient with respect to x into g.
template <class Raw>
CappedStats relax_core(std::vector<double>& x, Raw raw, const CappedParams& prm, Workspace& ws) {
    // The neighbour list is valid for one relaxation only: switched on here, off on every return.
    struct VerletScope {
        Workspace& w;
        explicit VerletScope(Workspace& ws_) : w(ws_) { w.verlet = true; w.verlet_valid = false; }
        ~VerletScope() { w.verlet = false; w.verlet_valid = false; }
    } verlet_scope(ws);
    const size_t N3 = x.size();
    const size_t m = static_cast<size_t>(std::max(1, prm.history));
    ws.g.resize(N3);
    ws.g_old.resize(N3);
    ws.d.resize(N3);
    ws.best = x;
    ws.hist_s.resize(m * N3);
    ws.hist_y.resize(m * N3);
    ws.rho.resize(m);
    ws.alpha.resize(m);

    // Gradient in the scaled variables (x, y, R theta).
    const auto eval = [&](double* g) {
        const CostInfo ci = raw(g);
        for (size_t i = 2; i < N3; i += 3) g[i] /= kR;
        return ci;
    };

    CappedStats st;
    CostInfo cur = eval(ws.g.data());
    CostInfo best = cur;
    if (within(cur, prm.tol)) {
        st.reached = true;
        st.info = cur;
        return st;
    }

    for (const RelaxStage& stage : prm.stages) {
        size_t k = 0, head = 0;  // history pairs held, ring position of the next one
        for (int it = 0; it < stage.iters; ++it) {
            // d = -H g, two-loop recursion over the ring, newest first.
            double* d = ws.d.data();
            for (size_t i = 0; i < N3; ++i) d[i] = -ws.g[i];
            for (size_t t = 0; t < k; ++t) {
                const size_t r = (head + m - 1 - t) % m;
                ws.alpha[r] = ws.rho[r] * dot(&ws.hist_s[r * N3], d, N3);
                for (size_t i = 0; i < N3; ++i) d[i] -= ws.alpha[r] * ws.hist_y[r * N3 + i];
            }
            if (k > 0) {
                const size_t r = (head + m - 1) % m;
                const double gamma = dot(&ws.hist_s[r * N3], &ws.hist_y[r * N3], N3) /
                                     dot(&ws.hist_y[r * N3], &ws.hist_y[r * N3], N3);
                for (size_t i = 0; i < N3; ++i) d[i] *= gamma;
            }
            for (size_t t = k; t-- > 0;) {
                const size_t r = (head + m - 1 - t) % m;
                const double beta = ws.rho[r] * dot(&ws.hist_y[r * N3], d, N3);
                for (size_t i = 0; i < N3; ++i) d[i] += (ws.alpha[r] - beta) * ws.hist_s[r * N3 + i];
            }
            if (k == 0 || !(dot(d, ws.g.data(), N3) < 0.0)) {  // no history, or uphill: steepest descent
                k = 0;
                for (size_t i = 0; i < N3; ++i) d[i] = -ws.g[i];
            }
            // The cap: no variable moves more than stage.cap. Steepest descent gets exactly the cap,
            // because its raw length is the gradient's, which is meaningless near a contact.
            double big = 0.0;
            for (size_t i = 0; i < N3; ++i) big = std::max(big, std::abs(d[i]));
            if (!(big > 0.0)) break;
            const double scale = (k == 0 || big > stage.cap) ? stage.cap / big : 1.0;
            for (size_t i = 0; i < N3; ++i) {
                d[i] *= scale;
                x[i] += (i % 3 == 2) ? d[i] / kR : d[i];
            }

            std::swap(ws.g, ws.g_old);
            const CostInfo next = eval(ws.g.data());
            ++st.iters;
            // History pair (s = d, y = g - g_old), kept only with positive curvature.
            double* hs = &ws.hist_s[head * N3];
            double* hy = &ws.hist_y[head * N3];
            for (size_t i = 0; i < N3; ++i) {
                hs[i] = d[i];
                hy[i] = ws.g[i] - ws.g_old[i];
            }
            const double sy = dot(hs, hy, N3);
            if (next.cost > cur.cost) {
                k = 0;  // the model broke: drop it, the cap keeps the next step sane
            } else if (sy > 1e-300) {
                ws.rho[head] = 1.0 / sy;
                head = (head + 1) % m;
                k = std::min(k + 1, m);
            }
            cur = next;
            if (cur.cost < best.cost) {
                best = cur;
                ws.best = x;
            }
            if (within(cur, prm.tol)) {
                st.reached = true;
                st.info = cur;
                return st;
            }
        }
    }
    x = ws.best;
    st.info = best;
    st.reached = within(best, prm.tol);
    return st;
}

CappedStats relax_capped(Packing& p, const CappedParams& prm, Workspace& ws) {
    return relax_core(p.xyt, [&](double* g) { return cost_grad(p, prm.penalty, ws, g); }, prm, ws);
}

CappedStats relax_symmetric(Packing& gen, const SymLayout& L, const CappedParams& prm, Workspace& ws) {
    if (L.sym == Sym::None) return relax_capped(gen, prm, ws);
    return relax_core(gen.xyt, [&](double* g) { return cost_symmetric(gen, L, prm.penalty, ws, g); }, prm, ws);
}

}  // namespace pack
