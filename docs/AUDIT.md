# Phase 1 audit — squares-in-squares engine

Audit date 2026-11-14. Read-only: nothing in the repo was changed except this file. Builds and
test runs went to a scratch directory (`make BUILD=<scratch>`), and probe runs used `--no-save`
with the CSV redirected, so `build/`, `data/` and `logs/` are untouched.

Line references are `file:line` as of this audit.

---

## 1. Build and run today

| Item | State |
|---|---|
| Build system | One hand-written `Makefile`, `clang++`, header dependencies generated automatically (`-MMD -MP`). |
| Flags (release) | `-std=c++20 -Wall -Wextra -Wpedantic -O2 -pthread` |
| Flags (debug) | `-O0 -g -fsanitize=address,undefined` (tests only) |
| Targets | `make all` = 7 test binaries + `squares` (main.cpp) + `squares_shrink` + `squares_compress`. `make check` builds and runs the tests. `make debug` = sanitizer tests. |
| Toolchain seen | Apple clang 21, arm64. A full `make -j8 all` takes about 8 s. |
| Tests | All 7 suites pass (`test_symmetry` about 13 s, `test_lnm` about 4 s, the rest about 1 s each). |
| Warnings | 9, all in `src/main.cpp` (unused `kick`, `best_side`, …): see bug B2. |

Running one instance today (it must be run from the repo root, because paths are relative to the working directory):

```sh
make shrink compress
./build/squares_shrink   --n 29 --seed 7 --kicks 40  [--cpu 60] [--no-save] [--csv out.csv] [--eps 1e-3] [--symmetrize w] [--propagate]
./build/squares_compress --n 29 --seed 7 --cycles 5  [--cpu 60] [--no-save] [--csv out.csv] [--T 0.01]  [--symmetrize w] [--propagate]
```

Each run loads `data/init_external/n%04d.cfg`, prints progress and a machine-readable
`result n=… seed=… cycles=… cpu=… start=… best_cycle=… final=… score=…` line, writes
`logs/moves_n%04d.csv` (even with `--no-save`), and unless `--no-save` is given also writes
`logs/n%04d.traj` and calls `save_record()` into `data/best/`.

Measured: `squares_shrink --n 29 --kicks 40` takes about 10 CPU s. `squares_compress --n 29 --cycles 5` takes about 8.5 CPU s.

**Apple-silicon suitability.** No x86 intrinsics or arch-specific code anywhere (checked with grep for
`immintrin`, `_mm_`, `__x86`). `-O2` without `-mcpu`, `-march` or LTO is portable but leaves
performance on the table. Suggested split, both settable from the Makefile:
- `PORTABLE` (default for results that must replay): `-O2 -ffp-contract=off`. clang on arm64
  fuses multiply-adds by default, while GCC on x86-64 without `-march` does not, so the same seed gives
  different floating-point results on the Mac and on Kaggle.
- `FAST`: `-O3 -mcpu=native -flto` for local throughput. Not bit-reproducible across machines.
- No `-ffast-math`. Keep it that way: the tolerance logic depends on IEEE behaviour.
- `CXX := clang++` is hard-set. Kaggle images have g++; `make CXX=g++` should work (C++20
  designated initializers are used, for example in `relax.hpp`'s `RelaxOptions relax{.hopeless_window=200}`) but this is untested.

## 2. Configuration today

There is no config file. It is three layers:
1. **`constexpr` constants at the top of each driver** (`main_shrink.cpp:39-64`,
   `main_compress.cpp:39-62`, `main.cpp:26-44`): kick probabilities, escalation, shake
   temperature and steps, compression `ds`, eps adaptation, restart patience. Changing one means editing and rebuilding.
2. **A small hand-written argv parser in each driver** (duplicated three times): `--n --seed --kicks/--cycles
   --cpu --no-save --csv --symmetrize --propagate --eps/--T`. Numbers are read with `strtod`, so
   the seed is parsed as a double (bug B9).
3. **Option structs with defaults** in the library (`RelaxOptions`, `CompressOptions`,
   `LnmOptions`, `RegionWeights`, `SymmetrizeOptions`, `PropagateOptions`, `StepSizes`,
   `Schedule`). These are already data-shaped and are the natural target for a schedule file.
   There are also file-scope constants not reachable from any struct: the MC move mix
   `kWTranslate/kWRotate/kWFlip/kWShove = .45/.45/.06/.04` and `kSweepRegions` (`mc.cpp:11-17`),
   `attempts = 100` in `large_neighbour_move` (`lnm.cpp:356`), the line-search constants in `relax()`.

## 3. Search components

### 3.1 Drivers (the "protocols")

| Driver | Loop | Already schedule-like |
|---|---|---|
| `main_shrink.cpp` (newest) | **Fixed-target ratchet**: target = best·(1−eps), rescale best into the target, kick, `relax()` once at the fixed box, accept if feasible, then `fit_box`. eps adapts ×1.2 on success and ×0.95 on failure, clamped to [1e-7, 0.05]. No restarts. | eps adaptation, large-move escalation (+0.1 per 25 failures, capped at 0.6) |
| `main_compress.cpp` | **Kick → free relax → compress cycle**: kick, `relax_free` (walls off, then refit the box), `anneal()` (full compress with ds=0.005), Metropolis on side with T = `side_t`·L/√n, restart to best after 200 stagnant cycles. | escalation, restart rule, side temperature |
| `main.cpp` | Older version of `main_compress`. **Truncated**: the cycle loop is gone (bug B2). | — |

Both working drivers begin identically: load the seed, **`mirror_symmetry()`**, relax, then `anneal(ds, 5000)`. See bug B1: this step destroys the seed.

### 3.2 Kicks (`kick()` in each driver)
- **Shake**, with probability 1−p_large: 100 MC `sweep()`s at T=0.05 with translate and rotate steps of 0.1.
- **Large move**, with probability p_large: `large_neighbour_move()` = `lnm_round` at T=1e9 (always accept), up to 100 attempts until one applies.

### 3.3 MC moves (`src/moves/mc.{hpp,cpp}`)
`sweep()` makes n attempts, each drawn from the mix above: `try_translate`, `try_rotate`, `try_flip`, `try_shove_region`. Acceptance is Metropolis on `particle_energy` (cell list, O(1) per move). `tune_step` adapts the step sizes, but the drivers build a fresh `StepSizes` every kick, so the tuning is discarded (B8). There are also group moves (`try_translate/rotate/slide/reflect_group`), `try_box` (NPT box move) and MC `anneal(Schedule)`. The drivers use none of these.

### 3.4 LNS moves (`src/moves/lnm.{hpp,cpp}`). Pattern: copy → perturb hard → `relax()` (300 iterations, tolerance `kFeasibleTol`=1e-5) → Metropolis on enthalpy
| Move | What | Weight | Params |
|---|---|---|---|
| `try_destroy_repair` | remove squares in a region (disc, band or corner), reinsert each as the best of `insert_tries` random poses (45° snap with probability 0.5) | 0.4 | `RegionWeights`, `insert_tries=12`, `snap_prob=.5` |
| `try_strip` | re-lay a band as a 45° tilted strip ("carve") or as an axis-aligned one ("straighten") | 0.3 | band width and direction weights |
| `try_reshape_corner` | re-lay a corner block as a random staircase | 0.2 | `corner_depth_lo/hi` |
| `try_loosen` | grow the box 3%, shake 20 sweeps (half the time with stiffness ×0.1), squeeze back | 0.1 | `loosen_*` |
| `try_symmetrize_region` (`symmetry.cpp`) | copy a region through a D4 element, repair the seams, own compress; arm chosen by a bandit | 0 (off) | `SymmetrizeOptions` |

### 3.5 Relaxation and compression (`src/moves/relax.{hpp,cpp}`)
- `relax()`: L-BFGS over (x, y, R·θ) with history `clamp(3·n_mobile, 10, 100)` and Armijo backtracking
  (up to 60 halvings, no Wolfe test). Stops on Feasible, Converged, Stalled, MaxIterations, Hopeless or Escaped.
  `hopeless_window` extrapolates the log-energy decay and gives up early.
- `compress()`: phase 1 grows the box ×1.05 until feasible. Phase 2 shrinks from `tight_side` by `step` and
  halves the step on failure, until `step ≤ max(min_shrink, √tol/side)` or `max_relaxations` is reached.
- `anneal(sys, ds, max_iter)`: wraps `compress` with step = ds/L. Unrelated to MC `anneal`, despite the name.
- `relax_free()`: sets `k_wall=0`, relaxes, then `fit_box()`.
- `compress_best_of()`: a threaded best-of over candidates. **Not called anywhere.**

### 3.6 Cross-N seeding (already exists, `symmetry.hpp:127-188`)
`symmetric_grow` (n+4 or n+1 from C4), `row_duplicate` and `row_remove`, `frame_grow` (L-shaped layer),
`delete_rattler`, `delete_for_compress` (n−1). `Pool` plus `propagate()` relax, compress and
validate these candidates in memory and dedupe them with a D4-canonical key. **The pool is lost when the process exits.**
This is most of Phase 3 item 7 already.

### 3.7 Symmetry utilities
D4 group elements, `detect_symmetry`, `canonical_key` (quantised and sorted, minimised over the 8 images),
`wrap90` and `angle_distance90`. These are useful for the Phase 3 diversity distance.

## 4. Solution representation, file format, RNG

- **In memory:** `System` holds `species` (polygon shapes; species 0 = unit square), `particles`
  (`Particle` = x, y, θ with cached cos and sin; `frozen`; species index), `L`, `k_pair`, `k_wall`,
  tracked `EnergyParts`, a persistent `CellList`, and `merge_of` (kept on purpose for later work).
  The box is `[0, L]²`. θ is in radians and is not wrapped.
- **On disk** (`io.hpp`): a line-based text format: a header (`n`, `L`, `k_pair`, `k_wall`, `feasible_tol`,
  `nspecies`), species blocks (`s i nv` plus `v x y` lines), then `p x y θ frozen species` lines,
  written with `%.17g`. `load_record()` reads both `data/init_external/*.cfg` (1000 seeds, source
  "submission.csv") and `data/best/*.log`. Trajectories add `f`/`m` lines.
- **`data/best/`**: `save_record()` re-validates (`why_not_a_packing`), writes the file and
  renames it to `<date>_<time>_n<N>_s<side %.6f>.log`, keeps the top `kKeepPerN=5` per N and
  deletes the rest. The ranking reads the side **from the filename** (6 decimals): B5.
- **RNG:** one `std::mt19937_64(seed)` per run, plus a second stream `seed ^ 0x9e37…` for the pool.
  Every move takes the rng by reference, and there is no hidden randomness (no `random_device` or `rand`).
- **Reproducibility (checked):** same binary, same machine, same seed, count budget (`--kicks`) → bit-identical `result` lines (ran it twice).
  It breaks when:
  - `--cpu` budget: stopping depends on the clock (`std::clock`, process CPU time).
  - Across platforms: `std::uniform_real_distribution` and `uniform_int_distribution` are
    implementation-defined, so libc++ (Mac) and libstdc++ (Linux/Kaggle) draw different numbers;
    FMA contraction also differs (see §1). `std::cos`/`sin` may differ in the last ulp across libms.
  - Fix for Phase 2: our own `uniform01(rng) = (rng() >> 11) * 0x1p-53` and integer draws
    by rejection; pin `-ffp-contract=off`. Then "same seed → same output" holds on one platform and
    gives "statistically equivalent" results across platforms.

## 5. Overlap and containment

- **Narrow phase** (`polygon.cpp`): SAT test `overlaps()` (touching = no overlap, strict
  `<= 0.0`); `overlap_area()` and `overlap_area_grad()` by convex polygon clipping, with an analytic
  gradient (boundary integral); `outside_area[_grad]()` for the box.
- **Broad phase** (`geometry/cells.{hpp,cpp}`): uniform cell list, cell width ≥ reach = 2·circumradius,
  at most (2⌈√n⌉+1)² cells, 3×3 neighbourhood, doubly linked lists in flat arrays. A persistent copy lives in
  `System::cells`, updated by `System::set_pose` (contract: every pose write goes through it or is
  followed by `recount()`). `energy_parts`, `max_overlap` and `energy_grad` **build a fresh grid on every
  call** (correct even after raw writes, but they allocate each time). Particles outside the box are
  clamped into the edge cells: still correct, but they bunch up in those cells (for example during `relax_free`).
- **Tolerances:**
  - `kCertifyTol = 1e-10`: largest pair-overlap **area** and outside-box area in a "valid" packing.
  - `kFeasibleTol = 1e-5`: LNS proposal tolerance.
  - `save_record`: claimed side must equal `tight_side()` to within 1e-9 relative; area floor √(nA).
  - An area tolerance is not a depth tolerance. A tilted corner poking into an edge has
    area ≈ d², so 1e-10 area allows d ≈ 1e-5 of penetration. `batch.out` reports sides as
    "± 1e-05" for this reason. See B4.
- **Smooth penalty:** yes, separate from the hard checks. `E = k_pair·Σ(area_ij)² + k_wall·Σ(out_i)²`,
  with k=100 and an analytic gradient (`energy_grad`, O(n) via the cell list). It has finite-difference tests
  (`test_relax.cpp:25`, `:174`) at 1e-6. Note that area² scales like d⁴ for corner contacts, so
  the energy is very flat as the overlap goes to zero. This is why the relaxations need `hopeless_window` and why
  certification stops at an area of 1e-10.

## 6. Legalization

**There is none.** "Valid" means max overlap area ≤ 1e-10 after L-BFGS. No step turns a
near-feasible configuration into a strictly non-overlapping one, and the stored sides carry about 1e-5
of uncertainty. Proposal for Phase 3 item 2: a double-precision pass that (1) uses the
SAT minimum-penetration axis per overlapping pair to push squares apart, iterated, (2) falls back to
uniform dilation of centres about the box centre by the smallest factor that clears all pairs, then
refits the box, and (3) passes the independent validator at depth tolerance 0 (see Q1).

## 7. Box side: fixed or free

It is fixed during every `relax()`. It changes only outside the relaxation: `compress()` rescales it between
relaxations, `main_shrink` sets a target, and `fit_box()` refits. The MC `try_box` (NPT, enthalpy
`E + P·L²`) makes L a free variable, but no driver uses it: the LNS calls pass `pressure=0`. So the
existing structure already matches the "fixed-box ratchet" (Phase 3 item 3); `main_shrink` is essentially
that ratchet for a population of one.

## 8. Threading and global state

- **No mutable globals.** Everything lives in `System`, the option structs, `MoveLog`, `LnmStats`, `Pool`
  and the rng, all per run. Many processes in parallel are safe, **except for shared files**: all runs write
  `logs/moves_n%04d.csv` and `logs/n%04d.traj` by N, so two runs of the same N clobber each other's files.
  `save_record` is careful about name races (claimed names, `_2` suffix), but it ranks and deletes across
  concurrent writers with no lock.
- **Single-threaded.** The only threads are in `compress_best_of` (unused). Many individuals per
  thread pool would work as long as each has its own `System`, rng and log. Costs:
  per-call allocations (grid arrays in every `energy_grad`/`total_energy`; L-BFGS `s_new`/`y_new` per
  iteration; `hist_s.erase(begin)`), which will contend on malloc under threads. `System` copies are heavy
  (species vectors, cell list, shared_ptrs).
- Paths are relative to the working directory and hard-coded (`data/init_external`, `data/best`, `logs/`).

## 9. Logging and telemetry today

| Exists | Where |
|---|---|
| Progress lines on new best | stdout |
| Per-kick tally (tries, accepted, improved, CPU, gain, gain per CPU s) | stdout table |
| `MoveLog` per move type (21 ids) | `logs/moves_n%04d.csv` |
| One `result …` key=value line | stdout (read by `scripts/benchmark_moves.py`) |
| Trajectory frames (rate-limited to 0.2 s) | `logs/n%04d.traj`, viewed by `src/viz/show.py --watch` |

Missing: wall-clock time (only process CPU is recorded); a time series of the best side (time-to-target is impossible to compute);
overlap energy and max overlap of the final and best states; `relax()` stop-reason counts (jam,
hopeless, stalled, escaped); escalation and restart events; run identity (engine version or git
hash, schedule hash, host, flags); a validity verdict from an independent checker; structured output
(JSON). Output goes to fixed paths keyed by N, not by run.

## 10. Bugs and fragile parts (recorded, not fixed)

| # | Severity | Finding |
|---|---|---|
| **B1** | **high** | Every driver calls `mirror_symmetry()` on the seed before anything else (`main_shrink.cpp:207`, `main_compress.cpp:200`, `main.cpp:168`). The seeds are already valid and tight (probe: n=29 L=5.933833, overlap 0), but after mirroring, the max overlap is 0.18–0.99 (n=10, 29, 37, 200, 552, 1000; n=100 unaffected). n=29 then compresses to 6.0005 and the run ends at about 5.988–6.000, worse than the seed. `data/best` holds such worse entries (n=29 at 6.000000, n=37 at 6.5995, n=552 at 24.12, against seeds at 5.9338, 6.5986, 23.946). **The search as configured cannot improve an asymmetric seed.** `benchmark_moves.py`'s docstring says "final can never be worse than the starting packing", which is false under this behaviour. Intentional (a symmetry experiment left on) or not? → Q3. |
| **B2** | high | `src/main.cpp` is truncated: after the initial compress it declares the loop state and returns 0, with no cycle loop and no `result` line. `build/squares` is still the default `--binary` of `scripts/benchmark_moves.py`, so that benchmark currently measures nothing. |
| B3 | med | `best_side = start.L` after the initial `anneal()` (`main_shrink.cpp:218`, `main_compress.cpp:213`), but `anneal` reports `tight_side()`, which can be below the box L. Effects: (a) `main_compress` compares a cycle's tight side against the box L, so the first cycle can "improve" by slack alone; (b) if nothing improves, `save_record(best, start.L)` fails its 1e-9 side check and the record is dropped. |
| B4 | med | Certification is by overlap **area** ≤ 1e-10, which allows ≈1e-5 penetration depth; stored sides are optimistic by up to about 1e-5. There is no strict legalization (§6). This matters when comparing against published records at 1e-6. |
| B5 | med | `data/best` ranking parses `%.6f` sides out of filenames: improvements below 1e-6 are invisible, and ties are broken by name. |
| B6 | med | Cross-platform replay is impossible (implementation-defined distributions, FMA contraction): §4. |
| B7 | low | `--cpu` uses `std::clock()` (process CPU, not wall clock). Under `--cpu` runs are not reproducible, and the budget would count every thread once threading exists. |
| B8 | low | The kick rebuilds `StepSizes` every time, so `sweep()`'s step tuning never carries over (may be intended). |
| B9 | low | The seed is parsed with `strtod`: seeds above 2⁵³ lose bits. `--n abc` silently becomes 0. |
| B10 | low | `--no-save` still writes `logs/moves_n%04d.csv` into the repo. The engine writes `data/best/` itself, which conflicts with the new rule that only the runner writes `best/`. |
| B11 | low | `CellList::bin` does `int(v/side)` on possibly NaN poses: undefined behaviour before any NaN check. |
| B12 | low | `make check` runs `./$(BUILD)/…`, so it breaks with an absolute `BUILD=`. |
| B13 | doc | `lnm.hpp` says a relaxation is "a few hundred O(n²) passes"; since the cell list it is O(n). |
| B14 | repo | Git has **no commits** and about 3,250 files staged: `build/` objects, `.DS_Store`, `data/` (44 MB), `logs/` (13 MB), 1000 PNGs in `src/viz/best/`, `batch.out`. No `.gitignore`. Needs cleaning before the "one branch per step" workflow (Q8). |
| B15 | info | `batch.out` was produced by a batch runner that is no longer in the repo; the `scripts/run_one`, `wins_and_losses` and `verdict` scripts referenced in older notes do not exist. Only `scripts/benchmark_moves.py` remains. |

---

## 11. Proposals

### (a) Schedule format

Authored as YAML; the runner converts it to JSON for the engine. The engine needs a small JSON
reader (vendor `nlohmann/json`, one MIT header: Q4). **Rule: every numeric field may be a scalar or
a curve**, `[[t0, v0], [t1, v1], …]`, piecewise-linear in `t` = elapsed fraction of the budget. This
one rule gives the compression-rate curve, temperature ramps and move-mix drift, without new syntax.

```yaml
schema: sqpack.schedule/1
name: shrink-baseline            # current main_shrink, expressed as data
start:   { from: seed, mirror: false }     # seed | best | file:<path>; mirror reproduces today's B1 behaviour
prepare: { relax: true, compress: { ds: 0.01, max_relax: 5000 } }
loop:
  mode: shrink                   # shrink (main_shrink) | cycle (main_compress)
  target: { eps: 1.0e-3, eps_min: 1.0e-7, eps_max: 0.05, on_success: 1.2, on_fail: 0.95 }
  kick:
    large_prob: { base: 0.1, max: 0.6, step: 0.1, every_fails: 25 }     # stall escalation rule
    shake: { sweeps: 100, T: 0.05, translate: 0.1, rotate: 0.1, keep_tuning: false,
             mix: { translate: 0.45, rotate: 0.45, flip: 0.06, shove: 0.04 } }
    lnm:   { weights: { destroy: 0.4, strip: 0.3, corner: 0.2, loosen: 0.1, symmetrize: 0.0 },
             relax_iterations: 300, insert_tries: 12, snap_prob: 0.5,
             regions: { disc: 0.2, band: 0.5, corner: 0.3 } }
  relax:   { max_iter: 5000, history: 0, hopeless_window: 0 }
  accept:  { side_T: 0.0 }                   # cycle mode: Metropolis on side, in L/sqrt(n) units
  restart: { after_stagnant: null, to: best }   # cycle mode: 200
  propagate: { enabled: false, every_best: 1 }
```

Unknown keys are errors, and every omitted key takes the current C++ default, so `{}` equals today's
default. The engine echoes the fully resolved schedule and its hash into telemetry. Acceptance
test (Phase 2 step 3): `shrink-baseline.yaml` and `compress-baseline.yaml` with `mirror: true` must
reproduce today's `result` lines bit for bit for a few N and seeds under `--kicks`/`--cycles` budgets.

### (b) Solution file format

JSON, one file per packing, readable without the engine:

```json
{ "format": "sqpack/1", "n": 29, "side": 5.933833462681862,
  "squares": [[x, y, theta], ...],
  "meta": { "engine": "<git sha>", "schedule": "<name>@<hash>", "seed": 7, "parent": "best/29.json@<sha256>",
            "wall_s": 61.2, "host": "mac-m2", "created": "2026-11-14T21:00:00Z" } }
```
Unit squares only (no species: merge records stay in the internal format), centres in `[0, side]²`,
θ in radians (any value; the validator wraps it), `%.17g` numbers. `side` is the claim. The validator
recomputes everything from `squares` alone. Converters both ways between this and the existing `.cfg`/`.log` format, so
`data/init_external` and `data/best` import as `best/N.json` (Q6).

**Validator** (`tools/validate.py`, shares no code with the engine): parse, check N, finite
values, compute the 4 vertices per square in float64, then **exact** checks on those vertices with
`fractions.Fraction`: containment in `[0, side]²` and SAT separation for candidate pairs from a
simple grid. Verdict: valid or invalid, plus diagnostics (max float penetration depth, worst pair, recomputed
tight side). Tolerance: Q1.

### (c) Telemetry record

One JSON object per run at `--out <file>.telemetry.json`, plus a fixed 2–3 line text summary
the runner concatenates:

```json
{ "run": {"id": "...", "n": 200, "seed": 7, "schedule": "shrink-baseline@ab12", "engine": "<sha>", "host": "...", "flags": "PORTABLE"},
  "budget": {"wall_s": 600, "cpu_s": 598, "stop": "budget|converged|error"},
  "side": {"start": 14.6569, "final": 14.6512, "t_best_s": 412.0,
           "trace": [[0.0, 14.6569], [37.2, 14.6551], ...]},          // every improvement, capped at 200
  "final_state": {"energy": 3.1e-21, "max_overlap_area": 2e-11, "validated": null},
  "moves": {"shake": {"tries": 173, "acc": 28, "impr": 28, "cpu_s": 52.4, "gain": 1.1e-2}, "destroy": {...}},
  "relax": {"calls": 197, "iters": 51234, "stops": {"feasible": 31, "hopeless": 0, "stalled": 12, "max_iter": 150, "escaped": 4}},
  "events": {"escalations": 6, "restarts": 0, "jams": 162, "eps_final": 5.7e-5} }
```
Summary line example: `n=200 s7 shrink-baseline 600s: 14.65685→14.65120 (−5.6e-3, best@412s) | kicks 197 acc 16% | shake 28/173 destroy 0/13 strip 2/7 | jams 162 | VALID`.

### (d) What the engine needs to relax a population cheaply

Already there: an analytic penalty and gradient (overlap area² + wall area²), an O(N) cell list, finite-difference tests,
L-BFGS, a fixed box during relaxation, D4 utilities, and N±1 candidate generators.

Gaps:
1. **A System-free kernel.** `energy_grad` needs a `System` (species vectors, persistent grid,
   merge pointers). For a population, use flat genomes `double xyt[3N]` plus `h`, and a function
   `cost_grad(const double* xyt, int n, double h, Workspace&, double* grad)` for unit squares. That function is also the
   shape any future GPU backend would mirror.
2. **A reusable `Workspace`** (grid arrays, gradient, L-BFGS ring buffer): no allocation per
   evaluation or iteration. Today every evaluation allocates a grid, and every iteration allocates `s_new`/`y_new`.
3. **A new relaxer** next to the existing one, not replacing it: L-BFGS with history ≈3, a hard per-variable step cap, no line search,
   stages (rough with a big cap, then fine with smaller caps). It must not change `relax()`.
4. **Backend interface**: `RelaxBackend::relax_batch(span<const Genome> in, const RelaxParams&, span<Genome> out, span<RelaxStats>)`,
   with a CPU implementation over a thread pool, one individual per task. A float32 CUDA backend would come later, only with your approval and after the micro-benchmark.
5. **A legalization pass** (§6), in double, on the CPU.
6. **A penalty choice to benchmark**: area² is d⁴-flat at corner contacts. A penetration-depth²
   penalty (SAT minimum overlap along the separating axes) may converge faster near feasibility. This would be a schedule
   field (`penalty: area2 | depth2`), gated by the benchmark.
7. **A thread pool**: none exists; `compress_best_of`'s static striding is the only precedent.
8. **Pose cost**: `Particle::set_pose` calls cos and sin on every write. That is fine for the kernel if the genome stores θ
   and the kernel computes cos and sin once per evaluation.

### Plan impact

`main_shrink` is already a population-of-one fixed-box ratchet, and the N±1 generators and
D4 canonicalisation exist, so Phase 3 items 3 and 7 are partly built. Phase 3 item 1 (cell list) is
**done**; what remains is gap (d)2: removing per-call allocations.

---

## 12. Open questions

1. **Validity standard.** Strict (zero penetration on the float vertices, exact arithmetic), or a
   tolerance (for example depth ≤ 1e-9)? This decides whether legalization is mandatory before anything enters `best/`.
2. **Which protocol is "my current best"**: `main_shrink` (newest, Nov 11 runs) or `main_compress`?
   Is there a specific command line or seed set whose results I should reproduce?
3. **B1, the mirror at startup**: deliberate? If not, the new CLI defaults to `mirror: false`, which
   changes behaviour. The acceptance test would run with `mirror: true` to match the old results.
4. OK to **vendor `nlohmann/json`** (single MIT header) for schedule parsing in C++? The alternative is
   a tiny hand-written parser.
5. **Budget semantics**: wall-clock (your spec: `--time <sec>`) or CPU? I propose wall-clock for the runner,
   and CPU plus an iteration count recorded for reproducibility.
6. **Migration of `data/best`**: import the best entry per N into `best/N.json` (after validation), and keep
   `data/` read-only from then on? Should `data/init_external` remain the seed source?
7. Should the new `engine` binary **replace** `squares_shrink`/`squares_compress` (sharing their loops through the
   schedule interpreter) or sit beside them? I suggest beside them, until the acceptance test passes.
8. **Git cleanup before Phase 2**: unstage `build/`, `logs/`, `.DS_Store`, `batch.out`, `src/viz/best/*.png`, and
   possibly `data/`, then add a `.gitignore` (`build/ runs/ logs/ .DS_Store`) and make an initial commit. Is `data/` (44 MB) meant to be in git?
9. Should `main.cpp` (B2) be deleted, or restored? `benchmark_moves.py` depends on it.
10. **Kaggle toolchain**: is g++ acceptable there, or should I plan for clang? It affects the replay guarantees.

---

## 13. Decisions (2026-11-14, after review)

- **Q1 validity:** legalization is required before anything enters `best/`; `best/` holds only
  validated configurations. Near-feasible results may go to `non_validated_best/`. Standard adopted:
  strict, i.e. no positive-area overlap and no part outside `[0, side]^2`, checked exactly on the float
  vertices (touching is allowed).
- **Q2 protocol:** there is no settled "best protocol". Both `main_shrink` and `main_compress` become
  baseline schedules. Better protocols are to be found, taking inspiration from jcottaar/packing.
- **Q3 mirror:** deliberate. Symmetric starts cut the search space, but they need serious work before they pay off.
  It stays a schedule field (`start.mirror`), on in the baseline schedules.
- **Q4 JSON:** vendor `nlohmann/json` (single header, MIT).
- **Q8 git:** visualisation output, venvs, build output, logs, scratch and the legacy `data/best/` are
  gitignored. `data/init_external/` (22 MB, never changes, cannot be regenerated) stays tracked.
