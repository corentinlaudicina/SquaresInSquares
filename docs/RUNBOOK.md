# Runbook

## Build and test
```sh
make -j8 all          # engine, old drivers, tests, bench_relax
make check            # C++ suites + Python tests (validator, runner, benchmark, engine paths)
make acceptance       # baseline schedules vs old drivers, bit for bit (~1 min)
make bench            # build/bench_relax: the population-relaxation micro-benchmark
python3 -m venv .venv && .venv/bin/pip install -r tools/requirements.txt   # PyYAML for queue.yaml
```
While a benchmark or the runner is using `build/engine`, do not rebuild it. Develop in a separate
build directory instead: `make -j2 BUILD=build-dev build-dev/engine` (gitignored).

## One run
```sh
build/engine --N 29 --seed 1 --time 300 --schedule schedules/ga.json --out /tmp/n29.json
python3 tools/validate.py /tmp/n29.json --n 29        # exit 0 = valid
build/engine --list-schedule-keys                      # the schedule vocabulary
build/engine --legalize <in.cfg|json> --out out.json   # legalize only, no search
```
- `--start <file>`: starting packing (default `data/init_external/nNNNN.cfg`; the runner passes `best/N.json`).
- `--kicks k`: count budget (kicks, cycles, or GA generations): reproducible, unlike `--time`.
- `--neighbors <dir>`: where GA seeding reads N+1 / N-1 (default `best/`).
- The engine writes only `--out` and its `.telemetry.json`. It never writes `best/`. Every output is legalized.

## Schedules (`schedules/`)
| file | mode | what |
|---|---|---|
| `shrink-baseline` | shrink | the old `squares_shrink`, bit for bit |
| `compress-baseline` | cycle | the old `squares_compress`, bit for bit |
| `shrink-legal` | shrink | + accept on the legal side |
| `shrink-capped` | shrink | + depth^2 penalty, capped L-BFGS, eps reset (passed the item 2 gate) |
| `ga` | ga | ratchet GA, selection by fitness (items 3-4, passed) |
| `ga-seeded` | ga | **default for N >= 82**: 6 islands seeded from N+1, N-1, lattices (item 7, passed) |
| `ga-div` | ga | `ga` with diversity selection (item 5, failed: for tuning only) |
| `ga-islands` | ga | 4 islands on a ring, with diversity (item 6, inconclusive) |
| `ga-c2`, `ga-c4` | ga | symmetric genotypes (item 8, failed); C4 only for N mod 4 in {0, 1} |
| `ga-seeded-adapt` | ga | production since 2026-11-15: ga-seeded + adaptive ratchet |
| `ga-explore` | ga | ga-seeded-adapt + large moves (destroy, strip, corner), explorer slots, resets from own source |
| `shrink-polish` | shrink | shrink-capped without the mirrored start |

The GA is explained in the headers of `src/pop/ga.hpp` (loop, ratchet, islands, seeding),
`src/pop/moves.hpp`, `src/pop/diversity.hpp` and `src/pop/symmetric.hpp`. The mirror start
(`start.mirror`) is on in the shrink and cycle baselines and off in the GA schedules.

### How the GA code fits together (`src/pop/`)
- `penalty`: a packing is a flat `Packing {h, xyt}`. `cost_grad` gives the penalty and its gradient,
  using a neighbour grid, or a Verlet list while a relaxation runs. A `Workspace` holds the buffers
  for one thread.
- `relax_capped`: L-BFGS in stages, with every step capped. `relax_core` holds the loop.
  `relax_capped` (plain) and `relax_symmetric` (generators) only choose the cost function.
- `backend`: `CpuBackend` is a thread pool with one `Workspace` per thread. `relax_batch` and
  `cost_batch` take an optional `SymLayout`; with C2/C4 they work on generators. It is a plain
  class: a GPU backend was measured and ruled out, so it has no interface.
- `ga.cpp` `run_ga`: legalize the start, build the seeding candidates (`seeding`), then seed each
  island with `seed_island` (rough relax, then fine). Each generation, for each island: make the
  offspring (`moves`), pre-filter them on the raw cost, relax them rough then fine, select
  (`diversity`; fitness only by default), then ratchet or count a stall. A run that has not
  improved for `island_patience` generations reseeds from the best legal packing.
- Boxes are sized from `out.best.L`. Once anything is legal, this equals `out.best_side`. If the
  start cannot be legalized, `out.best_side` stays infinite but `out.best.L` stays finite.
- Symmetric genotypes (`ga.symmetry`, item 8) and diversity (item 5) failed their gates. They are
  kept as opt-in features for later work, not dead code.

### Change log
- 2026-11-15, review cleanup, no change in results (same packing and `ga` telemetry from the old and
  new engine on `ga`, `ga-c2` and `ga-seeded`, N=29, `--kicks 15`):
  - Removed the `RelaxBackend` interface. It had one implementation, `CpuBackend`, which is now a
    plain class.
  - `run_ga`: the box for one island and for an island reset now comes from `out.best.L`. Before,
    it came from `out.best_side`, which is infinite when the start cannot be legalized; the box
    became infinite, the positions became NaN and the whole budget was wasted.
  - `run_ga`: symmetry telemetry is written once at the end (no save and restore around the
    `tel["ga"]` overwrite). The start-source pointer is set once. The `stop` variable, which never
    changed, is gone. `relax_all` takes `const SymLayout&` (every caller passed one).

## Unattended
```sh
cp queue.example.yaml queue.yaml    # edit
.venv/bin/python tools/runner.py --queue queue.yaml --workers 4
```
- `best/N.json`: validated, legal packings only. The runner replaces one only on a strict improvement.
  Gitignored: sync it between machines outside git.
- `non_validated_best/N.json`: smaller than `best/` but not valid. Kept for inspection, never trusted.
- `runs/ledger.jsonl`: one line per finished job (sides, validity, where it was stored, stop reason,
  per-move totals, summary line). What the planner reads, and the done marker.
- `runs/jobs/<id>/`: schedule, output, full telemetry, stdout. Kept only for improvements of `best/`,
  non-validated improvements and failures (each marked by a `keep` file with the reason), plus the
  newest 20 other jobs; older unmarked folders are deleted.
- `runs/summary.md`: short status for the planner. `runs/heartbeat`: touched every 30 s.
- Restart-safe: rerunning the same queue skips the jobs in the ledger and re-runs unfinished ones.
- GA schedules use `ga.threads` relaxation threads each: keep `workers x threads` near the core count.
- One-off import of the legacy records: `.venv/bin/python tools/runner.py --import-legacy`.

## The search machine (planner)
```sh
tools/setup_machine.sh --schedule   # once: build, test, best/, planner/, daily 07:00 job
claude                              # once: log the CLI in
tools/daily.sh                      # a cycle by hand: digest -> plan -> check -> runner until 06:40
```
Only N marked `open` in `data/known_optima_1_1000.csv` are ever queued. Details: `docs/PLANNER.md`.
The daily cycle also drafts a lab-notebook entry (`planner/notebook/<date>.md`, local);
`tools/publish_notebook.sh` appends the new drafts to `docs/NOTEBOOK.md`, commits and pushes (by hand:
this is what makes results public; `-n` previews).
Live view: `src/viz/.venv/bin/python tools/live_view.py` (every running job, one panel each; pass a job id
for its best and every island champion). It reads the snapshots the engine writes with `--live`
(`runs/jobs/<id>/live.json`, every ~2 s; the runner always passes it).
Watch it: `python3 tools/status.py` (one screen; `watch -n 60 python3 tools/status.py` to keep it open).

## Benchmark
```sh
python3 tools/benchmark.py --schedules schedules/shrink-baseline.json schedules/new.json \
    --set train --seeds 3 --time 120 --workers 4 --name new-vs-base
```
- Report in `runs/bench/<name>/report.md`. Do not tune on `--set holdout`.
- Targets are frozen from `best/` at start (or `--targets file.json`).
- **Starting from `best/` cannot show a search improving** (it returns the strong start). Use the
  recovery gate: `--starts 'runs/starts/recovery/n{N:04d}.json'` (legal, worse, different arrangements).
- `--combine name1 name2 --schedules ...` merges earlier runs into one paired report (same N, seeds, budget).
- Only `results.jsonl` and `report.md` are kept; `--keep` also keeps every run's packing and telemetry.
- GA runs use several threads: run them with fewer `--workers` (e.g. 2 x 4 threads on 8 cores).

## Micro-benchmark
```sh
build/bench_relax --N 200 --count 1000 --threads 1,2,4,8   # kernel cost, relax time, thread scaling
```

## Formats
- Solution: `sqpack/1` JSON, see `tools/validate.py`.
- Schedule: `sqpack.schedule/1` JSON, see `src/engine/schedule.hpp` and `schedules/*.json`.
- Telemetry: written by `src/engine/protocol.cpp` (`run_engine`), `src/pop/ga.cpp` (`run_ga`: per-move
  survival, per-generation trace, islands, seeds) and `src/main_engine.cpp`.
- Findings and gate results: `docs/AUDIT.md` (§13-15).
