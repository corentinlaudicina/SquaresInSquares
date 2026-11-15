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

The GA is explained in the headers of `src/pop/ga.hpp` (loop, ratchet, islands, seeding),
`src/pop/moves.hpp`, `src/pop/diversity.hpp` and `src/pop/symmetric.hpp`. The mirror start
(`start.mirror`) is on in the shrink and cycle baselines and off in the GA schedules.

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
