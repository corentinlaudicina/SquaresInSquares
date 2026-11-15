# Squares in Squares

An unattended search for dense packings of N unit squares (rotation allowed) in the smallest square
container, N = 1 to 1000. A C++20 engine does the search, an independent validator decides what counts,
a runner keeps only verified improvements, and a daily planning step (one Claude call) chooses what to
run next.

This repository holds code, scripts and documentation only: no packings. Starting configurations and
results (`data/init_external/`, `best/`) stay local.

## How it works

```
queue ──► runner ──► engine (one N, one schedule, one seed, a time budget)
            │           start from best/N.json → search → legalize → packing + telemetry
            │           ▼
            │        validate.py (independent, exact) ──► valid and smaller? ──► best/N.json
            ▼
      ledger + summary ──► daily planner (claude -p, no tools) ──► next queue
```

- **Engine** (`src/`, `build/engine`): three search modes, all described by schedule files.
  - `shrink`: one packing, a target box slightly below the best side, a kick (shake or large
    rearrangement), relaxation, keep if it fits.
  - `cycle`: one packing, kick, wall-free relaxation, full compression, Metropolis on the side.
  - `ga`: a population at a shared box side h, many offspring from local moves and region crossover,
    a raw-cost pre-filter, staged relaxation, selection, and a ratchet that shrinks h whenever the best
    individual fits. Optional islands, seeding from N+1 / N-1 / lattices, and C2/C4-symmetric genotypes.
- **Relaxation**: penetration-depth penalty with analytic gradients, staged L-BFGS (history 3, capped
  steps, no line search), cell lists and a neighbour list, a thread pool behind a batch interface.
- **Legalization**: every output is made strictly valid in double precision before it is written.
- **Validator** (`tools/validate.py`): checks containment and overlap exactly (rational arithmetic on the
  float corners) and shares no code with the engine. Touching is allowed; any positive overlap is not.
- **Runner** (`tools/runner.py`): parallel jobs with timeouts, validation, atomic updates of `best/` on
  strict improvement only, a one-line-per-job ledger, restart-safe.
- **Benchmark** (`tools/benchmark.py`): fixed training and held-out N, paired seeds, fixed wall clock.
  Every search method was kept only if it beat the baseline (results in `docs/AUDIT.md`, section 15).
- **Planner** (`tools/daily.sh`): a script digest of the ledger, one `claude -p` turn with no tools that
  answers with a JSON plan, a script that validates the plan (open N only, budget, limits) before
  anything runs. N whose optimum is proven (`data/known_optima_1_1000.csv`) are never worked on.

## Quick start

```sh
make -j8 all && make check                 # build (clang++ or g++, C++20) and run the tests
python3 -m venv .venv && .venv/bin/pip install -r tools/requirements.txt

build/engine --N 29 --seed 1 --time 60 --schedule schedules/ga-seeded.json --start my_start.json --out n29.json
python3 tools/validate.py n29.json --n 29  # exit 0 = valid
```

Unattended: `tools/setup_machine.sh --schedule`, then `tools/daily.sh`; watch with `python3 tools/status.py`.
Details in `docs/RUNBOOK.md` (operations), `docs/PLANNER.md` (the daily planner) and `docs/AUDIT.md`
(the original audit, design decisions and every benchmark gate).

Starting packings: the runner starts each job from `best/N.json`. A fresh checkout has none. Either
supply your own starts (legacy `.cfg` records in `data/init_external/nNNNN.cfg`, imported with
`tools/runner.py --import-legacy`), or start from scratch packings via `--start`.

## Credits

- **jcottaar/packing** (Santa 2025, 1st place; CC BY-SA 4.0): the population layer (fixed-box ratchet,
  overproduction with pre-filtering, rough/fine capped relaxation, the move set, diversity-aware
  selection, islands, seeding, symmetric genotypes) follows ideas from that solution. Everything here is
  an independent reimplementation for squares in C++; no code was copied.
- **nlohmann/json** (MIT, Niels Lohmann), vendored in `third_party/nlohmann/` with its license.
- The table of proven optima in `data/known_optima_1_1000.csv` cites its sources per row (for example
  Nagamochi 2005 for s(k^2 - 1) and s(k^2 - 2)).

## How this was built

The project, its goals, the architecture and every decision were Corentin Laudicina's: the phase plan,
the validity standard, keeping the symmetric start, which results to trust and what to merge. Most of the
code, tests and documentation were written by Claude (Anthropic's model, through Claude Code) under that
direction, step by step. Each step was checked with tests, and each search method was kept only if it
beat the baseline on a fixed benchmark. Every packing the system keeps is checked by a validator that
shares no code with the search engine. At run time the daily planner is also a Claude call; it has no
tools, and a script validates its plans. Commits are authored by Corentin and carry a `Co-Authored-By:
Claude` line.

## License

No license has been chosen yet: all rights reserved.
