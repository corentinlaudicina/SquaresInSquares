# Runbook

## Build and test
```sh
make -j8 all          # engine, old drivers, tests
make check            # C++ suites + Python tests (validator, runner, benchmark)
make acceptance       # baseline schedules vs old drivers, bit for bit (~1 min)
python3 -m venv .venv && .venv/bin/pip install -r tools/requirements.txt   # PyYAML for queue.yaml
```

## One run
```sh
build/engine --N 29 --seed 1 --time 300 --schedule schedules/shrink-baseline.json --out /tmp/n29.json
python3 tools/validate.py /tmp/n29.json --n 29        # exit 0 = valid
build/engine --list-schedule-keys                      # the schedule vocabulary
build/engine --legalize <in.cfg|json> --out out.json   # legalize only, no search
```
The engine writes only `--out` and its `.telemetry.json`. It never writes `best/`.

## Unattended
```sh
cp queue.example.yaml queue.yaml    # edit
.venv/bin/python tools/runner.py --queue queue.yaml --workers 4
```
- `best/N.json`: validated, legal packings only. The runner replaces one only on a strict improvement.
- `non_validated_best/N.json`: smaller than `best/` but not valid. Kept for inspection, never trusted.
- `runs/jobs/<id>/`: per-job schedule, output, telemetry, stdout, `result.json` (the done marker).
- `runs/summary.md`: short status for the planner. `runs/heartbeat`: touched every 30 s.
- Restart-safe: rerunning the same queue skips done jobs and re-runs unfinished ones.

## Benchmark
```sh
python3 tools/benchmark.py --schedules schedules/shrink-baseline.json schedules/new.json \
    --set train --seeds 3 --time 120 --workers 4 --name new-vs-base
```
Report in `runs/bench/<name>/report.md`. Do not tune on `--set holdout`.

## Formats
- Solution: `sqpack/1` JSON, see `tools/validate.py`.
- Schedule: `sqpack.schedule/1` JSON, see `src/engine/schedule.hpp` and `schedules/*.json`.
- Telemetry: written by `src/engine/protocol.cpp` (`run_engine`) and `src/main_engine.cpp`.
