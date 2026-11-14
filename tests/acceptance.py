#!/usr/bin/env python3
"""Acceptance test: the baseline schedules reproduce the old drivers exactly.

    make shrink compress engine && python3 tests/acceptance.py [--quick]

For each N and seed, runs build/squares_shrink (resp. squares_compress) with --no-save and
build/engine with schedules/shrink-baseline.json (resp. compress-baseline.json) on the same kick
budget, and requires the same best side to the last printed digit (12 decimals). Writes only to a
temporary directory.
"""
import argparse
import json
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FINAL = re.compile(r"final=([0-9.]+)")


def old(binary, budget_flag, n, seed, budget, tmp):
    out = subprocess.run([str(ROOT / "build" / binary), "--n", str(n), "--seed", str(seed), budget_flag,
                          str(budget), "--no-save", "--csv", f"{tmp}/{binary}_{n}_{seed}.csv"],
                         cwd=ROOT, capture_output=True, text=True, check=True).stdout
    return float(FINAL.search(out).group(1))


def new(schedule, n, seed, budget, tmp):
    out = f"{tmp}/engine_{schedule}_{n}_{seed}.json"
    subprocess.run([str(ROOT / "build/engine"), "--N", str(n), "--seed", str(seed), "--kicks", str(budget),
                    "--schedule", str(ROOT / "schedules" / f"{schedule}.json"), "--out", out],
                   cwd=ROOT, capture_output=True, text=True, check=True)
    tel = json.loads(Path(out.replace(".json", ".telemetry.json")).read_text())
    return tel["side"]["best"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="2 N x 1 seed")
    a = ap.parse_args()
    ns, seeds = ((10, 29), (1,)) if a.quick else ((10, 29, 37), (1, 2))
    cases = [(mode, n, s) for mode in ("shrink", "compress") for n in ns for s in seeds]
    with tempfile.TemporaryDirectory() as tmp, ThreadPoolExecutor(8) as pool:
        def run(case):
            mode, n, s = case
            if mode == "shrink":
                return old("squares_shrink", "--kicks", n, s, 20, tmp), new("shrink-baseline", n, s, 20, tmp)
            return old("squares_compress", "--cycles", n, s, 3, tmp), new("compress-baseline", n, s, 3, tmp)
        results = list(pool.map(run, cases))
    fails = 0
    for (mode, n, s), (o, e) in zip(cases, results):
        same = f"{o:.12f}" == f"{e:.12f}"
        fails += not same
        print(f"{mode:8s} n={n:3d} seed={s}  old {o:.12f}  engine {e:.12f}  {'same' if same else 'DIFFERENT'}")
    print(f"acceptance: {len(cases) - fails}/{len(cases)} identical")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
