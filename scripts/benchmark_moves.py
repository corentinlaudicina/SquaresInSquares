#!/usr/bin/env python3
"""Compare the baseline engine with the symmetry moves, on a fixed CPU budget and fixed seeds.

    make squares && python3 scripts/benchmark_moves.py

Three configurations of build/squares, for every n and seed: baseline, baseline + symmetrize_region
and baseline + expand/propagate. Nothing is written to data/best. Each seed runs all three
configurations, so they are compared seed by seed, and a configuration is only called better for an
n if it beats the baseline in most seeds and loses in none.

`final` can never be worse than the starting packing, so when the seeds are already at the best
known side every configuration ties there. `best_cycle` is the smallest side the search reached
by itself, which still separates them in that case.
"""

import argparse
import csv
import re
import statistics
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

CONFIGS = {"baseline": [], "symmetrize": None, "propagate": ["--propagate"]}
RESULT = re.compile(r"^result (.*)$", re.M)
TIE = 1e-9


def run_one(job):
    binary, n, seed, config, flags, budget, out_dir = job
    stem = f"n{n:04d}_s{seed}_{config}"
    csv_path = out_dir / f"{stem}.moves.csv"
    cmd = [binary, "--n", str(n), "--seed", str(seed), "--cpu", str(budget), "--no-save",
           "--csv", str(csv_path), *flags]
    done = subprocess.run(cmd, capture_output=True, text=True, check=False)
    (out_dir / f"{stem}.out").write_text(done.stdout + done.stderr)
    match = RESULT.search(done.stdout)
    if done.returncode != 0 or not match:
        sys.exit(f"run failed ({done.returncode}): {' '.join(cmd)}\n{done.stderr}")
    row = {"config": config}
    for field in match.group(1).split():
        key, value = field.split("=")
        row[key] = float(value)
    row["moves"] = list(csv.DictReader(csv_path.open()))
    return row


def wins_and_losses(rows, n, config, seeds, field):
    lookup = {(r["config"], int(r["seed"])): r[field] for r in rows if r["n"] == n}
    wins = losses = 0
    for seed in seeds:
        delta = lookup[(config, seed)] - lookup[("baseline", seed)]
        wins += delta < -TIE
        losses += delta > TIE
    return wins, losses


def verdict(wins, losses, seeds):
    if wins > len(seeds) / 2 and losses == 0:
        return "better"
    if losses > len(seeds) / 2 and wins == 0:
        return "worse"
    return "no difference shown"


def report(rows, ns, seeds, budget):
    budgets = ", ".join(f"n={n}: {budget_for(n, budget):g}" for n in ns)
    lines = [f"# Benchmark, CPU s per run ({budgets}), seeds {seeds}", ""]
    for field, title in (("final", "final L(n)"), ("best_cycle", "best side reached by a cycle")):
        lines += [f"## {title}", "",
                  "| n | config | mean | min | mean L/sqrt(n)-1 | wins | losses | verdict |",
                  "|---|---|---|---|---|---|---|---|"]
        for n in ns:
            for config in CONFIGS:
                values = [r[field] for r in rows if r["n"] == n and r["config"] == config]
                mean = statistics.mean(values)
                score = mean / n ** 0.5 - 1.0
                if config == "baseline":
                    tail = "- | - | -"
                else:
                    w, l = wins_and_losses(rows, n, config, seeds, field)
                    tail = f"{w} | {l} | {verdict(w, l, seeds)}"
                lines.append(f"| {n} | {config} | {mean:.6f} | {min(values):.6f} | {score:.6f} | {tail} |")
        lines.append("")

    lines += ["## Runs", "", "| n | config | mean cycles | mean cpu s | other n improved |",
              "|---|---|---|---|---|"]
    for n in ns:
        for config in CONFIGS:
            subset = [r for r in rows if r["n"] == n and r["config"] == config]
            lines.append(f"| {n} | {config} | {statistics.mean(r['cycles'] for r in subset):.1f} | "
                         f"{statistics.mean(r['cpu'] for r in subset):.1f} | "
                         f"{int(sum(r['other_n_improved'] for r in subset))} |")
    lines.append("")

    for config in CONFIGS:
        totals = {}
        for r in rows:
            if r["config"] != config:
                continue
            for m in r["moves"]:
                t = totals.setdefault(m["move"], [0, 0, 0, 0.0, 0.0])
                t[0] += int(m["attempts"])
                t[1] += int(m["accepted"])
                t[2] += int(m["improved"])
                t[3] += float(m["cpu_seconds"])
                t[4] += float(m["improvement"])
        lines += [f"## Moves: {config}", "",
                  "| move | attempts | accepted | improved | cpu s | improvement | improvement/cpu s |",
                  "|---|---|---|---|---|---|---|"]
        for name, t in totals.items():
            if t[0]:
                rate = t[4] / t[3] if t[3] > 0 else 0.0
                lines.append(f"| {name} | {t[0]} | {t[1]} | {t[2]} | {t[3]:.2f} | {t[4]:.3e} | {rate:.3e} |")
        lines.append("")
    return "\n".join(lines)


def budget_for(n, base):
    return base * max(1.0, n / 50.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--ns", type=int, nargs="+", default=[17, 29, 37, 50, 101, 200])
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4])
    ap.add_argument("--cpu", type=float, default=60.0,
                    help="CPU seconds per run up to n=50, then proportional to n: "
                         "a cycle at n=200 costs about 40 s")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--symmetrize", type=float, default=2.0,
                    help="w_symmetrize. A rejected attempt redraws, so this adds about "
                         "0.1 * w attempts per cycle and leaves the other moves' rates alone.")
    ap.add_argument("--binary", default="build/squares")
    ap.add_argument("--out", default="logs/benchmark")
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    CONFIGS["symmetrize"] = ["--symmetrize", str(args.symmetrize)]

    jobs = [(args.binary, n, seed, config, flags, budget_for(n, args.cpu), out_dir)
            for n in args.ns for seed in args.seeds for config, flags in CONFIGS.items()]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = list(pool.map(run_one, jobs))

    fields = ["n", "seed", "config", "cycles", "cpu", "start", "best_cycle", "final", "score",
              "other_n_improved", "other_n_gain"]
    with (out_dir / "runs.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    text = report(rows, args.ns, args.seeds, args.cpu)
    (out_dir / "report.md").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
