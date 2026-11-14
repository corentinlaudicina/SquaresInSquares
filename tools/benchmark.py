#!/usr/bin/env python3
"""Benchmark schedules on fixed N sets, paired seeds and a fixed wall-clock budget.

    python3 tools/benchmark.py --schedules schedules/shrink-baseline.json schedules/compress-baseline.json \\
        [--set train|holdout|all] [--seeds 3] [--time 120] [--workers 4] [--targets targets.json] [--name tag]

Every (N, seed) pair runs once per schedule, from the same start (best/N.json, else the seed file), so
schedules are compared on identical inputs. Writes runs/bench/<name>/{results.jsonl, report.md}; never
touches best/.

Per schedule:
  success   share of runs whose legal, validated side <= target * (1 + 1e-9)
  ttt       median time to reach the target (from the telemetry trace of the raw side; runs that
            never reach it are excluded, so read it together with success)
  mean side mean legal side, and mean (side - target)
Paired: wins / losses / ties on legal side against the first schedule, run for run.

Targets: best/N.json's side when the benchmark starts, unless --targets maps N -> side (for example
published records). shortcut: the ttt uses the pre-legalization side in the trace; legalization adds
~1e-5 (docs/AUDIT.md §14), so ttt is optimistic by that much until the search certifies tighter.
"""
import argparse
import json
import os
import statistics
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import validate  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
# Regimes: small | near k^2 | near k^2 - k | mid | large. Held-out N are never used for tuning.
TRAIN = [11, 17, 19, 26, 29, 37, 41, 50, 53, 56, 65, 82, 99, 110, 145, 200, 240, 333, 552, 800]
HOLDOUT = [18, 28, 39, 48, 72, 120, 168, 300, 420, 1000]
TOL = 1e-9


def start_file(n):
    p = ROOT / "best" / f"{n}.json"
    return p if p.exists() else ROOT / "data" / "init_external" / f"n{n:04d}.cfg"


def target_of(n, targets):
    if str(n) in targets:
        return float(targets[str(n)])
    p = ROOT / "best" / f"{n}.json"
    if p.exists():
        return json.loads(p.read_text())["side"]
    # no best yet: the seed's claimed side
    for line in (ROOT / "data" / "init_external" / f"n{n:04d}.cfg").read_text().splitlines():
        if line.startswith("L "):
            return float(line.split()[1])
    raise ValueError(f"no target for n={n}")


def run_one(job, outdir, engine, budget):
    sched, n, seed, target = job
    name = Path(sched).stem
    out = outdir / f"{name}_n{n:04d}_s{seed}.json"
    cmd = [str(engine), "--N", str(n), "--seed", str(seed), "--time", str(budget), "--schedule", str(ROOT / sched),
           "--start", str(start_file(n)), "--out", str(out)]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=budget * 1.5 + 180)
        status = "ok" if p.returncode == 0 else f"exit_{p.returncode}"
    except subprocess.TimeoutExpired:
        status = "timeout"
    r = {"schedule": name, "N": n, "seed": seed, "target": target, "status": status, "wall_s": round(time.time() - t0, 1),
         "valid": False, "side": None, "ttt": None}
    tel_path = Path(str(out).replace(".json", ".telemetry.json"))
    if out.exists():
        doc = validate.load(out)
        r["valid"] = validate.validate(doc, n)[0]
        r["side"] = doc["side"]
    if tel_path.exists():
        tel = json.loads(tel_path.read_text())
        start = tel.get("side", {}).get("start")
        trace = ([[0.0, start]] if start else []) + tel.get("side", {}).get("trace", [])
        hit = [t for t, s in trace if s <= target * (1 + TOL)]
        r["ttt"] = hit[0] if hit else None
        r["kicks"] = tel.get("budget", {}).get("kicks")
    r["success"] = bool(r["valid"] and r["side"] is not None and r["side"] <= target * (1 + TOL))
    return r


def report(results, schedules, ns, seeds, budget):
    names = [Path(s).stem for s in schedules]
    lines = [f"# Benchmark: {len(ns)} N x {len(seeds)} seeds x {budget:.0f}s, {len(results)} runs", "",
             "| schedule | success | median ttt s | mean legal side | mean side - target | invalid/failed |",
             "|---|---|---|---|---|---|"]
    by = {nm: [r for r in results if r["schedule"] == nm] for nm in names}
    for nm in names:
        rs = by[nm]
        ok = [r for r in rs if r["valid"]]
        ttt = [r["ttt"] for r in rs if r["ttt"] is not None and r["success"]]
        lines.append(f"| {nm} | {sum(r['success'] for r in rs)}/{len(rs)} | "
                     f"{statistics.median(ttt):.1f} | " if ttt else f"| {nm} | {sum(r['success'] for r in rs)}/{len(rs)} | - | ")
        lines[-1] += (f"{statistics.mean(r['side'] for r in ok):.6f} | "
                      f"{statistics.mean(r['side'] - r['target'] for r in ok):+.3e} | " if ok else "- | - | ")
        lines[-1] += f"{sum(not r['valid'] for r in rs)} |"
    if len(names) > 1:
        lines += ["", f"Paired against `{names[0]}` (legal side, lower is better; tie within {TOL:g} relative):", ""]
        ref = {(r["N"], r["seed"]): r for r in by[names[0]]}
        for nm in names[1:]:
            w = l = t = 0
            for r in by[nm]:
                a = ref.get((r["N"], r["seed"]))
                if not a or not a["valid"] or not r["valid"]:
                    continue
                if abs(r["side"] - a["side"]) <= TOL * a["side"]:
                    t += 1
                elif r["side"] < a["side"]:
                    w += 1
                else:
                    l += 1
            lines.append(f"- `{nm}`: {w} wins, {l} losses, {t} ties")
    lines += ["", "## Per N (mean legal side over seeds; * = all seeds reached the target)", "",
              "| N | target | " + " | ".join(names) + " |", "|---|---|" + "---|" * len(names)]
    for n in ns:
        row = [str(n), f"{next(r['target'] for r in results if r['N'] == n):.6f}"]
        for nm in names:
            rs = [r for r in by[nm] if r["N"] == n and r["valid"]]
            allhit = rs and all(r["success"] for r in [x for x in by[nm] if x["N"] == n])
            row.append(f"{statistics.mean(r['side'] for r in rs):.6f}{'*' if allhit else ''}" if rs else "-")
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--schedules", nargs="+", required=True)
    ap.add_argument("--set", choices=["train", "holdout", "all"], default="train")
    ap.add_argument("--n", type=int, nargs="*", help="explicit N list, overrides --set")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--time", type=float, default=120.0)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--targets", help="JSON {N: side}")
    ap.add_argument("--name", default=time.strftime("%Y%m%d_%H%M%S"))
    ap.add_argument("--combine", nargs="+", metavar="NAME",
                    help="no runs: merge runs/bench/<NAME>/results.jsonl into one report (same N, seeds, budget)")
    a = ap.parse_args(argv)
    if a.combine:
        keep = {Path(x).stem for x in a.schedules}
        results = [r for nm in a.combine for l in (ROOT / "runs" / "bench" / nm / "results.jsonl").open()
                   if (r := json.loads(l))["schedule"] in keep]
        ns = sorted({r["N"] for r in results})
        seeds = sorted({r["seed"] for r in results})
        text = report(results, a.schedules, ns, seeds, a.time)
        out = ROOT / "runs" / "bench" / a.name
        out.mkdir(parents=True, exist_ok=True)
        (out / "report.md").write_text(text)
        print(text)
        return 0
    engine = ROOT / "build" / "engine"
    if not engine.exists():
        sys.exit("benchmark: make engine first")
    ns = a.n or {"train": TRAIN, "holdout": HOLDOUT, "all": TRAIN + HOLDOUT}[a.set]
    seeds = list(range(1, a.seeds + 1))
    targets = json.loads(Path(a.targets).read_text()) if a.targets else {}
    tmap = {n: target_of(n, targets) for n in ns}  # frozen before any run
    outdir = ROOT / "runs" / "bench" / a.name
    outdir.mkdir(parents=True, exist_ok=True)
    # Large N first, so the long jobs do not trail at the end.
    jobs = sorted([(s, n, seed, tmap[n]) for s in a.schedules for n in ns for seed in seeds], key=lambda j: -j[1])
    print(f"benchmark {a.name}: {len(jobs)} runs, {a.workers} workers, ~{len(jobs) * a.time / a.workers / 60:.0f} min")
    results = []
    with ThreadPoolExecutor(a.workers) as pool, open(outdir / "results.jsonl", "w") as f:
        for r in pool.map(lambda j: run_one(j, outdir, engine, a.time), jobs):
            results.append(r)
            f.write(json.dumps(r) + "\n")
            f.flush()
    text = report(results, a.schedules, ns, seeds, a.time)
    (outdir / "report.md").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
