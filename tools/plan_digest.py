#!/usr/bin/env python3
"""The planner's daily digest: a short markdown summary of the ledger, best/ and the open N.

    python3 tools/plan_digest.py [--root .] > planner/digest.md

Only N whose status in data/known_optima_1_1000.csv is "open" are listed; every other N is a proven
optimum and never worked on. The planner (claude -p) reads this, never the raw logs.

Candidate ranking (a default the planner may override):
  waste(N) = side(N) - sqrt(N), the width lost to imperfect packing;
  slack(N) = waste(N) - median(waste of the open neighbours within 5), positive when N packs worse
             than its neighbours -- a bump in the curve, the likeliest place for an improvement;
  score    = slack / median |slack| + 0.5 if a neighbour improved in the last 3 days
             + 0.3 if never tried - 0.1 * ln(1 + tries).
"""
import argparse
import csv
import datetime as dt
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

RANGES = [(1, 81, "N<=81"), (82, 300, "82-300"), (301, 1000, "301-1000")]


def load_status(root):
    with (root / "data" / "known_optima_1_1000.csv").open() as f:
        return {int(r["N"]): r["status"].strip() for r in csv.DictReader(f)}


def open_ns(root):
    """The open N inside the machine's focus range (planner/machine.json n_min, n_max)."""
    m = machine(root)
    return sorted(n for n, s in load_status(root).items() if s == "open" and m["n_min"] <= n <= m["n_max"])


def load_ledger(root):
    out = []
    p = root / "runs" / "ledger.jsonl"
    if p.exists():
        for line in p.read_text().splitlines():
            try:
                out.append(json.loads(line))
            except ValueError:
                pass
    return out


def best_sides(root):
    sides = {}
    for p in (root / "best").glob("*.json"):
        try:
            sides[int(p.stem)] = json.loads(p.read_text())["side"]
        except (ValueError, KeyError, OSError):
            pass
    return sides


def machine(root):
    p = root / "planner" / "machine.json"
    m = {"workers": 2, "budget_hours": 23.3, "cores": 8, "n_min": 1, "n_max": 1000, "require_islands": False}
    if p.exists():
        m.update(json.loads(p.read_text()))
    return m


def when(r):
    try:
        return dt.datetime.strptime(r.get("finished", ""), "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc)
    except ValueError:
        return None


def gain(r):
    return (r["incumbent"] - r["side"]) if r.get("stored") == "best" and r.get("incumbent") not in (None, float("inf")) else 0.0


def rank(root, ledger, now):
    """Open N with side, waste, slack, tries, last tried, last improvement, neighbour flag, score."""
    ns, sides = open_ns(root), best_sides(root)
    waste = {n: sides[n] - math.sqrt(n) for n in ns if n in sides}
    tries, last_try, last_imp, recent_imp = defaultdict(int), {}, {}, set()
    for r in ledger:
        n, t = r.get("N"), when(r)
        tries[n] += 1
        if t:
            last_try[n] = max(last_try.get(n, t), t)
        if r.get("stored") == "best" and t:
            last_imp[n] = max(last_imp.get(n, t), t)
            if now - t < dt.timedelta(days=3):
                recent_imp.add(n)
    slack = {}
    for n in waste:
        nb = [waste[m] for m in range(n - 5, n + 6) if m != n and m in waste]
        slack[n] = waste[n] - statistics.median(nb) if nb else 0.0
    scale = statistics.median([abs(v) for v in slack.values()] or [1.0]) or 1.0
    rows = []
    for n in waste:
        score = slack[n] / scale + 0.5 * bool(recent_imp & {n - 1, n + 1}) + 0.3 * (tries[n] == 0) - 0.1 * math.log1p(tries[n])
        rows.append({"N": n, "side": sides[n], "waste": waste[n], "slack": slack[n], "tries": tries[n],
                     "last_try": last_try.get(n), "last_imp": last_imp.get(n),
                     "nb_imp": bool(recent_imp & {n - 1, n + 1}), "score": score})
    rows.sort(key=lambda r: -r["score"])
    return rows


def ranges(ns):
    """[1, 2, 3, 7, 9, 10] -> "1-3, 7, 9-10"."""
    out, i = [], 0
    while i < len(ns):
        j = i
        while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
            j += 1
        out.append(f"{ns[i]}-{ns[j]}" if j > i else str(ns[i]))
        i = j + 1
    return ", ".join(out)


def fmt_t(t, now):
    if not t:
        return "-"
    d = now - t
    return f"{d.days}d" if d.days else f"{int(d.total_seconds() // 3600)}h"


def digest(root, now=None, top=30):
    now = now or dt.datetime.now(dt.timezone.utc)
    ledger, m = load_ledger(root), machine(root)
    status = load_status(root)
    ns = open_ns(root)
    day = [r for r in ledger if (t := when(r)) and now - t < dt.timedelta(days=1)]
    L = [f"# Planner digest  {now:%Y-%m-%d %H:%M} UTC", ""]
    L += ["## Machine", f"- workers {m['workers']} (GA jobs use 4 threads each), cores {m['cores']}, "
          f"budget until the next plan: {m['budget_hours']:.1f} h of wall clock per worker",
          f"- focus: N from {m['n_min']} to {m['n_max']} only (N outside it are dropped)"
          + ("; exploit and explore jobs must use an island GA (mode ga, ga.islands >= 2); only H<id> "
             "experiments may use other schedules" if m["require_islands"] else ""), ""]
    tried = {r["N"] for r in ledger}
    improved = {r["N"] for r in ledger if r.get("stored") == "best"}
    L += ["## Coverage", f"- open N in focus: {len(ns)} (of {sum(s == 'open' for s in status.values())} open, "
          f"{len(status)} in all; the rest are proven optima or out of focus: never queue them)",
          f"- open N ever tried: {len(tried & set(ns))}, ever improved: {len(improved & set(ns))}",
          f"- the open N (only these may be queued): {ranges(ns)}",
          f"- last 24 h: {len(day)} jobs, {sum(r.get('stored') == 'best' for r in day)} improvements, "
          f"{sum(r.get('wall_s', 0) for r in day) / 3600:.1f} job-hours", ""]

    L += ["## Schedules by N range (all time | last 3 days)", "",
          "| schedule | range | jobs | improved | per job-hour | mean gain | jobs 3d | improved 3d |", "|---|---|---|---|---|---|---|---|"]
    agg = defaultdict(lambda: [0, 0, 0.0, 0.0, 0, 0])
    for r in ledger:
        rg = next((name for lo, hi, name in RANGES if lo <= r.get("N", 0) <= hi), "?")
        a = agg[(Path(r.get("schedule", "?")).stem, rg)]
        a[0] += 1
        a[1] += r.get("stored") == "best"
        a[2] += r.get("wall_s", 0) / 3600
        a[3] += gain(r)
        if (t := when(r)) and now - t < dt.timedelta(days=3):
            a[4] += 1
            a[5] += r.get("stored") == "best"
    for (s, rg), a in sorted(agg.items()):
        L.append(f"| {s} | {rg} | {a[0]} | {a[1]} | {a[1] / a[2] if a[2] else 0:.2f} | "
                 f"{a[3] / a[1] if a[1] else 0:.2e} | {a[4]} | {a[5]} |")
    if not agg:
        L.append("| (no jobs yet) | | | | | | | |")
    L.append("")

    tags = defaultdict(lambda: [0, 0, 0.0])
    for r in ledger:
        if r.get("tag"):
            t = tags[r["tag"]]
            t[0] += 1
            t[1] += r.get("stored") == "best"
            if r.get("valid") and r.get("start_side") and r.get("side"):
                t[2] += r["side"] - r["start_side"]
    L += ["## Results by tag (hypotheses are tagged H<id>)", "", "| tag | jobs | improved | mean side - start side |", "|---|---|---|---|"]
    for k, t in sorted(tags.items()):
        L.append(f"| {k} | {t[0]} | {t[1]} | {t[2] / t[0]:+.2e} |")
    if not tags:
        L.append("| (none yet) | | | |")
    L.append("")

    L += [f"## Candidates: top {top} open N by default score", "",
          "| N | side | waste | slack | tries | last tried | last improved | neighbour improved 3d | score |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in rank(root, ledger, now)[:top]:
        L.append(f"| {r['N']} | {r['side']:.6f} | {r['waste']:.4f} | {r['slack']:+.4f} | {r['tries']} | "
                 f"{fmt_t(r['last_try'], now)} | {fmt_t(r['last_imp'], now)} | {'yes' if r['nb_imp'] else ''} | {r['score']:.2f} |")
    L.append("")

    imps = sorted((r for r in ledger if r.get("stored") == "best"), key=lambda r: r.get("finished", ""))[-15:]
    L += ["## Recent improvements", ""]
    L += [f"- N={r['N']} {r.get('incumbent', 0):.9f} -> {r['side']:.9f} ({Path(r.get('schedule', '')).stem}, "
          f"tag {r.get('tag')}, {r.get('finished', '')})" for r in imps] or ["- none yet"]
    L.append("")

    bad = [r for r in day if r.get("status") != "ok" or not r.get("valid")]
    L += ["## Problems (last 24 h)", f"- {len(bad)} failed or invalid jobs" +
          (": " + ", ".join(f"N={r['N']} {r['status']}" for r in bad[:8]) if bad else ""), ""]

    props = sorted((root / "planner" / "proposals").glob("*.md"))
    L += ["## Move proposals", ""]
    for p in props:
        status = next((l for l in p.read_text().splitlines() if l.startswith("status:")), "status: ?")
        L.append(f"- `{p.stem}`: {status[8:]}")
    if not props:
        L.append("- none")
    L.append("")
    L += ["## Schedules available", ""]
    for p in sorted((root / "schedules").glob("*.json")) + sorted((root / "schedules" / "planner").glob("*.json")):
        try:
            desc = json.loads(p.read_text()).get("description", "")
        except ValueError:
            desc = "(unreadable)"
        L.append(f"- `{p.relative_to(root)}`: {desc[:160]}")
    return "\n".join(L) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent.parent))
    a = ap.parse_args(argv)
    sys.stdout.write(digest(Path(a.root)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
