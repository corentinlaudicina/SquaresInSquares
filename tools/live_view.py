#!/usr/bin/env python3
"""Live view of the running search, from the snapshots the engine writes with --live (every ~2 s).

    src/viz/.venv/bin/python tools/live_view.py                  # every running job, one panel each
    src/viz/.venv/bin/python tools/live_view.py <job-id | file>  # one job: its best and every island champion

The runner passes --live runs/jobs/<id>/live.json to every job. Squares are coloured by tilt (0 degrees
blue, 45 degrees orange); the box is the packing's side. A panel title gives N, the side and its change
since the job started, the generation and elapsed time. Close the window to stop.
"""
import json
import math
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

ROOT = Path(__file__).resolve().parent.parent
JOBS = ROOT / "runs" / "jobs"
FRESH = 30.0  # a snapshot older than this belongs to a job that is no longer running
CMAP = plt.get_cmap("coolwarm")


def corners(x, y, t):
    c, s = math.cos(t), math.sin(t)
    return [(x + c * u - s * v, y + s * u + c * v) for u, v in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))]


def tilt(t):  # 0 for axis-aligned, 1 for 45 degrees
    a = abs(math.degrees(t)) % 90.0
    return min(a, 90.0 - a) / 45.0


def draw(ax, squares, side, title):
    ax.clear()
    polys = [corners(*q) for q in squares]
    ax.add_collection(PolyCollection(polys, facecolors=[CMAP(0.15 + 0.7 * tilt(q[2])) for q in squares],
                                     edgecolors="#333333", linewidths=0.3))
    ax.add_patch(plt.Rectangle((0, 0), side, side, fill=False, edgecolor="black", linewidth=1.0))
    pad = 0.02 * side
    ax.set_xlim(-pad, side + pad)
    ax.set_ylim(-pad, side + pad)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=8)


def load(p):
    try:
        return json.loads(p.read_text())
    except (OSError, ValueError):
        return None  # being replaced right now: try again next tick


def running():
    now = time.time()
    return sorted((p for p in JOBS.glob("*/live.json") if now - p.stat().st_mtime < FRESH), key=lambda p: p.parent.name)


def job_title(name, s):
    d = s["best"]["side"] - s["start_side"]
    return (f"N={s['n']}  {name.split('_')[2] if name.count('_') >= 3 else name}\n"
            f"side {s['best']['side']:.9f} ({d:+.2e} since start)\n"
            f"generation {s['generation']}, {s['t']:.0f} s")


def view_all(fig):
    files = running()
    fig.clf()
    if not files:
        fig.text(0.5, 0.5, "no running job with a live snapshot\n(jobs started by the runner write runs/jobs/<id>/live.json)",
                 ha="center", va="center")
        return
    cols = min(len(files), 3)
    rows = math.ceil(len(files) / cols)
    for k, p in enumerate(files):
        s = load(p)
        if s:
            draw(fig.add_subplot(rows, cols, k + 1), s["best"]["squares"], s["best"]["side"], job_title(p.parent.name, s))


def view_one(fig, p):
    s = load(p)
    fig.clf()
    if not s:
        return
    panels = [("best (legal)", s["best"]["squares"], s["best"]["side"], "")] + [
        (f"island {k}: {i['source']}", i["squares"], i["h"],
         f"h {i['h']:.6f}, step {i['eps']:.0e}, ratchets {i['ratchets']}\npenalty {i['cost']:.1e}, depth {i['depth']:.1e}")
        for k, i in enumerate(s.get("islands", []))]
    cols = min(len(panels), 4)
    rows = math.ceil(len(panels) / cols)
    for k, (name, sq, side, extra) in enumerate(panels):
        title = (job_title(p.parent.name, s) if k == 0 else name + "\n" + extra)
        draw(fig.add_subplot(rows, cols, k + 1), sq, side, title)


def main():
    target = None
    if len(sys.argv) > 1:
        a = Path(sys.argv[1])
        target = a if a.is_file() else JOBS / sys.argv[1] / "live.json"
    fig = plt.figure("squares: live", figsize=(12, 8))
    seen = None
    while plt.fignum_exists(fig.number):
        files = [target] if target else running()
        stamp = tuple((str(p), p.stat().st_mtime) for p in files if p.exists())
        if stamp != seen:
            seen = stamp
            view_one(fig, target) if target else view_all(fig)
            fig.tight_layout()
            fig.canvas.draw_idle()
        plt.pause(2.0)


if __name__ == "__main__":
    main()
