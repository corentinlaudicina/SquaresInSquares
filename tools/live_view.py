#!/usr/bin/env python3
"""Live view of the running search, from the snapshots the engine writes with --live (every ~2 s).

    src/viz/.venv/bin/python tools/live_view.py                  # every running job, one panel each
    src/viz/.venv/bin/python tools/live_view.py <job-id | file>  # one job: its best and every island champion

The runner passes --live runs/jobs/<id>/live.json to every job. Squares are coloured by tilt (0 degrees
blue, 45 degrees orange); the box is the packing's side. In the one-job view, island champions and
explorers are coloured by distance from the best (grey = the same square, red = moved), and squares
that moved since the previous snapshot are outlined in black. A panel title gives N, the side and its change
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
DIFF = plt.get_cmap("Reds")
PREV = {}  # (job, panel) -> squares at the previous snapshot, to outline what moved since


def corners(x, y, t):
    c, s = math.cos(t), math.sin(t)
    return [(x + c * u - s * v, y + s * u + c * v) for u, v in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))]


def tilt(t):  # 0 for axis-aligned, 1 for 45 degrees
    a = abs(math.degrees(t)) % 90.0
    return min(a, 90.0 - a) / 45.0


def displacement(a, b):
    """For each square of a, the distance to the nearest centre in b (grid lookup, O(N))."""
    cell = {}
    for x, y, _ in b:
        cell.setdefault((int(x), int(y)), []).append((x, y))
    out = []
    for x, y, _ in a:
        c = [p for dx in (-1, 0, 1) for dy in (-1, 0, 1) for p in cell.get((int(x) + dx, int(y) + dy), [])]
        out.append(min(math.dist((x, y), p) for p in c) if c else 1.0)
    return out


def draw(ax, squares, side, title, ref=None, key=None):
    """Tilt colours; with `ref`, colour by distance from ref (grey = same, red = moved) and outline the
    squares that moved since this panel's previous snapshot."""
    ax.clear()
    polys = [corners(*q) for q in squares]
    if ref is None:
        faces = [CMAP(0.15 + 0.7 * tilt(q[2])) for q in squares]
    else:
        faces = [("#d9d9d6" if d < 0.02 else DIFF(0.3 + 0.7 * min(1.0, d / 0.6))) for d in displacement(squares, ref)]
    edges, widths = "#333333", 0.3
    if key is not None:
        prev = PREV.get(key)
        if prev is not None and len(prev) == len(squares):
            moved = [math.dist(q[:2], p[:2]) > 0.01 or abs(q[2] - p[2]) > 0.01 for q, p in zip(squares, prev)]
            edges = ["black" if m else "#333333" for m in moved]
            widths = [1.4 if m else 0.3 for m in moved]
        PREV[key] = squares
    ax.add_collection(PolyCollection(polys, facecolors=faces, edgecolors=edges, linewidths=widths))
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
    isl = s.get("islands", [])
    differ = sum(sum(x > 0.05 for x in displacement(i["squares"], s["best"]["squares"])) > 5 for i in isl)
    expl = sum("explorer" in i for i in isl)
    return (f"N={s['n']}  {name.split('_')[2] if name.count('_') >= 3 else name}\n"
            f"side {s['best']['side']:.9f} ({d:+.2e} since start)\n"
            f"generation {s['generation']}, {s['t']:.0f} s"
            + (f"\nislands differing from the best: {differ}/{len(isl)}, explorers {expl}" if isl else ""))


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
    best = s["best"]["squares"]
    panels = [("best", best, s["best"]["side"], job_title(p.parent.name, s), None)]
    for k, i in enumerate(s.get("islands", [])):
        panels.append((f"isl{k}", i["squares"], i["h"],
                       f"island {k}: {i['source']}  (champion)\nh {i['h']:.6f}, step {i['eps']:.0e}, ratchets {i['ratchets']}\n"
                       f"penalty {i['cost']:.1e}, depth {i['depth']:.1e}", best))
        if "explorer" in i:
            e = i["explorer"]
            fit = f", fitted {e['fit_h']:.5f}" if e.get("fit_h", math.inf) < math.inf else ", not fitted yet"
            eh = e.get("h", i["h"])  # an explorer has its own box; compare in the same scale or everything looks moved
            ref = [(x * eh / i["h"], y * eh / i["h"], t) for x, y, t in i["squares"]]
            panels.append((f"exp{k}", e["squares"], eh,
                           f"island {k}: explorer, own box {eh:.5f}{fit}\npenalty {e['cost']:.1e}", ref))
    cols = min(len(panels), 5)
    rows = math.ceil(len(panels) / cols)
    for k, (key, sq, side, title, ref) in enumerate(panels):
        draw(fig.add_subplot(rows, cols, k + 1), sq, side, title, ref=ref, key=(p.parent.name, key))


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
