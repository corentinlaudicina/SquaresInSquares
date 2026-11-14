#!/usr/bin/env python3
"""Open a packing in a matplotlib window.

    python show.py 11                  # data/init_external/n0011.cfg
    python show.py <file>              # any record, such as one from data/best
    python show.py --watch 3           # follow a run: logs/n0003.traj
    python show.py --watch <file.traj>

--watch can start before the run does, replays the frames already written, then keeps following.
A large move or symmetry operation is held on screen for a moment, then its relaxation plays.
Move the cursor over the packing to read the nearest square's details in the toolbar.
"""

import math
import sys
import time
from collections import deque
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[2]

SURFACE, INK, EDGE = "#fcfcfb", "#1a1a19", "#8a8a86"

STREAMING = {"", "relax", "relax free", "seam relax", "compress", "result"}  # played straight through
HOLD_TICKS = 10  # every other label is an event: stay on it this many ticks


def read_record(path):
    """Returns (L, species vertex lists, particles as (x, y, theta, species))."""
    L, species, particles = None, [], []
    for line in path.read_text().splitlines():
        f = line.split()
        if not f or f[0].startswith("#"):
            continue
        if f[0] == "L":
            L = float(f[1])
        elif f[0] == "s":
            species.append([])
        elif f[0] == "v":
            species[-1].append((float(f[1]), float(f[2])))
        elif f[0] == "p":
            particles.append((float(f[1]), float(f[2]), float(f[3]), int(f[5])))
    if L is None or not particles:
        sys.exit(f"{path}: no box side or no particles found")
    return L, species, particles


def read_new_frames(path, state):
    """Returns (frames finished since the last call, whether this is a new file).

    A frame is (step, L, energy, [x, y, theta, ...], label). `state` starts as {} and remembers the
    place. A last line with no newline yet is a frame still being written, so it waits.
    """
    try:
        inode = path.stat().st_ino
    except FileNotFoundError:
        return [], False
    restarted = inode != state.get("inode")
    if restarted:
        state.update(inode=inode, offset=0, label="")
    with path.open("rb") as f:
        f.seek(state["offset"])
        data = f.read()
    finished = data[: data.rfind(b"\n") + 1]
    state["offset"] += len(finished)

    frames = []
    for line in finished.splitlines():
        if line.startswith(b"m "):
            state["label"] = line[2:].decode()
        elif line.startswith(b"f "):
            w = line.split()
            frames.append((int(w[1]), float(w[2]), float(w[3]), [float(v) for v in w[4:]], state["label"]))
            state["label"] = ""
    return frames, restarted


def corners(particle, species):
    x, y, theta, sp = particle
    c, s = math.cos(theta), math.sin(theta)
    return [(x + vx * c - vy * s, y + vx * s + vy * c) for vx, vy in species[sp]]


def area(vertices):
    return 0.5 * sum(
        a[0] * b[1] - a[1] * b[0] for a, b in zip(vertices, vertices[1:] + vertices[:1])
    )


def tilt_degrees(theta):
    """Rotation away from the nearest axis, 0 to 45: a square looks the same every 90 degrees."""
    return abs(math.degrees(math.remainder(theta, math.pi / 2)))


def figure(name, L, species, particles):
    """Returns (fig, update); update(L, particles, note) redraws it with a new configuration."""
    particles = list(particles)
    covered = sum(area(species[p[3]]) for p in particles)

    fig, ax = plt.subplots(figsize=(8, 8.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    squares = PolyCollection(
        [corners(p, species) for p in particles],
        array=[tilt_degrees(p[2]) for p in particles],
        cmap="inferno",
        clim=(0.0, 45.0),
        edgecolors=EDGE,  # gray: the colormap's ends are near black and near the page's white
        linewidths=0.6,
    )
    ax.add_collection(squares)
    box = Rectangle((0, 0), L, L, fill=False, edgecolor=INK, linewidth=1.5)
    ax.add_patch(box)
    ax.set_aspect("equal")
    ax.set_axis_off()
    title = ax.set_title("", loc="left", color=INK)
    bar = fig.colorbar(squares, ax=ax, orientation="horizontal", shrink=0.6, pad=0.03)
    bar.set_label("tilt from the nearest axis (degrees)", color=INK)
    bar.set_ticks([0, 15, 30, 45])

    def update(L, new_particles, note=""):
        particles[:] = new_particles
        shapes = [corners(p, species) for p in particles]
        squares.set_verts(shapes)
        squares.set_array([tilt_degrees(p[2]) for p in particles])
        box.set_bounds(0, 0, L, L)
        # Span the squares as well as the box, so a square that left the box stays on screen.
        xs = [x for shape in shapes for x, _ in shape] + [0.0, L]
        ys = [y for shape in shapes for _, y in shape] + [0.0, L]
        pad = 0.02 * max(max(xs) - min(xs), max(ys) - min(ys))
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
        ax.set_ylim(min(ys) - pad, max(ys) + pad)
        title.set_text(
            f"n = {len(particles)}    side {L:.6f}    density {100 * covered / (L * L):.2f}%{note}"
        )
        fig.canvas.draw_idle()

    def describe(x, y):
        i = min(range(len(particles)), key=lambda k: math.dist((x, y), particles[k][:2]))
        p = particles[i]
        return f"square {i}: x={p[0]:.6f}, y={p[1]:.6f}, angle={math.degrees(p[2]) % 360:.3f} deg"

    ax.format_coord = describe
    fig.canvas.manager.set_window_title(name)
    update(L, particles)
    return fig, update


def follow(path):
    """Waits for the first frame. Returns (fig, tick); each tick shows at most one new frame."""
    state, pending = {}, deque()
    announced = False
    while not pending:
        pending.extend(read_new_frames(path, state)[0])
        if not pending:
            if not announced:
                print(f"waiting for {path} ...")
                announced = True
            time.sleep(0.2)

    # shortcut: a new run at this path is assumed to have the same shapes, restart the viewer if not.
    L, species, particles = read_record(path)
    fig, update = figure(path.name, L, species, particles)
    hold = 0

    def tick():
        nonlocal hold
        frames, restarted = read_new_frames(path, state)
        if restarted:
            pending.clear()
            hold = 0
        pending.extend(frames)
        if hold > 0:
            hold -= 1
        elif pending:
            step, L, energy, values, label = pending.popleft()
            hold = 0 if label in STREAMING else HOLD_TICKS
            if len(values) == 3 * len(particles):
                poses = [
                    (values[3 * i], values[3 * i + 1], values[3 * i + 2], p[3])
                    for i, p in enumerate(particles)
                ]
                lead = f"\n{label}" if label else ""
                update(L, poses, f"{lead}    step {step}    energy {energy:.4g}")

    return fig, tick


def main():
    args = sys.argv[1:]
    watch = args[:1] == ["--watch"]
    if watch:
        args = args[1:]
    if len(args) != 1:
        sys.exit(__doc__)

    arg = args[0]
    if arg.isdigit():
        number = int(arg)
        path = ROOT / "logs" / f"n{number:04d}.traj" if watch else ROOT / "data/init_external" / f"n{number:04d}.cfg"
    else:
        path = Path(arg)

    if watch:
        fig, tick = follow(path)
        timer = fig.canvas.new_timer(interval=100)
        timer.add_callback(tick)
        timer.start()
    else:
        if not path.is_file():
            sys.exit(f"no such file: {path}")
        figure(path.name, *read_record(path))
    plt.show()


if __name__ == "__main__":
    main()
