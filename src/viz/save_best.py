#!/usr/bin/env python3
"""Saves a PNG of the best packing held for every n, from data/best to src/viz/best/n0037.png.

    src/viz/.venv/bin/python src/viz/save_best.py
"""

import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

import show  # noqa: E402

RECORD = re.compile(r"_n(\d+)_s(\d+\.\d+)\.log$")


def best_records(directory):
    """Returns {n: path of the record with the smallest side}, read off the file names."""
    best = {}
    for path in directory.glob("*.log"):
        match = RECORD.search(path.name)
        if match:
            n, side = int(match[1]), float(match[2])
            if n not in best or side < best[n][0]:
                best[n] = (side, path)
    return {n: path for n, (_, path) in best.items()}


def main():
    out = show.ROOT / "src/viz/best"
    out.mkdir(parents=True, exist_ok=True)
    records = best_records(show.ROOT / "data/best")
    for n, path in sorted(records.items()):
        fig, _ = show.figure(path.name, *show.read_record(path))
        fig.savefig(out / f"n{n:04d}.png")
        plt.close(fig)
    print(f"saved {len(records)} images to {out}")


if __name__ == "__main__":
    main()
