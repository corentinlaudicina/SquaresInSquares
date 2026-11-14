#!/usr/bin/env python3
"""Independent validator for sqpack/1 solution files. Shares no code with the C++ engine.

    python3 tools/validate.py best/29.json [--n 29]      # exit 0 valid, 1 invalid, 2 unreadable

A file is valid when every unit square lies inside [0, side]^2 and no two squares overlap with
positive area (touching is allowed). The check is exact on the squares the file describes: each
square is centre + R * (+-1/2, +-1/2) where R = [[c, -s], [s, c]] and c, s are the float64 values
of cos(theta), sin(theta). Those floats are exact rationals, so the corners and the separating-axis
test are computed with fractions.Fraction and carry no rounding. (|det R| = c^2 + s^2 is 1 to within
~2e-16, the only approximation, and far below any tolerance a packing record cares about.)

File format (sqpack/1):
    {"format": "sqpack/1", "n": <int>, "side": <float>,
     "squares": [[x, y, theta], ...],        # centres, radians
     "meta": {...}}                           # free-form provenance, ignored here
"""
import argparse
import json
import math
import sys
from fractions import Fraction

FORMAT = "sqpack/1"
# Float pre-filter: a pair (or a wall) whose float gap exceeds this is certainly clear, because the
# float corner error is ~1e-15 at side 1000. Anything closer is decided exactly.
SURE_GAP = 1e-9
HALF = Fraction(1, 2)


def _corners_exact(x, y, t):
    c, s = Fraction(math.cos(t)), Fraction(math.sin(t))
    X, Y = Fraction(x), Fraction(y)
    return [(X + c * u - s * v, Y + s * u + c * v) for u, v in ((-HALF, -HALF), (HALF, -HALF), (HALF, HALF), (-HALF, HALF))]


def _corners_float(x, y, t):
    c, s = math.cos(t), math.sin(t)
    return [(x + c * u - s * v, y + s * u + c * v) for u, v in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))]


def _max_axis_gap(a, b):
    """Largest separation over the 4 edge normals of both quads (positive = separated on that axis).

    Works for exact (Fraction) and float corners alike. Normals are not normalised: only the sign
    matters for the exact test, and the float pre-filter uses unit squares, so |normal| ~ 1."""
    best = None
    for poly in (a, b):
        for k in range(2):  # a parallelogram has two edge directions
            ex, ey = poly[k + 1][0] - poly[k][0], poly[k + 1][1] - poly[k][1]
            nx, ny = -ey, ex
            pa = [p[0] * nx + p[1] * ny for p in a]
            pb = [p[0] * nx + p[1] * ny for p in b]
            gap = max(min(pb) - max(pa), min(pa) - max(pb))
            best = gap if best is None or gap > best else best
    return best


def _pairs_near(centres, reach):
    """Candidate pairs whose centres are within `reach`, via a uniform grid."""
    grid = {}
    for i, (x, y) in enumerate(centres):
        grid.setdefault((math.floor(x / reach), math.floor(y / reach)), []).append(i)
    for (gx, gy), members in grid.items():
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j in grid.get((gx + dx, gy + dy), ()):
                    for i in members:
                        if i < j and math.dist(centres[i], centres[j]) <= reach:
                            yield i, j


def validate(doc, expected_n=None):
    """Returns (ok, reasons, diagnostics). Never raises on bad content."""
    reasons, diag = [], {}
    if not isinstance(doc, dict):
        return False, ["top level is not an object"], diag
    if doc.get("format") != FORMAT:
        reasons.append(f"format is {doc.get('format')!r}, expected {FORMAT!r}")
    n, side, squares = doc.get("n"), doc.get("side"), doc.get("squares")
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        reasons.append("n is not a positive integer")
    if expected_n is not None and n != expected_n:
        reasons.append(f"n is {n}, expected {expected_n}")
    if not isinstance(side, (int, float)) or isinstance(side, bool) or not math.isfinite(side) or side <= 0:
        reasons.append("side is not a positive finite number")
    if not isinstance(squares, list):
        reasons.append("squares is not a list")
        return False, reasons, diag
    if isinstance(n, int) and len(squares) != n:
        reasons.append(f"{len(squares)} squares listed, n is {n}")
    for k, sq in enumerate(squares):
        if (not isinstance(sq, list) or len(sq) != 3
                or not all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in sq)):
            reasons.append(f"square {k} is not [x, y, theta]")
        elif not all(math.isfinite(v) for v in sq):
            reasons.append(f"square {k} has a non-finite value")
    if reasons:
        return False, reasons, diag

    m = len(squares)
    if side < math.sqrt(m) * (1 - 1e-12):
        reasons.append(f"side {side!r} is below the area bound sqrt(n) = {math.sqrt(m)!r}")

    fl = [_corners_float(*sq) for sq in squares]
    xs = [p[0] for q in fl for p in q]
    ys = [p[1] for q in fl for p in q]
    diag["tight_side"] = max(max(xs) - min(xs), max(ys) - min(ys))

    S = Fraction(side)
    for k, (q, sq) in enumerate(zip(fl, squares)):
        lo = min(min(p[0] for p in q), min(p[1] for p in q))
        hi = max(max(p[0] for p in q), max(p[1] for p in q))
        if lo > SURE_GAP and hi < side - SURE_GAP:
            continue
        ex = _corners_exact(*sq)
        if any(p[0] < 0 or p[1] < 0 or p[0] > S or p[1] > S for p in ex):
            reasons.append(f"square {k} is outside [0, side]^2")

    worst = math.inf  # smallest float gap seen over near pairs; negative = penetration depth
    overlaps = 0
    for i, j in _pairs_near([(sq[0], sq[1]) for sq in squares], math.sqrt(2) + 1e-6):
        g = _max_axis_gap(fl[i], fl[j])
        worst = min(worst, g)
        if g > SURE_GAP:
            continue
        if _max_axis_gap(_corners_exact(*squares[i]), _corners_exact(*squares[j])) < 0:
            overlaps += 1
            if overlaps <= 5:
                reasons.append(f"squares {i} and {j} overlap (float depth {-g:.3e})")
    if overlaps > 5:
        reasons.append(f"... {overlaps} overlapping pairs in total")
    diag["min_pair_gap"] = None if worst == math.inf else worst
    return not reasons, reasons, diag


def load(path):
    """json.load, but NaN/Infinity tokens are kept (as floats) so validate() can reject them."""
    with open(path) as f:
        return json.load(f)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("file")
    ap.add_argument("--n", type=int, default=None, help="expected square count")
    a = ap.parse_args(argv)
    try:
        doc = load(a.file)
    except (OSError, ValueError) as e:
        print(json.dumps({"file": a.file, "valid": False, "reasons": [f"unreadable: {e}"]}))
        return 2
    ok, reasons, diag = validate(doc, a.n)
    print(json.dumps({"file": a.file, "valid": ok, "reasons": reasons, **diag}))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
