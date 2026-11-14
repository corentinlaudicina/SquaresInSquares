#!/usr/bin/env python3
"""Converts the engine's text records (data/init_external/*.cfg, data/best/*.log) to sqpack/1 JSON.

    python3 tools/convert.py data/init_external/n0029.cfg -o n29.json

Unit squares only: a record with another species, or a frozen flag, is refused (exit 1).
"""
import argparse
import json
import sys

UNIT = {(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)}


def parse_record(text, source=""):
    n, side, verts, squares = None, None, [], []
    for line in text.splitlines():
        f = line.split()
        if not f or f[0].startswith("#"):
            continue
        if f[0] == "n":
            n = int(f[1])
        elif f[0] == "L":
            side = float(f[1])
        elif f[0] == "s" and f[1] != "0":
            raise ValueError(f"{source}: more than one species")
        elif f[0] == "v":
            verts.append((float(f[1]), float(f[2])))
        elif f[0] == "p":
            if f[4] != "0" or f[5] != "0":
                raise ValueError(f"{source}: frozen or non-base-species particle")
            squares.append([float(f[1]), float(f[2]), float(f[3])])
        elif f[0] == "f":
            raise ValueError(f"{source}: trajectory file, not a record")
    if set(verts) != UNIT:
        raise ValueError(f"{source}: species 0 is not the unit square")
    if n is None or side is None:
        raise ValueError(f"{source}: missing n or L")
    return {"format": "sqpack/1", "n": n, "side": side, "squares": squares,
            "meta": {"converted_from": source}}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("record")
    ap.add_argument("-o", "--out", help="output path (default: stdout)")
    a = ap.parse_args(argv)
    try:
        with open(a.record) as f:
            doc = parse_record(f.read(), a.record)
    except (OSError, ValueError, IndexError) as e:
        print(f"convert: {e}", file=sys.stderr)
        return 1
    text = json.dumps(doc, indent=None, separators=(",", ":"))
    if a.out:
        with open(a.out, "w") as f:
            f.write(text + "\n")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
