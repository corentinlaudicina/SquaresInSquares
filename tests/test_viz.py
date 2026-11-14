#!/usr/bin/env python3
"""Checks how show.py follows a growing trajectory. Run with the viewer's Python:

    src/viz/.venv/bin/python tests/test_viz.py
"""

import os
import sys
import tempfile
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src/viz"))
import show  # noqa: E402

with tempfile.TemporaryDirectory() as d:
    path = Path(d) / "run.traj"
    state = {}
    assert show.read_new_frames(path, state) == ([], False)  # not started yet

    path.write_bytes(b"# header\nL 3\n")
    assert show.read_new_frames(path, state) == ([], True)  # a header alone is no frame
    assert show.read_new_frames(path, state) == ([], False)

    with path.open("ab") as f:
        f.write(b"f 0 3 1.5 1 2 0.5\nf 10 3 1.")
        f.flush()
        frames, restarted = show.read_new_frames(path, state)
        assert frames == [(0, 3.0, 1.5, [1.0, 2.0, 0.5], "")] and not restarted  # second is unfinished

        f.write(b"25 4 5 0.75\n")
    frames, _ = show.read_new_frames(path, state)
    assert frames == [(10, 3.0, 1.25, [4.0, 5.0, 0.75], "")]
    assert show.read_new_frames(path, state) == ([], False)

    # A new run replaces the file, even with one larger than what has been read so far.
    path.unlink()
    path.write_bytes(b"L 3\nf 0 3 9 1 1 0\n" + b"f 5 3 8 1 1 0\n" * 20)
    frames, restarted = show.read_new_frames(path, state)
    assert restarted and len(frames) == 21 and frames[0][2] == 9.0

    # A label belongs to the next frame only, even if that frame arrives in a later read.
    with path.open("ab") as f:
        f.write(b"m symmetrize 2xdiagonal\n")
        f.flush()
        assert show.read_new_frames(path, state) == ([], False)
        f.write(b"f 6 3 7 1 1 0\nf 7 3 7 1 1 0\n")
    frames, _ = show.read_new_frames(path, state)
    assert [x[4] for x in frames] == ["symmetrize 2xdiagonal", ""]

assert show.tilt_degrees(0.0) == 0.0
print("viz: all tests passed")
