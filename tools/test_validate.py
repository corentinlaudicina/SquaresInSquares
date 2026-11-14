#!/usr/bin/env python3
"""Tests for tools/validate.py and tools/convert.py.   python3 -m unittest tools/test_validate.py -v"""
import json
import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))
import convert  # noqa: E402
import validate  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def doc(squares, side, n=None):
    return {"format": "sqpack/1", "n": len(squares) if n is None else n, "side": side, "squares": squares}


def grid(k):
    return [[i + 0.5, j + 0.5, 0.0] for i in range(k) for j in range(k)]


class Valid(unittest.TestCase):
    def test_touching_grid_is_valid(self):
        ok, why, d = validate.validate(doc(grid(3), 3.0))
        self.assertTrue(ok, why)
        self.assertEqual(d["tight_side"], 3.0)

    def test_n5_record_shape(self):
        # s(5) = 2 + 1/sqrt(2): four corner squares and one tilted 45 degrees in the middle,
        # shrunk by 1e-9 so the float representation is strictly clear.
        s = 2 + 1 / math.sqrt(2) + 1e-9
        sq = [[0.5, 0.5, 0], [s - 0.5, 0.5, 0], [0.5, s - 0.5, 0], [s - 0.5, s - 0.5, 0],
              [s / 2, s / 2, math.pi / 4]]
        ok, why, _ = validate.validate(doc(sq, s))
        self.assertTrue(ok, why)

    def test_quarter_turn_needs_slack(self):
        # theta = 2 pi is "the same square" in the reals, but sin(2 pi) is -2.4e-16 in float, so the
        # square the file describes is tilted and touches its neighbours with ~1e-16 overlap. Strict
        # means strict: invalid when touching, valid once there is slack. Legalization owns this.
        sq = grid(2)
        sq[0][2] = 2 * math.pi
        self.assertFalse(validate.validate(doc(sq, 2.0))[0])
        sq = [[x * (1 + 1e-12), y * (1 + 1e-12), t] for x, y, t in sq]
        ok, why, _ = validate.validate(doc(sq, 2.0 + 1e-11))
        self.assertTrue(ok, why)

    def test_expected_n_matches(self):
        self.assertTrue(validate.validate(doc(grid(2), 2.0), expected_n=4)[0])


class Invalid(unittest.TestCase):
    def assertInvalid(self, d, fragment, **kw):
        ok, why, _ = validate.validate(d, **kw)
        self.assertFalse(ok)
        self.assertTrue(any(fragment in w for w in why), why)

    def test_overlap(self):
        sq = grid(2)
        sq[1][1] -= 1e-12  # square 1 sits on square 0: 1e-12 penetration, below any float tolerance a sloppy checker would use
        self.assertInvalid(doc(sq, 2.0), "overlap")

    def test_tilted_corner_overlap(self):
        # a 45-degree square whose corner pokes 1e-10 into its neighbour's face
        d = 0.5 + math.sqrt(2) / 2 - 1e-10
        sq = [[0.5, 1.5, 0.0], [0.5 + d, 1.5, math.pi / 4]]
        self.assertInvalid(doc(sq, 3.0), "overlap")

    def test_out_of_bounds(self):
        sq = grid(2)
        sq[3][0] += 1e-12
        self.assertInvalid(doc(sq, 2.0), "outside")

    def test_negative_coordinate(self):
        sq = grid(2)
        sq[0][1] = 0.4999999
        self.assertInvalid(doc(sq, 2.0), "outside")

    def test_wrong_n(self):
        self.assertInvalid(doc(grid(2), 2.0, n=5), "squares listed")
        self.assertInvalid(doc(grid(2), 2.0), "expected 9", expected_n=9)

    def test_nan_and_inf(self):
        sq = grid(2)
        sq[2][2] = float("nan")
        self.assertInvalid(doc(sq, 2.0), "non-finite")
        self.assertInvalid(doc(grid(2), float("inf")), "side")

    def test_nan_token_in_file(self):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
            f.write('{"format":"sqpack/1","n":1,"side":1,"squares":[[0.5,NaN,0]]}')
        try:
            self.assertEqual(validate.main([f.name]), 1)
        finally:
            os.unlink(f.name)

    def test_malformed(self):
        self.assertInvalid(doc([[0.5, 0.5]], 1.0), "not [x, y, theta]")
        self.assertInvalid(doc([[0.5, 0.5, True]], 1.0), "not [x, y, theta]")
        self.assertInvalid({**doc(grid(1), 1.0), "format": "v0"}, "format")
        self.assertInvalid([], "not an object")

    def test_below_area_bound(self):
        self.assertInvalid(doc(grid(2), 1.9), "area bound")

    def test_identical_squares(self):
        self.assertInvalid(doc([[1, 1, 0], [1, 1, 0]], 2.0), "overlap")

    def test_unreadable_file(self):
        self.assertEqual(validate.main(["/nonexistent/x.json"]), 2)


class Convert(unittest.TestCase):
    def test_seed_roundtrip(self):
        path = os.path.join(ROOT, "data/init_external/n0005.cfg")
        if not os.path.exists(path):
            self.skipTest("no seeds")
        with open(path) as f:
            d = convert.parse_record(f.read(), path)
        self.assertEqual(d["n"], 5)
        json.loads(json.dumps(d))  # serialisable
        ok, why, _ = validate.validate(d, expected_n=5)
        print(f"\n  seed n5: valid={ok} side={d['side']} {why[:1]}")

    def test_refuses_other_species(self):
        with self.assertRaises(ValueError):
            convert.parse_record("n 1\nL 1\ns 0 3\nv 0 0\nv 1 0\nv 0 1\np .5 .5 0 0 0\n")


if __name__ == "__main__":
    unittest.main()
