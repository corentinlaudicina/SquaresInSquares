#!/usr/bin/env python3
"""benchmark.report() on synthetic results.   python3 -m unittest tools/test_benchmark.py"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))
import benchmark  # noqa: E402


def r(sched, n, seed, side, target=10.0, valid=True, ttt=None):
    return {"schedule": sched, "N": n, "seed": seed, "side": side, "target": target, "valid": valid, "ttt": ttt,
            "success": valid and side <= target * (1 + 1e-9)}


class Report(unittest.TestCase):
    def test_success_pairs_and_per_n(self):
        res = [r("a", 7, 1, 10.0, ttt=3.0), r("a", 7, 2, 10.5), r("b", 7, 1, 10.2), r("b", 7, 2, 10.1),
               r("a", 8, 1, 11.0, target=11.0, ttt=1.0), r("b", 8, 1, 11.0, target=11.0, ttt=5.0)]
        text = benchmark.report(res, ["x/a.json", "x/b.json"], [7, 8], [1, 2], 60)
        self.assertIn("| a | 2/3 | 2.0 |", text)          # median ttt of the successful runs (3.0, 1.0)
        self.assertIn("| b | 1/3 | 5.0 |", text)
        self.assertIn("`b`: 1 wins, 1 losses, 1 ties", text)
        self.assertIn("| 8 | 11.000000 | 11.000000* | 11.000000* |", text)

    def test_invalid_runs_are_counted_not_averaged(self):
        res = [r("a", 7, 1, 9.0, valid=False), r("a", 7, 2, 10.0)]
        text = benchmark.report(res, ["a.json"], [7], [1, 2], 60)
        self.assertIn("| a | 1/2 | - | 10.000000 | +0.000e+00 | 1 |", text)


if __name__ == "__main__":
    unittest.main()
