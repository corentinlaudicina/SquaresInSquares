#!/usr/bin/env python3
"""Engine-level checks of the Phase 3 item 2 paths.   python3 -m unittest tools/test_engine.py"""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
import validate  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
ENGINE = ROOT / "build" / "engine"


def run(schedule, n=10, kicks=6, seed=1):
    tmp = tempfile.mkdtemp()
    out = Path(tmp) / "o.json"
    subprocess.run([str(ENGINE), "--N", str(n), "--seed", str(seed), "--kicks", str(kicks), "--schedule",
                    str(ROOT / "schedules" / schedule), "--out", str(out)], cwd=ROOT, check=True, capture_output=True)
    return validate.load(out), json.loads((Path(tmp) / "o.telemetry.json").read_text())


@unittest.skipUnless(ENGINE.exists(), "make engine first")
class EnginePaths(unittest.TestCase):
    def test_accept_legal_reports_the_legal_side(self):
        doc, tel = run("shrink-legal.json")
        self.assertTrue(validate.validate(doc, 10)[0])
        # the search's best is already legal, so the output's side is that best (refit adds ulps at most)
        self.assertAlmostEqual(doc["side"], tel["side"]["best"], delta=1e-12)
        self.assertGreaterEqual(tel["legal"]["calls"], 1)  # the start, at least

    def test_capped_relaxer_runs_and_validates(self):
        doc, tel = run("shrink-capped.json")
        self.assertTrue(validate.validate(doc, 10)[0])
        self.assertEqual(sum(tel["relax"]["stops"].values()), tel["relax"]["calls"])
        self.assertTrue(set(tel["relax"]["stops"]) <= {"feasible", "not reached"})

    def test_same_seed_same_result(self):
        a, _ = run("shrink-capped.json", seed=3)
        b, _ = run("shrink-capped.json", seed=3)
        self.assertEqual(a["squares"], b["squares"])


if __name__ == "__main__":
    unittest.main()
