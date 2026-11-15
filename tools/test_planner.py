#!/usr/bin/env python3
"""Planner tools: digest and plan_check on a synthetic root.   python3 -m unittest tools/test_planner.py"""
import datetime as dt
import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
import plan_check  # noqa: E402
import plan_digest  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
ENGINE = REPO / "build" / "engine"
OPEN = [5, 6, 7, 8, 9, 90, 91]


class Planner(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        (self.root / "data").mkdir()
        rows = ["N,status,side_length,basis"] + [f"{n},{'open' if n in OPEN else 'A_classical'},," for n in range(1, 101)]
        (self.root / "data" / "known_optima_1_1000.csv").write_text("\n".join(rows) + "\n")
        (self.root / "best").mkdir()
        for n in range(1, 101):
            side = n ** 0.5 + (0.3 if n == 7 else 0.1)  # N = 7 packs worse than its neighbours
            (self.root / "best" / f"{n}.json").write_text(json.dumps({"side": side}))
        (self.root / "schedules" / "planner").mkdir(parents=True)
        for s in ("ga", "ga-seeded", "shrink-capped"):
            shutil.copy(REPO / "schedules" / f"{s}.json", self.root / "schedules")
        (self.root / "planner").mkdir()
        (self.root / "planner" / "machine.json").write_text(json.dumps({"workers": 2, "cores": 8, "budget_hours": 1.0}))
        (self.root / "runs").mkdir()
        now = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        ledger = [{"id": "a", "N": 6, "schedule": "schedules/ga.json", "status": "ok", "valid": True, "stored": "best",
                   "incumbent": 2.6, "side": 2.55, "wall_s": 600, "finished": now, "tag": "exploit", "start_side": 2.6}]
        (self.root / "runs" / "ledger.jsonl").write_text("".join(json.dumps(r) + "\n" for r in ledger))
        self.today = dt.date(2026, 11, 16)

    def tearDown(self):
        shutil.rmtree(self.root)

    def check(self, plan, engine=None):
        return plan_check.check_and_build(plan, self.root, engine, self.today)

    def test_digest_lists_only_open_n_and_ranks_the_bump(self):
        text = plan_digest.digest(self.root)
        cand = text.split("## Candidates")[1].split("## Recent")[0]
        listed = [int(line.split("|")[1]) for line in cand.splitlines() if line.startswith("| ") and line[2].isdigit()]
        self.assertTrue(listed and set(listed) <= set(OPEN))
        self.assertEqual(listed[0], 7)  # the bump in the waste curve ranks first
        self.assertIn("N=6 2.600000000 -> 2.550000000", text)
        self.assertIn("only these may be queued): 5-9, 90-91", text)
        self.assertIn("| exploit | 1 | 1 |", text)

    def test_good_plan_is_applied_with_fresh_seeds(self):
        plan = {"summary": "s", "jobs": [{"N": "5-7", "schedule": "schedules/shrink-capped.json", "time": 300, "seeds": 2, "tag": "exploit"}],
                "hypotheses": [{"id": "H1", "status": "open", "claim": "c"}],
                "move_proposals": [{"name": "row-shift", "idea": "i", "mechanism": "m", "why": "w"}, {"name": "Bad Name!"}]}
        queue, variants, errs = self.check(plan)
        self.assertEqual(errs, [])
        self.assertEqual(queue["jobs"][0]["seeds"], [202611160, 202611161])
        plan_check.apply(plan, queue, variants, self.root, self.today)
        self.assertTrue((self.root / "queue.json").exists())
        self.assertIn("H1", (self.root / "planner" / "hypotheses.jsonl").read_text())
        props = list((self.root / "planner" / "proposals").glob("*.md"))
        self.assertEqual([p.name for p in props], ["2026-11-16-row-shift.md"])  # the badly named one is dropped
        self.assertIn("2026-11-16-row-shift", plan_digest.digest(self.root))

    def test_rejections(self):
        def errs(job, **extra):
            return self.check({"jobs": [job], **extra})[2]
        base = {"N": 5, "schedule": "schedules/ga.json", "time": 300, "seeds": 1, "tag": "exploit"}
        q, _, e = self.check({"jobs": [{**base, "N": [4, 5]}]})  # closed N are dropped with a warning
        self.assertEqual((e, q["jobs"][0]["N"]), ([], 5))
        self.assertTrue(any("dropped N [4]" in w for w in q["warnings"]))
        self.assertTrue(any("no open N left" in e for e in errs({**base, "N": [1, 2]})))
        self.assertTrue(any("does not exist" in e for e in errs({**base, "schedule": "schedules/nope.json"})))
        self.assertTrue(any("budget" in e for e in errs({**base, "N": "5-9", "time": 3600, "seeds": 8})))
        self.assertTrue(any("time" in e for e in errs({**base, "time": 5})))
        self.assertTrue(any("no tag" in e for e in errs({**base, "tag": ""})))
        self.assertTrue(any("20%" in e for e in self.check({"jobs": [
            {**base, "time": 3000, "seeds": 1}, {**base, "N": 6, "time": 600}]})[2]))
        self.assertTrue(any("experiments" in e for e in self.check({"jobs": [
            {**base, "tag": "H1", "N": "5-9"}, {**base, "N": 90}]})[2]))

    @unittest.skipUnless(ENGINE.exists(), "make engine first")
    def test_variants(self):
        v = {"name": "ga-big", "base": "schedules/ga.json", "changes": {"ga.offspring": 240}, "why": "w"}
        job = {"N": 90, "schedule": "schedules/planner/ga-big.json", "time": 300, "seeds": 1, "tag": "H2"}
        queue, variants, errs = self.check({"jobs": [job, {**job, "schedule": "schedules/ga.json", "tag": "exploit", "N": "5-9"}],
                                            "new_schedules": [v]}, ENGINE)
        self.assertEqual(errs, [])
        self.assertEqual(variants["schedules/planner/ga-big.json"]["ga"]["offspring"], 240)
        bad = self.check({"jobs": [{**job, "tag": "exploit"}], "new_schedules": [v]}, ENGINE)[2]
        self.assertTrue(any("experiment" in e for e in bad))
        bad = self.check({"jobs": [job], "new_schedules": [{**v, "changes": {"ga.nope": 1}}]}, ENGINE)[2]
        self.assertTrue(any("not a schedule key" in e for e in bad))

    def test_extract_and_fallback(self):
        wrapped = json.dumps({"result": "Here is the plan:\n```json\n" + json.dumps({"jobs": [], "summary": "x"}) + "\n```"})
        self.assertEqual(plan_check.extract_plan(wrapped)["summary"], "x")
        with self.assertRaises(ValueError):
            plan_check.extract_plan("no plan here")
        q = plan_check.fallback(self.root, self.today)
        ns = [j["N"] for j in q["jobs"]]
        self.assertTrue(ns and set(ns) <= set(OPEN))
        for j in q["jobs"]:
            self.assertEqual(j["schedule"], "schedules/ga-seeded.json" if j["N"] >= 82 else "schedules/shrink-capped.json")


if __name__ == "__main__":
    unittest.main()
