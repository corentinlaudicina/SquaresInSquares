#!/usr/bin/env python3
"""Runner tests, in a temporary root with a fake engine (no C++ needed) plus one real-engine job.

    python3 -m unittest tools/test_runner.py -v
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
import runner  # noqa: E402

REPO = Path(__file__).resolve().parent.parent

# Fake engine: behaviour picked by --seed. 1 = valid 2x2 grid at side 2+1e-9*seed... see below.
FAKE = textwrap.dedent('''\
    #!/usr/bin/env python3
    import json, sys, time
    a = dict(zip(sys.argv[1::2], sys.argv[2::2]))
    seed, out = int(a["--seed"]), a["--out"]
    side = 2.0 + 0.001 * seed
    sq = [[0.5, 0.5, 0], [0.5, 1.5, 0], [1.5, 0.5, 0], [1.5, 1.5, 0]]
    if seed == 7:              # overlapping, and claims a smaller side
        sq[1][1] = 1.4; side = 1.999
    if seed == 8:              # hangs
        time.sleep(600)
    if seed == 9:              # crashes without output
        sys.exit(3)
    json.dump({"format": "sqpack/1", "n": 4, "side": side, "squares": sq, "meta": {}}, open(out, "w"))
    print(f"fake n=4 seed={seed} side={side}")
''')


class RunnerTest(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        (self.root / "schedules").mkdir()
        (self.root / "schedules" / "s.json").write_text('{"name": "s"}')
        (self.root / "data" / "init_external").mkdir(parents=True)
        self.engine = self.root / "fake_engine.py"
        self.engine.write_text(FAKE)
        self.engine.chmod(0o755)

    def tearDown(self):
        shutil.rmtree(self.root)

    def queue(self, jobs, name="queue.json"):
        (self.root / name).write_text(json.dumps({"defaults": {"schedule": "schedules/s.json", "time": 1}, "jobs": jobs}))
        return name

    def run_queue(self, name, workers=2):
        jobs = runner.expand(runner.load_structured(self.root / name), self.root)
        return runner.Runner(self.root, self.engine, workers).run(jobs)

    def best_side(self):
        return json.loads((self.root / "best" / "4.json").read_text())["side"]

    def test_improvement_order_and_restart(self):
        q = self.queue([{"N": 4, "seeds": [3, 2]}])
        rs = self.run_queue(q, workers=1)
        self.assertEqual([r["stored"] for r in rs], ["best", "best"])  # 2.003, then 2.002 beats it
        self.assertAlmostEqual(self.best_side(), 2.002)
        self.assertTrue(json.loads((self.root / "best" / "4.json").read_text())["meta"]["validated"])
        # a restart runs nothing; deleting one done-marker re-runs exactly that job
        self.assertEqual(self.run_queue(q), [])
        ledger = self.root / "runs" / "ledger.jsonl"
        lines = ledger.read_text().splitlines()
        self.assertEqual(len(lines), 2)
        ledger.write_text("\n".join(lines[1:]) + "\n{\"id\": \"cut sh")  # one job forgotten, a torn last line
        rs = self.run_queue(q)
        self.assertEqual(len(rs), 1)
        self.assertIsNone(rs[0]["stored"])  # no improvement: best/ untouched
        self.assertAlmostEqual(self.best_side(), 2.002)
        self.assertTrue((self.root / "runs" / "summary.md").exists())
        self.assertTrue((self.root / "runs" / "heartbeat").exists())

    def test_retention(self):
        # 25 jobs: 2.001 improves best (kept, full telemetry), seed 7 is invalid (kept), the other 23 change
        # nothing; of those only the newest KEEP_RECENT folders survive. Every job is in the ledger.
        rs = self.run_queue(self.queue([{"N": 4, "seeds": [1, 7] + list(range(10, 33))}]), workers=1)
        self.assertEqual(len(rs), 25)
        jobs = self.root / "runs" / "jobs"
        kept = {d.name: (d / "keep").read_text().strip() for d in jobs.iterdir() if (d / "keep").exists()}
        self.assertEqual(sorted(kept.values()), ["improvement", "non_validated"])
        unmarked = [d for d in jobs.iterdir() if not (d / "keep").exists()]
        self.assertEqual(len(unmarked), runner.KEEP_RECENT)
        imp = next(d for d in jobs.iterdir() if kept.get(d.name) == "improvement")
        self.assertTrue((imp / "out.json").exists() and (imp / "stdout.txt").exists())
        ledger = [json.loads(l) for l in (self.root / "runs" / "ledger.jsonl").read_text().splitlines()]
        self.assertEqual(len(ledger), 25)
        self.assertEqual(sum(r["kept"] == "recent" for r in ledger), 23)

    def test_invalid_never_reaches_best(self):
        rs = self.run_queue(self.queue([{"N": 4, "seeds": [7]}]))
        self.assertFalse(rs[0]["valid"])
        self.assertEqual(rs[0]["stored"], "non_validated")
        self.assertFalse((self.root / "best" / "4.json").exists())
        self.assertTrue((self.root / "non_validated_best" / "4.json").exists())

    def test_wrong_n_is_invalid(self):
        rs = self.run_queue(self.queue([{"N": 5, "seeds": [1]}]))
        self.assertFalse(rs[0]["valid"])
        self.assertFalse((self.root / "best" / "5.json").exists())

    def test_timeout_and_crash(self):
        orig = runner.Runner.run_job

        def short(self_, job):  # shrink the timeout for the test: time*1.5+120 is too long to wait
            job = {**job, "time": -79.0}  # -> timeout 1.5 s
            return orig(self_, job)
        runner.Runner.run_job = short
        try:
            rs = self.run_queue(self.queue([{"N": 4, "seeds": [8, 9]}]))
        finally:
            runner.Runner.run_job = orig
        status = sorted(r["status"] for r in rs)
        self.assertEqual(status, ["engine_exit_3", "timeout"])
        self.assertFalse((self.root / "best" / "4.json").exists())

    def test_queue_expansion(self):
        q = {"defaults": {"schedule": "schedules/s.json"}, "jobs": [{"N": "10-12", "seeds": 2}, {"N": [5, 6]}]}
        jobs = runner.expand(q, self.root)
        self.assertEqual(len(jobs), 3 * 2 + 2)
        self.assertEqual(len({j["id"] for j in jobs}), len(jobs))
        self.assertEqual(jobs[0]["id"], runner.expand(q, self.root)[0]["id"])  # ids are deterministic
        with self.assertRaises(ValueError):
            runner.expand({"jobs": [{"N": 4, "schedule": "schedules/s.json", "bogus": 1}]}, self.root)

    @unittest.skipUnless((REPO / "build" / "engine").exists(), "make engine first")
    def test_real_engine_job(self):
        shutil.copy(REPO / "data" / "init_external" / "n0005.cfg", self.root / "data" / "init_external")
        shutil.copy(REPO / "schedules" / "shrink-baseline.json", self.root / "schedules")
        self.engine = REPO / "build" / "engine"
        q = self.queue([{"N": 5, "seeds": [1], "kicks": 3, "time": 60, "schedule": "schedules/shrink-baseline.json"}])
        rs = self.run_queue(q)
        self.assertEqual(rs[0]["status"], "ok", rs[0])
        self.assertTrue(rs[0]["valid"], rs[0])
        self.assertEqual(rs[0]["stored"], "best")
        out = subprocess.run([sys.executable, str(REPO / "tools" / "validate.py"), str(self.root / "best" / "5.json"),
                              "--n", "5"], capture_output=True, text=True)
        self.assertEqual(out.returncode, 0, out.stdout)


if __name__ == "__main__":
    unittest.main()
