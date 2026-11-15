#!/usr/bin/env python3
"""Unattended runner: queue -> engine jobs in parallel -> validate -> best/ on improvement.

    python3 tools/runner.py [--queue queue.yaml] [--workers 4] [--root .]
    python3 tools/runner.py --import-legacy          # one-off: seeds + data/best -> best/

Queue (YAML, or JSON with a .json name):
    defaults: {schedule: schedules/shrink-baseline.json, time: 600, start: best}
    jobs:
      - {N: 29, seeds: [1, 2, 3]}
      - {N: "100-110", time: 300, schedule: schedules/compress-baseline.yaml}
      - {N: [200, 300], seeds: 4, kicks: 50}       # seeds: list, or a count (1..k)

Job fields: N (int, list, or "a-b"), seeds, schedule (.json or .yaml), time (wall seconds), kicks
(optional count budget, for reproducible runs), start (best | seed | <path>), tag (free text, copied
into the ledger: the planner tags exploit / explore / hypothesis jobs).

Each job has a deterministic id from its spec (N, seed, schedule content, budget, start) and runs in
runs/jobs/<id>/. Every finished job appends one line to runs/ledger.jsonl (what the planner reads, and
the done marker: a restart skips the ids it holds and re-runs the rest). Job folders are kept only
when they matter -- an improvement of best/ (full telemetry), a non-validated improvement, or a
failure (invalid, crash, timeout), each marked by a `keep` file -- plus the newest KEEP_RECENT others;
older unmarked folders are deleted. Every output is checked by tools/validate.py (no engine code involved). Only this
process writes best/ and non_validated_best/:
    best/N.json                 valid and strictly smaller side than what is there (atomic replace)
    non_validated_best/N.json   not valid, but smaller than both best/ and what is there
runs/summary.md is rewritten after every job; runs/heartbeat is touched every 30 s.
"""
import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import convert  # noqa: E402
import validate  # noqa: E402

IMPROVE = 1e-12  # relative: a side must beat the incumbent by more than this to replace it
KEEP_RECENT = 20  # unmarked job folders kept, newest first


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_structured(path):
    path = Path(path)
    if path.suffix in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError:
            sys.exit(f"runner: {path} is YAML but PyYAML is missing (pip install -r tools/requirements.txt)")
        return yaml.safe_load(path.read_text())
    return json.loads(path.read_text())


def expand_n(v):
    if isinstance(v, int):
        return [v]
    if isinstance(v, list):
        return [int(x) for x in v]
    lo, hi = str(v).split("-")
    return list(range(int(lo), int(hi) + 1))


def expand(queue, root):
    """Queue -> list of job dicts with ids. Raises ValueError on a bad queue."""
    defaults = {"time": 600, "start": "best", "seeds": [1], **(queue.get("defaults") or {})}
    jobs = []
    for spec in queue.get("jobs") or []:
        spec = {**defaults, **spec}
        unknown = set(spec) - {"N", "seeds", "schedule", "time", "kicks", "start", "tag"}
        if unknown:
            raise ValueError(f"unknown job keys {sorted(unknown)}")
        if "schedule" not in spec:
            raise ValueError("job without a schedule")
        sched_path = root / spec["schedule"]
        sched = load_structured(sched_path)  # YAML is converted here; the engine only reads JSON
        sched_text = json.dumps(sched, sort_keys=True)
        seeds = spec["seeds"] if isinstance(spec["seeds"], list) else list(range(1, int(spec["seeds"]) + 1))
        for n in expand_n(spec["N"]):
            for seed in seeds:
                key = json.dumps({"N": n, "seed": seed, "schedule": sched_text, "time": spec["time"],
                                  "kicks": spec.get("kicks"), "start": spec["start"], "tag": spec.get("tag")},
                                 sort_keys=True)
                jid = f"n{n:04d}_s{seed}_{sched.get('name', 'unnamed')}_{hashlib.sha256(key.encode()).hexdigest()[:10]}"
                jobs.append({"id": jid, "N": n, "seed": seed, "schedule": sched, "schedule_file": str(spec["schedule"]),
                             "time": float(spec["time"]), "kicks": spec.get("kicks"), "start": spec["start"],
                             "tag": spec.get("tag")})
    return jobs


def digest(tel_path):
    """The few telemetry fields the ledger keeps (the full file survives only for kept jobs)."""
    try:
        t = json.loads(Path(tel_path).read_text())
    except (OSError, ValueError):
        return {}
    side, budget, run = t.get("side", {}), t.get("budget", {}), t.get("run", {})
    out = {"schedule_name": run.get("schedule", {}).get("name"), "schedule_hash": run.get("schedule", {}).get("hash"),
           "start_side": side.get("start"), "best_raw": side.get("best"), "t_best_s": side.get("t_best_s"),
           "kicks": budget.get("kicks"), "stop": budget.get("stop")}
    if "ga" in t:
        g = t["ga"]
        out["ga"] = {k: g.get(k) for k in ("generations", "ratchets", "legal_improvements", "resets") if k in g}
        out["moves"] = {m: [v["produced"], v["survived"]] for m, v in g.get("moves", {}).items()}
    elif "kicks" in t:
        out["moves"] = {m: [v["tries"], v["impr"]] for m, v in t["kicks"].items()}
    return out


class Store:
    """best/ and non_validated_best/: compare-and-replace under a lock (thread and process safe)."""

    def __init__(self, root):
        self.best = root / "best"
        self.nonval = root / "non_validated_best"
        self.best.mkdir(exist_ok=True)
        self.nonval.mkdir(exist_ok=True)
        self.lock_path = root / "best" / ".lock"
        self.tlock = threading.Lock()

    @staticmethod
    def side_of(path):
        try:
            return json.loads(path.read_text())["side"]
        except (OSError, ValueError, KeyError):
            return float("inf")

    @staticmethod
    def atomic_write(path, doc):
        tmp = path.with_suffix(".json.tmp")
        with open(tmp, "w") as f:
            json.dump(doc, f, separators=(",", ":"))
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)

    def offer(self, n, doc, valid):
        """Returns 'best', 'non_validated' or None, with the incumbent side it beat."""
        with self.tlock, open(self.lock_path, "w") as lk:
            fcntl.flock(lk, fcntl.LOCK_EX)
            side = doc["side"]
            best_side = self.side_of(self.best / f"{n}.json")
            if valid:
                if side < best_side * (1 - IMPROVE):
                    self.atomic_write(self.best / f"{n}.json", doc)
                    return "best", best_side
                return None, best_side
            nv_side = self.side_of(self.nonval / f"{n}.json")
            if side < min(best_side, nv_side) * (1 - IMPROVE):
                self.atomic_write(self.nonval / f"{n}.json", doc)
                return "non_validated", min(best_side, nv_side)
            return None, best_side


class Runner:
    def __init__(self, root, engine, workers, until=None):
        self.root, self.engine, self.workers = root, engine, workers
        self.until = until  # epoch seconds: no job starts after it (it stays unfinished for the next queue)
        self.skipped = 0
        self.runs = root / "runs"
        (self.runs / "jobs").mkdir(parents=True, exist_ok=True)
        self.store = Store(root)
        self.lock = threading.Lock()
        self.running = 0
        self.pending = 0
        self.results = []  # result dicts finished in this session

    def job_dir(self, job):
        return self.runs / "jobs" / job["id"]

    def start_file(self, job):
        n, start = job["N"], job["start"]
        if start == "best":
            p = self.root / "best" / f"{n}.json"
            if p.exists():
                return p
            start = "seed"
        if start == "seed":
            return self.root / "data" / "init_external" / f"n{n:04d}.cfg"
        return self.root / start

    def run_job(self, job):
        d = self.job_dir(job)
        d.mkdir(parents=True, exist_ok=True)
        (d / "schedule.json").write_text(json.dumps(job["schedule"], indent=1))
        start = self.start_file(job)
        out = d / "out.json"
        for stale in (out, d / "out.telemetry.json"):
            stale.unlink(missing_ok=True)
        cmd = [str(self.engine), "--N", str(job["N"]), "--seed", str(job["seed"]), "--time", str(job["time"]),
               "--schedule", str(d / "schedule.json"), "--start", str(start), "--out", str(out)]
        if job["kicks"]:
            cmd += ["--kicks", str(job["kicks"])]
        # The engine checks the clock only between kicks, and legalization comes after: allow slack.
        timeout = job["time"] * 1.5 + 120
        res = {"id": job["id"], "tag": job.get("tag"), "N": job["N"], "seed": job["seed"], "schedule": job["schedule_file"],
               "start": str(start), "started": now(), "cmd": cmd}
        with self.lock:
            self.running += 1
        t0 = time.time()
        try:
            p = subprocess.run(cmd, cwd=self.root, capture_output=True, text=True, timeout=timeout)
            (d / "stdout.txt").write_text(p.stdout + p.stderr)
            res["status"] = "ok" if p.returncode == 0 else f"engine_exit_{p.returncode}"
            res["summary"] = (p.stdout.strip().splitlines() or [""])[-1]
        except subprocess.TimeoutExpired:
            res["status"] = "timeout"
            res["summary"] = f"killed after {timeout:.0f}s"
        res["wall_s"] = round(time.time() - t0, 1)

        if out.exists():
            try:
                doc = validate.load(out)
                ok, reasons, diag = validate.validate(doc, job["N"])
            except (OSError, ValueError) as e:
                doc, ok, reasons, diag = None, False, [f"unreadable: {e}"], {}
            res.update(valid=ok, reasons=reasons[:3], side=doc.get("side") if isinstance(doc, dict) else None)
            if isinstance(doc, dict) and isinstance(doc.get("side"), (int, float)):
                doc.setdefault("meta", {}).update(job=job["id"], validated=ok, validated_at=now())
                where, incumbent = self.store.offer(job["N"], doc, ok)
                res.update(stored=where, incumbent=incumbent)
        else:
            res.update(valid=False, side=None, stored=None)
        res["finished"] = now()
        res.update(digest(d / "out.telemetry.json"))
        keep = ("improvement" if res.get("stored") == "best" else
                "non_validated" if res.get("stored") == "non_validated" else
                "failure" if res["status"] != "ok" or not res.get("valid") else None)
        res["kept"] = keep or "recent"
        if keep:
            (d / "keep").write_text(keep + "\n")
        Store.atomic_write(d / "result.json", res)
        self.append_ledger(res)  # the done marker, written last
        with self.lock:
            self.running -= 1
            self.results.append(res)
            self.prune()
            self.write_summary()
        return res

    def append_ledger(self, res):
        line = json.dumps({k: v for k, v in res.items() if k != "cmd"}, separators=(",", ":")) + "\n"
        with self.lock, open(self.runs / "ledger.jsonl", "a") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            f.write(line)
            f.flush()
            os.fsync(f.fileno())

    def done_ids(self):
        ids = set()
        try:
            for line in (self.runs / "ledger.jsonl").read_text().splitlines():
                try:
                    ids.add(json.loads(line)["id"])
                except (ValueError, KeyError):
                    pass  # a line cut short by a crash: that job simply runs again
        except FileNotFoundError:
            pass
        return ids

    def prune(self):
        """Delete unmarked job folders beyond the newest KEEP_RECENT finished ones."""
        recent = [d for d in (self.runs / "jobs").iterdir()
                  if d.is_dir() and (d / "result.json").exists() and not (d / "keep").exists()]
        recent.sort(key=lambda d: (d / "result.json").stat().st_mtime, reverse=True)
        for d in recent[KEEP_RECENT:]:
            shutil.rmtree(d, ignore_errors=True)

    def write_summary(self):
        done = self.results
        lines = [f"# Runner summary  ({now()})", "",
                 f"session: {len(done)} finished, {self.running} running, {self.pending} queued | "
                 f"valid {sum(1 for r in done if r.get('valid'))} | stored best "
                 f"{sum(1 for r in done if r.get('stored') == 'best')} | failed "
                 f"{sum(1 for r in done if r['status'] != 'ok')}", ""]
        improved = [r for r in done if r.get("stored") == "best"]
        if improved:
            lines += ["## Improvements", ""]
            lines += [f"- N={r['N']}: {r['incumbent']:.9f} -> {r['side']:.9f} ({r['id']})" for r in improved[-20:]]
            lines.append("")
        lines += ["## Last jobs", "", "```"]
        for r in done[-15:]:
            tag = r.get("stored") or ("valid" if r.get("valid") else "INVALID")
            lines.append(f"[{r['status']}/{tag}] {r.get('summary', '')}"[:220])
        lines += ["```", ""]
        tmp = self.runs / "summary.md.tmp"
        tmp.write_text("\n".join(lines))
        os.replace(tmp, self.runs / "summary.md")

    def heartbeat(self, stop):
        while not stop.wait(30):
            self.touch()

    def touch(self):
        (self.runs / "heartbeat").write_text(f"{now()} running={self.running} pending={self.pending}\n")

    def run(self, jobs):
        done = self.done_ids()
        todo = [j for j in jobs if j["id"] not in done]
        print(f"runner: {len(jobs)} jobs, {len(jobs) - len(todo)} already done, {len(todo)} to run, "
              f"{self.workers} workers")
        self.pending = len(todo)
        stop = threading.Event()
        threading.Thread(target=self.heartbeat, args=(stop,), daemon=True).start()
        self.touch()

        def one(j):
            with self.lock:
                self.pending -= 1
                if self.until is not None and time.time() >= self.until:
                    self.skipped += 1
                    return {"id": j["id"], "status": "skipped_until"}
            r = self.run_job(j)
            print(f"  {r['status']:8s} {'VALID  ' if r.get('valid') else 'INVALID'} {r.get('stored') or '-':13s} "
                  f"{r.get('summary', '')[:150]}", flush=True)
            return r

        with ThreadPoolExecutor(self.workers) as pool:
            results = list(pool.map(one, todo))
        stop.set()
        self.touch()
        self.write_summary()
        if self.skipped:
            print(f"runner: --until reached, {self.skipped} jobs not started")
        return [r for r in results if r["status"] != "skipped_until"]


def import_legacy(root, engine):
    """Seeds and legacy records -> legalize (engine) -> validate -> best/ (or non_validated_best/)."""
    store = Store(root)
    cands = {}
    for p in sorted((root / "data" / "init_external").glob("*.cfg")) + sorted((root / "data" / "best").glob("*.log")):
        try:
            d = convert.parse_record(p.read_text(), str(p))
        except ValueError as e:
            print(f"  skip {p}: {e}")
            continue
        cands.setdefault(d["n"], []).append((d["side"], p))
    tmpdir = root / "runs" / "import"
    tmpdir.mkdir(parents=True, exist_ok=True)

    def one(n):
        best = None
        for side, p in sorted(cands[n])[:3]:  # the three smallest claims; legalization can reorder them
            out = tmpdir / f"{n}_{p.stem}.json"
            subprocess.run([str(engine), "--legalize", str(p), "--out", str(out)], cwd=root,
                           capture_output=True, text=True)
            if not out.exists():
                continue
            doc = validate.load(out)
            ok = validate.validate(doc, n)[0]
            if best is None or (ok, -doc["side"]) > (best[0], -best[1]["side"]):
                best = (ok, doc)
        if best:
            best[1].setdefault("meta", {}).update(validated=best[0], validated_at=now(), job="import-legacy")
            return n, best[0], store.offer(n, best[1], best[0])[0]
        return n, False, None

    with ThreadPoolExecutor(os.cpu_count() or 4) as pool:
        res = list(pool.map(one, sorted(cands)))
    print(f"import: {len(res)} N, valid {sum(r[1] for r in res)}, stored best "
          f"{sum(r[2] == 'best' for r in res)}, non-validated {sum(r[2] == 'non_validated' for r in res)}")
    return res


def parse_until(v):
    """HH:MM -> the next such local time as epoch seconds; a number is taken as epoch seconds."""
    if v is None:
        return None
    if ":" not in v:
        return float(v)
    hh, mm = (int(x) for x in v.split(":"))
    t = dt.datetime.now().replace(hour=hh, minute=mm, second=0, microsecond=0)
    if t <= dt.datetime.now():
        t += dt.timedelta(days=1)
    return t.timestamp()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--queue", default="queue.yaml")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument("--engine", default=None, help="default: <root>/build/engine")
    ap.add_argument("--import-legacy", action="store_true")
    ap.add_argument("--until", help="local time HH:MM (next occurrence) or epoch seconds: start no job after it")
    a = ap.parse_args(argv)
    root = Path(a.root).resolve()
    engine = Path(a.engine).resolve() if a.engine else root / "build" / "engine"
    if not engine.exists():
        sys.exit(f"runner: no engine at {engine} (make engine)")
    if a.import_legacy:
        import_legacy(root, engine)
        return 0
    try:
        jobs = expand(load_structured(root / a.queue), root)
    except (OSError, ValueError) as e:
        sys.exit(f"runner: bad queue: {e}")
    results = Runner(root, engine, a.workers, parse_until(a.until)).run(jobs)
    return 0 if all(r["status"] == "ok" for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
