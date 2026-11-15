#!/usr/bin/env python3
"""One screen on the running search.   python3 tools/status.py   (or: watch -n 60 python3 tools/status.py)

Runner alive and heartbeat age, jobs running (N, schedule, elapsed / budget), done and pending today,
improvements, the last finished jobs, and when the planner runs next.
"""
import datetime as dt
import json
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def alive(pidfile):
    try:
        os.kill(int(pidfile.read_text()), 0)
        return True
    except (OSError, ValueError):
        return False


def main():
    now = dt.datetime.now(dt.timezone.utc)
    pid = ROOT / "planner" / "runner.pid"
    hb = ROOT / "runs" / "heartbeat"
    print(f"runner: {'RUNNING' if alive(pid) else 'not running'}", end="")
    if hb.exists():
        age = now.timestamp() - hb.stat().st_mtime
        print(f"   heartbeat {int(age)} s ago: {hb.read_text().strip()}", end="")
    print()

    ps = subprocess.run(["ps", "-ax", "-o", "etime=,command="], capture_output=True, text=True).stdout
    print("\nrunning jobs:")
    for line in ps.splitlines():
        if "build/engine" in line and "--N" in line:
            m = re.search(r"--N (\d+).*--time ([\d.]+).*runs/jobs/([^/]+)/", line)
            if m:
                print(f"  N={m.group(1):>4}  elapsed {line.split()[0]:>8} of {float(m.group(2)):.0f} s   {m.group(3)}")

    ledger = []
    lp = ROOT / "runs" / "ledger.jsonl"
    if lp.exists():
        for line in lp.read_text().splitlines():
            try:
                ledger.append(json.loads(line))
            except ValueError:
                pass
    today = now.strftime("%Y-%m-%d")
    done = [r for r in ledger if r.get("finished", "").startswith(today)]
    imp = [r for r in ledger if r.get("stored") == "best"]
    q = ROOT / "queue.json"
    total = sum(len(j["N"]) if isinstance(j["N"], list) else 1 for j in json.loads(q.read_text())["jobs"]) if q.exists() else 0
    print(f"\nfinished today: {len(done)}   all time: {len(ledger)}   improvements all time: {len(imp)}"
          f"   failures today: {sum(r.get('status') != 'ok' or not r.get('valid') for r in done)}")
    for r in imp[-5:]:
        print(f"  IMPROVED N={r['N']}: {r.get('incumbent', 0):.9f} -> {r['side']:.9f}  ({r.get('schedule', '')}, {r.get('finished', '')})")
    print("\nlast finished:")
    for r in ledger[-6:]:
        print(f"  {r.get('finished', '')[11:16]}  N={r['N']:>4}  {Path(r.get('schedule', '')).stem:14s} {r.get('tag') or '':8s} "
              f"{r.get('status', '?'):8s} side {r.get('side') or 0:.9f}  {r.get('stored') or ''}")
    if not ledger:
        print("  (none yet)")

    plan = ROOT / "planner" / "plans"
    last = sorted(plan.glob("*.raw.json"))[-1:] if plan.exists() else []
    nxt = now.astimezone().replace(hour=7, minute=0, second=0, microsecond=0)
    if nxt <= now.astimezone():
        nxt += dt.timedelta(days=1)
    print(f"\nplanner: last plan {last[0].name[:10] if last else '-'}, next run {nxt:%a %H:%M}; log planner/daily.log")


if __name__ == "__main__":
    main()
