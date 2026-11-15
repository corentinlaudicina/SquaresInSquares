#!/usr/bin/env python3
"""Validates the planner's JSON plan and applies it; or writes the deterministic fallback plan.

    python3 tools/plan_check.py planner/plans/2026-11-16.raw.json    # claude -p --output-format json output
    python3 tools/plan_check.py --fallback                           # no usable plan: default queue

The planner has no tools: it answers with one JSON object (schema in planner_template/prompt.md), and
this script is the only thing that writes. Rules, all enforced here:
  * only N whose status is "open" in data/known_optima_1_1000.csv;
  * every schedule exists (schedules/ or a variant defined in the same plan) and the engine accepts it;
  * variants only change existing keys of an existing schedule, live in schedules/planner/, and may only
    run in experiment jobs (tag H<id>);
  * budget: sum of job times / workers <= budget_hours (+5%); no N over 20% of it; experiments <= 25%;
  * 30 <= time <= 3600 per job; seeds are a count, turned into fresh seed values here (job ids are
    deterministic, so reusing a seed would make the runner skip the job as already done).
On success it writes queue.json, any variants, planner/hypotheses.jsonl and planner/playbook.md, and keeps
the plan in planner/plans/. On failure it writes nothing and exits 1 with the reasons.
"""
import argparse
import copy
import datetime as dt
import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plan_digest  # noqa: E402
import runner  # noqa: E402

GA_MIN_N = 82  # the recovery gates: shrink-capped below, the seeded GA from here (docs/AUDIT.md section 15)
DEFAULT_GA = "schedules/ga-explore.json"  # islands + seeding + adaptive ratchet + large moves: production default


def extract_plan(text):
    """The JSON object in claude's answer (the --output-format json wrapper, or a bare object)."""
    try:
        wrapper = json.loads(text)
        if isinstance(wrapper, dict) and "result" in wrapper:
            text = wrapper["result"]
        elif isinstance(wrapper, dict) and "jobs" in wrapper:
            return wrapper
    except ValueError:
        pass
    start = text.find("{")
    while start >= 0:
        try:
            obj, _ = json.JSONDecoder().raw_decode(text[start:])
            if isinstance(obj, dict) and "jobs" in obj:
                return obj
        except ValueError:
            pass
        start = text.find("{", start + 1)
    raise ValueError("no JSON plan with a 'jobs' list in the answer")


def set_dotted(d, path, value):
    keys = path.split(".")
    for k in keys[:-1]:
        d = d.setdefault(k, {})
    d[keys[-1]] = value


def check_and_build(plan, root, engine, today):
    """Returns (queue, variants {path: doc}, errors); warnings go into queue["warnings"]."""
    errs, warnings = [], []
    status = plan_digest.load_status(root)
    m = plan_digest.machine(root)
    keys = set()
    if engine and Path(engine).exists():
        out = subprocess.run([str(engine), "--list-schedule-keys"], capture_output=True, text=True).stdout
        keys = {line.split(" ")[0] for line in out.splitlines()}

    variants = {}
    for v in plan.get("new_schedules", []) or []:
        name, base = str(v.get("name", "")), str(v.get("base", ""))
        if not re.fullmatch(r"[a-z0-9][a-z0-9-]{1,40}", name):
            errs.append(f"variant name {name!r}: lowercase letters, digits and dashes")
            continue
        bp = root / base
        if not bp.exists() or bp.parent.name not in ("schedules", "planner"):
            errs.append(f"variant {name}: base {base} is not an existing schedule")
            continue
        doc = copy.deepcopy(json.loads(bp.read_text()))
        for k, val in (v.get("changes") or {}).items():
            if keys and k not in keys:
                errs.append(f"variant {name}: {k} is not a schedule key")
            set_dotted(doc, k, val)
        doc["name"] = name
        doc["description"] = f"planner variant of {base}: {v.get('why', '')}"[:300]
        variants[f"schedules/planner/{name}.json"] = doc

    jobs, total, per_n, experiments = [], 0.0, {}, 0.0
    stamp = int(today.strftime("%Y%m%d"))
    for i, j in enumerate(plan.get("jobs") or []):
        tag, sched = str(j.get("tag", "")).strip(), str(j.get("schedule", ""))
        try:
            ns = runner.expand_n(j["N"])
        except (KeyError, ValueError):
            errs.append(f"job {i}: bad N {j.get('N')!r}")
            continue
        allowed = set(plan_digest.open_ns(root))  # open and inside the machine's focus range
        closed = [n for n in ns if n not in allowed]
        if closed:  # never run them; dropping N only shrinks the plan, so the rest of it stands
            warnings.append(f"job {i}: dropped N {closed[:8]} (not open, or outside the focus range)")
            ns = [n for n in ns if n in allowed]
            if not ns:
                errs.append(f"job {i}: no open N left")
                continue
        if not tag:
            errs.append(f"job {i}: no tag (exploit, explore or H<id>)")
        if sched in variants:
            if not re.fullmatch(r"H\d+", tag):
                errs.append(f"job {i}: variant {sched} may only run in an experiment (tag H<id>)")
        elif not (root / sched).exists():
            errs.append(f"job {i}: schedule {sched} does not exist")
        if m.get("require_islands") and not re.fullmatch(r"H\d+", tag):
            doc = variants.get(sched)
            if doc is None and (root / sched).exists():
                doc = json.loads((root / sched).read_text())
            if doc is not None and not (doc.get("mode") == "ga" and int(doc.get("ga", {}).get("islands", 1)) >= 2):
                errs.append(f"job {i}: {tag} jobs must use an island GA (mode ga, ga.islands >= 2); {sched} is not one")
        t = j.get("time")
        if not isinstance(t, (int, float)) or not 30 <= t <= 3600:
            errs.append(f"job {i}: time {t!r} not in [30, 3600] s")
            continue
        seeds = j.get("seeds", 1)
        if not isinstance(seeds, int) or not 1 <= seeds <= 8:
            errs.append(f"job {i}: seeds must be a count in [1, 8]")
            continue
        cost = float(t) * seeds
        total += cost * len(ns)
        for n in ns:
            per_n[n] = per_n.get(n, 0.0) + cost
        if re.fullmatch(r"H\d+", tag):
            experiments += cost * len(ns)
        jobs.append({"N": ns if len(ns) > 1 else ns[0], "schedule": sched, "time": t, "tag": tag,
                     "seeds": [stamp * 10 + k for k in range(seeds)]})
    if not jobs:
        errs.append("no jobs")
    budget = m["budget_hours"] * 3600 * m["workers"]
    if total > budget * 1.05:
        errs.append(f"budget: {total / 3600:.1f} job-hours asked, {budget / 3600:.1f} available")
    if per_n and total and max(per_n.values()) > 0.2 * max(total, budget * 0.5):
        n, c = max(per_n.items(), key=lambda kv: kv[1])
        errs.append(f"N={n} gets {c / 3600:.1f} job-hours, over 20% of the plan")
    if experiments > 0.25 * max(total, 1.0):
        errs.append(f"experiments take {experiments / total:.0%} of the plan, over 25%")
    return {"defaults": {"start": "best"}, "jobs": jobs, "warnings": warnings}, variants, errs


def fallback(root, today):
    """The playbook default: top-ranked open N in focus, the island GA (shrink-polish below GA_MIN_N unless
    islands are required), 1800 s, one seed."""
    m = plan_digest.machine(root)
    slots = int(m["budget_hours"] * 3600 * m["workers"] // 1800)
    rows = plan_digest.rank(root, plan_digest.load_ledger(root), dt.datetime.now(dt.timezone.utc))
    stamp = int(today.strftime("%Y%m%d"))
    jobs = [{"N": r["N"], "time": 1800, "tag": "fallback", "seeds": [stamp * 10],
             "schedule": DEFAULT_GA if r["N"] >= GA_MIN_N or m.get("require_islands") else "schedules/shrink-polish.json"}
            for r in rows[:max(1, slots)]]
    return {"defaults": {"start": "best"}, "jobs": jobs}


def write_json(path, doc):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(doc, indent=1) + "\n")
    tmp.replace(path)


def apply(plan, queue, variants, root, today):
    pdir = root / "planner"
    for rel, doc in variants.items():
        write_json(root / rel, doc)
    write_json(root / "queue.json", {k: v for k, v in queue.items() if k != "warnings"})
    hyp_path = pdir / "hypotheses.jsonl"
    hyps = {}
    if hyp_path.exists():
        for line in hyp_path.read_text().splitlines():
            try:
                h = json.loads(line)
                hyps[h["id"]] = h
            except (ValueError, KeyError):
                pass
    for h in plan.get("hypotheses") or []:
        if isinstance(h, dict) and re.fullmatch(r"H\d+", str(h.get("id", ""))):
            hyps[h["id"]] = {**hyps.get(h["id"], {}), **h, "updated": today.isoformat()}
    hyp_path.write_text("".join(json.dumps(h) + "\n" for h in hyps.values()))
    for mp in plan.get("move_proposals") or []:
        name = str(mp.get("name", "")) if isinstance(mp, dict) else ""
        if re.fullmatch(r"[a-z0-9][a-z0-9-]{1,40}", name):
            prop = pdir / "proposals" / f"{today.isoformat()}-{name}.md"
            prop.parent.mkdir(parents=True, exist_ok=True)
            prop.write_text(f"# Move proposal: {name}\n\nstatus: proposed ({today.isoformat()})\n\n"
                            + "".join(f"**{k}**: {mp.get(k, '')}\n\n" for k in ("idea", "mechanism", "why")))
    adds = [a for a in plan.get("playbook_additions") or [] if isinstance(a, str) and re.search(r"\bH\d+\b", a)]
    if adds:
        with open(pdir / "playbook.md", "a") as f:
            f.write("".join(f"- {today.isoformat()}: {a}\n" for a in adds))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("plan", nargs="?")
    ap.add_argument("--fallback", action="store_true")
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent.parent))
    a = ap.parse_args(argv)
    root = Path(a.root).resolve()
    today = dt.date.today()
    if a.fallback:
        write_json(root / "queue.json", fallback(root, today))
        print(f"plan_check: fallback queue written ({len(fallback(root, today)['jobs'])} jobs)")
        return 0
    try:
        plan = extract_plan(Path(a.plan).read_text())
    except (OSError, ValueError) as e:
        print(f"plan_check: {e}")
        return 1
    queue, variants, errs = check_and_build(plan, root, root / "build" / "engine", today)
    if not errs:  # the engine itself must accept every variant
        for rel, doc in variants.items():
            tmp = root / "planner" / "check.tmp.json"
            write_json(tmp, doc)
            p = subprocess.run([str(root / "build" / "engine"), "--check-schedule", str(tmp)], capture_output=True, text=True)
            tmp.unlink(missing_ok=True)
            if p.returncode != 0:
                errs.append(f"variant {rel}: engine refuses it: {p.stderr.strip()}")
    if errs:
        print("plan_check: plan rejected:\n  " + "\n  ".join(errs))
        return 1
    for w in queue["warnings"]:
        print(f"plan_check: warning: {w}")
    apply(plan, queue, variants, root, today)
    print(f"plan_check: plan accepted: {len(queue['jobs'])} job specs, {len(variants)} variants. {plan.get('summary', '')[:200]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
