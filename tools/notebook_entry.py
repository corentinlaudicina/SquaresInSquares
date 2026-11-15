#!/usr/bin/env python3
"""Drafts one lab-notebook entry for the last 24 hours: planner/notebook/<date>.md (local, unpublished).

    python3 tools/notebook_entry.py            # run by tools/daily.sh after each plan
Published only by tools/publish_notebook.sh, by hand. An entry reports counts, which N improved and
their side lengths, the hypotheses and the new plan; never a packing.
"""
import datetime as dt
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plan_check  # noqa: E402
import plan_digest  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


def entry(root, now=None):
    now = now or dt.datetime.now(dt.timezone.utc)
    day = [r for r in plan_digest.load_ledger(root) if (t := plan_digest.when(r)) and now - t < dt.timedelta(days=1)]
    L = [f"## {now.astimezone():%Y-%m-%d}: daily entry", ""]
    by = Counter((Path(r.get("schedule", "")).stem, r.get("tag") or "-") for r in day)
    L.append(f"- {len(day)} runs finished, {sum(r.get('wall_s', 0) for r in day) / 3600:.1f} job-hours, "
             f"{sum(r.get('status') != 'ok' or not r.get('valid') for r in day)} failed or invalid.")
    if by:
        L.append("- By schedule and tag: " + ", ".join(f"{s} / {t}: {n}" for (s, t), n in sorted(by.items())) + ".")
    imps = [r for r in day if r.get("stored") == "best"]
    if imps:
        L.append(f"- **{len(imps)} improvements**:")
        L += [f"  - N = {r['N']}: {r['incumbent']:.9f} -> {r['side']:.9f} ({Path(r.get('schedule', '')).stem}, {r.get('tag')})"
              for r in imps]
    else:
        L.append("- No improvement.")
    hp = root / "planner" / "hypotheses.jsonl"
    if hp.exists():
        hyps = [json.loads(l) for l in hp.read_text().splitlines() if l.strip()]
        if hyps:
            L.append("- Hypotheses: " + "; ".join(f"**{h['id']}** ({h.get('status', '?')}): {h.get('claim', '')}"
                                                  for h in hyps))
    raw = root / "planner" / "plans" / f"{now.astimezone():%Y-%m-%d}.raw.json"
    try:
        plan = plan_check.extract_plan(raw.read_text())
        L.append(f"- Next plan: {plan.get('summary', '').strip()}")
    except (OSError, ValueError):
        L.append("- Next plan: the fallback plan (no planner answer today).")
    return "\n".join(L) + "\n"


def main():
    out = ROOT / "planner" / "notebook" / f"{dt.date.today():%Y-%m-%d}.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(entry(ROOT))
    print(f"notebook entry drafted: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
