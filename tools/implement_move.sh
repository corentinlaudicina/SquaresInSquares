#!/usr/bin/env bash
# Implements one planner move proposal on its own branch, then gates it. Started by hand:
#   tools/implement_move.sh planner/proposals/2026-11-16-row-shift.md
# 1. git worktree ../<repo>-move-<name> on a new branch move/<name>, from main
# 2. claude -p with edit and shell tools inside that worktree (planner_template/implement_move.md)
# 3. make check, then a recovery gate: schedules/planner/ga-<name>.json vs schedules/ga.json
# 4. the gate report is appended to the proposal; you review the branch and merge it (or not)
# Nothing here touches main, best/ or the running search.
set -eu
PROP="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
NAME="$(basename "$PROP" .md | sed -E 's/^[0-9]{4}-[0-9]{2}-[0-9]{2}-//')"
WT="$ROOT/../$(basename "$ROOT")-move-$NAME"
BUDGET="${CLAUDE_BUDGET_USD:-5}"
cd "$ROOT"
git worktree add -b "move/$NAME" "$WT" main
cd "$WT"
mkdir -p best runs && cp -R "$ROOT/best/." best/ && cp -R "$ROOT/runs/starts" runs/  # gate inputs (gitignored)

{ cat planner_template/implement_move.md; echo; cat "$PROP"; } |
    claude -p --permission-mode acceptEdits --allowedTools "Read" "Edit" "Write" "Bash" \
        --no-session-persistence --max-budget-usd "$BUDGET" > "$WT/implement.log" 2>&1 || true

STATUS="implemented on branch move/$NAME"
if ! make -j8 all >/dev/null 2>&1 || ! make check > "$WT/check.log" 2>&1; then
    STATUS="FAILED: build or tests (see $WT/check.log)"
elif [ ! -f "schedules/planner/ga-$NAME.json" ]; then
    STATUS="FAILED: no schedules/planner/ga-$NAME.json"
else
    python3 tools/benchmark.py --schedules schedules/ga.json "schedules/planner/ga-$NAME.json" --set train \
        --n 82 99 110 145 200 240 --seeds 2 --time 45 --workers 2 \
        --starts 'runs/starts/recovery/n{N:04d}.json' --name "move-$NAME" > "$WT/gate.log" 2>&1 || true
    STATUS="gated: $(grep -E '^- `ga-' "runs/bench/move-$NAME/report.md" 2>/dev/null || echo 'gate failed to run')"
fi
{ echo; echo "status: $STATUS ($(date +%Y-%m-%d))"; echo; echo "branch: move/$NAME, worktree: $WT"; } >> "$PROP"
python3 - "$PROP" "$STATUS" <<'PY'
import sys
p, status = sys.argv[1], sys.argv[2]
lines = open(p).read().splitlines()
i = next((k for k, l in enumerate(lines) if l.startswith("status: ")), None)
if i is not None:
    lines[i] = "status: " + status  # the first status line is what the digest shows
open(p, "w").write("\n".join(lines) + "\n")
PY
echo "$STATUS"
echo "review:  git -C $ROOT diff main...move/$NAME"
echo "merge:   git -C $ROOT merge move/$NAME     cleanup: git -C $ROOT worktree remove $WT"
