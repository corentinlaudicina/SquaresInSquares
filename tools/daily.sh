#!/usr/bin/env bash
# The daily cycle on the search machine (run by launchd or cron, once a day):
#   1. stop the runner if it is still busy (it was told to start no job after planner/machine.json "until")
#   2. digest -> planner/digest.md
#   3. claude -p, no tools, context on stdin -> planner/plans/<date>.raw.json
#   4. tools/plan_check.py: validate and apply the plan, or fall back to the default queue
#   5. start the runner on queue.json in the background, until tomorrow's "until"
# Everything is logged to planner/daily.log. Safe to run by hand: `tools/daily.sh`.
set -u
export PATH="$HOME/.local/bin:/usr/local/bin:/opt/homebrew/bin:$PATH"  # cron and launchd give a bare PATH
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 1
PY="$ROOT/.venv/bin/python"; [ -x "$PY" ] || PY=python3
P="$ROOT/planner"
DAY="$(date +%Y-%m-%d)"
mkdir -p "$P/plans"
exec >>"$P/daily.log" 2>&1
echo "=== $(date '+%Y-%m-%d %H:%M:%S') daily cycle"

cfg() { "$PY" -c "import json,sys; print(json.load(open('$P/machine.json')).get('$1', '$2'))"; }
WORKERS="$(cfg workers 2)"; UNTIL="$(cfg until 06:40)"; BUDGET="$(cfg claude_budget_usd 1.0)"

# 1. the runner should be idle by now; give a running job up to 20 minutes, then stop it
if [ -f "$P/runner.pid" ] && kill -0 "$(cat "$P/runner.pid")" 2>/dev/null; then
    for _ in $(seq 1 40); do kill -0 "$(cat "$P/runner.pid")" 2>/dev/null || break; sleep 30; done
    if kill -0 "$(cat "$P/runner.pid")" 2>/dev/null; then
        echo "runner still busy: stopping it (its unfinished jobs simply run again later)"
        kill "$(cat "$P/runner.pid")"; sleep 5
    fi
fi

# 2. digest
"$PY" tools/plan_digest.py > "$P/digest.md" || { echo "digest failed"; exit 1; }

# 3. plan: the model sees prompt + playbook + hypotheses + digest, and has no tools
RAW="$P/plans/$DAY.raw.json"
if command -v claude >/dev/null 2>&1; then
    {
        cat "$P/prompt.md"; echo; echo "# Playbook"; cat "$P/playbook.md"
        echo; echo "# Hypothesis log (JSON lines)"; cat "$P/hypotheses.jsonl" 2>/dev/null
        echo; cat "$P/digest.md"
    } | claude -p --tools "" --output-format json --no-session-persistence --max-budget-usd "$BUDGET" > "$RAW"
    echo "claude exit $? -> $RAW"
else
    echo "claude CLI not found: using the fallback plan"
    : > "$RAW"
fi

# 4. validate and apply, or fall back
if ! "$PY" tools/plan_check.py "$RAW"; then
    "$PY" tools/plan_check.py --fallback
fi
cp queue.json "$P/plans/$DAY.queue.json"

# 5. run until tomorrow's deadline
nohup "$PY" tools/runner.py --queue queue.json --workers "$WORKERS" --until "$UNTIL" >> "$P/runner.log" 2>&1 &
echo $! > "$P/runner.pid"
# macOS: no idle sleep while the runner lives (a closed laptop lid still sleeps unless on power + display)
command -v caffeinate >/dev/null && nohup caffeinate -i -w "$(cat "$P/runner.pid")" >/dev/null 2>&1 &
echo "runner started (pid $(cat "$P/runner.pid"), $WORKERS workers, no new job after $UNTIL)"
