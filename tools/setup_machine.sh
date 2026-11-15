#!/usr/bin/env bash
# One-time setup of a search machine (macOS or Linux), from a fresh clone of the repo:
#   tools/setup_machine.sh                 check, build, test, create best/ and planner/
#   tools/setup_machine.sh --schedule      ... and install the daily 07:00 run (launchd or cron)
# Re-running it is safe: nothing that exists (best/, planner state) is overwritten.
set -eu
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
SCHEDULE=0; [ "${1:-}" = "--schedule" ] && SCHEDULE=1
OS="$(uname -s)"
say() { printf '\n== %s\n' "$*"; }
die() { printf 'setup: %s\n' "$*" >&2; exit 1; }

say "1. tools"
command -v git >/dev/null || die "git missing"
command -v make >/dev/null || die "make missing"
if command -v clang++ >/dev/null; then CXX=clang++; elif command -v g++ >/dev/null; then CXX=g++; else die "no C++ compiler (clang++ or g++ >= 10)"; fi
command -v python3 >/dev/null || die "python3 missing"
python3 -c 'import sys; assert sys.version_info >= (3, 9)' || die "python3 >= 3.9 needed"
if [ "$OS" = Darwin ]; then CORES="$(sysctl -n hw.ncpu)"; else CORES="$(nproc)"; fi
echo "compiler $CXX, $CORES cores, $OS"
command -v claude >/dev/null && echo "claude CLI: $(claude --version 2>/dev/null)" || echo "claude CLI not found: the planner will use its fallback plan until you install it"

say "2. build and test"
make -j"$CORES" CXX="$CXX" all
make CXX="$CXX" check

say "3. python environment"
[ -x .venv/bin/python ] || python3 -m venv .venv
.venv/bin/pip install -q -r tools/requirements.txt

say "4. best/"
COUNT="$(ls best/*.json 2>/dev/null | wc -l | tr -d ' ')"
if [ "$COUNT" -lt 1000 ]; then
    echo "best/ holds $COUNT packings: importing the seeds (legalized, validated; ~20 min)."
    echo "(To continue from another machine's best/ instead, stop here and copy that folder over.)"
    .venv/bin/python tools/runner.py --import-legacy
else
    echo "best/ already holds $COUNT packings: kept"
fi

say "5. planner state"
mkdir -p planner/plans
[ -f planner/prompt.md ] || cp planner_template/prompt.md planner/prompt.md
[ -f planner/playbook.md ] || cp planner_template/playbook.md planner/playbook.md
[ -f planner/hypotheses.jsonl ] || : > planner/hypotheses.jsonl
if [ ! -f planner/machine.json ]; then
    WORKERS=$(( CORES / 4 )); [ "$WORKERS" -ge 1 ] || WORKERS=1
    printf '{\n  "workers": %d,\n  "cores": %d,\n  "until": "06:40",\n  "plan_at": "07:00",\n  "budget_hours": 23.3,\n  "claude_budget_usd": 1.0\n}\n' "$WORKERS" "$CORES" > planner/machine.json
fi
cat planner/machine.json

if [ "$SCHEDULE" = 1 ]; then
    say "6. daily schedule (07:00)"
    if [ "$OS" = Darwin ]; then
        PLIST="$HOME/Library/LaunchAgents/local.squares.planner.plist"
        sed -e "s#@ROOT@#$ROOT#g" -e "s#@HOME@#$HOME#g" tools/planner.launchd.plist > "$PLIST"
        launchctl unload "$PLIST" 2>/dev/null || true
        launchctl load "$PLIST"
        echo "launchd job installed: $PLIST"
    else
        ( crontab -l 2>/dev/null | grep -v 'tools/daily.sh'; echo "0 7 * * * $ROOT/tools/daily.sh" ) | crontab -
        echo "cron entry installed: 0 7 * * * $ROOT/tools/daily.sh"
    fi
fi

say "done"
cat <<EOF
Next:
  - log the claude CLI in once, interactively:   claude        (then /exit)
  - first cycle by hand (plans, then starts the runner in the background):   tools/daily.sh
  - watch:  tail -f planner/daily.log   |   runs/summary.md   |   runs/ledger.jsonl
  - stop:   kill \$(cat planner/runner.pid)
EOF
