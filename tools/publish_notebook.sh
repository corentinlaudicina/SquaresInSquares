#!/usr/bin/env bash
# Appends the drafted daily entries (planner/notebook/*.md) not yet published to docs/NOTEBOOK.md,
# commits that file alone and pushes. By hand only: this is what makes results public.
#   tools/publish_notebook.sh            (add -n to print what would be appended, without committing)
set -eu
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
DONE=planner/notebook/.published; touch "$DONE"
NEW=$(ls planner/notebook/*.md 2>/dev/null | sort | grep -vxF -f "$DONE" || true)
[ -n "$NEW" ] || { echo "nothing new to publish"; exit 0; }
if [ "${1:-}" = "-n" ]; then cat $NEW; exit 0; fi
for f in $NEW; do { echo; cat "$f"; } >> docs/NOTEBOOK.md; echo "$f" >> "$DONE"; done
git add docs/NOTEBOOK.md
git commit -q -m "Notebook: $(echo "$NEW" | xargs -n1 basename | sed 's/\.md$//' | tr '\n' ' ')"
git push -q && echo "published: $(echo "$NEW" | wc -l | tr -d ' ') entries"
