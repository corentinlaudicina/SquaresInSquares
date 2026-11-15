# The daily planner

One `claude -p` turn a day decides what the search machine runs next. Scripts do everything else:
the model never reads raw logs, never touches code or `best/`, and has **no tools**. It reads a
short context on stdin and answers with one JSON plan, which a script validates before anything is
written.

## The cycle (`tools/daily.sh`, 07:00 by launchd or cron)
1. The runner was started with `--until 06:40`, so it starts no job after that; a job still running
   gets up to 20 minutes, then is stopped (it simply runs again later).
2. `tools/plan_digest.py` writes `planner/digest.md` (~80 lines): machine and budget, coverage of the
   open N, improvements per job-hour by schedule and N range, results by tag (hypotheses), the top 30
   candidate N, recent improvements, failures, and the schedules available.
3. `claude -p --tools "" --output-format json --no-session-persistence --max-budget-usd <cap>` gets
   `planner/prompt.md` + `planner/playbook.md` + `planner/hypotheses.jsonl` + the digest and answers
   with the plan (schema in the prompt).
4. `tools/plan_check.py` validates and applies it; if it is rejected (or claude is unavailable) the
   deterministic fallback runs instead: the top-ranked open N, `shrink-capped` below N = 82, `ga` above,
   600 s each, filling the budget.
5. The runner starts on `queue.json` in the background until the next 06:40.

## What counts as work
- Only N whose status in `data/known_optima_1_1000.csv` is `open` (822 of 1000). Every other N is a
  proven optimum: the digest never lists it and `plan_check` rejects any job for it.
- There are no record values for open N: priorities come from our own data. The default score favours
  N whose packing wastes more width than its neighbours' (`slack`: a bump in side - sqrt(N)), N whose
  neighbours just improved (seeding passes gains to N +- 1), and N never tried.

## Rules `plan_check.py` enforces
- open N only; schedules exist and the engine accepts them (`engine --check-schedule`);
- variants change only existing keys of an existing schedule, live in `schedules/planner/`, and run
  only in experiment jobs tagged `H<id>`;
- total job time / workers <= budget (+5%); no N over 20% of the plan; experiments <= 25%;
  30 s <= time <= 3600 s; seeds is a count (1-8), turned into fresh seed values per day.

## New schedules and new moves
- **Schedules** are data: the planner designs them freely in `new_schedules` (a `base` file plus any
  existing keys in `changes`, including the mode, islands, seeding, symmetry, moves, relaxation,
  ratchet). They run only in experiment jobs (`H<id>`) until a confirmed hypothesis puts them in the
  playbook; from then on the planner uses them in exploit jobs.
- **Moves** are code, so the planner only *designs* them (`move_proposals`: idea, mechanism, why);
  `plan_check` saves each to `planner/proposals/<date>-<name>.md` and the digest lists their status.
  Building one is a separate, hand-started step:
  ```sh
  tools/implement_move.sh planner/proposals/2026-11-16-row-shift.md
  ```
  It creates a git worktree on branch `move/<name>`, runs `claude -p` there with edit and shell tools
  (`planner_template/implement_move.md`: the move goes in with weight 0, plus keys, a test and
  `schedules/planner/ga-<name>.json`), runs `make check`, then a recovery gate of that schedule against
  `schedules/ga.json` on 6 N x 2 seeds. The result is written back into the proposal. You review the
  branch and merge it; after that the planner can weight `ga.moves.<name>` in a new schedule.
  The implementer has an unrestricted shell inside its worktree (spending cap `CLAUDE_BUDGET_USD`,
  default $5): it is started by you, not by the daily cycle, for that reason.

## State (per machine, gitignored, under `planner/`)
| file | what |
|---|---|
| `machine.json` | workers, cores, `until`, budget hours, claude spending cap per turn |
| `prompt.md`, `playbook.md` | copied from `planner_template/` at setup; the playbook grows only with confirmed hypotheses |
| `hypotheses.jsonl` | the planner's hypothesis log, updated by `plan_check` |
| `proposals/<date>-<name>.md` | move proposals and their status (proposed, implemented, gated) |
| `digest.md`, `plans/<date>.raw.json`, `plans/<date>.queue.json` | each day's input, answer and queue |
| `daily.log`, `runner.log`, `runner.pid` | logs; `kill $(cat planner/runner.pid)` stops the search |

## Setting up the spare machine
```sh
git clone <this repo> && cd SquaresInSquares
tools/setup_machine.sh --schedule   # tools, build, tests, venv, best/ (import ~20 min), planner/, 07:00 job
claude                              # log the CLI in once, then /exit
tools/daily.sh                      # first cycle now, by hand
```
To continue from this Mac's `best/` instead of re-importing, copy the folder over before running the
setup (it keeps an existing `best/`). Setup is safe to re-run. On Linux it uses cron and g++ if clang++ is
absent (the g++ build is untested so far).
