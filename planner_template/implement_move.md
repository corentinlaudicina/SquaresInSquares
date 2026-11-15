You are implementing one new GA move in a C++20 square-packing engine, in a git worktree on its own
branch. A human reviews and merges your branch; nothing you do reaches main by itself. The move
proposal is at the end of this prompt.

Read first: src/pop/moves.hpp, src/pop/moves.cpp, src/pop/ga.hpp, src/engine/schedule.cpp (the
`fields()` table), tests/test_ga.cpp, and docs/RUNBOOK.md.

Do exactly this:
1. Add the move to `GaMove` in src/pop/moves.hpp, as the last entry before `Count`, and its name in
   `ga_move_name`. Append its weight to `MoveParams::weight` with the value 0.0, so every existing
   schedule behaves exactly as before. Add any parameters it needs to `MoveParams` with defaults.
2. Implement it in `apply_move` (src/pop/moves.cpp). Rules: draw randomness only through `u01`,
   `uniform`, `uniform_int`, `normal` (portable across platforms); keep the square count; the move may
   leave overlaps (relaxation repairs them) but must leave every value finite; if it cannot apply,
   fall back to a small jiggle as Crossover does.
3. Add the schedule keys in src/engine/schedule.cpp: `ga.moves.<name>` for the weight (curve-able) and
   `ga.move.<param>` for each parameter.
4. Add a test to tests/test_ga.cpp: the move keeps N and finiteness over many draws, and with its weight
   0 a GA run gives the same result as before (compare against the same run on main's behaviour by
   checking the move is never drawn).
5. Write schedules/planner/ga-<name>.json: schedules/ga.json with the new weight set to about 0.10 and
   the other weights scaled down to keep the total, `name` and `description` set.
6. Run `make -j8 all` and `make check`; fix until both pass. Do not change any other behaviour, file or
   default. Commit on the current branch with a message that cites the proposal.
Then stop; the script that started you runs the benchmark gate.
