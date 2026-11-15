# Lab notebook

A dated log of what was tried, what the measurements said, and what was decided. Newest entries at the
bottom. Benchmark numbers are "wins / losses / ties" on the legal (strictly valid) side, paired by N and
seed, unless stated otherwise. Details of every gate are in `docs/AUDIT.md`, section 15. No packings are
published here, only side lengths and counts.

---

## 2026-11-14: audit of the existing engine

Starting point: a C++ engine with Monte Carlo moves, large-neighbourhood moves (destroy-and-repair,
strips, corner reshapes, loosen), an L-BFGS relaxation of an overlap-area penalty, and two drivers
(`squares_shrink`, `squares_compress`). The audit (`docs/AUDIT.md`) found, among 15 issues:

- Every driver mirrors its start through the diagonal before searching. The 1000 starting packings are
  already valid and tight; mirroring breaks them (overlap up to almost a whole square), so the old
  search could not improve an asymmetric start. **Decision:** the mirror is deliberate (symmetric
  packings are common among good solutions); keep it as a schedule option.
- Validity was judged by overlap *area* <= 1e-10, which allows ~1e-5 of overlap *depth* at a tilted
  corner. Stored sides were optimistic by up to ~1e-5. **Decision:** results must be strictly valid;
  a legalization step is required before anything is kept.

## 2026-11-14: foundation

- A validator independent of the engine (exact rational arithmetic on the float corners). On the 1000
  seeds: 687 strictly valid, the other 313 out of the box by at most 1.8e-15 (float rounding of corners
  on the wall).
- An engine CLI with JSON schedules (the two old drivers reproduced bit for bit: 12/12), JSON telemetry,
  and a legalizer (inflated squares, cumulative box growth, bisection).
- Finding: legalizing the old engine's certified packings costs real side, +8.6e-6 at n = 29 and +1.1e-4 at
  n = 1000. It is not the legalizer being clumsy: the search measured its side on overlapping packings.
- A runner that keeps only validated strict improvements, and a benchmark with paired seeds.

## 2026-11-14: relaxation (Phase 3 item 2)

- New penalty: squared penetration depth (smallest separating-axis overlap) instead of squared overlap
  area, which grows like depth^4 at a corner contact. With a little slack the new relaxation reaches
  exactly zero overlap; the area penalty stalls at ~2e-8 of depth.
- Staged L-BFGS (history 3, hard step cap, no line search). At N = 200, 15 ms per relaxation on one
  thread, 247 per second on 8, against 908 ms for the old relaxation.
- Gate (20 N x 2 seeds x 45 s, starts = best known): `shrink-capped` (new penalty and relaxer, accept on
  the legal side) beat the old baseline **25 / 3 / 12** (26 / 2 / 12 in the later recovery gate).
  Accepting on the legal side alone did not help (9 / 13 / 18).

## 2026-11-14: a gate that measured nothing

The first GA gate started every method from the best known packings. The GA "won" 36 / 0 / 4, but it
ratcheted **zero times in 80 runs**: it simply returned the strong start, while the baseline's mirror
degraded its own. Lesson: starting from the best cannot show a search improving anything.
**Decision:** gate on *recovery*: every method starts from the same legal but worse and different
arrangements, scored against the best known.

## 2026-11-15: population search (Phase 3 items 3-8), recovery gates

Mean change of the side versus the start: old baseline +0.26% (worse), shrink-capped +0.08%, plain GA
-0.19%, seeded GA -0.50%.

| item | comparison | result | decision |
|---|---|---|---|
| 3-4: ratchet GA | vs old baseline | 32 / 4 / 4 | **kept**; better than shrink-capped for N >= 82 (13 / 3 / 2), worse for N <= 65 (6 / 14 / 2) |
| 5: diversity-aware selection | vs plain GA | 5 / 15 / 20 | off by default; selection cost and untuned thresholds |
| 6: islands (with diversity) | vs plain GA | 12 / 11 / 17 | inconclusive; islands without diversity untested |
| 7: seeding from N+1, N-1, lattices | vs plain GA | 22 / 4 / 14 | **kept**; vs shrink-capped 16 / 0 / 2 for N >= 82 |
| 8: C2 / C4 symmetric genotypes | vs plain GA | 3 / 22 / 15, 2 / 16 / 8 | off by default |

Kernel optimisations (gradient only for the deepest axis; a neighbour list reused during a relaxation):
relaxation 1.33x faster on one thread, 1.45x on eight.

**Defaults:** shrink-capped for N <= 81, seeded GA for N >= 82.

## 2026-11-15: unattended search starts

- Scope: only N whose status in `data/known_optima_1_1000.csv` is `open` (822 of 1000); proven optima are
  never worked on. There are no record values for the open N, so priorities come from our own data: N
  whose best side wastes more width than its neighbours' (side - sqrt(N), against the median of N +- 5).
- The daily planner (one `claude -p` turn, no tools, plan validated by a script) wrote its first plan:
  106 runs, ~21 h on 2 workers; 60% exploit, 24% explore, 15% experiments. The first plan was rejected
  because it guessed six non-open N; the digest now lists the open N, and non-open N are dropped with a
  warning instead of rejecting the plan.
- Open hypotheses: **H1** at 600 s the seeded GA beats shrink-capped on small high-slack N (the 45 s gate
  may have been too short for the GA); **H2** at N >= 600, 1800 s jobs give as many improvements per hour
  as 3600 s jobs.

---

Daily entries from here on are drafted by `tools/notebook_entry.py` and published with
`tools/publish_notebook.sh`.

## 2026-11-15: the first live protocol could not produce a record

After 22 runs (4.7 job-hours) no improvement, and the ledger showed why, by construction:
`shrink-capped` inherits the mirrored start, so it rebuilt every record from a broken copy (all 18 runs
ended worse, median +1.1e-3); `ga-seeded` aimed its first ratchet 0.02% below the record, which a
near-optimal packing never fits (all runs ended exactly at the record). The recovery gates had started
about 1% worse than the best, where neither problem shows.

- Fix: an **adaptive ratchet** per island (start 1e-6 below the last fit, x1.5 after a fit, halve and
  loosen after 8 generations without one, floor 1e-9). On a 5 x 5 grid with 1e-4 of slack the fixed
  ratchet stays at 5.000400, the adaptive one reaches 5.0000044 (optimum 5).
- Polish test (start = record scaled up by 1e-4): N = 491, 600 s: closed 96% of the gap. N = 300, 180 s:
  adaptive island GA 90%, single-packing `shrink-polish` (no mirror) 99%.
- **Decision** (Corentin): island GA on essentially all N, focus on N > 300 for now. Production schedule
  `ga-seeded-adapt`; the machine enforces an island GA for exploit and explore jobs.
