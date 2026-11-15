# Playbook

Evidence-backed rules. The planner appends a line only for a confirmed hypothesis (citing its id).
Initial rules from the Phase 3 gates (docs/AUDIT.md section 15; 20 N x 2 seeds x 45 s per run):

- **Production schedule: `schedules/ga-seeded-adapt.json` for every N** (6 seeded islands on a ring and
  the adaptive ratchet). The machine enforces an island GA for exploit and explore jobs, and focuses on
  N > 300 for now (see the digest).
- Never use `ga-seeded` (fixed ratchet) or `shrink-capped` (mirrored start) from best/: from a
  near-optimal packing the fixed ratchet's first target (0.02% below) never fits and the mirror destroys
  the start. The first 22 live runs ended at best/ exactly (ga-seeded) or worse (shrink-capped, median
  +1.1e-3).
- Polish test (start = best/N scaled up by 1e-4): at N = 491 / 600 s ga-seeded-adapt closed 96% of the
  gap (36 legal improvements); at N = 300 / 180 s it closed 90%, while `shrink-polish` (no mirror) closed
  99%. Comparing the two at N > 300 with longer jobs is an open question for an experiment.
- Diversity-aware selection (`schedules/ga-div.json`) lost to plain fitness selection 5/15/20: do not
  use it except to test a retuned threshold.
- Islands (`schedules/ga-islands.json`) were neutral against the plain GA (12/11/17), tested with
  diversity on; islands without diversity are untested.
- Symmetric genotypes (`ga-c2`, `ga-c4`) lost to the plain GA (3/22/15 and 2/16/8): the symmetric
  subspace is too narrow as is. Only worth testing with a retuned schedule; C4 needs N mod 4 in {0, 1}.
- The best/ packings are strong: expect few improvements per job. 600 s jobs are a reasonable unit; the
  GA does ~15 generations at N = 800 in 45 s, so large N need longer jobs (1800-3600 s).
- Seeding needs good neighbours: an improvement at N makes N-1 and N+1 worth a seeded run.
- Old shrink/cycle baselines mirror their start and usually end worse than it: do not queue them.
