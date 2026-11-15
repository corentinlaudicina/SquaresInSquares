# Playbook

Evidence-backed rules. The planner appends a line only for a confirmed hypothesis (citing its id).
Initial rules from the Phase 3 gates (docs/AUDIT.md section 15; 20 N x 2 seeds x 45 s per run):

- N >= 82: `schedules/ga-seeded.json` (6 islands seeded from N+1 minus a square, N-1 plus a square and
  lattices). It beat the plain GA 22/4/14 and shrink-capped 16/0/2 there; mean -0.50% vs its start.
- N <= 81: `schedules/shrink-capped.json` (beat the plain GA 14/6; tied with ga-seeded 10/10/2).
- Seeding needs good neighbours: an improvement at N makes N-1 and N+1 worth a seeded run.
- Diversity-aware selection (`schedules/ga-div.json`) lost to plain fitness selection 5/15/20: do not
  use it except to test a retuned threshold.
- Islands (`schedules/ga-islands.json`) were neutral against the plain GA (12/11/17), tested with
  diversity on; islands without diversity are untested.
- Symmetric genotypes (`ga-c2`, `ga-c4`) lost to the plain GA (3/22/15 and 2/16/8): the symmetric
  subspace is too narrow as is. Only worth testing with a retuned schedule; C4 needs N mod 4 in {0, 1}.
- The best/ packings are strong: expect few improvements per job. 600 s jobs are a reasonable unit; the
  GA does ~15 generations at N = 800 in 45 s, so large N need longer jobs (1800-3600 s).
- Old shrink/cycle baselines mirror their start and usually end worse than it: do not queue them.
