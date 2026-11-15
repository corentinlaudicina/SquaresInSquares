# Playbook

Evidence-backed rules. The planner appends a line only for a confirmed hypothesis (citing its id).
Initial rules from the Phase 3 gates (docs/AUDIT.md section 15; 20 N x 2 seeds x 45 s per run):

- N <= 81: `schedules/shrink-capped.json` (beat the GA 14/6 there). N >= 82: `schedules/ga.json`
  (beat shrink-capped 13/3, and was the only family that improved most starts: -0.17% mean).
- Diversity-aware selection (`schedules/ga-div.json`) lost to plain fitness selection 5/15/20: do not
  use it except to test a retuned threshold.
- Islands (`schedules/ga-islands.json`) were neutral against the plain GA (12/11/17), tested with
  diversity on; islands without diversity are untested.
- Seeding from N+1 / N-1 / lattices (`ga-seeded`) and symmetric genotypes (`ga-c2`, `ga-c4`): see the
  items 7-8 gate in docs/AUDIT.md; C4 only applies to N mod 4 in {0, 1}.
- The best/ packings are strong: expect few improvements per job. 600 s jobs are a reasonable unit; the
  GA does ~15 generations at N = 800 in 45 s, so large N need longer jobs (1800-3600 s).
- Old shrink/cycle baselines mirror their start and usually end worse than it: do not queue them.
