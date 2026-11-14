#pragma once

#include "moves/move_log.hpp"
#include "moves/regions.hpp"
#include "system.hpp"

#include <functional>
#include <random>
#include <vector>

namespace pack {

struct Schedule {
    double t_hot = 1.0;
    double t_cold = 1e-10;
    long sweeps = 100000;  // one sweep = particles.size() attempted moves
};

struct StepSizes {
    double translate = 0.1;
    double rotate = 0.1;  // radians
    double box = 0.02;    // change in box side L, fixed-pressure runs only
};

struct Stats {
    long attempted = 0;
    long accepted = 0;
};

bool metropolis(double delta, double temperature, std::mt19937_64& rng);

bool try_translate(System& sys, int i, const StepSizes& steps, double temperature,
                   std::mt19937_64& rng);
bool try_rotate(System& sys, int i, const StepSizes& steps, double temperature,
                std::mt19937_64& rng);

bool try_flip(System& sys, int i, double temperature, std::mt19937_64& rng);

void tune_step(double& step, const Stats& trial, double lo, double hi);

bool try_box(System& sys, const StepSizes& steps, double pressure, double temperature,
             std::mt19937_64& rng);

// `log`, if given, gets one row entry per attempted move. It never touches the random stream.
Stats sweep(System& sys, StepSizes& steps, double temperature, std::mt19937_64& rng,
            MoveLog* log = nullptr);

bool try_translate_group(System& sys, const std::vector<int>& group, const StepSizes& steps,
                         double temperature, std::mt19937_64& rng);
bool try_rotate_group(System& sys, const std::vector<int>& group, const StepSizes& steps,
                      double temperature, std::mt19937_64& rng);

bool try_rotate_group(System& sys, const std::vector<int>& group, Vec2 pivot,
                      const StepSizes& steps, double temperature, std::mt19937_64& rng);

bool try_slide_group(System& sys, const std::vector<int>& group, Vec2 dir,
                     const StepSizes& steps, double temperature, std::mt19937_64& rng);

bool try_reflect_group(System& sys, const std::vector<int>& group, Vec2 point, Vec2 axis,
                       double temperature, std::mt19937_64& rng);

bool try_shove_region(System& sys, const RegionWeights& regions, const StepSizes& steps,
                      double temperature, std::mt19937_64& rng);

// `watch`, if given, is called after every sweep with that sweep's index.
Stats anneal(System& sys, const Schedule& sched, StepSizes& steps, std::mt19937_64& rng,
             const std::function<void(long, const System&)>& watch = {});


}  // namespace pack
