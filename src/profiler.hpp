#pragma once

#define PROFILE_FUNCTIONS

#ifdef PROFILE_FUNCTIONS

#include <igl/Timer.h>
#include <iostream>
#include <vector>

enum ProfiledPoint {
    DETECTING_COLLISIONS,
    COMPUTING_CONSTRAINTS,
    COMPUTING_GRADIENT,
    COMPUTING_HESSIAN,
    UPDATE_SOLVE,
    SUMMING_HESSIAN,
    _COUNT
};

static const char* ProfiledPointNames[]
    = { "Detecting collisions", "Computing constraint values",
          "Computing constraint gradients", "Computing constraint hessians",
          "Solving for an update", "Summing hessians" };

extern double time_spent_at_profiled_points[];
extern long number_of_evals_profiled_points[];

void reset_profiler();
void print_profile(double total_time);

#define PROFILE(op, point)                                                     \
    {                                                                          \
        number_of_evals_profiled_points[int(point)]++;                         \
        igl::Timer timer;                                                      \
        timer.start();                                                         \
        op;                                                                    \
        timer.stop();                                                          \
        time_spent_at_profiled_points[int(point)] += timer.getElapsedTime();   \
    }

#else

#define PROFILE(op, point)                                                     \
    {                                                                          \
        op;                                                                    \
    }

#endif
