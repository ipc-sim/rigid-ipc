#pragma once

#ifdef PROFILE_FUNCTIONS

#include <igl/Timer.h>
#include <iostream>
#include <vector>

// To add new profile points: add a new value to the `ProfiledPoint` enum, and
// then add a name for the point in `ProfiledPointNames`.

enum ProfiledPoint {
    DETECTING_COLLISIONS_BROAD_PHASE,
    DETECTING_COLLISIONS_NARROW_PHASE,
    COMPUTING_CONSTRAINTS,
    COMPUTING_GRADIENT,
    COMPUTING_HESSIAN,
    UPDATE_SOLVE,
    SUMMING_HESSIAN,
    // LINE_SEARCH,
    _COUNT
};

static const char* ProfiledPointNames[] = {
    "Broad-phase collision detection", "Narrow-phase collision detection",
    "Computing constraint values", "Computing constraint gradients",
    "Computing constraint hessians", "Solving for an update",
    "Summing hessians",
    // "Line search"
};

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
