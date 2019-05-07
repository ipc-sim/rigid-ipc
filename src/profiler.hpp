#pragma once

#define PROFILE_FUNCTIONS

#ifdef PROFILE_FUNCTIONS

#include <igl/Timer.h>
#include <iostream>

// Computing the constraints
extern long number_of_constraint_calls;
extern double time_spent_computing_constraint;
extern long number_of_gradient_calls;
extern double time_spent_computing_gradient;
extern long number_of_hessian_calls;
extern double time_spent_computing_hessian;

// CCD
extern long number_of_collision_detection_calls;
extern double time_spent_detecting_collisions;

void reset_profiler();
void print_profile(double total_time);

#endif
