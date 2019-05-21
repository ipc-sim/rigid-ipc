#include <profiler.hpp>

#ifdef PROFILE_FUNCTIONS

void reset_profiler()
{
    // Computing the constraints
    number_of_constraint_calls = 0;
    time_spent_computing_constraint = 0;
    number_of_gradient_calls = 0;
    time_spent_computing_gradient = 0;
    number_of_hessian_calls = 0;
    time_spent_computing_hessian = 0;

    // CCD
    number_of_collision_detection_calls = 0;
    time_spent_detecting_collisions = 0;

    // Optimization
    number_of_update_solves = 0;
    time_spent_solving_for_update = 0;
    number_of_hessian_summations = 0;
    time_spent_summing_hessians = 0;
}

void print_profile(double total_time)
{
    std::cout << std::endl
              << "Total time spent optimizing displacements: " << total_time
              << " seconds\n"

              << "Detecting collisions:\n"
              << "\tNumber of calls: " << number_of_collision_detection_calls
              << "\n\tTime: " << time_spent_detecting_collisions << " seconds ("
              << time_spent_detecting_collisions / total_time * 100 << "%)\n"
              << "\tAverage: "
              << time_spent_detecting_collisions
            / number_of_collision_detection_calls
              << " seconds per call\n"

              << "Computing constraints:\n"
              << "\tValue:\n"
              << "\t\tNumber of calls: " << number_of_constraint_calls
              << "\n\t\tTime: " << time_spent_computing_constraint
              << " seconds ("
              << time_spent_computing_constraint / total_time * 100 << "%)\n"
              << "\t\tAverage: "
              << time_spent_computing_constraint / number_of_constraint_calls
              << " seconds per call\n"

              << "\tGradient:\n"
              << "\t\tNumber of calls: " << number_of_gradient_calls
              << "\n\t\tTime: " << time_spent_computing_gradient << " seconds ("
              << time_spent_computing_gradient / total_time * 100 << "%)\n"
              << "\t\tAverage: "
              << time_spent_computing_gradient / number_of_gradient_calls
              << " seconds per call\n"

              << "\tHessian:\n"
              << "\t\tNumber of calls: " << number_of_hessian_calls
              << "\n\t\tTime: " << time_spent_computing_hessian << " seconds ("
              << time_spent_computing_hessian / total_time * 100 << "%)\n"
              << "\t\tAverage: "
              << time_spent_computing_hessian / number_of_hessian_calls
              << " seconds per call\n"

              << "Solving for an update:\n"
              << "\tNumber of solves: " << number_of_update_solves
              << "\n\tTime: " << time_spent_solving_for_update << " seconds ("
              << time_spent_solving_for_update / total_time * 100 << "%)\n"
              << "\tAverage: "
              << time_spent_solving_for_update / number_of_update_solves
              << " seconds per call\n"

              << "Summing hessians:\n"
              << "\tNumber of summations: " << number_of_hessian_summations
              << "\n\tTime: " << time_spent_summing_hessians << " seconds ("
              << time_spent_summing_hessians / total_time * 100 << "%)\n"
              << "\tAverage: "
              << time_spent_summing_hessians / number_of_hessian_summations
              << " seconds per call\n"
              << std::flush;
}

#endif
