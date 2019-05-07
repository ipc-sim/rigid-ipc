#include <profiler.hpp>

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
              << " seconds per call"

              << std::endl;
}
