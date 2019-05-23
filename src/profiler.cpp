#include <profiler.hpp>

#ifdef PROFILE_FUNCTIONS

double time_spent_at_profiled_points[int(ProfiledPoint::_COUNT)];
long number_of_evals_profiled_points[int(ProfiledPoint::_COUNT)];

void reset_profiler()
{
    for (int i = 0; i < int(ProfiledPoint::_COUNT); i++) {
        time_spent_at_profiled_points[i] = 0;
        number_of_evals_profiled_points[i] = 0;
    }
}

void print_profile(double total_time)
{
    std::cout << "Total time spent optimizing displacements: " << total_time
              << " seconds" << std::endl;

    for (int i = 0; i < int(ProfiledPoint::_COUNT); i++) {
        std::cout << ProfiledPointNames[i] << ":\n"
                  << "\tNumber of evals: " << number_of_evals_profiled_points[i]
                  << "\n\tTime: " << time_spent_at_profiled_points[i]
                  << " seconds ("
                  << (time_spent_at_profiled_points[i] / total_time * 100)
                  << "%%)\n\tAverage: "
                  << (time_spent_at_profiled_points[i]
                         / number_of_evals_profiled_points[i])
                  << " seconds per eval" << std::endl;
    }
    std::cout << std::endl;
}

#endif
