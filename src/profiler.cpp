#include <profiler.hpp>

#ifdef PROFILE_FUNCTIONS

double time_spent_at_profiled_points[int(ProfiledPoint::_COUNT)];
long number_of_evals_profiled_points[int(ProfiledPoint::_COUNT)];

void reset_profile_point(ProfiledPoint point)
{
    time_spent_at_profiled_points[point] = 0;
    number_of_evals_profiled_points[point] = 0;
}

void reset_profiler()
{
    for (int i = 0; i < int(ProfiledPoint::_COUNT); i++) {
        reset_profile_point(ProfiledPoint(i));
    }
}

void print_profile_point(ProfiledPoint point, double total_time)
{
    std::cout << ProfiledPointNames[int(point)] << ":\n"
              << "\tNumber of evals: "
              << number_of_evals_profiled_points[int(point)]
              << "\n\tTime: " << time_spent_at_profiled_points[int(point)]
              << " seconds ("
              << (time_spent_at_profiled_points[int(point)] / total_time * 100)
              << "%)\n\tAverage: "
              << (time_spent_at_profiled_points[int(point)]
                     / number_of_evals_profiled_points[int(point)])
              << " seconds per eval" << std::endl;
}

void print_profile(double total_time)
{
    std::cout << "Total time spent optimizing displacements: " << total_time
              << " seconds" << std::endl;

    for (int i = 0; i < int(ProfiledPoint::_COUNT); i++) {
        print_profile_point(ProfiledPoint(i), total_time);
    }
    std::cout << std::endl;
}

#endif
