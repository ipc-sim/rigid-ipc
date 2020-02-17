// Compare timings of two functions

#include <Eigen/Geometry>
#include <igl/timer.h>

#include <ccd/hash_grid.hpp>
#include <logger.hpp>

inline bool
foo(const Eigen::AlignedBox<double, 3>& bbox1,
    const Eigen::AlignedBox<double, 3>& bbox2)
{
    return bbox1.intersects(bbox2);
}

inline bool bar(const ccd::AABB& bbox1, const ccd::AABB& bbox2)
{
    return ccd::AABB::are_overlaping(bbox1, bbox2);
}

int main(int argc, char* argv[])
{
    const int num_calls = 1000000;
    igl::Timer timer;

    Eigen::Vector3d min1, max1;
    min1 = Eigen::Vector3d::Random();
    max1 = min1.array() + 10;
    Eigen::AlignedBox<double, 3> eigen_bbox1(min1, max1);
    ccd::AABB my_bbox1(min1, max1);

    Eigen::Vector3d min2, max2;
    min2 = Eigen::Vector3d::Random();
    max2 = min2.array() + 10;
    Eigen::AlignedBox<double, 3> eigen_bbox2(min2, max2);
    ccd::AABB my_bbox2(min2, max2);

    int count = 0;
    timer.start();
    for (int i = 0; i < num_calls; i++) {
        count += foo(eigen_bbox1, eigen_bbox2);
    }
    timer.stop();
    spdlog::info(
        "function1: num_calls={:d} total_time={:g}s avg_time={:g}s", num_calls,
        timer.getElapsedTime(), timer.getElapsedTime() / num_calls);

    timer.start();
    for (int i = 0; i < num_calls; i++) {
        count += bar(my_bbox1, my_bbox2);
    }
    timer.stop();
    spdlog::info(
        "function2: num_calls={:d} total_time={:g}s avg_time={:g}s", num_calls,
        timer.getElapsedTime(), timer.getElapsedTime() / num_calls);
    std::cout << count << std::endl;
}
