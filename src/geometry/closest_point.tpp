#pragma once

#include "distance.hpp"

#include <utils/not_implemented_error.hpp>

namespace ccd {

template <typename T>
T point_segment_closest_point(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& point,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment_vertex0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment_vertex1)
{
    if (point.size() == 2) {

    } else {
        assert(point.size() == 3);
        throw NotImplementedError(
            "point_segment_closest_point() not implmented!");
    }
}

template <typename T>
T segment_segment_closest_points(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment0_vertex0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment0_vertex1,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment1_vertex0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment1_vertex1,
    Eigen::Matrix<T, Eigen::Dynamic, 1>& segment0_point,
    Eigen::Matrix<T, Eigen::Dynamic, 1>& segment1_point)
{
    throw NotImplementedError(
        "segment_segment_closest_point() not implmented!");
}

template <typename T>
T point_triangle_closest_point(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& point,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& triangle_vertex0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& triangle_vertex1,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& triangle_vertex2)
{
    throw NotImplementedError("point_triangle_closest_point() not implmented!");
}

} // namespace ccd
