#pragma once

#include "distance_barrier_constraint.hpp"
#include <opt/barrier.hpp>

namespace ccd {
namespace opt {

    template <typename T>
    T DistanceBarrierConstraint::distance_barrier(
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& c)
    {
        T distance = sqrt(point_to_edge_sq_distance<T>(a, b, c));
        return opt::spline_barrier<T>(distance, m_barrier_epsilon);
    }

} // namespace opt
} // namespace ccd
