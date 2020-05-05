#pragma once
#include "distance_barrier_constraint.hpp"

#include <barrier/barrier.hpp>
#include <geometry/distance.hpp>
#include <logger.hpp>

namespace ccd {
namespace opt {

    template <typename T>
    T DistanceBarrierConstraint::distance_barrier(
        const T& distance, const double eps) const
    {
#ifdef USE_DISTANCE_SQUARED
        return barrier(distance - min_distance, eps * eps, barrier_type);
#else
        return barrier(distance - min_distance, eps, barrier_type);
#endif
    }

} // namespace opt
} // namespace ccd
