#pragma once
#include "distance_barrier_constraint.hpp"

#include <barrier/barrier.hpp>
#include <geometry/distance.hpp>
#include <logger.hpp>

namespace ccd {
namespace opt {

    template <typename T>
    T DistanceBarrierConstraint::distance_barrier(
        const T& distance, const double dhat) const
    {
#ifdef USE_DISTANCE_SQUARED
        return barrier(distance, dhat * dhat, barrier_type);
#else
        return barrier(distance, dhat, barrier_type);
#endif
    }

} // namespace opt
} // namespace ccd
