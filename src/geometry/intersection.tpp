#pragma once
#include "intersection.hpp"

#include <ccd/interval.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace geometry {

    template <typename T>
    inline void point_segment_intersection(
        const Eigen::VectorX3<T>& point,
        const Eigen::VectorX3<T>& segment_start,
        const Eigen::VectorX3<T>& segment_end,
        T& alpha)
    {
        // Project the point onto the edge by computing its scalar projection
        Eigen::VectorX3<T> segment_vec = segment_end - segment_start;
        alpha = (point - segment_start).dot(segment_vec)
            / segment_vec.squaredNorm();
    }

    template <typename T>
    inline bool segment_segment_intersection(
        const Eigen::VectorX3<T>& segment0_start,
        const Eigen::VectorX3<T>& segment0_end,
        const Eigen::VectorX3<T>& segment1_start,
        const Eigen::VectorX3<T>& segment1_end,
        T& alpha0,
        T& alpha1)
    {
        throw NotImplementedError(
            "segment_segment_intersection not implemented!");
    }

} // namespace geometry
} // namespace ccd

#include "intersection.tpp"
