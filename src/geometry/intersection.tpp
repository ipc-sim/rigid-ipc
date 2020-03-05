#pragma once
#include "intersection.hpp"

#include <ccd/interval.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace geometry {

    template <typename T, int dim, int max_dim>
    inline T point_segment_intersection(
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& point,
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& segment_start,
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& segment_end)
    {
        // https://zalo.github.io/blog/closest-point-between-segments/
        Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim> segment_dir =
            segment_end - segment_start;
        T segment_length_sqr = segment_dir.squaredNorm();
        if (segment_length_sqr == 0.0) {
            // Segment is degenerate so return a point
            return T(0); // Either point will do
        }
        return (point - segment_start).dot(segment_dir) / segment_length_sqr;
    }

    template <typename T>
    inline bool segment_segment_intersection(
        const Eigen::Vector3<T>& segment0_start,
        const Eigen::Vector3<T>& segment0_end,
        const Eigen::Vector3<T>& segment1_start,
        const Eigen::Vector3<T>& segment1_end,
        T& alpha0,
        T& alpha1)
    {
        throw NotImplementedError(
            "segment_segment_intersection not implemented!");
    }

} // namespace geometry
} // namespace ccd

#include "intersection.tpp"
