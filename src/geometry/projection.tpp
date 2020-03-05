#pragma once
#include "projection.tpp"

#include <utils/not_implemented_error.hpp>

/// NOTE: Naming Convention
/// Point: Either a 2D or 3D point in space.
/// Segment: A line segment in either 2D or 3D defined by its endpoints.
/// Triangle: A triangle in 3D

namespace ccd {
namespace geometry {

    template <typename T, int dim, int max_dim>
    inline T project_point_to_line(
        const Eigen::Vector<T, dim, max_dim>& point,
        const Eigen::Vector<T, dim, max_dim>& line_point,
        const Eigen::Vector<T, dim, max_dim>& line_dir)
    {
        // https://zalo.github.io/blog/closest-point-between-segments/
        T dir_length_sqr = line_dir.squaredNorm();
        if (dir_length_sqr == 0.0) {
            // Segment is degenerate so return a point
            return T(0); // Either point will do
        }
        return (point - line_point).dot(line_dir) / dir_length_sqr;
    }

    template <typename T>
    inline Eigen::Vector3<T> project_segment_to_segment(
        const Eigen::Vector3<T>& segment0_start,
        const Eigen::Vector3<T>& segment0_end,
        const Eigen::Vector3<T>& segment1_start,
        const Eigen::Vector3<T>& segment1_end,
        T& alpha0,
        T& alpha1)
    {
        throw NotImplementedError(
            "project_segment_to_segment not implemented!");
    }

    template <typename T>
    inline Eigen::Vector3<T> project_point_to_plane(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& plane_point,
        const Eigen::Vector3<T>& normal)
    {
        return project_point_to_plane(
            point, plane_point, normal, normal.squaredNorm());
    }

    template <typename T>
    inline Eigen::Vector3<T> project_point_to_plane(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& plane_point,
        const Eigen::Vector3<T>& normal,
        const T& normal_sqrnorm)
    {
        return point
            - (point - plane_point).dot(normal) / normal_sqrnorm * normal;
    }

} // namespace geometry
} // namespace ccd
