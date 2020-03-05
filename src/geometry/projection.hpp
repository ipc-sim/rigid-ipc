#pragma once

#include <Eigen/Core>

#include <utils/eigen_ext.hpp>

/// NOTE: Naming Convention
/// Point: Either a 2D or 3D point in space.
/// Segment: A line segment in either 2D or 3D defined by its endpoints.
/// Triangle: A triangle in 3D

namespace ccd {
namespace geometry {

    /// Project the point to the line returning the parameter along the line.
    template <typename T, int dim, int max_dim>
    inline T project_point_to_line(
        const Eigen::Vector<T, dim, max_dim>& point,
        const Eigen::Vector<T, dim, max_dim>& line_point,
        const Eigen::Vector<T, dim, max_dim>& line_dir);

    template <typename T>
    inline Eigen::Vector3<T> project_segment_to_segment(
        const Eigen::Vector3<T>& segment0_start,
        const Eigen::Vector3<T>& segment0_end,
        const Eigen::Vector3<T>& segment1_start,
        const Eigen::Vector3<T>& segment1_end,
        T& alpha0,
        T& alpha1);

    /// Project the point to the plane returning the projected point.
    /// NOTE: The normal need not be unit length.
    template <typename T>
    inline Eigen::Vector3<T> project_point_to_plane(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& plane_point,
        const Eigen::Vector3<T>& normal);

    template <typename T>
    inline Eigen::Vector3<T> project_point_to_plane(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& plane_point,
        const Eigen::Vector3<T>& normal,
        const T& normal_sqrnorm);

} // namespace geometry
} // namespace ccd

#include "projection.tpp"
