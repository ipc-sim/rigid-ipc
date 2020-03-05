#pragma once

#include <Eigen/Core>

#include <ccd/interval.hpp>
#include <constants.hpp>
#include <utils/eigen_ext.hpp>

/// NOTE: Naming Convention
/// Point: Either a 2D or 3D point in space.
/// Segment: A line segment in either 2D or 3D defined by its endpoints.
/// Triangle: A triangle in 3D

namespace ccd {
namespace geometry {

    template <typename T, int dim, int max_dim>
    inline T point_segment_intersection(
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& point,
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& segment_start,
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& segment_end);

    bool is_point_along_segment(
        const Eigen::VectorX3<Interval>& point,
        const Eigen::VectorX3<Interval>& segment_start,
        const Eigen::VectorX3<Interval>& segment_end);

    template <typename T>
    inline bool segment_segment_intersection(
        const Eigen::Vector3<T>& segment0_start,
        const Eigen::Vector3<T>& segment0_end,
        const Eigen::Vector3<T>& segment1_start,
        const Eigen::Vector3<T>& segment1_end,
        T& alpha0,
        T& alpha1);

    bool is_point_inside_triangle(
        const Eigen::Vector3<Interval>& point,
        const Eigen::Vector3<Interval>& triangle_vertex0,
        const Eigen::Vector3<Interval>& triangle_vertex1,
        const Eigen::Vector3<Interval>& triangle_vertex2);

} // namespace geometry
} // namespace ccd

#include "intersection.tpp"
