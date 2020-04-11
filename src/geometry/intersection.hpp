#pragma once

#include <ECCD.hpp>
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

    bool is_point_along_segment(
        const Eigen::VectorX3I& point,
        const Eigen::VectorX3I& segment_start,
        const Eigen::VectorX3I& segment_end);

    bool are_segments_intersecting(
        const Eigen::Vector3I& segment0_start,
        const Eigen::Vector3I& segment0_end,
        const Eigen::Vector3I& segment1_start,
        const Eigen::Vector3I& segment1_end);

    bool is_point_inside_triangle(
        const Eigen::Vector3I& point,
        const Eigen::Vector3I& triangle_vertex0,
        const Eigen::Vector3I& triangle_vertex1,
        const Eigen::Vector3I& triangle_vertex2);

    bool segment_triangle_intersect(
        const Eigen::Vector3d& segment_vertex0,
        const Eigen::Vector3d& segment_vertex1,
        const Eigen::Vector3d& triangle_vertex0,
        const Eigen::Vector3d& triangle_vertex1,
        const Eigen::Vector3d& triangle_vertex2);

} // namespace geometry
} // namespace ccd
