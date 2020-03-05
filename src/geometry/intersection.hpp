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

    bool is_point_along_segment(
        const Eigen::VectorX3<Interval>& point,
        const Eigen::VectorX3<Interval>& segment_start,
        const Eigen::VectorX3<Interval>& segment_end);

    bool are_segments_intersecting(
        const Eigen::Vector3<Interval>& segment0_start,
        const Eigen::Vector3<Interval>& segment0_end,
        const Eigen::Vector3<Interval>& segment1_start,
        const Eigen::Vector3<Interval>& segment1_end);

    bool is_point_inside_triangle(
        const Eigen::Vector3<Interval>& point,
        const Eigen::Vector3<Interval>& triangle_vertex0,
        const Eigen::Vector3<Interval>& triangle_vertex1,
        const Eigen::Vector3<Interval>& triangle_vertex2);

} // namespace geometry
} // namespace ccd
