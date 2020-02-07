#pragma once
#include "distance.hpp"

#include <utils/not_implemented_error.hpp>

/// NOTE: Naming Convention
/// Point: Either a 2D or 3D point in space.
/// Segment: A line segment in either 2D or 3D defined by its endpoints.
///     Segment Vertex 0: The first endpoint of a segment (α = 0).
///     Segment Vertex 1: The second endpoint of a segment (α = 1).
/// Triangle: A triangle in 3D

namespace ccd {

//-----------------------------------------------------------------------------
// Unsigned Distances
//-----------------------------------------------------------------------------

template <typename T>
T point_point_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& point0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& point1)
{
    return (point1 - point0).norm();
}

template <typename T>
T point_segment_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& point,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment_vertex0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment_vertex1)
{
    return point_point_distance(
        point,
        point_segment_closest_point(point, segment_vertex0, segment_vertex1));
}

template <typename T>
T segment_segment_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment0_vertex0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment0_vertex1,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment1_vertex0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment1_vertex1)
{
    Eigen::Matrix<T, Eigen::Dynamic, 1> segment0_point, segment1_point;
    segment_segment_closest_points(
        segment0_vertex0, segment0_vertex1, segment1_vertex0, segment1_vertex1,
        segment0_point, segment1_point);
    return point_point_distance(segment0_point, segment1_point);
}

template <typename T>
T point_triangle_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& point,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& triangle_vertex0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& triangle_vertex1,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& triangle_vertex2)
{
    return point_triangle_closest_point(
        point,
        point_segment_closest_point(
            point, triangle_vertex0, triangle_vertex1, triangle_vertex2));
}

template <typename T>
T point_plane_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& point,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& plane_point,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& plane_normal)
{
    return abs(point_plane_signed_distance(point, plane_point, plane_normal));
}

//-----------------------------------------------------------------------------
// Signed Distances
//-----------------------------------------------------------------------------

template <typename T>
T point_plane_signed_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& p,  // Vertex
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& p0, // Point on the plane.
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& n)  // Unit normal of the plane.
{
    return (p - p0).dot(n);
}

} // namespace ccd
