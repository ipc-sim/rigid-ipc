#pragma once
#include "distance.hpp"

namespace ipc::rigid {

//-------------------------------------------------------------------------
// Signed Distances
//-------------------------------------------------------------------------

/// Compute the signed distance between a point and a line
template <typename T>
inline T point_line_signed_distance(
    const Vector2<T>& point,
    const Vector2<T>& line_point0,
    const Vector2<T>& line_point1)
{
    Vector2<T> line_dir = line_point1 - line_point0;
    Vector2<T> normal(-line_dir.y(), line_dir.x());
    return (point - line_point0).dot(normal);
}

// Compute the signed distance between two lines
template <typename T>
inline T line_line_signed_distance(
    const Vector3<T>& line0_point0,
    const Vector3<T>& line0_point1,
    const Vector3<T>& line1_point0,
    const Vector3<T>& line1_point1)
{
    Vector3<T> line0_dir = line0_point1 - line0_point0;
    Vector3<T> line1_dir = line1_point1 - line1_point0;

    Vector3<T> normal = (line0_dir).cross(line1_dir);

    // WARNING: Parallel edges will result in a distance of 0.
    // WARNING: This is not Euclidean distance
    return (line0_point0 - line1_point0).dot(normal);
}

/// Compute the distance between a point and a plane.
template <typename T>
inline T point_plane_signed_distance(
    const Vector3<T>& point,
    const Vector3<T>& triangle_vertex0,
    const Vector3<T>& triangle_vertex1,
    const Vector3<T>& triangle_vertex2)
{
    Vector3<T> normal = (triangle_vertex1 - triangle_vertex0)
                            .cross(triangle_vertex2 - triangle_vertex0);
    return (point - triangle_vertex0).dot(normal);
}

} // namespace ipc::rigid
