#pragma once

#include <Eigen/Core>

namespace ipc::rigid {

//-------------------------------------------------------------------------
// Signed Distances (useful for CCD)
//-------------------------------------------------------------------------

/// Compute the signed distance between a point and a line
/// WARNING: Produces the same sign as euclidean distance but may be
/// different scales.
template <typename T>
inline T point_line_signed_distance(
    const Vector2<T>& point,
    const Vector2<T>& line_point0,
    const Vector2<T>& line_point1);

/// Compute the signed distance between two lines
/// WARNING: Produces the same sign as euclidean distance but may be
/// different scales.
/// WARNING: Parallel edges results in zero distance
template <typename T>
inline T line_line_signed_distance(
    const Vector3<T>& line0_point0,
    const Vector3<T>& line0_point1,
    const Vector3<T>& line1_point0,
    const Vector3<T>& line1_point1);

/// Compute the distance between a point and a plane.
/// Normal is assumed to be unit length.
/// WARNING: Produces the same sign as euclidean distance but may be
/// different scales.
template <typename T>
inline T point_plane_signed_distance(
    const Vector3<T>& point,
    const Vector3<T>& triangle_vertex0,
    const Vector3<T>& triangle_vertex1,
    const Vector3<T>& triangle_vertex2);

} // namespace ipc::rigid

#include "distance.tpp"
