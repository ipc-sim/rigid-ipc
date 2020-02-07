#pragma once

#include <Eigen/Core>

namespace ccd {

//-----------------------------------------------------------------------------
// Unsigned Distances
//-----------------------------------------------------------------------------

template <typename T>
T point_point_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& point0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& point1);

template <typename T>
T point_segment_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& point,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment_vertex0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment_vertex1);

template <typename T>
T segment_segment_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment0_vertex0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment0_vertex1,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment1_vertex0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& segment1_vertex1);

template <typename T>
T point_triangle_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& point,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& triangle_vertex0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& triangle_vertex1,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& triangle_vertex2);

template <typename T>
T point_plane_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& point,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& plane_point,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& plane_normal);

//-----------------------------------------------------------------------------
// Signed Distances
//-----------------------------------------------------------------------------

/// Compute the distance between a point and a plane.
/// Normal is assumed to be unit length.
template <typename T>
T point_plane_signed_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& point,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& plane_point,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& plane_normal);

} // namespace ccd

#include "distance.tpp"
