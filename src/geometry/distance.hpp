#pragma once

#include <Eigen/Core>

#include <utils/eigen_ext.hpp>

/// NOTE: Naming Convention
/// Point: Either a 2D or 3D point in space.
/// Segment: A line segment in either 2D or 3D defined by its endpoints.
/// Triangle: A triangle in 3D

// Uncomment this line to switch to distance squared
//#define USE_DISTANCE_SQUARED

namespace ccd {
namespace geometry {

    //-------------------------------------------------------------------------
    // Unsigned Distances
    //-------------------------------------------------------------------------

    template <typename T, int dim, int max_dim = dim>
    inline T point_point_distance(
        const Eigen::Vector<T, dim, max_dim>& point0,
        const Eigen::Vector<T, dim, max_dim>& point1);

    template <typename T, int dim, int max_dim>
    inline T point_segment_distance(
        const Eigen::Vector<T, dim, max_dim>& point,
        const Eigen::Vector<T, dim, max_dim>& segment_start,
        const Eigen::Vector<T, dim, max_dim>& segment_end);

    template <typename T>
    inline T segment_segment_distance(
        const Eigen::Vector3<T>& segment0_start,
        const Eigen::Vector3<T>& segment0_end,
        const Eigen::Vector3<T>& segment1_start,
        const Eigen::Vector3<T>& segment1_end);

    template <typename T>
    inline T point_triangle_distance(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& triangle_vertex0,
        const Eigen::Vector3<T>& triangle_vertex1,
        const Eigen::Vector3<T>& triangle_vertex2);

    //-------------------------------------------------------------------------
    // Signed Distances (useful for CCD)
    //-------------------------------------------------------------------------

    /// Compute the signed distance between a point and a line
    /// WARNING: Produces the same sign as euclidean distance but may be
    /// different scales.
    template <typename T>
    inline T point_line_signed_distance(
        const Eigen::Vector2<T>& point,
        const Eigen::Vector2<T>& line_point0,
        const Eigen::Vector2<T>& line_point1);

    /// Compute the signed distance between two lines
    /// WARNING: Produces the same sign as euclidean distance but may be
    /// different scales.
    /// WARNING: Parallel edges results in zero distance
    template <typename T>
    inline T line_line_signed_distance(
        const Eigen::Vector3<T>& line0_point0,
        const Eigen::Vector3<T>& line0_point1,
        const Eigen::Vector3<T>& line1_point0,
        const Eigen::Vector3<T>& line1_point1);

    /// Compute the distance between a point and a plane.
    /// Normal is assumed to be unit length.
    /// WARNING: Produces the same sign as euclidean distance but may be
    /// different scales.
    template <typename T>
    inline T point_plane_signed_distance(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& triangle_vertex0,
        const Eigen::Vector3<T>& triangle_vertex1,
        const Eigen::Vector3<T>& triangle_vertex2);

} // namespace geometry
} // namespace ccd

#include "distance.tpp"
