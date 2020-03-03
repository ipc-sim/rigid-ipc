#pragma once
#include "distance.hpp"

#include <Eigen/Geometry>
#include <constants.hpp>
#include <geometry/normal.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace geometry {

    //-------------------------------------------------------------------------
    // Unsigned Distances
    //-------------------------------------------------------------------------

    template <typename T>
    inline T point_segment_distance_2D(
        const Eigen::Vector2<T>& point,
        const Eigen::Vector2<T>& segment_start,
        const Eigen::Vector2<T>& segment_end)
    {
        assert(point.size() == 2);
        assert(segment_start.size() == 2);
        assert(segment_end.size() == 2);

        Eigen::Vector2<T> ab = segment_end - segment_start;
        Eigen::Vector2<T> ac = point - segment_start;
        Eigen::Vector2<T> bc = point - segment_end;

        // Handle cases where point projects outside ab
        if (ac.dot(ab) <= T(0.0)) {
            return ac.norm();
        }
        if (bc.dot(ab) >= T(0.0)) {
            return bc.norm();
        }
        Eigen::Vector2<T> abperp =
            segment_normal(segment_start, segment_end, /*normalized=*/false);

        T g = abperp.dot(ac);
        if (g < 0) {
            g *= -1; // Avoid using abs for autodiff
        }
        T e = ab.norm();
        return g / e;
    }

    template <typename T>
    inline T point_segment_distance_3D(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& segment_start,
        const Eigen::Vector3<T>& segment_end)
    {
        throw NotImplementedError(
            "point_segment_distance() not implemented in 3D!");
    }

    template <typename T>
    inline T segment_segment_distance(
        const Eigen::Vector3<T>& segment0_start,
        const Eigen::Vector3<T>& segment0_end,
        const Eigen::Vector3<T>& segment1_start,
        const Eigen::Vector3<T>& segment1_end)
    {
        throw NotImplementedError(
            "segment_segment_distance() not implmeneted!");
    }

    template <typename T>
    inline T point_triangle_distance(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& triangle_vertex0,
        const Eigen::Vector3<T>& triangle_vertex1,
        const Eigen::Vector3<T>& triangle_vertex2)
    {
        throw NotImplementedError("point_triangle_distance() not implmeneted!");
    }

    //-------------------------------------------------------------------------
    // Signed Distances
    //-------------------------------------------------------------------------

    /// Compute the signed distance between a point and a line
    template <typename T>
    inline T point_line_signed_distance(
        const Eigen::Vector2<T>& point,
        const Eigen::Vector2<T>& line_point1,
        const Eigen::Vector2<T>& line_point0)
    {
        Eigen::Vector2<T> normal =
            segment_normal(line_point0, line_point1, /*normalized=*/false);
        return (point - line_point0).dot(normal);
    }

    // Compute the signed distance between two lines
    template <typename T>
    inline T line_line_signed_distance(
        const Eigen::Vector3<T>& line0_point0,
        const Eigen::Vector3<T>& line0_point1,
        const Eigen::Vector3<T>& line1_point0,
        const Eigen::Vector3<T>& line1_point1)
    {
        Eigen::Vector3<T> line0_vec = line0_point1 - line0_point0;
        Eigen::Vector3<T> line1_vec = line1_point1 - line1_point0;

        Eigen::Vector3<T> normal = (line0_vec).cross(line1_vec);

        // WARNING: Parallel edges will result in a distance of 0.
        // WARNING: This is not Euclidean distance
        return (line0_point0 - line1_point0).dot(normal);
    }

    /// Compute the distance between a point and a plane.
    /// Normal is assumed to be unit length.
    template <typename T>
    inline T point_plane_signed_distance(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& triangle_vertex0,
        const Eigen::Vector3<T>& triangle_vertex1,
        const Eigen::Vector3<T>& triangle_vertex2)
    {
        Eigen::Vector3<T> normal = geometry::triangle_normal(
            triangle_vertex0, triangle_vertex1, triangle_vertex2,
            /*normalized=*/false);
        return (point - triangle_vertex0).dot(normal);
    }

} // namespace geometry
} // namespace ccd
