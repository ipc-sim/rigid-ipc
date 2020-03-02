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
    inline T point_point_distance(
        const Eigen::VectorX3<T>& point0, const Eigen::VectorX3<T>& point1)
    {
        return (point1 - point0).norm();
    }

    template <typename T>
    inline T point_segment_distance(
        const Eigen::VectorX3<T>& point,
        const Eigen::VectorX3<T>& segment_start,
        const Eigen::VectorX3<T>& segment_end)
    {
        if (point.size() == 2) {
            return point_segment_distance_2D(point, segment_start, segment_end);
        }
        throw NotImplementedError(
            "point_segment_distance() not implemented in 3D!");
    }

    template <typename T>
    inline T point_segment_distance_2D(
        const Eigen::VectorX3<T>& point,
        const Eigen::VectorX3<T>& segment_start,
        const Eigen::VectorX3<T>& segment_end)
    {
        assert(point.size() == 2);
        assert(segment_start.size() == 2);
        assert(segment_end.size() == 2);

        Eigen::VectorX3<T> ab = segment_end - segment_start;
        Eigen::VectorX3<T> ac = point - segment_start;
        Eigen::VectorX3<T> bc = point - segment_end;

        // Handle cases where point projects outside ab
        if (ac.dot(ab) <= T(0.0)) {
            return ac.norm();
        }
        if (bc.dot(ab) >= T(0.0)) {
            return bc.norm();
        }
        Eigen::VectorX3<T> abperp =
            segment_normal(segment_start, segment_end, /*normalized=*/false);

        T g = abperp.dot(ac);
        if (g < 0) {
            g *= -1; // Avoid using abs for autodiff
        }
        T e = ab.norm();
        return g / e;
    }

    template <typename T>
    inline T point_line_distance(
        const Eigen::VectorX3<T>& point,
        const Eigen::VectorX3<T>& line_point0,
        const Eigen::VectorX3<T>& line_point1)
    {
        switch (point.size()) {
        case 2:
            return point_line_distance_2D(point, line_point0, line_point1);
        case 3:
            return point_line_distance_3D(point, line_point0, line_point1);
        default:
            throw NotImplementedError(
                "point_line_distance() only implemented for 2D and 3D!");
        }
    }

    template <typename T>
    inline T point_line_distance_2D(
        const Eigen::VectorX3<T>& point,
        const Eigen::VectorX3<T>& line_point0,
        const Eigen::VectorX3<T>& line_point1)
    {
        assert(point.size() == 2);
        assert(line_point0.size() == 2);
        assert(line_point1.size() == 2);

        T distance =
            point_line_signed_distance(point, line_point0, line_point1);
        if (distance < 0) {
            distance *= -1; // Avoid using abs for autodiff
        }
        return distance;
    }

    template <typename T>
    inline T point_line_distance_3D(
        const Eigen::VectorX3<T>& point,
        const Eigen::VectorX3<T>& line_point0,
        const Eigen::VectorX3<T>& line_point1)
    {
        assert(point.size() == 3);
        assert(line_point0.size() == 3);
        assert(line_point1.size() == 3);

        // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        T segment_length = (line_point1 - line_point0).norm();
        if (segment_length == 0) {
            return point_point_distance(point, line_point0);
        }
        Eigen::Vector3<T> line_to_point0 = point - line_point0;
        Eigen::Vector3<T> line_to_point1 = point - line_point1;
        return (line_to_point0).cross(line_to_point1).norm() / segment_length;
    }

    template <typename T>
    inline T segment_segment_distance(
        const Eigen::VectorX3<T>& segment0_start,
        const Eigen::VectorX3<T>& segment0_end,
        const Eigen::VectorX3<T>& segment1_start,
        const Eigen::VectorX3<T>& segment1_end)
    {
        throw NotImplementedError(
            "segment_segment_distance() not implmeneted!");
    }

    template <typename T>
    inline T point_triangle_distance(
        const Eigen::VectorX3<T>& point,
        const Eigen::VectorX3<T>& triangle_vertex0,
        const Eigen::VectorX3<T>& triangle_vertex1,
        const Eigen::VectorX3<T>& triangle_vertex2)
    {
        throw NotImplementedError("point_triangle_distance() not implmeneted!");
    }

    //-------------------------------------------------------------------------
    // Signed Distances
    //-------------------------------------------------------------------------

    /// Compute the signed distance between a point and a line
    template <typename T>
    inline T point_line_signed_distance(
        const Eigen::VectorX3<T>& point,
        const Eigen::VectorX3<T>& line_point1,
        const Eigen::VectorX3<T>& line_point0)
    {
        if (point.size() != 2) {
            throw NotImplementedError(
                "point_line_signed_distance() not implmeneted in 3D!");
        }
        Eigen::VectorX3<T> normal = segment_normal(line_point0, line_point1);
        return (point - line_point0).dot(normal);
    }

    // Compute the signed distance between two lines
    template <typename T>
    inline T line_line_signed_distance(
        const Eigen::VectorX3<T>& line0_point0,
        const Eigen::VectorX3<T>& line0_point1,
        const Eigen::VectorX3<T>& line1_point0,
        const Eigen::VectorX3<T>& line1_point1)
    {
        if (line0_point0.size() != 3) {
            throw NotImplementedError(
                "line_line_signed_distance() not implmeneted in 2D!");
        }
        Eigen::Vector3<T> line0_vec = line0_point1 - line0_point0;
        Eigen::Vector3<T> line1_vec = line1_point1 - line1_point0;
        Eigen::Vector3<T> normal = (line0_vec).cross(line1_vec);
        T normal_norm = normal.norm();
        if (normal_norm <= 1e-30) { // parallel lines
            // TODO: Define a signed version of the point line distance
            return point_line_distance(
                line0_point0, line1_point0, line1_point1);
        }
        return (line0_point0 - line1_point0).dot(normal) / normal_norm;
    }

    /// Compute the distance between a point and a plane.
    /// Normal is assumed to be unit length.
    template <typename T>
    inline T point_plane_signed_distance(
        const Eigen::VectorX3<T>& point,
        const Eigen::VectorX3<T>& plane_point,
        const Eigen::VectorX3<T>& plane_normal)
    {
        return (point - plane_point).dot(plane_normal);
    }

} // namespace geometry
} // namespace ccd
