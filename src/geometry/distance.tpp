#pragma once
#include "distance.hpp"

#include <Eigen/Geometry>
#include <constants.hpp>
#include <geometry/barycentric_coordinates.hpp>
#include <geometry/intersection.hpp>
#include <geometry/normal.hpp>
#include <geometry/projection.hpp>
#include <logger.hpp>
#include <utils/clamp.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace geometry {

    //-------------------------------------------------------------------------
    // Unsigned Distances
    //-------------------------------------------------------------------------

    template <typename T, int dim, int max_dim>
    inline T point_point_distance(
        const Eigen::Vector<T, dim, max_dim>& point0,
        const Eigen::Vector<T, dim, max_dim>& point1)
    {
#ifdef USE_DISTANCE_SQUARED
        return (point1 - point0).squaredNorm();
#else
        return (point1 - point0).norm();
#endif
    }

    // Find the closest point on the segment to the point.
    template <typename T, int dim, int max_dim>
    Eigen::Vector<T, dim, max_dim> point_segment_closest_point(
        const Eigen::Vector<T, dim, max_dim>& point,
        const Eigen::Vector<T, dim, max_dim>& segment_start,
        const Eigen::Vector<T, dim, max_dim>& segment_end)
    {
        Eigen::Vector<T, dim, max_dim> segment_dir =
            segment_end - segment_start;
        T alpha = clamp_to_01(
            project_point_to_line(point, segment_start, segment_dir));
        return segment_dir * alpha + segment_start;
    }

    template <typename T, int dim, int max_dim>
    inline T point_segment_distance(
        const Eigen::Vector<T, dim, max_dim>& point,
        const Eigen::Vector<T, dim, max_dim>& segment_start,
        const Eigen::Vector<T, dim, max_dim>& segment_end)
    {
        return point_point_distance(
            point,
            point_segment_closest_point(point, segment_start, segment_end));
    }

    template <typename T>
    inline T segment_segment_distance(
        const Eigen::Vector3<T>& segment0_start,
        const Eigen::Vector3<T>& segment0_end,
        const Eigen::Vector3<T>& segment1_start,
        const Eigen::Vector3<T>& segment1_end)
    {
        T alpha0, alpha1;
        project_segment_to_segment(
            segment0_start, segment0_end, segment1_start, segment1_end, alpha0,
            alpha1);

        // Compute the closesty point from seg0 to seg 1
        Eigen::Vector3<T> seg0_to_seg1 =
            (segment0_end - segment0_start) * alpha0 + segment0_start;
        // Compute the closesty point from seg1 to seg 0
        Eigen::Vector3<T> seg1_to_seg0 =
            (segment1_end - segment1_start) * alpha1 + segment1_start;

        // Return the distance between the two closest points
        return point_point_distance(seg0_to_seg1, seg1_to_seg0);
    }

    template <typename T>
    inline T point_triangle_distance(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& triangle_vertex0,
        const Eigen::Vector3<T>& triangle_vertex1,
        const Eigen::Vector3<T>& triangle_vertex2)
    {
        // Compute the barycentric coordinates of the projected_point
        const Eigen::Vector3<T> normal = triangle_normal(
            triangle_vertex0, triangle_vertex1, triangle_vertex2,
            /*normalized=*/false);
        Eigen::Vector3<T> projected_point =
            project_point_to_plane(point, triangle_vertex0, normal);
        T u, v, w;
        barycentric_coordinates(
            point, triangle_vertex0, triangle_vertex1, triangle_vertex2, //
            u, v, w);

        // Find the closest point using the barycentric coordinates
        // https://math.stackexchange.com/a/589362

        // Is closest point in the plane inside the trianlge?
        if (u >= 0 && v >= 0 && w >= 0) {
            return point_point_distance(point, projected_point);
        }

        // Check if a vertex is the closest point on the triangle
        if (u >= 0 && v < 0 && w < 0) {
            // vertex 0 is the closest
            return point_point_distance(point, triangle_vertex0);
        }
        if (u < 0 && v >= 0 && w < 0) {
            // vertex 0 is the closest
            return point_point_distance(point, triangle_vertex1);
        }
        if (u < 0 && v < 0 && w >= 0) {
            // vertex 0 is the closest
            return point_point_distance(point, triangle_vertex2);
        }

        // Check if an edge is the closest point on the triangle
        if (u >= 0 && v >= 0 && w < 0) {
            // vertex 0 is the closest
            return point_segment_distance(
                point, triangle_vertex0, triangle_vertex1);
        }
        if (u >= 0 && v < 0 && w >= 0) {
            // vertex 0 is the closest
            return point_segment_distance(
                point, triangle_vertex2, triangle_vertex0);
        }
        if (u < 0 && v >= 0 && w >= 0) {
            // vertex 0 is the closest
            return point_segment_distance(
                point, triangle_vertex1, triangle_vertex2);
        }

        // This should never happen.because u + v + w = 1.
        throw NotImplementedError(
            "point_triangle_distance is not implemented correctly!");
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
        Eigen::Vector3<T> line0_dir = line0_point1 - line0_point0;
        Eigen::Vector3<T> line1_dir = line1_point1 - line1_point0;

        Eigen::Vector3<T> normal = (line0_dir).cross(line1_dir);

        // WARNING: Parallel edges will result in a distance of 0.
        // WARNING: This is not Euclidean distance
        return (line1_point0 - line0_point0).dot(normal);
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
