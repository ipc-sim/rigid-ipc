#pragma once
#include "distance.hpp"

#include <Eigen/Geometry>
#include <constants.hpp>
#include <geometry/barycentric_coordinates.hpp>
#include <geometry/intersection.hpp>
#include <geometry/normal.hpp>
#include <logger.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace geometry {

    //-------------------------------------------------------------------------
    // Unsigned Distances
    //-------------------------------------------------------------------------

    template <typename T, int dim, int max_dim>
    inline T point_point_distance(
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& point0,
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& point1)
    {
        return (point1 - point0).norm();
    }

    template <typename T> inline T clamp_to_01(const T& x)
    {
        return x > T(1.0) ? T(1.0) : (x < T(0.0) ? T(0.0) : x);
    }

    // Find the closest point on the segment to the point.
    template <typename T, int dim, int max_dim>
    Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>
    point_segment_closest_point(
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& point,
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& segment_start,
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& segment_end)
    {
        T alpha = clamp_to_01(
            point_segment_intersection(point, segment_start, segment_end));
        return (segment_end - segment_start) * alpha + segment_start;
    }

    template <typename T, int dim, int max_dim>
    inline T point_segment_distance(
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& point,
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& segment_start,
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& segment_end)
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
        // https://zalo.github.io/blog/closest-point-between-segments/
        // Project the points of segment 0 to the plane orthogonal to segment 1
        Eigen::Vector3<T> normal = segment1_end - segment1_start;
        T normal_sqrnorm = normal.squaredNorm();
        if (normal_sqrnorm == 0) {
            // The second segment is degenerate
            return point_segment_distance(
                segment1_start, segment0_start, segment0_end);
        }

        auto project_to_plane = [&](const Eigen::Vector3<T>& p) {
            T alpha = (p - segment1_start).dot(normal) / normal_sqrnorm;
            return p - (alpha * normal);
        };

        Eigen::Vector3<T> s00_plane = project_to_plane(segment0_start);
        Eigen::Vector3<T> s01_plane = project_to_plane(segment0_end);

        // Compute the segment 0 direction in the plane
        Eigen::Vector3<T> s0_plane_dir = s01_plane - s00_plane;
        T s0_plane_len_sqr = s0_plane_dir.squaredNorm();

        // Compute the point-segment distance in the plane
        T alpha = (s0_plane_len_sqr == 0.0)
            ? T(0.0) // Zero if parallel or first segment is degenerate
            : (segment1_start - s00_plane).dot(s0_plane_dir) / s0_plane_len_sqr;
        // Clamp the parameter to the endpoints of segment 0
        alpha = clamp_to_01(alpha);

        // This is the closest point between segment 0 and the line through
        // segment 1.
        Eigen::Vector3<T> seg0_to_line1 =
            (segment0_end - segment0_start) * alpha + segment0_start;

        // Compute the closesty point from seg1 to seg 0
        Eigen::Vector3<T> seg1_to_seg0 = point_segment_closest_point(
            seg0_to_line1, segment1_start, segment1_end);
        // Compute the closesty point from seg0 to seg 1
        Eigen::Vector3<T> seg0_to_seg1 = point_segment_closest_point(
            seg1_to_seg0, segment0_start, segment0_end);

        // Return the distance between the two closest points
        return point_point_distance(seg1_to_seg0, seg0_to_seg1);
    }

    template <typename T>
    inline T point_triangle_distance(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& triangle_vertex0,
        const Eigen::Vector3<T>& triangle_vertex1,
        const Eigen::Vector3<T>& triangle_vertex2)
    {
        // Project the point to the plane without the need for a sqrt
        const Eigen::Vector3<T> normal = triangle_normal(
            triangle_vertex0, triangle_vertex1, triangle_vertex2,
            /*normalized=*/false);
        Eigen::Vector3<T> projected_point = point
            - ((point - triangle_vertex0).dot(normal) / normal.squaredNorm())
                * normal;

        // Compute the barycentric coordinates of the projected_point
        T u, v, w;
        barycentric_coordinates(
            projected_point, triangle_vertex0, triangle_vertex1,
            triangle_vertex2, u, v, w);

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
        if (u < 0 && v <= 0 && w < 0) {
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
