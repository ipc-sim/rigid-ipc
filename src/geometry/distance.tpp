#pragma once
#include "distance.hpp"

#include <Eigen/Geometry>
#include <constants.hpp>
#include <geometry/barycentric_coordinates.hpp>
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

    template <typename T, int dim, int max_dim>
    inline T point_segment_distance(
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& point,
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& segment_start,
        const Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim>& segment_end)
    {
        // https://zalo.github.io/blog/closest-point-between-segments/
        Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim> segment_vec =
            segment_end - segment_start;
        T segment_length_sqr = segment_vec.squaredNorm();
        // if (segment_length_sqr == T(0.0)) {
        //     return point_point_distance(point, segment_start);
        // }
        T alpha = (point - segment_start).dot(segment_vec) / segment_length_sqr;
        if (alpha > T(1.0)) {
            alpha = T(1.0);
        } else if (alpha < T(0.0)) {
            alpha = T(0.0);
        }
        return point_point_distance(
            point, (segment_start + alpha * segment_vec).eval());
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
        const Eigen::Vector3<T> normal = triangle_normal(
            triangle_vertex0, triangle_vertex1, triangle_vertex2);
        const Eigen::Vector3<T> projected_point =
            point - ((point - triangle_vertex0).dot(normal)) * normal;
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
