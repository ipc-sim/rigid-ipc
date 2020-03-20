#pragma once
#include "projection.tpp"

#include <logger.hpp>
#include <utils/clamp.hpp>
#include <utils/is_zero.hpp>
#include <utils/not_implemented_error.hpp>

/// NOTE: Naming Convention
/// Point: Either a 2D or 3D point in space.
/// Segment: A line segment in either 2D or 3D defined by its endpoints.
/// Triangle: A triangle in 3D

namespace ccd {
namespace geometry {

    template <typename T, int dim, int max_dim>
    inline T project_point_to_line(
        const Eigen::Vector<T, dim, max_dim>& point,
        const Eigen::Vector<T, dim, max_dim>& line_point,
        const Eigen::Vector<T, dim, max_dim>& line_dir)
    {
        // https://zalo.github.io/blog/closest-point-between-segments/
        T dir_length_sqr = line_dir.squaredNorm();
        if (is_zero(dir_length_sqr)) {
            // Segment is degenerate so return a point
            // Halfway so the gradient of an endpoints does not become zero.
            return T(0.5);
        }
        return (point - line_point).dot(line_dir) / dir_length_sqr;
    }

    template <typename T, int dim, int max_dim>
    inline T project_point_to_segment(
        const Eigen::Vector<T, dim, max_dim>& point,
        const Eigen::Vector<T, dim, max_dim>& segment_start,
        const Eigen::Vector<T, dim, max_dim>& segment_dir)
    {
        return clamp_to_01(
            project_point_to_line(point, segment_start, segment_dir));
    }

    template <typename T>
    inline void project_segment_to_segment(
        const Eigen::Vector3<T>& segment0_start,
        const Eigen::Vector3<T>& segment0_end,
        const Eigen::Vector3<T>& segment1_start,
        const Eigen::Vector3<T>& segment1_end,
        T& alpha0,
        T& alpha1)
    {
        // https://zalo.github.io/blog/closest-point-between-segments/
        Eigen::Vector3<T> s0_dir = segment0_end - segment0_start;
        Eigen::Vector3<T> s1_dir = segment1_end - segment1_start;

        // Project the points of segment 0 to the plane orthogonal to segment 1
        T s1_len_sqr = s1_dir.squaredNorm();
        if (is_zero(s1_len_sqr)) {
            // Halfway so the gradient of an endpoints does not become zero.
            alpha1 = T(0.5);
            // The second segment is degenerate
            Eigen::Vector3<T> s1_point = segment1_start + alpha1 * s1_dir;
            alpha0 = project_point_to_segment(s1_point, segment0_start, s0_dir);
            return;
        }

        Eigen::Vector3<T> s00_plane = project_point_to_plane(
            segment0_start, segment1_start, s1_dir, s1_len_sqr);
        Eigen::Vector3<T> s01_plane = project_point_to_plane(
            segment0_end, segment1_start, s1_dir, s1_len_sqr);

        // Compute the point-segment projection in the plane
        Eigen::Vector3<T> s0_plane_dir = s01_plane - s00_plane;
        alpha0 =
            project_point_to_segment(segment1_start, s00_plane, s0_plane_dir);

        // This is the closest point between segment 0 and the line through
        // segment 1.
        Eigen::Vector3<T> seg0_to_line1 = segment0_start + alpha0 * s0_dir;

        // Compute the closesty point from seg1 to seg 0
        alpha1 =
            project_point_to_segment(seg0_to_line1, segment1_start, s1_dir);
        Eigen::Vector3<T> seg1_to_seg0 = s1_dir * alpha1 + segment1_start;
        // Compute the closest point from seg0 to seg 1
        alpha0 = project_point_to_segment(seg1_to_seg0, segment0_start, s0_dir);
    }

    template <typename T>
    inline Eigen::Vector3<T> project_point_to_plane(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& plane_point,
        const Eigen::Vector3<T>& normal)
    {
        return project_point_to_plane(
            point, plane_point, normal, normal.squaredNorm());
    }

    template <typename T>
    inline Eigen::Vector3<T> project_point_to_plane(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& plane_point,
        const Eigen::Vector3<T>& normal,
        const T& normal_sqrnorm)
    {
        return point
            - (point - plane_point).dot(normal) / normal_sqrnorm * normal;
    }

} // namespace geometry
} // namespace ccd
