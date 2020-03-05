#include "intersection.hpp"

#include <ccd/interval.hpp>
#include <geometry/distance.hpp>
#include <geometry/normal.hpp>
#include <geometry/projection.hpp>
#include <utils/is_zero.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace geometry {

    bool is_point_along_segment(
        const Eigen::VectorX3<Interval>& point,
        const Eigen::VectorX3<Interval>& segment_start,
        const Eigen::VectorX3<Interval>& segment_end)
    {
        Eigen::VectorX3<Interval> segment_dir = segment_end - segment_start;
        Interval alpha =
            project_point_to_line(point, segment_start, segment_dir);
        // Check this in case empty intervals are not allowed
        if (!overlap(alpha, Interval(0, 1))) {
            return false;
        }
        // Check the distance to the closest point is small
        Interval valid_alpha = boost::numeric::intersect(alpha, Interval(0, 1));
        Eigen::VectorX3<Interval> segment_to_point =
            segment_start + valid_alpha * segment_dir - point;

        return is_zero(segment_to_point);
    }

    bool are_segments_intersecting(
        const Eigen::Vector3<Interval>& segment0_start,
        const Eigen::Vector3<Interval>& segment0_end,
        const Eigen::Vector3<Interval>& segment1_start,
        const Eigen::Vector3<Interval>& segment1_end)
    {
        // TODO: This can be made more efficient
        Interval distance = segment_segment_distance(
            segment0_start, segment0_end, segment1_start, segment1_end);
        return is_zero(distance);
    }

    bool is_point_inside_triangle(
        const Eigen::Vector3<Interval>& point,
        const Eigen::Vector3<Interval>& triangle_vertex0,
        const Eigen::Vector3<Interval>& triangle_vertex1,
        const Eigen::Vector3<Interval>& triangle_vertex2)
    {
        Eigen::Vector3<Interval> normal0 =
            triangle_normal(triangle_vertex0, triangle_vertex1, point);
        Eigen::Vector3<Interval> normal1 =
            triangle_normal(triangle_vertex0, point, triangle_vertex2);
        Eigen::Vector3<Interval> normal2 =
            triangle_normal(triangle_vertex1, triangle_vertex2, point);
        for (int i = 0; i < normal0.size(); i++) {
            if (!overlap(intersect(normal0(i), normal1(i)), normal2(i))) {
                return false;
            }
        }
        return true;
    }

} // namespace geometry
} // namespace ccd
