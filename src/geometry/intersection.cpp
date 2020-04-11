#include "intersection.hpp"

#include <igl/predicates/predicates.h>

#include <ccd/interval.hpp>
#include <geometry/distance.hpp>
#include <geometry/normal.hpp>
#include <geometry/projection.hpp>
#include <utils/is_zero.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace geometry {

    bool is_point_along_segment(
        const Eigen::VectorX3I& point,
        const Eigen::VectorX3I& segment_start,
        const Eigen::VectorX3I& segment_end)
    {
        Eigen::VectorX3I segment_dir = segment_end - segment_start;
        Interval alpha =
            project_point_to_line(point, segment_start, segment_dir);
        // Check this in case empty intervals are not allowed
        if (!overlap(alpha, Interval(0, 1))) {
            return false;
        }
        // Check the distance to the closest point is small
        Interval valid_alpha = boost::numeric::intersect(alpha, Interval(0, 1));
        Eigen::VectorX3I segment_to_point =
            segment_start + valid_alpha * segment_dir - point;

        return is_zero(segment_to_point);
    }

    bool are_segments_intersecting(
        const Eigen::Vector3I& segment0_start,
        const Eigen::Vector3I& segment0_end,
        const Eigen::Vector3I& segment1_start,
        const Eigen::Vector3I& segment1_end)
    {
        // TODO: This can be made more efficient
        Interval distance = segment_segment_distance(
            segment0_start, segment0_end, segment1_start, segment1_end);
        return is_zero(distance);
    }

    bool is_point_inside_triangle(
        const Eigen::Vector3I& point,
        const Eigen::Vector3I& triangle_vertex0,
        const Eigen::Vector3I& triangle_vertex1,
        const Eigen::Vector3I& triangle_vertex2)
    {
        Eigen::Vector3I normal0 =
            triangle_normal(triangle_vertex0, triangle_vertex1, point);
        Eigen::Vector3I normal1 =
            triangle_normal(triangle_vertex0, point, triangle_vertex2);
        Eigen::Vector3I normal2 =
            triangle_normal(triangle_vertex1, triangle_vertex2, point);
        for (int i = 0; i < normal0.size(); i++) {
            if (!overlap(intersect(normal0(i), normal1(i)), normal2(i))) {
                return false;
            }
        }
        return true;
    }

    bool segment_triangle_intersect(
        const Eigen::Vector3d& segment_vertex0,
        const Eigen::Vector3d& segment_vertex1,
        const Eigen::Vector3d& triangle_vertex0,
        const Eigen::Vector3d& triangle_vertex1,
        const Eigen::Vector3d& triangle_vertex2)
    {
        igl::predicates::exactinit();
        const auto ori1 = igl::predicates::orient3d(
            triangle_vertex0, triangle_vertex1, triangle_vertex2,
            segment_vertex0);
        const auto ori2 = igl::predicates::orient3d(
            triangle_vertex0, triangle_vertex1, triangle_vertex2,
            segment_vertex1);

        if (ori1 == igl::predicates::Orientation::COPLANAR
            || ori2 == igl::predicates::Orientation::COPLANAR) {
            // coplanar, we can detect it by d(EE)=0 or d(PT)=0
            return false;
        }

        if (ori1 == ori2) {
            // edge is on one side of the plane that triangle is in
            return false;
        }

        int res = eccd::segment_triangle_inter(
            segment_vertex0, segment_vertex1, //
            triangle_vertex0, triangle_vertex1, triangle_vertex2);
        return res == 1;
    }

} // namespace geometry
} // namespace ccd
