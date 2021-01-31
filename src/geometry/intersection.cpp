#include "intersection.hpp"

#include <Eigen/Geometry>

#ifdef RIGID_IPC_USE_EXACT_INTERSECTION
#include <ECCD.hpp>
#endif
#include <igl/predicates/predicates.h>
#include <ipc/friction/closest_point.hpp>

#include <utils/is_zero.hpp>

namespace ipc::rigid {

bool is_point_along_edge(
    const Eigen::VectorX3I& p,
    const Eigen::VectorX3I& e0,
    const Eigen::VectorX3I& e1)
{
    Eigen::VectorX3I e = e1 - e0;

    Interval alpha = point_edge_closest_point(p, e0, e1);
    // Check this in case empty intervals are not allowed
    if (!overlap(alpha, Interval(0, 1))) {
        return false;
    }
    Interval valid_alpha = boost::numeric::intersect(alpha, Interval(0, 1));

    // Check the distance to the closest point is small
    Eigen::VectorX3I edge_to_point = e0 + valid_alpha * e - p;

    return is_zero(edge_to_point);
}

bool are_edges_intersecting(
    const Eigen::Vector3I& ea0,
    const Eigen::Vector3I& ea1,
    const Eigen::Vector3I& eb0,
    const Eigen::Vector3I& eb1)
{
    // Check if the origin is withing the 3D difference of edge points
    // WARNING: This is a very converative estimate.
    Eigen::Vector3I ea_alpha = (ea1 - ea0) * Interval(0, 1) + ea0;
    Eigen::Vector3I eb_alpha = (eb1 - eb0) * Interval(0, 1) + eb0;
    return is_zero(Eigen::Vector3I(ea_alpha - eb_alpha));
}

inline bool are_points_on_same_side_of_edge(
    const Eigen::Vector3I& p1,
    const Eigen::Vector3I& p2,
    const Eigen::Vector3I& a,
    const Eigen::Vector3I& b)
{
    Eigen::Vector3I cp1 = (b - a).cross(p1 - a);
    Eigen::Vector3I cp2 = (b - a).cross(p2 - a);
    return cp1.dot(cp2).upper() >= 0;
}

bool is_point_inside_triangle(
    const Eigen::Vector3I& point,
    const Eigen::Vector3I& triangle_vertex0,
    const Eigen::Vector3I& triangle_vertex1,
    const Eigen::Vector3I& triangle_vertex2)
{
    return are_points_on_same_side_of_edge(
               point, triangle_vertex0, triangle_vertex1, triangle_vertex2)
        && are_points_on_same_side_of_edge(
               point, triangle_vertex1, triangle_vertex0, triangle_vertex2)
        && are_points_on_same_side_of_edge(
               point, triangle_vertex2, triangle_vertex0, triangle_vertex1);
}

bool are_edge_triangle_intersecting(
    const Eigen::Vector3d& e0,
    const Eigen::Vector3d& e1,
    const Eigen::Vector3d& t0,
    const Eigen::Vector3d& t1,
    const Eigen::Vector3d& t2)
{
    igl::predicates::exactinit();
    const auto ori1 = igl::predicates::orient3d(t0, t1, t2, e0);
    const auto ori2 = igl::predicates::orient3d(t0, t1, t2, e1);

    if (ori1 != igl::predicates::Orientation::COPLANAR
        && ori2 != igl::predicates::Orientation::COPLANAR && ori1 == ori2) {
        // edge is on one side of the plane that triangle is in
        return false;
    }

#ifdef RIGID_IPC_USE_EXACT_INTERSECTION
    int res = eccd::segment_triangle_inter(e0, e1, t0, t1, t2);
    return res == 1;
#else
    Eigen::Matrix3d M;
    M.col(0) = t1 - t0;
    M.col(1) = t2 - t0;
    M.col(2) = e0 - e1;
    Eigen::Vector3d uvt = M.fullPivLu().solve(e0 - t0);
    if (uvt[0] >= 0.0 && uvt[1] >= 0.0 && uvt[0] + uvt[1] <= 1.0
        && uvt[2] >= 0.0 && uvt[2] <= 1.0) {
        return true;
    } else {
        return false;
    }
#endif
}

} // namespace ipc::rigid
