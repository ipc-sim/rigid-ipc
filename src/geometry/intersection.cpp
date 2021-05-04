#include "intersection.hpp"

#include <Eigen/Geometry>

#include <igl/predicates/predicates.h>
#include <ipc/friction/closest_point.hpp>

#include <utils/is_zero.hpp>

namespace ipc::rigid {

bool is_point_along_edge(
    const VectorMax3I& p, const VectorMax3I& e0, const VectorMax3I& e1)
{
    VectorMax3I e = e1 - e0;

    Interval alpha = point_edge_closest_point(p, e0, e1);
    // Check this in case empty intervals are not allowed
    if (!overlap(alpha, Interval(0, 1))) {
        return false;
    }
    Interval valid_alpha = boost::numeric::intersect(alpha, Interval(0, 1));

    // Check the distance to the closest point is small
    VectorMax3I edge_to_point = e0 + valid_alpha * e - p;

    return is_zero(edge_to_point);
}

bool are_edges_intersecting(
    const Vector3I& ea0,
    const Vector3I& ea1,
    const Vector3I& eb0,
    const Vector3I& eb1)
{
    // Check if the origin is withing the 3D difference of edge points
    // WARNING: This is a very converative estimate.
    Vector3I ea_alpha = (ea1 - ea0) * Interval(0, 1) + ea0;
    Vector3I eb_alpha = (eb1 - eb0) * Interval(0, 1) + eb0;
    return is_zero(Vector3I(ea_alpha - eb_alpha));
}

inline bool are_points_on_same_side_of_edge(
    const Vector3I& p1,
    const Vector3I& p2,
    const Vector3I& a,
    const Vector3I& b)
{
    Vector3I cp1 = (b - a).cross(p1 - a);
    Vector3I cp2 = (b - a).cross(p2 - a);
    return cp1.dot(cp2).upper() >= 0;
}

bool is_point_inside_triangle(
    const Vector3I& point,
    const Vector3I& triangle_vertex0,
    const Vector3I& triangle_vertex1,
    const Vector3I& triangle_vertex2)
{
    return are_points_on_same_side_of_edge(
               point, triangle_vertex0, triangle_vertex1, triangle_vertex2)
        && are_points_on_same_side_of_edge(
               point, triangle_vertex1, triangle_vertex0, triangle_vertex2)
        && are_points_on_same_side_of_edge(
               point, triangle_vertex2, triangle_vertex0, triangle_vertex1);
}

} // namespace ipc::rigid
