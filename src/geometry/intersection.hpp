#pragma once

#include <Eigen/Core>

#include <interval/interval.hpp>

/// NOTE: Naming Convention
/// Point: Either a 2D or 3D point in space.
/// Edge: A line segment in either 2D or 3D defined by its endpoints.
/// Triangle: A triangle in 3D

namespace ipc::rigid {

bool is_point_along_edge(
    const VectorMax3I& p, const VectorMax3I& e0, const VectorMax3I& e1);

bool are_edges_intersecting(
    const Vector3I& ea0,
    const Vector3I& ea1,
    const Vector3I& eb0,
    const Vector3I& eb1);

bool is_point_inside_triangle(
    const Vector3I& p,
    const Vector3I& t0,
    const Vector3I& t1,
    const Vector3I& t2);

} // namespace ipc::rigid
