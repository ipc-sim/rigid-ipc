#pragma once

#include <Eigen/Core>

#include <interval/interval.hpp>

/// NOTE: Naming Convention
/// Point: Either a 2D or 3D point in space.
/// Edge: A line segment in either 2D or 3D defined by its endpoints.
/// Triangle: A triangle in 3D

namespace ipc::rigid {

bool is_point_along_edge(
    const Eigen::VectorX3I& p,
    const Eigen::VectorX3I& e0,
    const Eigen::VectorX3I& e1);

bool are_edges_intersecting(
    const Eigen::Vector3I& ea0,
    const Eigen::Vector3I& ea1,
    const Eigen::Vector3I& eb0,
    const Eigen::Vector3I& eb1);

bool is_point_inside_triangle(
    const Eigen::Vector3I& p,
    const Eigen::Vector3I& t0,
    const Eigen::Vector3I& t1,
    const Eigen::Vector3I& t2);

bool are_edge_triangle_intersecting(
    const Eigen::Vector3d& e0,
    const Eigen::Vector3d& e1,
    const Eigen::Vector3d& t0,
    const Eigen::Vector3d& t1,
    const Eigen::Vector3d& t2);

} // namespace ipc::rigid
