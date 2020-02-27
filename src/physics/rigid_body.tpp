#pragma once
#include "rigid_body.hpp"

#include <Eigen/Geometry>

#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace physics {

    template <typename T>
    Eigen::MatrixX<T> RigidBody::world_vertices(const Pose<T>& _pose) const
    {
        Eigen::MatrixXX3<T> R = _pose.construct_rotation_matrix();
        return (vertices.cast<T>() * R.transpose()).rowwise()
            + _pose.position.transpose();
    }

    template <typename T>
    Eigen::VectorX3<T>
    RigidBody::world_vertex(const Pose<T>& _pose, const int vertex_idx) const
    {
        // compute X[i] = R(θ) * rᵢ + X
        Eigen::MatrixXX3<T> R = _pose.construct_rotation_matrix();
        return (vertices.row(vertex_idx).cast<T>() * R.transpose())
            + _pose.position.transpose();
    }

} // namespace physics
} // namespace ccd
