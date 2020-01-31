#pragma once
#include "rigid_body.hpp"

#include <Eigen/Geometry>

#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace physics {

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> RigidBody::world_vertices(
        const Pose<T>& _pose) const
    {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> R
            = _pose.construct_rotation_matrix();
        return (vertices.cast<T>() * R.transpose()).rowwise()
            + _pose.position.transpose();
    }

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, 1> RigidBody::world_vertex(
        const Pose<T>& _pose, const int vertex_idx) const
    {
        // compute X[i] = R(theta) * r_i + X
        if (_pose.dim() == 2) {
            Eigen::Matrix<T, 2, 2> R = _pose.construct_rotation_matrix();
            return (vertices.row(vertex_idx).cast<T>() * R.transpose())
                + _pose.position.transpose();
        } else if (_pose.dim() == 3) {
            Eigen::Matrix<T, 3, 3> R = _pose.construct_rotation_matrix();
            return (vertices.row(vertex_idx).cast<T>() * R.transpose())
                + _pose.position.transpose();
        }
        throw NotImplementedError(
            "RigidBody::world_vertex<T>() does not work for ND!");
    }

} // namespace physics
} // namespace ccd
