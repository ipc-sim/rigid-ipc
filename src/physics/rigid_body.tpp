#pragma once
#include "rigid_body.hpp"

#include <Eigen/Geometry>

#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace physics {

    template <typename T>
    Eigen::MatrixX<T> RigidBody::world_vertices(
        const Eigen::MatrixXX3<T>& R, const Eigen::VectorX3<T>& p) const
    {
        return (vertices.cast<T>() * R.transpose()).rowwise() + p.transpose();
    }

    template <typename T>
    Eigen::VectorX3<T> RigidBody::world_vertex(
        const Eigen::MatrixXX3<T>& R,
        const Eigen::VectorX3<T>& p,
        const int vertex_idx) const
    {
        // compute X[i] = R(θ) * rᵢ + X
        return (vertices.row(vertex_idx).cast<T>() * R.transpose())
            + p.transpose();
    }

} // namespace physics
} // namespace ccd
