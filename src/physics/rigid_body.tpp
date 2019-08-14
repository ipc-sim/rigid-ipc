#pragma once

#include <Eigen/Geometry>

namespace ccd {
namespace physics {

    template <typename T>
    Eigen::Matrix<T, 2, 1> RigidBody::world_vertex(
        const Eigen::Matrix<T, 3, 1>& _position, const int vertex_idx) const
    {
        typedef Eigen::Matrix<T, 2, 2> Matrix2T;

        // compute X[i] = R(theta) * r_i + X
        Matrix2T R = Eigen::Rotation2D<T>(_position.z()).toRotationMatrix();
        return (vertices.row(vertex_idx).cast<T>() * R.transpose())
            + _position.head(2).transpose();
    }

} // namespace opt
} // namespace ccd
