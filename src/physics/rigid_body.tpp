#pragma once
#include "rigid_body.hpp"

#include <Eigen/Geometry>

#include <autodiff/autodiff_types.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace physics {

    template <typename T>
    Eigen::MatrixX<T> RigidBody::world_vertices(
        const Eigen::MatrixXX3<T>& R, const Eigen::VectorX3<T>& p) const
    {
        return (vertices * R.transpose()).rowwise() + p.transpose();
    }

    template <typename T>
    Eigen::VectorX3<T> RigidBody::world_vertex(
        const Eigen::MatrixXX3<T>& R,
        const Eigen::VectorX3<T>& p,
        const int vertex_idx) const
    {
        // compute X[i] = R(θ) * rᵢ + X
        return (vertices.row(vertex_idx) * R.transpose()) + p.transpose();
    }

    template <typename DScalar>
    Eigen::MatrixXd RigidBody::world_vertices_diff(
        const Pose<double>& pose,
        long rb_v0_i,
        Eigen::MatrixXd& V,
        Eigen::MatrixXd& jac,
        std::vector<Eigen::MatrixXd>& hess) const
    {
        // We will only use auto diff to compute the derivatives of the rotation
        // matrix.
        typedef AutodiffType<Eigen::Dynamic, /*maxN=*/3> Diff;
        // Activate autodiff with the correct number of variables.
        Diff::activate(rot_ndof());

        assert(rb_v0_i >= 0 && rb_v0_i <= V.rows() - vertices.rows());
        assert(V.cols() == dim());
        assert(rb_v0_i <= jac.rows() - vertices.size());
        assert(jac.cols() == ndof());
        bool compute_hess = hess.size() >= vertices.size()
            && std::is_base_of<Diff::DDouble2, DScalar>();
        assert(!compute_hess || rb_v0_i <= hess.size() - vertices.size());

        auto R = construct_rotation_matrix(
            Eigen::VectorX3<DScalar>(Diff::dTvars<DScalar>(0, pose.rotation)));
        Eigen::MatrixX<DScalar> V_diff = vertices * R.transpose();

        for (int i = 0; i < V_diff.rows(); i++) {
            for (int j = 0; j < V_diff.cols(); j++) {
                V(rb_v0_i + i, j) = get_value(V_diff(i, j)) + pose.position(j);

                // Fill in gradient of V(i, j) (∈ R⁶ for 3D)
                int vij_flat = (rb_v0_i + i) * V.cols() + j;
                jac(vij_flat, j) = 1; // ∇p V = I
                jac.row(vij_flat).tail(rot_ndof()) = get_gradient(V_diff(i, j));

                if (compute_hess) {
                    // Fill in hessian of V(i, j) (∈ R⁶ˣ⁶ for 3D)
                    // Hessian of position is zero
                    // ∇²_p V = ∇_p∇_r V = ∇_r∇_p V = 0
                    assert(hess[vij_flat].rows() == ndof());
                    assert(hess[vij_flat].cols() == ndof());
                    hess[vij_flat].bottomRightCorner(rot_ndof(), rot_ndof()) =
                        get_hessian(V_diff(i, j));
                }
            }
        }

        return V;
    }

} // namespace physics
} // namespace ccd
