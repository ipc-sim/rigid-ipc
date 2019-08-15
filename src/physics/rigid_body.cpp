﻿#include <physics/rigid_body.hpp>

#include <Eigen/Geometry>

#include <physics/mass_matrix.hpp>
#include <physics/moment_of_inertia.hpp>

#include <autodiff/autodiff_types.hpp>
#include <autodiff/finitediff.hpp>
#include <logger.hpp>
#include <utils/flatten.hpp>

namespace ccd {

namespace physics {

    RigidBody RigidBody::from_velocity(const Eigen::MatrixXd& vertices,
        const Eigen::MatrixX2i& edges,
        const Eigen::Vector3d& position,
        const Eigen::Vector3d& velocity,
        const Eigen::Vector3b& is_dof_fixed)
    {
        // move vertices so they center of mass is at 0,0
        Eigen::MatrixXd vertices_ = vertices;
        vertices_.rowwise() += position.head(2).transpose();

        Eigen::RowVector2d x = center_of_mass(vertices_, edges);
        Eigen::MatrixX2d centered_vertices = vertices_.rowwise() - x;

        // set position so current vertices match input
        Eigen::Vector3d position_t1 = Eigen::Vector3d::Zero();
        position_t1.segment(0, 2) += x;
        position_t1(2) = position(2);

        // set previous_step position to:
        Eigen::Vector3d position_t0 = position_t1;

        return RigidBody(centered_vertices, edges, velocity, is_dof_fixed,
            position_t1, position_t0);
    }

    RigidBody RigidBody::from_displacement(const Eigen::MatrixXd& vertices,
        const Eigen::MatrixX2i& edges,
        const Eigen::Vector3d& displacement,
        const Eigen::Vector3b& is_dof_fixed)
    {
        // move vertices so they center of mass is at 0,0
        Eigen::RowVector2d x = center_of_mass(vertices, edges);
        Eigen::MatrixX2d centered_vertices = vertices.rowwise() - x;

        // set position so current vertices match input
        Eigen::Vector3d position_t0 = Eigen::Vector3d::Zero();
        position_t0.segment(0, 2) = x;

        // set previous_step position to:
        Eigen::Vector3d position_t1 = position_t0 + displacement;

        // set velocity to zero
        Eigen::Vector3d velocity = Eigen::Vector3d::Zero();

        return RigidBody(centered_vertices, edges, velocity, is_dof_fixed,
            position_t1, position_t0);
    }

    RigidBody::RigidBody(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2i& edges,
        const Eigen::Vector3d& v,
        const Eigen::Vector3b& is_dof_fixed,
        const Eigen::Vector3d& x,
        const Eigen::Vector3d& x_prev)
        : vertices(vertices)
        , edges(edges)
        , is_dof_fixed(is_dof_fixed)
        , velocity(v)
        , position(x)
        , position_prev(x_prev)
    {
        Eigen::VectorXd masses;
        mass_vector(vertices, edges, masses);
        mass = masses.sum();
        moment_of_inertia = physics::moment_of_inertia(vertices, masses);

        mass_matrix << mass, 0, 0, 0, mass, 0, 0, 0, moment_of_inertia;
        inv_mass_matrix << 1.0 / mass, 0, 0, 0, 1.0 / mass, 0, 0, 0,
            1.0 / moment_of_inertia;
    }

    Eigen::MatrixXd RigidBody::world_vertices(const Step step) const
    {
        return world_vertices(step == PREVIOUS_STEP ? position_prev : position);
    }

    Eigen::MatrixXd RigidBody::world_vertices(
        const Eigen::Vector3d& _position) const
    {
        // compute X[i] = R(theta) * r_i + X
        Eigen::Matrix2d R
            = Eigen::Rotation2D<double>(_position.z()).toRotationMatrix();
        return  (vertices * R.transpose()).rowwise()
            + _position.head(2).transpose();
    }

    Eigen::MatrixXd RigidBody::world_velocities() const
    {
        // compute X[i] = dR(theta)/d\theta * r_i * d\theta/dt + dX/dt
        Eigen::Matrix2d grad_theta_ = grad_theta(position(2));
        Eigen::MatrixXd grad_theta_r_i = vertices * grad_theta_.transpose();

        return (grad_theta_r_i * velocity(2)).rowwise()
            + velocity.head(2).transpose();
    }

    Eigen::Matrix2d RigidBody::grad_theta(const double theta) const
    {
        Eigen::Matrix2d gradtheta_R;
        gradtheta_R << -sin(theta), -cos(theta), cos(theta), -sin(theta);
        return gradtheta_R;
    }

    Eigen::MatrixXd RigidBody::world_vertices_gradient_exact(
        const Eigen::Vector3d& _position) const
    {
        /// For a single vertice ri, gradient is
        /// [[1,0,−ri_x * sin(\theta) − ri_v * cos(\theta)],
        ///  [0,1, ri_x * cos(\theta) − ri_v * sin(\theta)]]
        ///
        ///
        /// gradient has shape vertices.size() by 3
        /// order of rows is x-positions, y-positions

        // gradient of r wrt position.x()
        Eigen::MatrixXd gradx_U
            = Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols());
        gradx_U.col(0).setOnes();
        flatten<double>(gradx_U);

        // gradient of r wrt position.y()
        Eigen::MatrixXd grady_U
            = Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols());
        grady_U.col(1).setOnes();
        flatten<double>(grady_U);

        // gradient of r wrt \theta
        double theta = _position(2);
        Eigen::Matrix2d gradtheta_R = grad_theta(theta);
        Eigen::MatrixXd gradtheta_U = vertices * gradtheta_R.transpose();
        flatten<double>(gradtheta_U);

        Eigen::MatrixXd gradient(vertices.size(), 3);
        gradient << gradx_U, grady_U, gradtheta_U;
        return gradient;
    }

    std::vector<Eigen::Matrix3d> RigidBody::world_vertices_hessian_exact(
        const Eigen::Vector3d& _position) const
    {
        /// For a single vertice ri, hessian is
        // [[[0, 0, 0], [0, 0, 0], [0, 0, -ri_x*cos(\theta) +
        // ri_y*sin(\theta)]],
        //  [[0, 0, 0], [0, 0, 0], [0, 0, -ri_x*sin(\theta) -
        //  ri_y*cos(\theta)]]]
        ///
        ///
        /// each hessian has shape 3x3, we return a list of vertice.size()
        /// order of list is x-positions, y-positions
        ///

        Eigen::Matrix2d hesstheta_R
            = -Eigen::Rotation2D<double>(_position(2)).toRotationMatrix();
        Eigen::MatrixXd gradtheta_gradtheta_U
            = vertices * hesstheta_R.transpose();
        flatten<double>(gradtheta_gradtheta_U);

        std::vector<Eigen::Matrix3d> hessian(
            size_t(vertices.size()), Eigen::Matrix3d::Zero());
        for (long i = 0; i < vertices.size(); i++) {
            hessian[size_t(i)](2, 2) = gradtheta_gradtheta_U(i);
        }
        return hessian;
    }

} // namespace physics
} // namespace ccd
