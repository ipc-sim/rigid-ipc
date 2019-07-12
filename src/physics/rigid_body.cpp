#include <physics/rigid_body.hpp>

#include <Eigen/Geometry>

#include <physics/mass_matrix.hpp>
#include <physics/moment_of_inertia.hpp>

#include <autodiff/finitediff.hpp>
#include <logger.hpp>
#include <utils/flatten.hpp>

namespace ccd {

namespace physics {

    RigidBody RigidBody::from_velocity(const Eigen::MatrixXd& vertices,
        const Eigen::MatrixX2i& edges,
        const Eigen::Vector3d& velocity,
        const Eigen::Vector3b& is_dof_fixed)
    {
        // move vertices so they center of mass is at 0,0
        Eigen::RowVector2d x = center_of_mass(vertices, edges);
        Eigen::MatrixX2d centered_vertices = vertices.rowwise() - x;

        // set position so current vertices match input
        Eigen::Vector3d position_t1 = Eigen::Vector3d::Zero();
        position_t1.segment(0, 2) = x;

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

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> RigidBody::world_vertices(
        const Eigen::Matrix<T, 3, 1>& _position) const
    {
        typedef Eigen::Matrix<T, 2, 2> Matrix2T;

        // compute X[i] = R(theta) * r_i + X
        Matrix2T R = Eigen::Rotation2D<T>(_position.z()).toRotationMatrix();
        return (vertices.cast<T>() * R.transpose()).rowwise()
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

    // AUTODIFF DERIVATIVES
    // ---------------------------------------------------------------
    Eigen::MatrixXd RigidBody::world_vertices_gradient(
        const Eigen::Vector3d& _position) const
    {
        RBDiff::activate();
        RBDiff::D1Vector3 dpos = RBDiff::d1vars(0, _position);
        RBDiff::D1MatrixXd dx = world_vertices<RBDiff::DScalar1>(dpos);

        flatten<RBDiff::DScalar1>(dx);
        return RBDiff::get_gradient(dx);
    }

    std::vector<Eigen::Matrix3d> RigidBody::world_vertices_hessian(
        const Eigen::Vector3d& _position) const
    {
        RBDiff::activate();
        RBDiff::D2Vector3 dpos = RBDiff::d2vars(0, _position);
        RBDiff::D2MatrixXd dx = world_vertices<RBDiff::DScalar2>(dpos);

        flatten<RBDiff::DScalar2>(dx);
        return RBDiff::get_hessian(dx);
    }

    // FINITE DERIVATIVES
    // ---------------------------------------------------------------
    Eigen::MatrixXd RigidBody::world_vertices_gradient_finite(
        const Eigen::Vector3d& _position) const
    {
        auto f = [&](const Eigen::VectorXd& v) -> Eigen::VectorXd {
            assert(v.rows() == 3);
            Eigen::Vector3d v3 = v;
            Eigen::MatrixXd d = world_vertices<double>(v3);
            flatten<double>(d);
            return d;
        };

        Eigen::MatrixXd grad;
        ccd::finite_jacobian(_position, f, grad);
        return grad;
    }

    // EXACT DERIVATIVES
    // ---------------------------------------------------------------
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
