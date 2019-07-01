#include <physics/rigid_body.hpp>

#include <Eigen/Geometry>

#include <autodiff/finitediff.hpp>
#include <logger.hpp>
#include <utils/flatten.hpp>

namespace ccd {

namespace physics {

    RigidBody::RigidBody(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2i& edges, const Eigen::Vector2d& position,
        const Eigen::Vector3d& velocity)
        : vertices(vertices)
        , edges(edges)
        , position(position)
        , velocity(velocity)
    {
        spdlog::trace("rigid_body velocity={},{},{}", velocity[0], velocity[1],
            velocity[2]);
    }

    Eigen::MatrixXd RigidBody::world_displacements() const
    {
        return world_displacements(velocity);
    }

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    RigidBody::world_displacements(const Eigen::Matrix<T, 3, 1>& vel) const
    {
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXT;
        typedef Eigen::Matrix<T, 2, 2> Matrix2T;

        // Create a 2D 2x2 rotation matrix
        Matrix2T R = Eigen::Rotation2D<T>(vel(2)).toRotationMatrix();

        MatrixXT verticesT = vertices.cast<T>();

        MatrixXT transformed_vertices = (verticesT * R.transpose()).rowwise()
            + vel.topRows(2).transpose();

        MatrixXT displacements = transformed_vertices - verticesT;

        return displacements;
    }

    Eigen::MatrixX2d RigidBody::world_vertices() const
    {
        return vertices.rowwise() + position.transpose();
    }

    Eigen::MatrixXd RigidBody::world_displacements_gradient(
        const Eigen::Vector3d& vel) const
    {
        RBDiff::activate();
        RBDiff::D1Vector3 dvel = RBDiff::d1vars(0, vel);
        RBDiff::D1MatrixXd dx = world_displacements<RBDiff::DScalar1>(dvel);

        flatten<RBDiff::DScalar1>(dx);
        return RBDiff::get_gradient(dx);
    }

    std::vector<Eigen::Matrix3d> RigidBody::world_displacements_hessian(
        const Eigen::Vector3d& vel) const
    {
        RBDiff::activate();
        RBDiff::D2Vector3 dvel = RBDiff::d2vars(0, vel);
        RBDiff::D2MatrixXd dx = world_displacements<RBDiff::DScalar2>(dvel);

        flatten<RBDiff::DScalar2>(dx);
        return RBDiff::get_hessian(dx);
    }

    Eigen::MatrixXd RigidBody::world_displacements_gradient_finite(
        const Eigen::Vector3d& vel) const
    {
        auto f = [&](const Eigen::VectorXd& v) -> Eigen::VectorXd {
            assert(v.rows() == 3);
            Eigen::Vector3d v3 = v;
            Eigen::MatrixXd d = world_displacements(v3);

            flatten<double>(d);
            return d;
        };

        Eigen::MatrixXd grad;
        ccd::finite_jacobian(vel, f, grad);
        return grad;
    }

    // Return the gradient of the world displacments with respect to the
    // velocities. ∇U(σ) ∈ R^{2n × 3} where n = |V|.
    Eigen::MatrixXd RigidBody::world_displacements_gradient_exact(
        const Eigen::Vector3d& vel) const
    {
        // With respect to vel.x()
        Eigen::MatrixXd gradx_U
            = Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols());
        gradx_U.col(0).setOnes();
        flatten<double>(gradx_U);

        // With respect to vel.y()
        Eigen::MatrixXd grady_U
            = Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols());
        grady_U.col(1).setOnes();
        flatten<double>(grady_U);

        // With respect to theta = vel(2)
        double theta = vel(2);
        Eigen::Matrix2d gradtheta_R;
        // clang-format off
        gradtheta_R << -sin(theta), -cos(theta),
                    cos(theta), -sin(theta);
        // clang-format on
        Eigen::MatrixXd gradtheta_U = vertices * gradtheta_R.transpose();
        flatten<double>(gradtheta_U);

        Eigen::MatrixXd gradient(vertices.size(), 3);
        gradient << gradx_U, grady_U, gradtheta_U;
        return gradient;
    }

    std::vector<Eigen::Matrix3d> RigidBody::world_displacements_hessian_exact(
        const Eigen::Vector3d& vel) const
    {
        // With respect to vel.x()
        Eigen::VectorXd gradx_gradx_U, gradx_grady_U, gradx_gradtheta_U;
        gradx_gradx_U = gradx_grady_U = gradx_gradtheta_U
            = Eigen::VectorXd::Zero(vertices.size());

        // With respect to vel.y()
        Eigen::VectorXd grady_gradx_U, grady_grady_U, grady_gradtheta_U;
        grady_gradx_U = grady_grady_U = grady_gradtheta_U
            = Eigen::VectorXd::Zero(vertices.size());

        // With respect to theta = vel(2)
        // -cos(theta)  sin(theta)
        // -sin(theta) -cos(theta)
        Eigen::Matrix2d hesstheta_R
            = -Eigen::Rotation2D<double>(vel(2)).toRotationMatrix();
        Eigen::VectorXd gradtheta_gradx_U, gradtheta_grady_U;
        gradtheta_gradx_U = gradtheta_grady_U
            = Eigen::VectorXd::Zero(vertices.size());
        Eigen::MatrixXd gradtheta_gradtheta_U
            = vertices * hesstheta_R.transpose();
        flatten<double>(gradtheta_gradtheta_U);

        std::vector<Eigen::Matrix3d> hessian(
            static_cast<unsigned long>(vertices.size()), Eigen::Matrix3d());
        for (long i = 0; i < vertices.size(); i++) {
            // clang-format off
            hessian[static_cast<unsigned long>(i)] <<
                gradx_gradx_U(i), gradx_grady_U(i), gradx_gradtheta_U(i),
                grady_gradx_U(i), grady_grady_U(i), grady_gradtheta_U(i),
                gradtheta_gradx_U(i), gradtheta_grady_U(i), gradtheta_gradtheta_U(i);
            // clang-format on
        }
        return hessian;
    }

} // namespace physics
} // namespace ccd
