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
        typedef Eigen::Translation<T, 2> Translation2d;
        typedef Eigen::Rotation2D<T> Rotation2Dd;
        typedef Eigen::Matrix<T, 3, 3> Matrix3d;
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
        typedef Eigen::Matrix<T, Eigen::Dynamic, 3> MatrixX3d;

        // create transform matrix matrix
        Matrix3d Tm = Eigen::Transform<T, 2, Eigen::Affine>(
            Translation2d(vel.x(), vel.y())
            * Translation2d(T(position.x()), T(position.y()))
            * Rotation2Dd(vel(2))
            * Translation2d(T(position.x()), -T(position.y())))
                          .matrix();

        MatrixXd v = world_vertices().cast<T>();

        MatrixX3d hvertices = v.rowwise().homogeneous();
        MatrixX3d tvertices = hvertices * Tm.transpose();

        MatrixXd displacements;
        displacements = tvertices.rowwise().hnormalized();
        displacements -= v;

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

    std::vector<Eigen::MatrixX2d> RigidBody::world_displacements_gradient_exact(
        const Eigen::Vector3d& vel) const
    {
        typedef Eigen::Transform<double, 2, Eigen::Affine> AffineTransform;
        Eigen::Matrix3d T_c
            = AffineTransform(Eigen::Translation2d(position.x(), position.y()))
                  .matrix();
        Eigen::Matrix3d T_negc = AffineTransform(
            Eigen::Translation2d(-position.x(), -position.y()))
                                     .matrix();
        Eigen::Matrix3d T_xy
            = AffineTransform(Eigen::Translation2d(vel.x(), vel.y())).matrix();
        double θ = vel(2);
        Eigen::Matrix3d R_θ = AffineTransform(Eigen::Rotation2Dd(θ)).matrix();

        Eigen::Matrix3d gradx_T
            = AffineTransform(Eigen::Translation2d(1, 0)).matrix() * T_c * R_θ
            * T_negc;
        Eigen::Matrix3d grady_T = Eigen::Transform<double, 2, Eigen::Affine>(
                                      Eigen::Translation2d(0, 1))
                                      .matrix()
            * T_c * R_θ * T_negc;
        Eigen::Matrix3d gradθ_R;
        // clang-format off
    gradθ_R <<
        -sin(θ), -cos(θ), 0,
         cos(θ), -sin(θ), 0,
              0,       0, 1;
        // clang-format on
        Eigen::Matrix3d gradθ_T = T_xy * T_c * gradθ_R * T_negc;

        const auto& this_vertices = world_vertices();
        auto transform_vertices
            = [&this_vertices](const Eigen::Matrix3d& T) -> Eigen::MatrixX2d {
            return (this_vertices.rowwise().homogeneous() * T.transpose())
                .rowwise()
                .hnormalized();
        };

        std::vector<Eigen::MatrixX2d> gradient;
        gradient.reserve(3);
        gradient.push_back(transform_vertices(gradx_T));
        gradient.push_back(transform_vertices(grady_T));
        gradient.push_back(transform_vertices(gradθ_T));

        return gradient;
    }

    std::vector<std::vector<Eigen::MatrixX2d>>
    RigidBody::compute_world_displacements_hessian() const
    {
        Eigen::Matrix3d T_c = Eigen::Transform<double, 2, Eigen::Affine>(
            Eigen::Translation2d(position.x(), position.y()))
                                  .matrix();
        Eigen::Matrix3d T_negc = Eigen::Transform<double, 2, Eigen::Affine>(
            Eigen::Translation2d(-position.x(), -position.y()))
                                     .matrix();
        Eigen::Matrix3d T_xy = Eigen::Transform<double, 2, Eigen::Affine>(
            Eigen::Translation2d(velocity.x(), velocity.y()))
                                   .matrix();
        double θ = velocity(2);
        Eigen::Matrix3d R_θ
            = Eigen::Transform<double, 2, Eigen::Affine>(Eigen::Rotation2Dd(θ))
                  .matrix();

        Eigen::Matrix3d gradx_T_xy = Eigen::Transform<double, 2, Eigen::Affine>(
            Eigen::Translation2d(1, 0))
                                         .matrix();
        Eigen::Matrix3d grady_T_xy = Eigen::Transform<double, 2, Eigen::Affine>(
            Eigen::Translation2d(0, 1))
                                         .matrix();
        Eigen::Matrix3d gradx_T = gradx_T_xy * T_c * R_θ * T_negc;
        Eigen::Matrix3d grady_T = grady_T_xy * T_c * R_θ * T_negc;
        Eigen::Matrix3d gradθ_R;
        // clang-format off
    // R_θ <<
    //     cos(θ), -sin(θ), 0,
    //     sin(θ),  cos(θ), 0,
    //          0,       0, 0;
    gradθ_R <<
        -sin(θ), -cos(θ), 0,
         cos(θ), -sin(θ), 0,
              0,       0, 0;
    // gradθ_gradθ_R <<
    //     -cos(θ),  sin(θ), 0,
    //     -sin(θ), -cos(θ), 0,
    //           0,       0, 0;
        // clang-format on
        Eigen::Matrix3d gradθ_T = T_xy * T_c * gradθ_R * T_negc;

        Eigen::Matrix3d gradx_gradx_T = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d gradx_grady_T = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d gradx_gradθ_T = gradx_T_xy * T_c * gradθ_R * T_negc;

        Eigen::Matrix3d grady_gradx_T = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d grady_grady_T = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d grady_gradθ_T = grady_T_xy * T_c * gradθ_R * T_negc;

        Eigen::Matrix3d gradθ_gradx_T = gradx_gradθ_T;
        Eigen::Matrix3d gradθ_grady_T = grady_gradθ_T;
        Eigen::Matrix3d gradθ_gradθ_T = T_xy * T_c * -R_θ * T_negc;

        Eigen::MatrixX3d homogeneous_vertices(vertices.rows(), 3);
        homogeneous_vertices.leftCols(2) = world_vertices();
        ;
        homogeneous_vertices.col(2).setOnes();

        std::vector<std::vector<Eigen::MatrixX2d>> hessian(
            3, std::vector<Eigen::MatrixX2d>());

        const auto& this_vertices = world_vertices();
        auto transform_vertices
            = [&this_vertices](const Eigen::Matrix3d& T) -> Eigen::MatrixX2d {
            return (this_vertices.rowwise().homogeneous() * T.transpose())
                .rowwise()
                .hnormalized();
        };

        // ∇_x∇U
        hessian[0].reserve(3);
        hessian[0].push_back(transform_vertices(gradx_gradx_T));
        hessian[0].push_back(transform_vertices(gradx_grady_T));
        hessian[0].push_back(transform_vertices(gradx_gradθ_T));

        // ∇_y∇U
        hessian[1].reserve(3);
        hessian[1].push_back(transform_vertices(grady_gradx_T));
        hessian[1].push_back(transform_vertices(grady_grady_T));
        hessian[1].push_back(transform_vertices(grady_gradθ_T));

        // ∇_θ∇U
        hessian[2].reserve(3);
        hessian[2].push_back(transform_vertices(gradθ_gradx_T));
        hessian[2].push_back(transform_vertices(gradθ_grady_T));
        hessian[2].push_back(transform_vertices(gradθ_gradθ_T));

        return hessian;
    }

} // namespace physics
} // namespace ccd
