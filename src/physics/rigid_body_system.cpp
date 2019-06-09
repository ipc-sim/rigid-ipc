#include <physics/rigid_body_system.hpp>

#include <Eigen/Geometry>
#include <logger.hpp>
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

    Eigen::MatrixX2d RigidBody::world_vertices() const
    {
        return vertices.rowwise() + position.transpose();
    }

    Eigen::MatrixXd RigidBody::world_displacements() const
    {
        return world_displacements(velocity);
    }

    Eigen::MatrixXd RigidBody::world_displacements(
        const Eigen::Vector3d& vel) const
    {
        typedef Eigen::Translation<double, 2> Translation2d;
        typedef Eigen::Rotation2D<double> Rotation2Dd;

        // create transform matrix matrix
        Eigen::Matrix<double, 3, 3> Tm
            = Eigen::Transform<double, 2, Eigen::Affine>(
                Translation2d(vel.x(), vel.y())
                * Translation2d(position.x(), position.y())
                * Rotation2Dd(vel(2))
                * Translation2d(position.x(), -position.y()))
                  .matrix();

        Eigen::MatrixXd v = world_vertices();
        Eigen::Matrix<double, Eigen::Dynamic, 3> hvertices
            = v.rowwise().homogeneous();
        Eigen::Matrix<double, Eigen::Dynamic, 3> tvertices
            = hvertices * Tm.transpose();

        Eigen::MatrixXd displacements;
        displacements = tvertices.rowwise().hnormalized();
        displacements -= v;

        return displacements;
    }

    std::vector<Eigen::MatrixX2d>
    RigidBody::compute_world_displacements_gradient() const
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

        Eigen::Matrix3d gradx_T = Eigen::Transform<double, 2, Eigen::Affine>(
                                      Eigen::Translation2d(1, 0))
                                      .matrix()
            * T_c * R_θ * T_negc;
        Eigen::Matrix3d grady_T = Eigen::Transform<double, 2, Eigen::Affine>(
                                      Eigen::Translation2d(0, 1))
                                      .matrix()
            * T_c * R_θ * T_negc;
        Eigen::Matrix3d gradθ_R;
        // clang-format off
    gradθ_R <<
        -sin(θ), -cos(θ), 0,
         cos(θ), -sin(θ), 0,
              0,       0, 0;
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

    void RigidBodySystem::clear()
    {
        rigid_bodies.clear();
        acc_vertex_id.clear();
        acc_edge_id.clear();
    }
    void RigidBodySystem::assemble()
    {
        size_t num_bodies = rigid_bodies.size();

        velocities.resize(int(3 * num_bodies));

        acc_vertex_id.clear();
        acc_edge_id.clear();
        acc_vertex_id.resize(num_bodies + 1);
        acc_edge_id.resize(num_bodies + 1);
        acc_vertex_id[0] = 0;
        acc_edge_id[0] = 0;
        for (size_t i = 0; i < num_bodies; ++i) {
            auto& rb = rigid_bodies[i];
            velocities.segment(int(3 * i), 3) = rb.velocity;
            acc_vertex_id[i + 1] = acc_vertex_id[i] + rb.vertices.rows();
            acc_edge_id[i + 1] = acc_edge_id[i] + rb.edges.rows();
        }

        vertex_to_body_map.resize(acc_vertex_id.back());
        vertices.resize(acc_vertex_id.back(), 2);
        edges.resize(acc_edge_id.back(), 2);
        for (size_t i = 0; i < num_bodies; ++i) {
            auto& rb = rigid_bodies[i];
            vertices.block(acc_vertex_id[i], 0, rb.vertices.rows(), 2)
                = rb.world_vertices();
            vertex_to_body_map.segment(acc_vertex_id[i], rb.vertices.rows())
                .setConstant(int(i));
            edges.block(acc_edge_id[i], 0, rb.edges.rows(), 2)
                = (rb.edges.array() + int(acc_vertex_id[i]));
        }

        assemble_displacements();
    }

    void RigidBodySystem::set_velocity(
        const size_t rb_id, const Eigen::Vector3d vel)
    {
        rigid_bodies[rb_id].velocity = vel;
        velocities.segment(int(3 * rb_id), 3) = vel;
    }

    void RigidBodySystem::assemble_displacements()
    {
        assemble_displacements(velocities, displacements);
    }

    void RigidBodySystem::assemble_displacements(
        const Eigen::VectorXd& v, Eigen::MatrixXd& u)
    {
        u.resize(acc_vertex_id.back(), 2);
        for (size_t i = 0; i < rigid_bodies.size(); ++i) {
            auto& rb = rigid_bodies[i];
            u.block(acc_vertex_id[i], 0, rb.vertices.rows(), 2)
                = rb.world_displacements(v.segment(int(3 * i), 3));
        }
    }

} // namespace physics
} // namespace ccd
