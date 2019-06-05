#include <rigid_bodies/rigid_body.hpp>

RigidBody::RigidBody(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2i& edges, const Eigen::Vector3d& velocity)
    : vertices(vertices)
    , edges(edges)
    , velocity(velocity)
{
    this->update_center_of_mass();
}

void RigidBody::update_center_of_mass()
{
    Eigen::VectorXd vertex_masses = Eigen::VectorXd::Zero(vertices.rows());
    double total_edge_length = 0;
    for (long i = 0; i < edges.rows(); i++) {
        double edge_length
            = (vertices.row(edges(i, 1)) - vertices.row(edges(i, 0))).norm();
        vertex_masses(edges(i, 0)) += edge_length / 2;
        vertex_masses(edges(i, 1)) += edge_length / 2;
        total_edge_length += edge_length;
    }
    center_of_mass = (vertex_masses.asDiagonal() * vertices).colwise().sum()
        / total_edge_length;
}

Eigen::MatrixX2d RigidBody::compute_particle_displacements() const
{
    // T = T_xy * T_c * R_θ * T_-c
    Eigen::Matrix3d T = Eigen::Transform<double, 2, Eigen::Affine>(
        Eigen::Translation2d(velocity.x(), velocity.y())
        * Eigen::Translation2d(center_of_mass.x(), center_of_mass.y())
        * Eigen::Rotation2Dd(velocity(2))
        * Eigen::Translation2d(-center_of_mass.x(), -center_of_mass.y()))
                            .matrix();
    // V' = V + U = V * T^T
    Eigen::MatrixX2d transformed_vertices
        = (vertices.rowwise().homogeneous() * T.transpose())
              .rowwise()
              .hnormalized();
    // U = V' - V
    return transformed_vertices - vertices;
}

std::vector<Eigen::MatrixX2d>
RigidBody::compute_particle_displacements_gradient() const
{
    Eigen::Matrix3d T_c = Eigen::Transform<double, 2, Eigen::Affine>(
        Eigen::Translation2d(center_of_mass.x(), center_of_mass.y()))
                              .matrix();
    Eigen::Matrix3d T_negc = Eigen::Transform<double, 2, Eigen::Affine>(
        Eigen::Translation2d(-center_of_mass.x(), -center_of_mass.y()))
                                 .matrix();
    Eigen::Matrix3d T_xy = Eigen::Transform<double, 2, Eigen::Affine>(
        Eigen::Translation2d(velocity.x(), velocity.y()))
                               .matrix();
    double θ = velocity(2);
    Eigen::Matrix3d R_θ
        = Eigen::Transform<double, 2, Eigen::Affine>(Eigen::Rotation2Dd(θ))
              .matrix();

    Eigen::Matrix3d gradx_T
        = Eigen::Transform<double, 2, Eigen::Affine>(Eigen::Translation2d(1, 0))
              .matrix()
        * T_c * R_θ * T_negc;
    Eigen::Matrix3d grady_T
        = Eigen::Transform<double, 2, Eigen::Affine>(Eigen::Translation2d(0, 1))
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

    const auto& this_vertices = vertices;
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
RigidBody::compute_particle_displacements_hessian() const
{
    Eigen::Matrix3d T_c = Eigen::Transform<double, 2, Eigen::Affine>(
        Eigen::Translation2d(center_of_mass.x(), center_of_mass.y()))
                              .matrix();
    Eigen::Matrix3d T_negc = Eigen::Transform<double, 2, Eigen::Affine>(
        Eigen::Translation2d(-center_of_mass.x(), -center_of_mass.y()))
                                 .matrix();
    Eigen::Matrix3d T_xy = Eigen::Transform<double, 2, Eigen::Affine>(
        Eigen::Translation2d(velocity.x(), velocity.y()))
                               .matrix();
    double θ = velocity(2);
    Eigen::Matrix3d R_θ
        = Eigen::Transform<double, 2, Eigen::Affine>(Eigen::Rotation2Dd(θ))
              .matrix();

    Eigen::Matrix3d gradx_T_xy
        = Eigen::Transform<double, 2, Eigen::Affine>(Eigen::Translation2d(1, 0))
              .matrix();
    Eigen::Matrix3d grady_T_xy
        = Eigen::Transform<double, 2, Eigen::Affine>(Eigen::Translation2d(0, 1))
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
    homogeneous_vertices.leftCols(2) = vertices;
    homogeneous_vertices.col(2).setOnes();

    std::vector<std::vector<Eigen::MatrixX2d>> hessian(
        3, std::vector<Eigen::MatrixX2d>());

    const auto& this_vertices = vertices;
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
