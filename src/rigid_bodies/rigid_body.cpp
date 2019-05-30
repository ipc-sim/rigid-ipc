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

void RigidBody::compute_particle_displacements(
    Eigen::MatrixX2d& displacements) const
{
    Eigen::Matrix3d T = Eigen::Transform<double, 2, Eigen::Affine>(
        Eigen::Translation2d(velocity.x(), velocity.y())
        * Eigen::Translation2d(center_of_mass.x(), center_of_mass.y())
        * Eigen::Rotation2Dd(velocity(2))
        * Eigen::Translation2d(-center_of_mass.x(), -center_of_mass.y()))
                            .matrix();
    // V' = (R * (V - c).T).T + c + t
    Eigen::MatrixX3d transformed_vertices
        = Eigen::MatrixX3d::Ones(vertices.rows(), 3);
    transformed_vertices.block(0, 0, vertices.rows(), vertices.cols())
        = vertices;
    transformed_vertices = transformed_vertices * T.transpose();
    displacements
        = transformed_vertices.block(0, 0, vertices.rows(), vertices.cols())
        - vertices;
}
