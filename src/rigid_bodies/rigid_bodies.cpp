#include <rigid_bodies/rigid_bodies.hpp>

void compute_particle_displacments(const Eigen::MatrixX2d& vertices,
    const Eigen::Vector3d& body_velocity, Eigen::MatrixX2d& displacements)
{
    Eigen::Vector2d center_of_mass = vertices.colwise().sum() / vertices.rows();
    Eigen::Matrix3d T = Eigen::Transform<double, 2, Eigen::Affine>(
        Eigen::Translation2d(body_velocity.x(), body_velocity.y())
        * Eigen::Translation2d(center_of_mass.x(), center_of_mass.y())
        * Eigen::Rotation2Dd(body_velocity(2))
        * Eigen::Translation2d(-center_of_mass.x(), -center_of_mass.y()))
                            .matrix();
    // V' = (R * (V - c).T).T + c + t
    Eigen::MatrixX3d new_vertices = Eigen::MatrixX3d::Ones(vertices.rows(), 3);
    new_vertices.block(0, 0, vertices.rows(), vertices.cols()) = vertices;
    new_vertices = new_vertices * T.transpose();
    displacements
        = new_vertices.block(0, 0, vertices.rows(), vertices.cols()) - vertices;
}
