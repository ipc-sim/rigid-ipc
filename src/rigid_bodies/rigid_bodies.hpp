#include <Eigen/Core>
#include <Eigen/Geometry>

void compute_particle_displacments(const Eigen::MatrixX2d& vertices,
    const Eigen::Vector3d& body_velocity, Eigen::MatrixX2d& displacements);
