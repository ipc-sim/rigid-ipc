#include <Eigen/Core>
#include <Eigen/Geometry>

class RigidBody {
public:
    Eigen::MatrixX2d vertices;
    Eigen::MatrixX2i edges;
    Eigen::Vector3d velocity;

    Eigen::Vector2d center_of_mass;

    RigidBody();
    RigidBody(const Eigen::MatrixX2d& vertices, const Eigen::MatrixX2i& edges,
        const Eigen::Vector3d& velocity);

    void update_center_of_mass();

    void compute_particle_displacements(Eigen::MatrixX2d& displacements) const;
};
