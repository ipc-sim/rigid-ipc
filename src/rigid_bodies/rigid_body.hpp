#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

namespace ccd {
namespace opt {

    class RigidBody {
    public:
        RigidBody(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::Vector3d& velocity);

        void compute_particle_displacements(
            Eigen::MatrixXd& displacements) const;

        void compute_particle_displacements(
            const Eigen::Vector3d& velocity,
            Eigen::MatrixXd& displacements) const;


    Eigen::MatrixX2d compute_particle_displacements() const;
    std::vector<Eigen::MatrixX2d>
    compute_particle_displacements_gradient() const;
    std::vector<std::vector<Eigen::MatrixX2d>>
    compute_particle_displacements_hessian() const;
        const Eigen::MatrixX2d& vertices() const { return _vertices; }
        const Eigen::MatrixX2i& edges() const { return _edges; }

        Eigen::Vector3d velocity;

    private:
        // attributes are privite to ensure consistency of center of mass
        Eigen::MatrixX2d _vertices;
        Eigen::MatrixX2i _edges;

        Eigen::Vector2d center_of_mass;
        void compute_center_of_mass(Eigen::Vector2d& cm) const;
    };
} // namespace opt
} // namespace ccd
>>>>>>> 5ba3f2a... moved RigidBody class to namespace ccd::opt
