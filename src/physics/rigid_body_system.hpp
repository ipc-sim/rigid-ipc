#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <physics/center_of_mass.hpp>

namespace ccd {
namespace physics {

    class RigidBody {
    public:
        RigidBody(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::Vector2d& position,
            const Eigen::Vector3d& velocity);

        Eigen::MatrixXd world_displacements() const;

        Eigen::MatrixXd world_displacements(
            const Eigen::Vector3d& velocity) const;

        Eigen::MatrixX2d world_vertices() const
        {
            return vertices.rowwise() + position.transpose();
        }

        std::vector<Eigen::MatrixX2d>
        compute_world_displacements_gradient() const;
        std::vector<std::vector<Eigen::MatrixX2d>>
        compute_world_displacements_hessian() const;

        /// \brief vertices as distances from the center of mass
        Eigen::MatrixX2d vertices;
        Eigen::MatrixX2i edges;

        /// \brief position: position of the object in world space
        Eigen::Vector2d position;

        /// \brief velocity: linear and angular velicity
        Eigen::Vector3d velocity;

        ///
        /// \brief Centered: Factory method to create a RB with
        /// with position equal to current center of mass.
        /// \param vertices
        /// \param edges
        /// \param velocity
        /// \return
        ///
        static inline RigidBody Centered(const Eigen::MatrixXd& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::Vector3d& velocity)
        {
            Eigen::RowVector2d x = center_of_mass(vertices, edges);
            Eigen::MatrixX2d centered_vertices = vertices.rowwise() - x;
            return RigidBody(centered_vertices, edges, x, velocity);
        }
    };

    class RigidBodySystem {
    public:
        RigidBodySystem() {}
        ~RigidBodySystem() {}
        void assemble();
        void assemble_displacements(
            const Eigen::VectorXd& v, Eigen::MatrixXd& u);
        void add_rigid_body(RigidBody rb) { rigid_bodies.push_back(rb); }

        std::vector<RigidBody> rigid_bodies;

        ///> velocities: flatten array with all RB velocities
        Eigen::VectorXd velocities;
        std::vector<long> acc_vertex_id;
        std::vector<long> acc_edge_id;

        // world-space units of the RB, used for ccd
        Eigen::MatrixXd vertices;
        Eigen::MatrixXd displacements;
        Eigen::MatrixX2i edges;
        Eigen::SparseMatrix<double> mass_matrix;
    };
} // namespace physics
} // namespace ccd
