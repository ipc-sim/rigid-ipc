#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <physics/rigid_body.hpp>

#include <autodiff/autodiff.h>

namespace ccd {
namespace physics {

    class RigidBodySystem {
    public:
        RigidBodySystem() {}
        ~RigidBodySystem() {}

        void clear();
        void assemble();

        ///
        /// \brief compute_displacements: computes particles displacements for the CURRENT velocities
        /// updates `displacements` attribute.
        ///
        void assemble_displacements();

        ///
        /// \brief compute_displacements: computes particles displacements for the GIVEN velocities
        /// \param v[in]:   rigid bodies velocities  Bx3 by 1
        /// \param u[out]:  particles displacements N by 2
        ///
        void compute_displacements(
            const Eigen::VectorXd& v, Eigen::MatrixXd& u);

        void compute_displacements_gradient(
            const Eigen::VectorXd& v, Eigen::SparseMatrix<double>& grad_u);

        void compute_displacements_hessian(
            const Eigen::VectorXd& v, std::vector<Eigen::SparseMatrix<double>>& hess_u);

        void add_rigid_body(RigidBody rb) { rigid_bodies.push_back(rb); }
        void set_velocity(const size_t rb_id, const Eigen::Vector3d vel);
        const Eigen::Vector3d& get_velocity(const size_t rb_id)
        {
            return rigid_bodies[rb_id].velocity;
        }

        ///> velocities: flatten array with all RB velocities
        Eigen::VectorXd velocities;
        std::vector<long> acc_vertex_id;
        std::vector<long> acc_edge_id;
        Eigen::VectorXi vertex_to_body_map;

        // world-space units of the RB, used for ccd
        Eigen::MatrixXd vertices;
        Eigen::MatrixXd displacements;
        Eigen::MatrixX2i edges;
        Eigen::SparseMatrix<double> mass_matrix;

    private:
        std::vector<RigidBody> rigid_bodies;
    };
} // namespace physics
} // namespace ccd
