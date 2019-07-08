#pragma once

#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <physics/rigid_body.hpp>

#include <autodiff/autodiff.h>

namespace ccd {
namespace physics {

    class RigidBodyAssembler {
    public:
        RigidBodyAssembler() {}
        ~RigidBodyAssembler() {}

        /// \brief inits assembler to use this set of rigid-bodies
        void init(const std::vector<RigidBody>& rbs);

        // World Vertices Functions
        // ----------------------------------------------------------------------
        /// \brief assemble world coordinates of the rigid-body vertices
        Eigen::MatrixXd world_vertices(const Eigen::VectorXd& positions);
        Eigen::MatrixXd world_vertices(const bool previous = false);
        void world_vertices_gradient(const Eigen::VectorXd& positions,
            Eigen::SparseMatrix<double>& grad);
        void world_vertices_hessian(const Eigen::VectorXd& positions,
                    std::vector<Eigen::SparseMatrix<double>>& hess);

        // Ridig Body CM Functions
        // ----------------------------------------------------------------------
        /// \brief assemble rigid body positions to a single vector
        Eigen::VectorXd rb_positions(const bool previous = false);
        /// \brief set rigid body positions
        void set_rb_positions(const Eigen::VectorXd&);

        long num_vertices() { return m_body_vertex_id.back(); }
        long num_edges() { return m_body_edge_id.back(); }


        std::vector<RigidBody> m_rbs;

        /// @brief indexes for body
        std::vector<long> m_body_vertex_id;
        std::vector<long> m_body_edge_id;

        /// @brief indexes for vertices
        Eigen::VectorXi m_vertex_to_body_map;

        /// @brief  re-indexed edges of the whole system
        Eigen::MatrixXi m_edges;

        /// @brief  mass_matrix of the rigid bodies
        Eigen::SparseMatrix<double> m_mass_matrix;
    };

    class RigidBodySystem {
    public:
        RigidBodySystem() {}
        ~RigidBodySystem() {}

        // ------------------------------------------------------------------------
        // State Functions
        // ------------------------------------------------------------------------

        void clear();
        void add_rigid_body(RigidBody rb) { rigid_bodies.push_back(rb); }
        void add_rigid_body(std::vector<RigidBody> rbs)
        {
            rigid_bodies.insert(rigid_bodies.end(), rbs.begin(), rbs.end());
        }
        unsigned long get_rigid_body_count() { return rigid_bodies.size(); }

        // ------------------------------------------------------------------------
        // CCD Functions
        // ------------------------------------------------------------------------
        void assemble();

        /// @brief compute_displacements: computes particles displacements for
        /// the CURRENT velocities updates `displacements` attribute.
        void assemble_displacements();

        /// @brief compute_displacements: computes particles displacements for
        /// the GIVEN velocities
        ///
        /// @param v[in]:   rigid bodies velocities  Bx3 by 1
        /// @param u[out]:  particles displacements N by 2
        void compute_displacements(
            const Eigen::VectorXd& v, Eigen::MatrixXd& u);

        void compute_displacements_gradient(
            const Eigen::VectorXd& v, Eigen::SparseMatrix<double>& grad_u);

        void compute_displacements_hessian(const Eigen::VectorXd& v,
            std::vector<Eigen::SparseMatrix<double>>& hess_u);

        void set_velocity(const size_t rb_id, const Eigen::Vector3d vel);
        const Eigen::Vector3d& get_velocity(const size_t rb_id)
        {
            return rigid_bodies[rb_id].velocity;
        }

        // ------------------------------------------------------------------------
        // CCD structures
        // ------------------------------------------------------------------------

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

        std::vector<RigidBody> rigid_bodies;
    };
} // namespace physics
} // namespace ccd
