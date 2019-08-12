#pragma once

#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <utils/eigen_ext.hpp>

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
        Eigen::MatrixXd world_vertices(
            const std::vector<Eigen::Vector3d>& positions) const;
        Eigen::MatrixXd world_vertices(const Eigen::VectorXd& positions) const;
        Eigen::MatrixXd world_vertices(
            const RigidBody::Step step = RigidBody::CURRENT_STEP) const;
        Eigen::MatrixXd world_vertices_t0() const
        {
            return world_vertices(RigidBody::PREVIOUS_STEP);
        }
        Eigen::MatrixXd world_vertices_t1() const
        {
            return world_vertices(RigidBody::CURRENT_STEP);
        }
        void world_vertices_gradient(const Eigen::VectorXd& positions,
            Eigen::SparseMatrix<double>& grad) const;
        void world_vertices_hessian(const Eigen::VectorXd& positions,
            std::vector<Eigen::SparseMatrix<double>>& hess) const;

        Eigen::MatrixXd world_velocities() const;

        void global_to_local(const int global_vertex_id,
            int& rigid_body_id,
            int& local_vertex_id);

        // Ridig Body CM Functions
        // ----------------------------------------------------------------------
        /// \brief assemble rigid body positions to a single vector
        Eigen::VectorXd rb_positions(const bool previous = false) const;
        Eigen::VectorXd rb_positions_t0() const { return rb_positions(true); }
        Eigen::VectorXd rb_positions_t1() const { return rb_positions(false); }
        /// \brief set rigid body positions
        void set_rb_positions(const Eigen::VectorXd&);

        // ----------------------------------------------------------------------

        long num_vertices() const { return m_body_vertex_id.back(); }
        long num_edges() const { return m_body_edge_id.back(); }

        std::vector<RigidBody> m_rbs;

        /// @brief indexes for body
        std::vector<long> m_body_vertex_id;
        std::vector<long> m_body_edge_id;

        /// @brief indexes for vertices
        Eigen::VectorXi m_vertex_to_body_map;

        /// @brief  re-indexed edges of the whole system
        Eigen::MatrixXi m_edges;

        /// @brief  mass_matrix of the rigid bodies
        Eigen::SparseMatrix<double> m_rb_mass_matrix;

        /// @brief  mass_matrix of the vertices
        Eigen::SparseMatrix<double> m_mass_matrix;

        /// @brief inverted mass_matrix of the vertices
        Eigen::SparseMatrix<double> m_inv_mass_matrix;

        /// @brief flag for rb degrees of freedom
        Eigen::VectorXb is_rb_dof_fixed;

        /// @brief flag for vertices degrees of freedom (used for visualization)
        Eigen::MatrixXb is_dof_fixed;
    };

} // namespace physics
} // namespace ccd
