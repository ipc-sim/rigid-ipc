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

        /// @brief inits assembler to use this set of rigid-bodies
        void init(const std::vector<RigidBody>& rbs);

        // World Vertices Functions
        // ----------------------------------------------------------------------

        template <typename T>
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
        world_vertices(const std::vector<Pose<T>>& v) const;
        template <typename T>
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
        world_vertices(const Eigen::VectorX<T>& dof) const
        {
            return world_vertices(Pose<T>::dofs_to_poses(dof, dim()));
        }

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

        Eigen::MatrixXd world_velocities() const;

        void world_vertices_gradient(
            const std::vector<Pose<double>>& poses,
            Eigen::SparseMatrix<double>& grad) const;
        void world_vertices_gradient(
            const Eigen::VectorXd& dof, Eigen::SparseMatrix<double>& grad) const
        {
            return world_vertices_gradient(
                Pose<double>::dofs_to_poses(dof, dim()), grad);
        }

        void global_to_local(
            const int global_vertex_id,
            int& rigid_body_id,
            int& local_vertex_id);

        // Ridig Body CM Functions
        // ----------------------------------------------------------------------
        /// @brief assemble rigid body poses to a single vector
        std::vector<Pose<double>> rb_poses(const bool previous = false) const;
        std::vector<Pose<double>> rb_poses_t0() const { return rb_poses(true); }
        std::vector<Pose<double>> rb_poses_t1() const
        {
            return rb_poses(false);
        }
        /// @brief set rigid body poses
        void set_rb_poses(const std::vector<Pose<double>>& poses);

        // ----------------------------------------------------------------------

        long num_vertices() const { return m_body_vertex_id.back(); }
        long num_edges() const { return m_body_edge_id.back(); }
        long num_faces() const { return m_body_face_id.back(); }
        long num_bodies() const { return m_rbs.size(); }
        int dim() const { return m_rbs.size() ? m_rbs[0].dim() : 0; }

        std::vector<RigidBody> m_rbs;

        /// @brief indexes for body
        std::vector<long> m_body_vertex_id;
        std::vector<long> m_body_face_id;
        std::vector<long> m_body_edge_id;

        /// @brief indexes for vertices
        Eigen::VectorXi m_vertex_to_body_map;

        /// @brief  re-indexed edges of the whole system
        Eigen::MatrixXi m_edges;

        /// @brief  re-indexed faces of the whole system
        Eigen::MatrixXi m_faces;

        /// @brief  mass_matrix of the rigid bodies dof
        Eigen::SparseMatrix<double> m_rb_mass_matrix;

        /// @brief transforms dof to poses
        Eigen::SparseMatrix<double> m_dof_to_pose;

        /// @brief transforms poses to dof
        Eigen::SparseMatrix<double> m_pose_to_dof;

        /// @brief  mass_matrix of the vertices
        Eigen::SparseMatrix<double> m_mass_matrix;

        /// @brief inverted mass_matrix of the vertices
        Eigen::SparseMatrix<double> m_inv_mass_matrix;

        /// @brief flags for rb degrees of freedom
        Eigen::VectorXb is_rb_dof_fixed;

        /// @brief flag for vertices degrees of freedom (used for visualization)
        Eigen::MatrixXb is_dof_fixed;
    };

} // namespace physics
} // namespace ccd

#include "rigid_body_assembler.tpp"
