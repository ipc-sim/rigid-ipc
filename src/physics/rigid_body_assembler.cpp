#include <physics/rigid_body_assembler.hpp>

#include <Eigen/Geometry>

#include <autodiff/finitediff.hpp>

#include <physics/mass_matrix.hpp>

#include <logger.hpp>
#include <utils/eigen_ext.hpp>
#include <utils/flatten.hpp>

namespace ccd {

namespace physics {

    void RigidBodyAssembler::init(const std::vector<RigidBody>& rigid_bodies)
    {
        m_rbs = rigid_bodies;

        size_t num_bodies = rigid_bodies.size();
        m_body_vertex_id.resize(num_bodies + 1);
        m_body_edge_id.resize(num_bodies + 1);

        for (size_t i = 0; i < num_bodies; ++i) {
            auto& rb = rigid_bodies[i];
            m_body_vertex_id[i + 1] = m_body_vertex_id[i] + rb.vertices.rows();
            m_body_edge_id[i + 1] = m_body_edge_id[i] + rb.edges.rows();
        }

        // global edges
        m_edges.resize(m_body_edge_id.back(), 2);
        for (size_t i = 0; i < num_bodies; ++i) {
            auto& rb = rigid_bodies[i];
            m_edges.block(m_body_edge_id[i], 0, rb.edges.rows(), 2)
                = (rb.edges.array() + int(m_body_vertex_id[i]));
        }
        // vertex to body map
        m_vertex_to_body_map.resize(m_body_vertex_id.back());
        for (size_t i = 0; i < num_bodies; ++i) {
            auto& rb = rigid_bodies[i];
            m_vertex_to_body_map
                .segment(m_body_vertex_id[i], rb.vertices.rows())
                .setConstant(int(i));
        }

        // rigid body mass-matrix
        Eigen::VectorXd mass_vector(int(num_bodies) * 3);
        for (int i = 0; i < int(num_bodies); ++i) {
            auto& rb = rigid_bodies[size_t(i)];
            mass_vector[3 * i + 0] = rb.mass;
            mass_vector[3 * i + 1] = rb.mass;
            mass_vector[3 * i + 2] = rb.moment_of_inertia;
        }
        m_rb_mass_matrix = Eigen::SparseDiagonal<double>(mass_vector);

        // particles mass-matrix
        Eigen::VectorXd vertex_masses;
        physics::mass_vector(world_vertices_t0(), m_edges, vertex_masses);
        // Repeat the mass vector to make a mass matrix per dof
        m_mass_matrix
            = Eigen::SparseDiagonal<double>(vertex_masses.replicate(2, 1));
        m_inv_mass_matrix = m_mass_matrix.cwiseInverse();

        // rigid_body dof_fixed flag
        is_rb_dof_fixed.resize(int(num_bodies) * 3);
        for (int i = 0; i < int(num_bodies); ++i) {
            auto& rb = rigid_bodies[size_t(i)];
            is_rb_dof_fixed[3 * i + 0] = rb.is_dof_fixed[0];
            is_rb_dof_fixed[3 * i + 1] = rb.is_dof_fixed[1];
            is_rb_dof_fixed[3 * i + 2] = rb.is_dof_fixed[2];
        }

        is_dof_fixed.resize(m_body_vertex_id.back(), 3);
        for (size_t i = 0; i < num_bodies; ++i) {
            auto& rb = rigid_bodies[i];
            is_dof_fixed.block(m_body_vertex_id[i], 0, rb.vertices.rows(), 3)
                = rb.is_dof_fixed.transpose().replicate(rb.vertices.rows(), 1);
        }
    }

    void RigidBodyAssembler::global_to_local(
        const int global_vertex_id, int& rigid_body_id, int& local_vertex_id)
    {
        rigid_body_id = m_vertex_to_body_map(global_vertex_id);
        local_vertex_id
            = global_vertex_id - m_body_vertex_id[size_t(rigid_body_id)];
    }

    Eigen::VectorXd RigidBodyAssembler::rb_positions(const bool previous) const
    {
        Eigen::VectorXd X(int(m_rbs.size()) * 3);
        for (size_t i = 0; i < m_rbs.size(); ++i) {
            auto& rb = m_rbs[i];
            X.segment(int(3 * i), 3)
                = previous ? rb.position_prev : rb.position;
        }
        return X;
    }

    void RigidBodyAssembler::set_rb_positions(const Eigen::VectorXd& x)
    {
        for (size_t i = 0; i < m_rbs.size(); ++i) {
            auto& rb = m_rbs[i];
            rb.position = x.segment(int(3 * i), 3);
        }
    }

    Eigen::MatrixXd RigidBodyAssembler::world_vertices(
        const RigidBody::Step step) const
    {
        Eigen::MatrixXd V(num_vertices(), 2);
        for (size_t i = 0; i < m_rbs.size(); ++i) {
            auto& rb = m_rbs[i];
            V.block(m_body_vertex_id[i], 0, rb.vertices.rows(), 2)
                = rb.world_vertices(step);
        }
        return V;
    }

    Eigen::MatrixXd RigidBodyAssembler::world_vertices(
        const Eigen::VectorXd& positions) const
    {
        Eigen::MatrixXd V(num_vertices(), 2);
        for (size_t i = 0; i < m_rbs.size(); ++i) {
            auto& rb = m_rbs[i];
            Eigen::Vector3d p_i = positions.segment(3 * int(i), 3);
            V.block(m_body_vertex_id[i], 0, rb.vertices.rows(), 2)
                = rb.world_vertices(p_i);
        }
        return V;
    }

    Eigen::MatrixXd RigidBodyAssembler::world_velocities() const
    {
        Eigen::MatrixXd V(num_vertices(), 2);
        for (size_t i = 0; i < m_rbs.size(); ++i) {
            auto& rb = m_rbs[i];
            V.block(m_body_vertex_id[i], 0, rb.vertices.rows(), 2)
                = rb.world_velocities();
        }
        return V;
    }

} // namespace physics
} // namespace ccd
