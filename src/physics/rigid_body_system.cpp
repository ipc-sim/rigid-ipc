#include <physics/rigid_body_system.hpp>

#include <Eigen/Geometry>

#include <logger.hpp>
#include <utils/eigen_ext.hpp>

#include <autodiff/finitediff.hpp>
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
            m_body_vertex_id[i + 1]
                = m_body_vertex_id[i] + rb.vertices.rows();
            m_body_edge_id[i + 1] = m_body_edge_id[i] + rb.edges.rows();
        }

        // global edges
        m_edges.resize(m_body_edge_id.back(), 2);
        for (size_t i = 0; i < num_bodies; ++i) {
            auto& rb = rigid_bodies[i];
            m_edges.block(m_body_edge_id[i], 0, rb.edges.rows(), 2)
                = (rb.edges.array() + int(m_body_vertex_id[i]));
        }

        // mass-matrix
        Eigen::VectorXd mass_vector(int(num_bodies) * 3);
        for (int i = 0; i < int(num_bodies); ++i) {
            auto& rb = rigid_bodies[size_t(i)];
            mass_vector[3 * i + 0] = rb.mass;
            mass_vector[3 * i + 1] = rb.mass;
            mass_vector[3 * i + 2] = rb.moment_of_inertia;
        }
        m_mass_matrix = Eigen::SparseDiagonal<double>(mass_vector);
    }

    Eigen::VectorXd RigidBodyAssembler::rb_positions(const bool previous)
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

    Eigen::MatrixXd RigidBodyAssembler::world_vertices(const bool previous)
    {
        Eigen::MatrixXd V(num_vertices(), 2);
        for (size_t i = 0; i < m_rbs.size(); ++i) {
            auto& rb = m_rbs[i];
            V.block(m_body_vertex_id[i], 0, rb.vertices.rows(), 2)
                = rb.world_vertices(previous);
        }
        return V;
    }

    Eigen::MatrixXd RigidBodyAssembler::world_vertices(
        const Eigen::VectorXd& positions)
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

    void RigidBodyAssembler::world_vertices_gradient(
        const Eigen::VectorXd& positions, Eigen::SparseMatrix<double>& grad_u)
    {
        typedef Eigen::Triplet<double> M;
        std::vector<M> triplets;

        long num_vertices = m_body_vertex_id.back();
        triplets.reserve(size_t(num_vertices) * 2);

        grad_u.resize(int(num_vertices * 2), int(m_rbs.size()) * 3);
        for (size_t i = 0; i < m_rbs.size(); ++i) {
            auto& rb = m_rbs[i];
            Eigen::Vector3d p_i = positions.segment(3 * int(i), 3);
            Eigen::MatrixXd el_grad = rb.world_vertices_gradient(p_i);

            long d = el_grad.rows() / 2;
            // x-axis entries
            for (int j = 0; j < d; ++j) {
                for (int k = 0; k < el_grad.cols(); ++k) {
                    triplets.push_back(M(int(m_body_vertex_id[i]) + j,
                        int(3 * i) + k, el_grad(j, k)));
                }
            }
            // y-axis entries
            for (int j = 0; j < d; ++j) {
                for (int k = 0; k < el_grad.cols(); ++k) {
                    triplets.push_back(
                        M(int(m_body_vertex_id[i] + num_vertices) + j,
                            int(3 * i) + k, el_grad(d + j, k)));
                }
            }
        }

        grad_u.setFromTriplets(triplets.begin(), triplets.end());
    }

    void RigidBodyAssembler::world_vertices_hessian(const Eigen::VectorXd& v,
        std::vector<Eigen::SparseMatrix<double>>& hess_u)
    {
        // for local entry hessian
        typedef Eigen::Triplet<double> M;
        std::array<M, 9> triplets;

        long num_vertices = m_body_vertex_id.back();
        hess_u.resize(size_t(num_vertices) * 2);

        for (size_t i = 0; i < m_rbs.size(); ++i) {
            auto& rb = m_rbs[i];
            std::vector<Eigen::Matrix3d> el_hessians
                = rb.world_displacements_hessian(v.segment(int(3 * i), 3));
            for (size_t ii = 0; ii < el_hessians.size(); ++ii) {
                auto& el_hess = el_hessians[ii];
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        triplets[size_t(3 * j + k)]
                            = M(int(3 * i) + j, int(3 * i) + k, el_hess(j, k));
                    }
                }
                Eigen::SparseMatrix<double> global_el_hessian(
                    int(m_rbs.size()) * 3, int(m_rbs.size()) * 3);
                global_el_hessian.setFromTriplets(
                    triplets.begin(), triplets.end());

                int d = int(el_hessians.size() / 2);
                int padding
                    = int(ii) >= d ? int(ii) - d + int(num_vertices) : int(ii);
                hess_u[size_t(m_body_vertex_id[i] + padding)]
                    = global_el_hessian;
            }
        }
    }

    // ------------------------------------------------------------------------
    // State Functions
    // ------------------------------------------------------------------------
    void RigidBodySystem::clear()
    {
        rigid_bodies.clear();
        acc_vertex_id.clear();
        acc_edge_id.clear();
    }

    // ------------------------------------------------------------------------
    // CCD Functions
    // ------------------------------------------------------------------------
    void RigidBodySystem::assemble()
    {
        size_t num_bodies = rigid_bodies.size();

        velocities.resize(int(3 * num_bodies));

        acc_vertex_id.clear();
        acc_edge_id.clear();
        acc_vertex_id.resize(num_bodies + 1);
        acc_edge_id.resize(num_bodies + 1);
        acc_vertex_id[0] = 0;
        acc_edge_id[0] = 0;
        for (size_t i = 0; i < num_bodies; ++i) {
            auto& rb = rigid_bodies[i];
            velocities.segment(int(3 * i), 3) = rb.velocity;
            acc_vertex_id[i + 1] = acc_vertex_id[i] + rb.vertices.rows();
            acc_edge_id[i + 1] = acc_edge_id[i] + rb.edges.rows();
        }

        vertex_to_body_map.resize(acc_vertex_id.back());
        vertices.resize(acc_vertex_id.back(), 2);
        edges.resize(acc_edge_id.back(), 2);
        for (size_t i = 0; i < num_bodies; ++i) {
            auto& rb = rigid_bodies[i];
            vertices.block(acc_vertex_id[i], 0, rb.vertices.rows(), 2)
                = rb.world_vertices();
            vertex_to_body_map.segment(acc_vertex_id[i], rb.vertices.rows())
                .setConstant(int(i));
            edges.block(acc_edge_id[i], 0, rb.edges.rows(), 2)
                = (rb.edges.array() + int(acc_vertex_id[i]));
        }

        assemble_displacements();
    }

    void RigidBodySystem::set_velocity(
        const size_t rb_id, const Eigen::Vector3d vel)
    {
        rigid_bodies[rb_id].velocity = vel;
        velocities.segment(int(3 * rb_id), 3) = vel;
    }

    void RigidBodySystem::assemble_displacements()
    {
        compute_displacements(velocities, displacements);
    }

    void RigidBodySystem::compute_displacements(
        const Eigen::VectorXd& v, Eigen::MatrixXd& u)
    {
        u.resize(acc_vertex_id.back(), 2);
        for (size_t i = 0; i < rigid_bodies.size(); ++i) {
            auto& rb = rigid_bodies[i];
            u.block(acc_vertex_id[i], 0, rb.vertices.rows(), 2)
                = rb.world_displacements<double>(v.segment(int(3 * i), 3));
        }
    }

    void RigidBodySystem::compute_displacements_gradient(
        const Eigen::VectorXd& v, Eigen::SparseMatrix<double>& grad_u)
    {
        std::vector<Eigen::Triplet<double>> triplets;

        long num_vertices = acc_vertex_id.back();
        triplets.reserve(size_t(num_vertices) * 2);

        for (size_t i = 0; i < rigid_bodies.size(); ++i) {
            auto& rb = rigid_bodies[i];
            Eigen::MatrixXd el_grad
                = rb.world_displacements_gradient(v.segment(int(3 * i), 3));
#ifdef WITH_DERIVATIVE_CHECK
            auto foo = [&](const Eigen::Vector3d& x) -> Eigen::VectorXd {
                auto U = rb.world_displacements(x);
                flatten(U);
                return U;
            };
            Eigen::MatrixXd approx_el_grad;
            finite_jacobian(v.segment(int(3 * i), 3), foo, approx_el_grad);
            assert(compare_jacobian(el_grad, approx_el_grad));
#endif

            long d = el_grad.rows() / 2;
            // x-axis entries
            for (int j = 0; j < d; ++j) {
                for (int k = 0; k < el_grad.cols(); ++k) {
                    triplets.emplace_back(int(acc_vertex_id[i]) + j,
                        int(3 * i) + k, el_grad(j, k));
                }
            }
            // y-axis entries
            for (int j = 0; j < d; ++j) {
                for (int k = 0; k < el_grad.cols(); ++k) {
                    triplets.emplace_back(
                        int(acc_vertex_id[i] + num_vertices) + j,
                        int(3 * i) + k, el_grad(d + j, k));
                }
            }
        }
        grad_u.resize(int(num_vertices * 2), int(rigid_bodies.size()) * 3);
        grad_u.setFromTriplets(triplets.begin(), triplets.end());
    }

    void RigidBodySystem::compute_displacements_hessian(
        const Eigen::VectorXd& v,
        std::vector<Eigen::SparseMatrix<double>>& hess_u)
    {
        // for local entry hessian
        typedef Eigen::Triplet<double> M;
        std::array<M, 9> triplets;

        long num_vertices = acc_vertex_id.back();
        hess_u.resize(size_t(num_vertices) * 2);

        for (size_t i = 0; i < rigid_bodies.size(); ++i) {
            auto& rb = rigid_bodies[i];
            std::vector<Eigen::Matrix3d> el_hessians
                = rb.world_displacements_hessian(v.segment(int(3 * i), 3));
            for (size_t ii = 0; ii < el_hessians.size(); ++ii) {
                auto& el_hess = el_hessians[ii];
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        triplets[size_t(3 * j + k)]
                            = M(int(3 * i) + j, int(3 * i) + k, el_hess(j, k));
                    }
                }
                Eigen::SparseMatrix<double> global_el_hessian(
                    int(rigid_bodies.size()) * 3, int(rigid_bodies.size()) * 3);
                global_el_hessian.setFromTriplets(
                    triplets.begin(), triplets.end());

                int d = int(el_hessians.size() / 2);
                int padding
                    = int(ii) >= d ? int(ii) - d + int(num_vertices) : int(ii);
                hess_u[size_t(acc_vertex_id[i] + padding)] = global_el_hessian;
            }
        }
    }
} // namespace physics
} // namespace ccd
