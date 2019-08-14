#include "volume_rb_problem.hpp"

#include <autogen/collision_volume.hpp>
#include <ccd/time_of_impact.hpp>

#include <logger.hpp>

namespace ccd {

namespace opt {
    VolumeRBProblem::VolumeRBProblem(const std::string& name)
        : RigidBodyProblem(name)
    {
    }

    void VolumeRBProblem::settings(const nlohmann::json& params)
    {
        constraint_.settings(params["volume_constraint"]);
        opt_solver_.settings(params["ncp_solver"]);
        opt_solver_.set_problem(*this);
        RigidBodyProblem::settings(params["rigid_body_problem"]);
    }

    Eigen::VectorXd VolumeRBProblem::eval_g(const Eigen::VectorXd& sigma)
    {
        Eigen::MatrixXd xk = m_assembler.world_vertices(sigma);
        Eigen::MatrixXd uk = xk - vertices_t0;
        constraint_.update_collision_set(uk);

        Eigen::VectorXd g_uk;
        constraint_.compute_constraints(uk, g_uk);
        return g_uk;
    }

    void VolumeRBProblem::eval_g(const Eigen::VectorXd& sigma,
        Eigen::VectorXd& g_uk,
        Eigen::SparseMatrix<double>& g_uk_jacobian,
        Eigen::VectorXi& g_uk_active)
    {
        Eigen::MatrixXd xk = m_assembler.world_vertices(sigma);
        Eigen::MatrixXd uk = xk - vertices_t0;
        constraint_.update_collision_set(uk);

        constraint_.compute_constraints(uk, g_uk);
        constraint_.dense_indices(g_uk_active);
        eval_jac_g_core(sigma, g_uk_jacobian);

    }

    void VolumeRBProblem::eval_jac_g_core(
        const Eigen::VectorXd& sigma, Eigen::SparseMatrix<double>& jac_gx)
    {
        typedef AutodiffType<6> Diff;
        Diff::activate();

        int num_edges = m_assembler.num_edges();

        typedef Eigen::Triplet<double> M;
        std::vector<M> triplets;

        for (size_t i = 0; i < constraint_.ee_impacts.size(); ++i) {
            auto& ee_impact = constraint_.ee_impacts[i];

            Eigen::Vector2i e_ij
                = m_assembler.m_edges.row(ee_impact.impacted_edge_index);
            Eigen::Vector2i e_kl
                = m_assembler.m_edges.row(ee_impact.impacting_edge_index);
            int body_ij_id, body_kl_id, local_i_id, local_j_id, local_k_id,
                local_l_id;

            m_assembler.global_to_local(e_ij(0), body_ij_id, local_i_id);
            m_assembler.global_to_local(e_ij(1), body_ij_id, local_j_id);
            m_assembler.global_to_local(e_kl(0), body_kl_id, local_k_id);
            m_assembler.global_to_local(e_kl(1), body_kl_id, local_l_id);

            Diff::D1Vector3d position_ij, position_kl;
            position_ij = Diff::d1vars(0, sigma.segment(3 * body_ij_id, 3));
            position_kl = Diff::d1vars(3, sigma.segment(3 * body_kl_id, 3));

            const auto& rbs = m_assembler.m_rbs;
            // _differentiable_ final positions
            Diff::D1Vector2d v1i, v1j, v1k, v1l;
            v1i = rbs[size_t(body_ij_id)].world_vertex<Diff::DDouble1>(
                position_ij, local_i_id);
            v1j = rbs[size_t(body_ij_id)].world_vertex<Diff::DDouble1>(
                position_ij, local_j_id);
            v1k = rbs[size_t(body_kl_id)].world_vertex<Diff::DDouble1>(
                position_kl, local_k_id);
            v1l = rbs[size_t(body_kl_id)].world_vertex<Diff::DDouble1>(
                position_kl, local_l_id);

            // initial positions
            Eigen::Vector2d v0i, v0j, v0k, v0l;
            v0i = vertices_t0.row(e_ij(0));
            v0j = vertices_t0.row(e_ij(1));
            v0k = vertices_t0.row(e_kl(0));
            v0l = vertices_t0.row(e_kl(1));

            // displacements
            Diff::D1Vector2d ui, uj, uk, ul;
            ui = v1i - v0i.cast<Diff::DDouble1>();
            uj = v1j - v0j.cast<Diff::DDouble1>();
            uk = v1k - v0k.cast<Diff::DDouble1>();
            ul = v1l - v0l.cast<Diff::DDouble1>();

            // impacting vertice
            Eigen::Vector2d v0;
            Diff::D1Vector2d u;
            if (ee_impact.impacting_node() == 0) {
                v0 = v0k;
                u = uk;
            } else {
                v0 = v0l;
                u = ul;
            }
            Diff::DDouble1 toi, alpha_ij, alpha_kl;
            alpha_kl = Diff::DDouble1(ee_impact.impacting_node());

            // get toi and alpha
            bool success;
            success = ccd::autodiff::compute_edge_vertex_time_of_impact<
                Diff::DDouble1>(v0i, v0j, v0, ui, uj, u, toi);
            success = success
                && ccd::autodiff::temporal_parameterization_to_spatial<
                       Diff::DDouble1>(v0i, v0j, v0, ui, uj, u, toi, alpha_ij);
            Diff::DDouble1 v_ij(0), v_kl(0);
            if (success) {
                v_ij
                    = ccd::autogen::space_time_collision_volume<Diff::DDouble1>(
                        v0i, v0j, ui, uj, toi, alpha_ij,
                        constraint_.volume_epsilon);
                v_kl
                    = ccd::autogen::space_time_collision_volume<Diff::DDouble1>(
                        v0k, v0l, uk, ul, toi, alpha_kl,
                        constraint_.volume_epsilon);
            }
            long c_ij
                = get_constraint_index(ee_impact, /*impacted=*/true, num_edges);
            long c_kl = get_constraint_index(
                ee_impact, /*impacted=*/false, num_edges);

            Eigen::VectorXd gradient_ij = v_ij.getGradient();
            Eigen::VectorXd gradient_kl = v_kl.getGradient();

            // 3 * rbc.vertex_body_id + dim
            for (size_t dim = 0; dim < 3; dim++) {
                triplets.push_back(
                    M(int(c_ij), 3 * body_ij_id + dim, gradient_ij(dim)));
                triplets.push_back(M(
                    int(c_ij), 3 * body_kl_id + dim, gradient_ij(3 + dim)));

                triplets.push_back(
                    M(int(c_kl), 3 * body_ij_id + dim, gradient_kl(dim)));
                triplets.push_back(M(
                    int(c_kl), 3 * body_kl_id + dim, gradient_kl(3 + dim)));
            }
        }
        const long num_constr = get_constraints_size(num_edges);
        jac_gx.resize(int(num_constr), int(sigma.size()));
        jac_gx.setFromTriplets(triplets.begin(), triplets.end());
    }

} // namespace opt
} // namespace ccd
