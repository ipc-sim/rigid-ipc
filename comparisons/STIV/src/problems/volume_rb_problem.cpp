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
        Eigen::VectorXd qk = sigma;
        Eigen::MatrixXd uk = m_assembler.world_vertices(qk) - vertices_t0;

        Eigen::VectorXd g_uk;
        constraint_.compute_constraints(uk, g_uk);
        return g_uk;
    }

    void VolumeRBProblem::eval_g(
        const Eigen::VectorXd& sigma,
        Eigen::VectorXd& g_uk,
        Eigen::MatrixXd& g_uk_jacobian)
    {
        Eigen::VectorXd qk = sigma;
        Eigen::MatrixXd uk = m_assembler.world_vertices(qk) - vertices_t0;
        auto ee_impacts = constraint_.get_ee_collision_set(uk);

        constraint_.compute_constraints(uk, ee_impacts, g_uk);
        eval_jac_g_core(sigma, ee_impacts, g_uk_jacobian);
    }

    void VolumeRBProblem::eval_g_normal(
        const Eigen::VectorXd& sigma,
        Eigen::VectorXd& g_uk,
        Eigen::MatrixXd& g_uk_jacobian)
    {
        Eigen::VectorXd qk = sigma;
        Eigen::MatrixXd uk = m_assembler.world_vertices(qk) - vertices_t0;
        auto ee_impacts = constraint_.get_ee_collision_set(uk);

        constraint_.compute_constraints(uk, ee_impacts, g_uk);

        Eigen::SparseMatrix<double> jac_xk_sigma;
        m_assembler.world_vertices_gradient(sigma, jac_xk_sigma);

        Eigen::MatrixXd jac_g_uk;
        constraint_.compute_constraints_normals(uk, ee_impacts, jac_g_uk);

        g_uk_jacobian = jac_g_uk * jac_xk_sigma;
    }

    void VolumeRBProblem::eval_jac_g_core(
        const Eigen::VectorXd& sigma,
        const EdgeEdgeImpacts& ee_impacts,
        Eigen::MatrixXd& jac_gx)
    {
        typedef AutodiffType<6> Diff;
        Diff::activate();

        jac_gx.resize(int(ee_impacts.size()) * 2, int(sigma.size()));
        jac_gx.setZero();

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

            Eigen::Vector2i e_ij =
                m_assembler.m_edges.row(ee_impact.impacted_edge_index);
            Eigen::Vector2i e_kl =
                m_assembler.m_edges.row(ee_impact.impacting_edge_index);
            long body_ij_id, body_kl_id, local_i_id, local_j_id, local_k_id,
                local_l_id;

            m_assembler.global_to_local_vertex(e_ij(0), body_ij_id, local_i_id);
            m_assembler.global_to_local_vertex(e_ij(1), body_ij_id, local_j_id);
            m_assembler.global_to_local_vertex(e_kl(0), body_kl_id, local_k_id);
            m_assembler.global_to_local_vertex(e_kl(1), body_kl_id, local_l_id);

            Diff::D1Vector3d sigma_ij, sigma_kl;
            sigma_ij = Diff::d1vars(0, sigma.segment(3 * body_ij_id, 3));
            sigma_kl = Diff::d1vars(3, sigma.segment(3 * body_kl_id, 3));

            Diff::D1Vector3d pose_ij, pose_kl;

            pose_ij = sigma_ij;

            pose_kl = sigma_kl;

            const auto& rbs = m_assembler.m_rbs;
            // _differentiable_ final poses
            Diff::D1Vector2d v1i, v1j, v1k, v1l;
            v1i = rbs[size_t(body_ij_id)].world_vertex<Diff::DDouble1>(
                pose_ij, local_i_id);
            v1j = rbs[size_t(body_ij_id)].world_vertex<Diff::DDouble1>(
                pose_ij, local_j_id);
            v1k = rbs[size_t(body_kl_id)].world_vertex<Diff::DDouble1>(
                pose_kl, local_k_id);
            v1l = rbs[size_t(body_kl_id)].world_vertex<Diff::DDouble1>(
                pose_kl, local_l_id);

            // initial poses
            Eigen::Vector2d v_i, v_j, v_k, v_l;
            v_i = vertices_t0.row(e_ij(0));
            v_j = vertices_t0.row(e_ij(1));
            v_k = vertices_t0.row(e_kl(0));
            v_l = vertices_t0.row(e_kl(1));

            // displacements
            Diff::D1Vector2d u_i, u_j, u_k, u_l;
            u_i = v1i - v_i;
            u_j = v1j - v_j;
            u_k = v1k - v_k;
            u_l = v1l - v_l;

            // impacting vertice
            Eigen::Vector2d v_c;
            Diff::D1Vector2d u_c;
            if (ee_impact.impacting_node() == 0) {
                v_c = v_k;
                u_c = u_k;
            } else {
                v_c = v_l;
                u_c = u_l;
            }
            Diff::DDouble1 toi, alpha_ij, alpha_kl;
            alpha_kl = Diff::DDouble1(ee_impact.impacting_node());

            // get toi and alpha
            bool success;
            success = ccd::autodiff::compute_edge_vertex_time_of_impact<
                Diff::DDouble1>(v_i, v_j, v_c, u_i, u_j, u_c, toi);
            success = success
                && ccd::autodiff::temporal_parameterization_to_spatial<
                          Diff::DDouble1>(
                          v_i, v_j, v_c, u_i, u_j, u_c, toi, alpha_ij);

            Diff::D1Vector2d avg_u = (u_i + alpha_ij * (u_j - u_i) - u_k);
            Diff::DDouble1 avg_d = (avg_u).norm();

            // toi =  toi - Diff::DDouble1(constraint_.time_epsilon);
            toi = toi * (1.0 - constraint_.time_epsilon / avg_d);
            if (toi < 0)
                toi = Diff::DDouble1(0);

            Diff::DDouble1 vol_ij(0), vol_kl(0);
            if (success) {
                vol_ij =
                    ccd::autogen::space_time_collision_volume<Diff::DDouble1>(
                        v_i, v_j, u_i, u_j, toi, alpha_ij,
                        constraint_.volume_epsilon);
                vol_kl =
                    ccd::autogen::space_time_collision_volume<Diff::DDouble1>(
                        v_k, v_l, u_k, u_l, toi, alpha_kl,
                        constraint_.volume_epsilon);
            }
            long c_ij = 2 * int(i) + 0;
            long c_kl = 2 * int(i) + 1;

            Eigen::VectorXd gradient_ij = vol_ij.getGradient();
            Eigen::VectorXd gradient_kl = vol_kl.getGradient();

            // 3 * rbc.vertex_body_id + dim
            for (int dim = 0; dim < 3; dim++) {
                jac_gx(int(c_ij), 3 * body_ij_id + dim) = gradient_ij(dim);
                jac_gx(int(c_ij), 3 * body_kl_id + dim) = gradient_ij(3 + dim);

                jac_gx(int(c_kl), 3 * body_ij_id + dim) = gradient_kl(dim);
                jac_gx(int(c_kl), 3 * body_kl_id + dim) = gradient_kl(3 + dim);
            }
        }
    }

} // namespace opt
} // namespace ccd
