#include "volume_constraint.hpp"

#include <iostream>

#include <autogen/collision_volume.hpp>
#include <ccd/prune_impacts.hpp>
#include <ccd/time_of_impact.hpp>

#include <logger.hpp>

namespace ccd {

namespace opt {

    VolumeConstraint::VolumeConstraint()
        : VolumeConstraint("volume_constraint")
    {
    }
    VolumeConstraint::VolumeConstraint(const std::string& name)
        : CollisionConstraint(name)
        , volume_epsilon(1E-3)
    {
    }

    void VolumeConstraint::settings(const nlohmann::json& json)
    {
        CollisionConstraint::settings(json);
        volume_epsilon = json["volume_epsilon"].get<double>();
    }

    nlohmann::json VolumeConstraint::settings() const
    {
        nlohmann::json json = CollisionConstraint::settings();
        json["volume_epsilon"] = volume_epsilon;
        return json;
    }

    void VolumeConstraint::initialize(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2i& edges,
        const Eigen::VectorXi& group_ids,
        const Eigen::MatrixXd& Uk)
    {
        m_ee_impacts.clear();
        m_edge_impact_map.resize(edges.rows());
        m_edge_impact_map.setZero();

        CollisionConstraint::initialize(vertices, edges, group_ids, Uk);
        num_constraints = get_constraints_size(edges.rows());
    }

    void VolumeConstraint::update_collision_set(const Eigen::MatrixXd& Uk)
    {
        CollisionConstraint::update_collision_set(Uk);
        ccd::convert_edge_vertex_to_edge_edge_impacts(
            edges, m_ev_impacts, m_ee_impacts);
        prune_impacts(m_ee_impacts, m_edge_impact_map);
    }

    int VolumeConstraint::number_of_constraints() { return num_constraints; }

    void VolumeConstraint::compute_constraints(
        const Eigen::MatrixXd& Uk, Eigen::VectorXd& g_uk)
    {
        g_uk.resize(num_constraints);
        g_uk.setZero();

        for (size_t i = 0; i < m_ee_impacts.size(); ++i) {
            auto& ee_impact = m_ee_impacts[i];

            Eigen::Vector2i e_ij = edges.row(ee_impact.impacted_edge_index);
            Eigen::Vector2i e_kl = edges.row(ee_impact.impacting_edge_index);

            int node_i, node_j, node_k, node_l;
            node_i = e_ij(0);
            node_j = e_ij(1);
            node_k = e_kl(0);
            node_l = e_kl(1);


            Eigen::VectorXd v_i, v_j, v_k, v_l;
            v_i = vertices.row(node_i);
            v_j = vertices.row(node_j);
            v_k = vertices.row(node_k);
            v_l = vertices.row(node_l);

            Eigen::VectorXd u_i, u_j, u_k, u_l;
            u_i = Uk.row(node_i);
            u_j = Uk.row(node_j);
            u_k = Uk.row(node_k);
            u_l = Uk.row(node_l);

            Eigen::VectorXd v_c, u_c;
            if (ee_impact.impacting_node() == 0) {
                v_c = v_k;
                u_c = u_k;
            } else {
                v_c = v_l;
                u_c = u_l;
            }

            double toi, alpha_ij, alpha_kl;
            alpha_kl = ee_impact.impacting_node();
            // get toi and alpha
            bool success;
            success = ccd::autodiff::compute_edge_vertex_time_of_impact<double>(
                v_i, v_j, v_c, u_i, u_j, u_c, toi);
            success = success
                && ccd::autodiff::temporal_parameterization_to_spatial<double>(
                       v_i, v_j, v_c, u_i, u_j, u_c, toi, alpha_ij);

            double vol_ij(0), vol_kl(0);
            if (success) {
                vol_ij = ccd::autogen::space_time_collision_volume<double>(
                    v_i, v_j, u_i, u_j, toi, alpha_ij, volume_epsilon);
                vol_kl = ccd::autogen::space_time_collision_volume<double>(
                    v_k, v_l, u_k, u_l, toi, alpha_kl, volume_epsilon);
            }
            long c_ij = get_constraint_index(
                ee_impact, /*impacted=*/true, edges.rows());
            long c_kl = get_constraint_index(
                ee_impact, /*impacted=*/false, edges.rows());

            g_uk(c_ij) = vol_ij;
            g_uk(c_kl) = vol_kl;
        }
    }

    void VolumeConstraint::compute_constraints_jacobian(
        const Eigen::MatrixXd& Uk, Eigen::SparseMatrix<double>& jac_uk)
    {

        typedef Eigen::Triplet<double> M;
        std::vector<M> triplets;

        typedef AutodiffType<8> Diff;
        Diff::activate();

        for (size_t i = 0; i < m_ee_impacts.size(); ++i) {
            auto& ee_impact = m_ee_impacts[i];

            Eigen::Vector2i e_ij = edges.row(ee_impact.impacted_edge_index);
            Eigen::Vector2i e_kl = edges.row(ee_impact.impacting_edge_index);

            int node_i, node_j, node_k, node_l;
            node_i = e_ij(0);
            node_j = e_ij(1);
            node_k = e_kl(0);
            node_l = e_kl(1);

            Eigen::VectorXd v_i, v_j, v_k, v_l;
            v_i = vertices.row(node_i);
            v_j = vertices.row(node_j);
            v_k = vertices.row(node_k);
            v_l = vertices.row(node_l);

            Diff::D1Vector2d u_i, u_j, u_k, u_l;
            u_i = Diff::d1vars(0, Uk.row(node_i));
            u_j = Diff::d1vars(2, Uk.row(node_j));
            u_k = Diff::d1vars(4, Uk.row(node_k));
            u_l = Diff::d1vars(6, Uk.row(node_l));

            Eigen::VectorXd v_c;
            Diff::D1Vector2d u_c;
            if (ee_impact.impacting_node() == 0) {
                v_c = v_k;
                u_c = u_k;
            } else {
                v_c = v_l;
                u_c = u_l;
            }

            Diff::DDouble1 toi, alpha_ij, alpha_kl;
            alpha_kl = ee_impact.impacting_node();
            // get toi and alpha
            bool success;
            success = ccd::autodiff::compute_edge_vertex_time_of_impact<
                Diff::DDouble1>(v_i, v_j, v_c, u_i, u_j, u_c, toi);
            success = success
                && ccd::autodiff::temporal_parameterization_to_spatial<
                       Diff::DDouble1>(
                       v_i, v_j, v_c, u_i, u_j, u_c, toi, alpha_ij);

            Diff::DDouble1 vol_ij(0), vol_kl(0);
            if (success) {
                vol_ij
                    = ccd::autogen::space_time_collision_volume<Diff::DDouble1>(
                        v_i, v_j, u_i, u_j, toi, alpha_ij, volume_epsilon);
                vol_kl
                    = ccd::autogen::space_time_collision_volume<Diff::DDouble1>(
                        v_k, v_l, u_k, u_l, toi, alpha_kl, volume_epsilon);
            }
            long c_ij = get_constraint_index(
                ee_impact, /*impacted=*/true, edges.rows());
            long c_kl = get_constraint_index(
                ee_impact, /*impacted=*/false, edges.rows());

            // local gradients
            Eigen::VectorXd grad_vol_ij = vol_ij.getGradient();
            Eigen::VectorXd grad_vol_kl = vol_kl.getGradient();

            // x entries:
            int nodes[4] = { node_i, node_j, node_k, node_l };
            for (int n_id = 0; n_id < 4; ++n_id) {
                // x and y entries
                triplets.emplace_back(
                    int(c_ij), nodes[n_id], grad_vol_ij(2 * n_id));
                triplets.emplace_back(int(c_ij), nodes[n_id] + vertices.rows(),
                    grad_vol_ij(2 * n_id + 1));

                // x and y entries
                triplets.emplace_back(
                    int(c_kl), nodes[n_id], grad_vol_kl(2 * n_id));
                triplets.emplace_back(int(c_kl), nodes[n_id] + vertices.rows(),
                    grad_vol_kl(2 * n_id + 1));
            }
        }

        jac_uk.resize(int(num_constraints), int(vertices.size()));
        jac_uk.setFromTriplets(triplets.begin(), triplets.end());
    }

    void VolumeConstraint::compute_constraints(const Eigen::MatrixXd& Uk,
        Eigen::VectorXd& g_uk,
        Eigen::SparseMatrix<double>& g_uk_jacobian,
        Eigen::VectorXi& g_uk_active)
    {
        compute_constraints(Uk, g_uk);
        compute_constraints_jacobian(Uk, g_uk_jacobian);
        dense_indices(g_uk_active);
    }

    void VolumeConstraint::dense_indices(Eigen::VectorXi& dense_indices)
    {
        dense_indices.resize(int(m_ee_impacts.size()) * 2);

        const int num_edges = int(edges.rows());

        for (size_t ee = 0; ee < m_ee_impacts.size(); ++ee) {
            auto& ee_impact = m_ee_impacts[ee];

            long c_ij
                = get_constraint_index(ee_impact, /*impacted=*/true, num_edges);
            long c_kl = get_constraint_index(
                ee_impact, /*impacted=*/false, num_edges);

            dense_indices(int(2 * ee) + 0) = int(c_ij);
            dense_indices(int(2 * ee) + 1) = int(c_kl);
        }
    }

    long get_constraint_index(
        const EdgeEdgeImpact& impact, const bool impacted, const int num_edges)
    {
        long e1 = impact.impacted_edge_index;
        long e2 = impact.impacting_edge_index;
        int p = impact.impacting_alpha > 0.5 ? 1 : 0;
        int q = impacted ? 0 : 1;

        // unravel index Q * P * E2 * e1 + Q * P * e2 + Q * p + q
        const int Q = 2, P = 2, E2 = num_edges;
        return Q * P * E2 * e1 + Q * P * e2 + Q * p + q;
    }

    long get_constraints_size(const int num_edges)
    {
        const int Q = 2, P = 2, E2 = num_edges, E1 = num_edges;
        return Q * P * E2 * E1;
    }
} // namespace opt
} // namespace ccd
