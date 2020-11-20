#include "volume_constraint.hpp"

#include <iostream>

#include <autogen/collision_volume.hpp>
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
        , time_epsilon(1e-16)
    {
    }

    void VolumeConstraint::settings(const nlohmann::json& json)
    {
        CollisionConstraint::settings(json);
        volume_epsilon = json["volume_epsilon"].get<double>();
        time_epsilon = json["time_epsilon"].get<double>();
    }

    nlohmann::json VolumeConstraint::settings() const
    {
        nlohmann::json json = CollisionConstraint::settings();
        json["volume_epsilon"] = volume_epsilon;
        json["time_epsilon"] = time_epsilon;
        return json;
    }

    std::vector<EdgeVertexImpact> VolumeConstraint::initialize(
        const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2i& edges,
        const Eigen::VectorXi& group_ids,
        const Eigen::MatrixXd& Uk)
    {
        m_edge_impact_map.resize(edges.rows());
        m_edge_impact_map.setZero();

        num_constraints = get_constraints_size(edges.rows());
        return CollisionConstraint::initialize();
    }

    /**
     * Prune the impacts to only include the earliest impact for each edge.
     *
     * @param all_impacts All of the impacts to be pruned.
     * @param pruned_impact_indices A vector of lengeth the numbe of edges.
     * Stores the index of the earliest impact in this vector.
     * @return The number of edges that have an impact.
     */
    int prune_impacts(
        const std::vector<EdgeEdgeImpact>& all_impacts,
        Eigen::VectorXi& pruned_impact_indices)
    {
        // An index value of -1 indicates no impact for that edge
        pruned_impact_indices.setConstant(-1);

        // Loop over all impacts and add them to the prunced impacts if they are
        // the earliest.
        int num_pruned_impacts = 0;
        for (size_t i = 0; i < all_impacts.size(); i++) {
            EdgeEdgeImpact ee_impact = all_impacts[i];
            for (int index : { ee_impact.impacted_edge_index,
                               ee_impact.impacting_edge_index }) {
                // If the edge is not in the map or the current impact happens
                // later replace it in the pruned impacts.
                if (pruned_impact_indices[index] == -1
                    || all_impacts[size_t(pruned_impact_indices[index])].time
                        > ee_impact.time) {
                    num_pruned_impacts +=
                        int(pruned_impact_indices[index] == -1);
                    pruned_impact_indices[index] = int(i);
                }
            }
        }

        return num_pruned_impacts;
    }

    std::vector<EdgeEdgeImpact>
    VolumeConstraint::get_ee_collision_set(const Eigen::MatrixXd& Uk)
    {
        auto ev_impacts = get_collision_set(Uk);

        std::vector<EdgeEdgeImpact> ee_impacts;
        ccd::convert_edge_vertex_to_edge_edge_impacts(
            edges, ev_impacts, ee_impacts);
        prune_impacts(ee_impacts, m_edge_impact_map);
        return ee_impacts;
    }

    int VolumeConstraint::number_of_constraints() { return num_constraints; }

    void VolumeConstraint::compute_constraints(
        const Eigen::MatrixXd& Uk, Eigen::VectorXd& g_uk)
    {
        std::vector<EdgeEdgeImpact> ee_impacts = get_ee_collision_set(Uk);
        compute_constraints(Uk, ee_impacts, g_uk);
    }

    void VolumeConstraint::compute_constraints(
        const Eigen::MatrixXd& Uk,
        Eigen::VectorXd& g_uk,
        Eigen::MatrixXd& g_uk_jacobian)
    {
        std::vector<EdgeEdgeImpact> ee_impacts = get_ee_collision_set(Uk);
        compute_constraints(Uk, ee_impacts, g_uk, g_uk_jacobian);
    }

    void VolumeConstraint::compute_constraints(
        const Eigen::MatrixXd& Uk,
        const std::vector<EdgeEdgeImpact>& ee_impacts,
        Eigen::VectorXd& g_uk)
    {
        //        g_uk.resize(num_constraints);
        g_uk.resize(int(ee_impacts.size()) * 2);
        g_uk.setZero();

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

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

            double toi, alpha_ij, alpha_kl = ee_impact.impacting_node();
            alpha_kl = ee_impact.impacting_node();
            // get toi and alpha
            bool success;
            success = ccd::autodiff::compute_edge_vertex_time_of_impact<double>(
                v_i, v_j, v_c, u_i, u_j, u_c, toi);
            success = success
                && ccd::autodiff::temporal_parameterization_to_spatial<double>(
                          v_i, v_j, v_c, u_i, u_j, u_c, toi, alpha_ij);

            // move toi to add some fake min-distance
            Eigen::Vector2d avg_u = (u_i + alpha_ij * (u_j - u_i) - u_k);
            double avg_d = (avg_u).norm();
            ;
            toi = toi * (1.0 - time_epsilon / avg_d);
            toi = std::max(0.0, toi);

            double vol_ij(0), vol_kl(0);
            if (success) {
                vol_ij = ccd::autogen::space_time_collision_volume<double>(
                    v_i, v_j, u_i, u_j, toi, alpha_ij, volume_epsilon);
                vol_kl = ccd::autogen::space_time_collision_volume<double>(
                    v_k, v_l, u_k, u_l, toi, alpha_kl, volume_epsilon);
            }

            // long c_ij = get_constraint_index(
            //     ee_impact, /*impacted=*/true, edges.rows());
            // long c_kl = get_constraint_index(
            //     ee_impact, /*impacted=*/false, edges.rows());

            int c_ij = 2 * int(i) + 0;
            int c_kl = 2 * int(i) + 1;
            g_uk(c_ij) = vol_ij;
            g_uk(c_kl) = vol_kl;
        }
    }

    void VolumeConstraint::compute_constraints_jacobian(
        const Eigen::MatrixXd& Uk,
        const std::vector<EdgeEdgeImpact>& ee_impacts,
        Eigen::MatrixXd& jac_uk)
    {
        jac_uk.resize(int(ee_impacts.size()) * 2, int(vertices.size()));

        typedef AutodiffType<8> Diff;
        Diff::activate();

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

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

            Diff::D1Vector2d avg_u = (u_i + alpha_ij * (u_j - u_i) - u_k);
            Diff::DDouble1 avg_d = (avg_u).norm();
            toi = toi * (1.0 - time_epsilon / avg_d);
            if (toi < 0)
                toi = Diff::DDouble1(0);

            Diff::DDouble1 vol_ij(0), vol_kl(0);
            if (success) {
                vol_ij =
                    ccd::autogen::space_time_collision_volume<Diff::DDouble1>(
                        v_i, v_j, u_i, u_j, toi, alpha_ij, volume_epsilon);
                vol_kl =
                    ccd::autogen::space_time_collision_volume<Diff::DDouble1>(
                        v_k, v_l, u_k, u_l, toi, alpha_kl, volume_epsilon);
            }

            // local gradients
            Eigen::VectorXd grad_vol_ij = vol_ij.getGradient();
            Eigen::VectorXd grad_vol_kl = vol_kl.getGradient();

            // x entries:
            int nodes[4] = { node_i, node_j, node_k, node_l };
            int c_ij = 2 * int(i) + 0;
            int c_kl = 2 * int(i) + 1;
            for (int n_id = 0; n_id < 4; ++n_id) {
                // x and y entries
                jac_uk(int(c_ij), nodes[n_id]) = grad_vol_ij(2 * n_id);
                jac_uk(int(c_ij), nodes[n_id] + vertices.rows()) =
                    grad_vol_ij(2 * n_id + 1);

                // x and y entries
                jac_uk(int(c_kl), nodes[n_id]) = grad_vol_kl(2 * n_id);
                jac_uk(int(c_kl), nodes[n_id] + vertices.rows()) =
                    grad_vol_kl(2 * n_id + 1);
            }
        }
    }

    void VolumeConstraint::compute_constraints_normals(
        const Eigen::MatrixXd& Uk,
        const std::vector<EdgeEdgeImpact>& ee_impacts,
        Eigen::MatrixXd& jac_uk)
    {
        jac_uk.resize(int(ee_impacts.size()) * 2, int(vertices.size()));
        jac_uk.setZero();

        typedef AutodiffType<8> Diff;
        Diff::activate();

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

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

            double toi = 0, alpha_ij = 0, alpha_kl = ee_impact.impacting_node();
            alpha_kl = ee_impact.impacting_node();
            // get toi and alpha
            bool success;
            success = ccd::autodiff::compute_edge_vertex_time_of_impact<double>(
                v_i, v_j, v_c, u_i, u_j, u_c, toi);
            success = success
                && ccd::autodiff::temporal_parameterization_to_spatial<double>(
                          v_i, v_j, v_c, u_i, u_j, u_c, toi, alpha_ij);

            // move toi to add some fake min-distance
            Eigen::Vector2d avg_u = (u_i + alpha_ij * (u_j - u_i) - u_k);
            double avg_d = (avg_u).norm();
            ;
            toi = toi * (1.0 - time_epsilon / avg_d);
            toi = std::max(0.0, toi);

            Eigen::VectorXd grad_vol_ij = Eigen::VectorXd::Zero(8);
            Eigen::VectorXd grad_vol_kl = Eigen::VectorXd::Zero(8);
            if (success) {
                Eigen::Vector2d e((v_j + u_j * toi) - (v_i + u_i * toi));
                // clockwise rotation works for solid objects!!!
                Eigen::Vector2d n(e[1], -e[0]);
                n.normalize();

                Eigen::VectorXd u_e = u_i + alpha_ij * (u_j - u_i);
                if ((u_c - u_e).dot(n) > 0) {
                    n *= -1;
                }
                grad_vol_ij.segment(0, 2) = -alpha_ij * n;         // wrt i
                grad_vol_ij.segment(2, 2) = -(1.0 - alpha_ij) * n; // wrt j

                grad_vol_kl.segment(0, 2) = -alpha_ij * n;         // wrt i
                grad_vol_kl.segment(2, 2) = -(1.0 - alpha_ij) * n; // wrt j

                if (ee_impact.impacting_node() == 0) {
                    grad_vol_kl.segment(4, 2) = n; // wrt k
                    grad_vol_ij.segment(4, 2) = n; // wrt k
                } else {
                    grad_vol_kl.segment(6, 2) = n; // wrt k
                    grad_vol_ij.segment(6, 2) = n; // wrt k
                }
            }

            // x entries:
            int nodes[4] = { node_i, node_j, node_k, node_l };
            int c_ij = 2 * int(i) + 0;
            int c_kl = 2 * int(i) + 1;
            for (int n_id = 0; n_id < 4; ++n_id) {
                // x and y entries
                jac_uk(int(c_ij), nodes[n_id]) = grad_vol_ij(2 * n_id);
                jac_uk(int(c_ij), nodes[n_id] + vertices.rows()) =
                    grad_vol_ij(2 * n_id + 1);

                // x and y entries
                jac_uk(int(c_kl), nodes[n_id]) = grad_vol_kl(2 * n_id);
                jac_uk(int(c_kl), nodes[n_id] + vertices.rows()) =
                    grad_vol_kl(2 * n_id + 1);
            }
        }
    }

    void VolumeConstraint::compute_constraints(
        const Eigen::MatrixXd& Uk,
        const std::vector<EdgeEdgeImpact>& ee_impacts,
        Eigen::VectorXd& g_uk,
        Eigen::MatrixXd& g_uk_jacobian)
    {
        compute_constraints(Uk, ee_impacts, g_uk);
        compute_constraints_jacobian(Uk, ee_impacts, g_uk_jacobian);
    }

    void VolumeConstraint::dense_indices(
        const std::vector<EdgeEdgeImpact>& ee_impacts,
        Eigen::VectorXi& dense_indices)
    {
        dense_indices.resize(int(ee_impacts.size()) * 2);

        const int num_edges = int(edges.rows());

        for (size_t ee = 0; ee < ee_impacts.size(); ++ee) {
            auto& ee_impact = ee_impacts[ee];

            long c_ij =
                get_constraint_index(ee_impact, /*impacted=*/true, num_edges);
            long c_kl =
                get_constraint_index(ee_impact, /*impacted=*/false, num_edges);

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
