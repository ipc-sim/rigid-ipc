#include "volume_constraint.hpp"

#include <iostream>

#include <autogen/collision_volume.hpp>
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

    int VolumeConstraint::number_of_constraints()
    {
        return int(get_constraints_size(int(edges->rows())));
    }

    void VolumeConstraint::compute_constraints(
        const Eigen::MatrixXd& displacements, Eigen::VectorXd& volumes)
    {
        std::vector<double> impact_volumes;
        compute_constraints_per_impact(displacements, impact_volumes);
        assemble_constraints_all(impact_volumes, volumes);
    }

    void VolumeConstraint::compute_constraints_jacobian(
        const Eigen::MatrixXd& Uk, Eigen::MatrixXd& g_uk_jacobian)
    {
        std::vector<DScalar> impact_volumes;
        compute_constraints_per_impact(Uk, impact_volumes);
        assemble_jacobian_all(impact_volumes, g_uk_jacobian);
    }

    void VolumeConstraint::compute_constraints_hessian(
        const Eigen::MatrixXd& Uk,
        std::vector<Eigen::SparseMatrix<double>>& volumes_hessian)
    {
        std::vector<DScalar> impact_volumes;
        compute_constraints_per_impact(Uk, impact_volumes);
        assemble_hessian(impact_volumes, volumes_hessian);
    }

    void VolumeConstraint::compute_constraints_jacobian(
        const Eigen::MatrixXd& Uk, Eigen::SparseMatrix<double>& g_uk_jacobian)
    {
        std::vector<DScalar> impact_volumes;
        compute_constraints_per_impact(Uk, impact_volumes);
        assemble_sparse_jacobian_all(impact_volumes, g_uk_jacobian);
    }

    void VolumeConstraint::compute_constraints(const Eigen::MatrixXd& Uk,
        Eigen::VectorXd& g_uk,
        Eigen::SparseMatrix<double>& g_uk_jacobian,
        Eigen::VectorXi& g_uk_active)
    {
        std::vector<DScalar> impact_volumes;
        compute_constraints_per_impact(Uk, impact_volumes);
        assemble_constraints_all(impact_volumes, g_uk);
        assemble_sparse_jacobian_all(impact_volumes, g_uk_jacobian);
        dense_indices(g_uk_active);
    }

    void VolumeConstraint::compute_active_constraints(const Eigen::MatrixXd& Uk,
        Eigen::VectorXd& g_uk,
        Eigen::MatrixXd& g_uk_jacobian)
    {
        std::vector<DScalar> impact_volumes;
        compute_constraints_per_impact(Uk, impact_volumes);
        assemble_constraints(impact_volumes, g_uk);
        assemble_jacobian(impact_volumes, g_uk_jacobian);
    }

    Eigen::VectorXd getGradient(const double& x) { return Eigen::VectorXd(); }
    Eigen::VectorXd getGradient(const DScalar& x) { return x.getGradient(); }
    double getValue(const double& x) { return 0.0; }
    double getValue(const DScalar& x) { return x.getValue(); }
    template <typename T>
    void VolumeConstraint::compute_constraints_per_impact(
        const Eigen::MatrixXd& displacements, std::vector<T>& constraints)
    {
        constraints.clear();
        constraints.reserve(2 * ee_impacts.size());

        if (differentiable<T>()) {
            // !!! All definitions using DScalar must be done after this !!!
            DiffScalarBase::setVariableCount(8);
        }

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

            ImpactTData<T> data = get_impact_data<T>(displacements, ee_impact);

            T toi, alpha_ij, alpha_kl;
            T v_ij(0), v_kl(0);

            if (compute_toi_alpha(data, toi, alpha_ij, alpha_kl)) {

                v_ij = ccd::autogen::space_time_collision_volume<T>(data.v[0],
                    data.v[1], data.u[0], data.u[1], toi, alpha_ij,
                    volume_epsilon);

                v_kl = ccd::autogen::space_time_collision_volume<T>(data.v[2],
                    data.v[3], data.u[2], data.u[3], toi, alpha_kl,
                    volume_epsilon);

//                if (differentiable<T>()) {
//                    int vi = edges->coeffRef(ee_impact.impacted_edge_index, 0);
//                    int vj = edges->coeffRef(ee_impact.impacted_edge_index, 1);

//                    int vk = edges->coeffRef(ee_impact.impacting_edge_index, 0);
//                    int vl = edges->coeffRef(ee_impact.impacting_edge_index, 1);

//                    spdlog::trace(
//                        "constraints vi={} vj={} vk={} vl={} a={} "
//                        "\ntoi_ij   {}"
//                        "\nalpha_ij {} "
//                        "\nvol_ij   {} "
//                        "\nvol_kl   {}"
//                        "\ntoi_ij={:10e} alpha_ij={:10e}  vol_ij={:10e} vol_kl={:10e}",
//                        vi, vj, vk, vl, ee_impact.impacting_alpha,
//                        log::fmt_eigen(getGradient(toi)),
//                        log::fmt_eigen(getGradient(alpha_ij)),
//                        log::fmt_eigen(getGradient(v_ij)),
//                        log::fmt_eigen(getGradient(v_kl)), getValue(toi),
//                        getValue(alpha_ij), getValue(v_ij), getValue(v_kl));
//                }
            }
            constraints.push_back(v_ij);
            constraints.push_back(v_kl);
        }
    }

    void VolumeConstraint::assemble_constraints_all(
        const std::vector<double>& impact_volumes,
        Eigen::VectorXd& dense_volumes)
    {

        const int num_edges = int(edges->rows());
        const long num_constr = get_constraints_size(num_edges);
        dense_volumes.resize(num_constr);
        dense_volumes.setZero();

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

            long c_ij
                = get_constraint_index(ee_impact, /*impacted=*/true, num_edges);
            long c_kl = get_constraint_index(
                ee_impact, /*impacted=*/false, num_edges);

            dense_volumes[c_ij] = impact_volumes[2 * i + 0];
            dense_volumes[c_kl] = impact_volumes[2 * i + 1];
        }
    }

    void VolumeConstraint::assemble_constraints_all(
        const std::vector<DScalar>& impact_volumes,
        Eigen::VectorXd& dense_volumes)
    {
        Eigen::VectorXd volumes;
        assemble_constraints(impact_volumes, volumes);

        const int num_edges = int(edges->rows());
        const long num_constr = get_constraints_size(num_edges);
        dense_volumes.resize(num_constr);
        dense_volumes.setZero();

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

            long c_ij
                = get_constraint_index(ee_impact, /*impacted=*/true, num_edges);
            long c_kl = get_constraint_index(
                ee_impact, /*impacted=*/false, num_edges);

            dense_volumes[c_ij] = volumes(2 * int(i) + 0);
            dense_volumes[c_kl] = volumes(2 * int(i) + 1);
        }
    }

    [[deprecated("Replaced by sparse jabobian")]] void
    VolumeConstraint::assemble_jacobian_all(
        const std::vector<DScalar>& impact_volumes,
        Eigen::MatrixXd& volumes_jac)
    {

        Eigen::MatrixXd impact_volumes_jac;
        assemble_jacobian(impact_volumes, impact_volumes_jac);

        assert(impact_volumes_jac.cols() == vertices->size());
        assert(impact_volumes_jac.rows() == int(impact_volumes.size()));

        const int num_edges = int(edges->rows());
        const long num_constr = get_constraints_size(num_edges);
        volumes_jac.resize(num_constr, impact_volumes_jac.cols());
        volumes_jac.setZero();

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

            long c_ij
                = get_constraint_index(ee_impact, /*impacted=*/true, num_edges);
            long c_kl = get_constraint_index(
                ee_impact, /*impacted=*/false, num_edges);

            volumes_jac.row(c_ij) = impact_volumes_jac.row(2 * int(i) + 0);
            volumes_jac.row(c_kl) = impact_volumes_jac.row(2 * int(i) + 1);
        }
    }

    void VolumeConstraint::assemble_sparse_jacobian_all(
        const std::vector<DScalar>& constraints,
        Eigen::SparseMatrix<double>& jacobian)
    {

        std::vector<DoubleTriplet> tripletList;
        assemble_jacobian_triplets(constraints, tripletList);

        const int num_edges = int(edges->rows());
        const long num_constr = get_constraints_size(num_edges);

        size_t counter = 0;
        for (size_t ee = 0; ee < ee_impacts.size(); ++ee) {
            auto& ee_impact = ee_impacts[ee];

            long c_ij
                = get_constraint_index(ee_impact, /*impacted=*/true, num_edges);
            long c_kl = get_constraint_index(
                ee_impact, /*impacted=*/false, num_edges);
            int rows[2] = { int(c_ij), int(c_kl) };

            for (size_t k = 0; k < 2; k++) {
                for (int i = 0; i < 8; i++) {
                    tripletList[counter++].set_row(rows[k]);
                }
            }
        }

        jacobian.resize(int(num_constr), int(vertices->size()));
        jacobian.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    void VolumeConstraint::dense_indices(Eigen::VectorXi& dense_indices)
    {
        dense_indices.resize(int(ee_impacts.size()) * 2);

        const int num_edges = int(edges->rows());

        for (size_t ee = 0; ee < ee_impacts.size(); ++ee) {
            auto& ee_impact = ee_impacts[ee];

            long c_ij
                = get_constraint_index(ee_impact, /*impacted=*/true, num_edges);
            long c_kl = get_constraint_index(
                ee_impact, /*impacted=*/false, num_edges);

            dense_indices(int(2 * ee) + 0) = int(c_ij);
            dense_indices(int(2 * ee) + 1) = int(c_kl);
        }
    }

    template void VolumeConstraint::compute_constraints_per_impact<double>(
        const Eigen::MatrixXd& displacements, std::vector<double>& constraints);

    template void VolumeConstraint::compute_constraints_per_impact<DScalar>(
        const Eigen::MatrixXd& displacements,
        std::vector<DScalar>& constraints);

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
