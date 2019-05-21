#include "volume_constraint.hpp"

#include <iostream>

#include <autogen/collision_volume.hpp>
#include <ccd/collision_constraint_diff.hpp>
#include <ccd/collision_volume_diff.hpp>

namespace ccd {

namespace opt {

    VolumeConstraint::VolumeConstraint()
        : volume_epsilon(1E-3)
    {
    }


    int VolumeConstraint::number_of_constraints()
    {
        return int(ccd::autodiff::get_constraints_size(int(edges->rows())));
    }

    void VolumeConstraint::compute_constraints(
        const Eigen::MatrixXd& displacements, Eigen::VectorXd& volumes)
    {
        std::vector<double> impact_volumes;
        compute_constraints_per_impact(displacements, impact_volumes);
        assemble_dense_constraints(impact_volumes, volumes);
    }

    void VolumeConstraint::compute_constraints_jacobian(

        const Eigen::MatrixXd& Uk, Eigen::MatrixXd& volumes_jac)
    {
        std::vector<DScalar> impact_volumes;
        compute_constraints_per_impact(Uk, impact_volumes);
        assemble_dense_jacobian(impact_volumes, volumes_jac);
    }

    void VolumeConstraint::compute_constraints_hessian(
        const Eigen::MatrixXd& Uk,
        std::vector<Eigen::SparseMatrix<double>>& volumes_hessian)
    {
        std::vector<DScalar> impact_volumes;
        compute_constraints_per_impact(Uk, impact_volumes);
        assemble_hessian(impact_volumes, volumes_hessian);
    }

    void VolumeConstraint::compute_constraints_and_derivatives(const Eigen::MatrixXd& Uk,
        Eigen::VectorXd& g_uk, Eigen::MatrixXd& g_uk_jacobian,
        std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian)
    {
        std::vector<DScalar> impact_volumes;
        compute_constraints_per_impact(Uk, impact_volumes);
        assemble_dense_constraints(impact_volumes, g_uk);
        assemble_dense_jacobian(impact_volumes, g_uk_jacobian);
        assemble_hessian(impact_volumes, g_uk_hessian);

    }

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
            }
            constraints.push_back(v_ij);
            constraints.push_back(v_kl);
        }
    }

    void VolumeConstraint::assemble_dense_constraints(
        const std::vector<double>& impact_volumes,
        Eigen::VectorXd& dense_volumes)
    {

        const int num_edges = int(edges->rows());
        const long num_constr = ccd::autodiff::get_constraints_size(num_edges);
        dense_volumes.resize(num_constr);
        dense_volumes.setZero();

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

            long c_ij = ccd::autodiff::get_constraint_index(
                ee_impact, /*impacted=*/true, num_edges);
            long c_kl = ccd::autodiff::get_constraint_index(
                ee_impact, /*impacted=*/false, num_edges);

            dense_volumes[c_ij] = impact_volumes[2 * i + 0];
            dense_volumes[c_kl] = impact_volumes[2 * i + 1];
        }
    }

    void VolumeConstraint::assemble_dense_constraints(
        const std::vector<DScalar>& impact_volumes,
            Eigen::VectorXd& dense_volumes){
        const int num_edges = int(edges->rows());
        const long num_constr = ccd::autodiff::get_constraints_size(num_edges);
        dense_volumes.resize(num_constr);
        dense_volumes.setZero();

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

            long c_ij = ccd::autodiff::get_constraint_index(
                ee_impact, /*impacted=*/true, num_edges);
            long c_kl = ccd::autodiff::get_constraint_index(
                ee_impact, /*impacted=*/false, num_edges);

            dense_volumes[c_ij] = impact_volumes[2 * i + 0].getValue();
            dense_volumes[c_kl] = impact_volumes[2 * i + 1].getValue();
        }
    }

    void VolumeConstraint::assemble_dense_jacobian(
        const std::vector<DScalar>& impact_volumes,
        Eigen::MatrixXd& volumes_jac)
    {

        Eigen::MatrixXd impact_volumes_jac;
        assemble_jacobian(impact_volumes, impact_volumes_jac);

        assert(impact_volumes_jac.cols() == vertices->size());
        assert(impact_volumes_jac.rows() == int(impact_volumes.size()));

        const int num_edges = int(edges->rows());
        const long num_constr = ccd::autodiff::get_constraints_size(num_edges);
        volumes_jac.resize(num_constr, impact_volumes_jac.cols());
        volumes_jac.setZero();

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

            long c_ij = ccd::autodiff::get_constraint_index(
                ee_impact, /*impacted=*/true, num_edges);
            long c_kl = ccd::autodiff::get_constraint_index(
                ee_impact, /*impacted=*/false, num_edges);

            volumes_jac.row(c_ij) = impact_volumes_jac.row(2 * int(i) + 0);
            volumes_jac.row(c_kl) = impact_volumes_jac.row(2 * int(i) + 1);
        }
    }

    template void VolumeConstraint::compute_constraints_per_impact<double>(
        const Eigen::MatrixXd& displacements, std::vector<double>& constraints);

    template void VolumeConstraint::compute_constraints_per_impact<DScalar>(
        const Eigen::MatrixXd& displacements,
        std::vector<DScalar>& constraints);

} // namespace opt
} // namespace ccd
