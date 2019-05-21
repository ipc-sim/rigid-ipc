#pragma once

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>
#include <opt/collision_constraint.hpp>

namespace ccd {

namespace opt {

    class VolumeConstraint : public CollisionConstraint {
    public:
        VolumeConstraint();

        void compute_constraints(
            const Eigen::MatrixXd& Uk, Eigen::VectorXd& g_uk) override;
        void compute_constraints_jacobian(
            const Eigen::MatrixXd& Uk, Eigen::MatrixXd& g_uk_jacobian) override;
        void compute_constraints_hessian(const Eigen::MatrixXd& Uk,
            std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian) override;
        void compute_constraints_and_derivatives(const Eigen::MatrixXd& Uk,
            Eigen::VectorXd& g_uk, Eigen::MatrixXd& g_uk_jacobian,
            std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian) override;

        void compute_constraints(const Eigen::MatrixXd& Uk,
            Eigen::VectorXd& g_uk, Eigen::MatrixXd& g_uk_jacobian);

        int number_of_constraints() override;

        template <typename T>
        void compute_constraints_per_impact(
            const Eigen::MatrixXd& displacements, std::vector<T>& constraints);

        // Settings
        // ----------
        double volume_epsilon;

        void assemble_dense_constraints(
            const std::vector<double>& impact_volumes,
            Eigen::VectorXd& dense_volumes);
        void assemble_dense_constraints(
            const std::vector<DScalar>& impact_volumes,
            Eigen::VectorXd& dense_volumes);
        void assemble_dense_jacobian(const std::vector<DScalar>& impact_volumes,
            Eigen::MatrixXd& volumes_jac);
    };

} // namespace opt
} // namespace ccd
