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

        void compute_constraints_jacobian(const Eigen::MatrixXd& Uk,
            Eigen::SparseMatrix<double>& g_uk_jacobian) override;

        void compute_constraints_hessian(const Eigen::MatrixXd& Uk,
            std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian) override;

        void compute_constraints(const Eigen::MatrixXd& Uk,
            Eigen::VectorXd& g_uk, Eigen::SparseMatrix<double>& g_uk_jacobian,
            Eigen::VectorXi& g_uk_active) override;

        /// \brief compute_constraints sparsly including only active
        void compute_active_constraints(const Eigen::MatrixXd& Uk,
            Eigen::VectorXd& g_uk, Eigen::MatrixXd& g_uk_jacobian);

        int number_of_constraints() override;

        template <typename T>
        void compute_constraints_per_impact(
            const Eigen::MatrixXd& displacements, std::vector<T>& constraints);

        // Settings
        // ----------
        double volume_epsilon;

        void dense_indices(Eigen::VectorXi& dense_indices);
        void assemble_constraints_all(const std::vector<double>& impact_volumes,
            Eigen::VectorXd& dense_volumes);
        void assemble_constraints_all(
            const std::vector<DScalar>& impact_volumes,
            Eigen::VectorXd& dense_volumes);
        void assemble_jacobian_all(const std::vector<DScalar>& impact_volumes,
            Eigen::MatrixXd& volumes_jac);
        void assemble_sparse_jacobian_all(
            const std::vector<DScalar>& impact_volumes,
            Eigen::SparseMatrix<double>& volumes_jac);
    };

    long get_constraint_index(
        const EdgeEdgeImpact& impact, const bool impacted, const int num_edges);

    long get_constraints_size(const int num_edges);

} // namespace opt
} // namespace ccd
