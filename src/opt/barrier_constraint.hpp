#pragma once

#include <Eigen/Core>

#include <opt/collision_constraint.hpp>

#include <ccd/collision_detection.hpp>

namespace ccd {
namespace opt {

    enum InitialBarrierEpsilon { MIN_TOI, MAX_TOI, CUSTOM };

    class BarrierConstraint : public CollisionConstraint {
    public:
        BarrierConstraint();
        BarrierConstraint(const std::string& name);

        void settings(const nlohmann::json& json) override;
        nlohmann::json settings() const override;

        void resetBarrierEpsilon();

        bool is_barrier() override { return true; }
        double get_barrier_epsilon() override { return barrier_epsilon; }
        void set_barrier_epsilon(const double eps) override { barrier_epsilon = eps; }

        void initialize(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges,
            const Eigen::MatrixXd& Uk) override;

        void detectCollisions(const Eigen::MatrixXd& Uk) override;

        int number_of_constraints() override;

        void compute_constraints(
            const Eigen::MatrixXd& Uk, Eigen::VectorXd& barriers) override;

        void compute_constraints_jacobian(const Eigen::MatrixXd& Uk,
            Eigen::MatrixXd& barriers_jacobian) override;

        void compute_constraints_hessian(const Eigen::MatrixXd& Uk,
            std::vector<Eigen::SparseMatrix<double>>& barriers_hessian)
            override;

        void compute_constraints_and_derivatives(const Eigen::MatrixXd& Uk,
            Eigen::VectorXd& g_uk,
            Eigen::MatrixXd& g_uk_jacobian,
            std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian) override;

        template <typename T>
        void compute_constraints_per_impact(
            const Eigen::MatrixXd& displacements, std::vector<T>& constraints);

        // Settings
        // ----------
        double barrier_epsilon;
        double custom_inital_epsilon;
        InitialBarrierEpsilon initial_epsilon;
    };

} // namespace opt
} // namespace ccd
