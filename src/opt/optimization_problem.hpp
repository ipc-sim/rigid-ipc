#pragma once

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <utils/eigen_ext.hpp>

namespace ccd {
namespace opt {

    class OptimizationProblem {
    public:
        virtual ~OptimizationProblem() = default;

        /// Compute the objective function f(x)
        virtual void compute_objective(
            const Eigen::VectorXd& x,
            double& fx,
            Eigen::VectorXd& grad_fx,
            Eigen::SparseMatrix<double>& hess_fx,
            bool compute_grad = true,
            bool compute_hess = true) = 0;

        // --------------------------------------------------------------------
        // Convience functions
        // --------------------------------------------------------------------

        virtual void
        compute_objective(const Eigen::VectorXd& x, double& fx) final
        {
            Eigen::VectorXd grad_fx;
            Eigen::SparseMatrix<double> hess_fx;
            return compute_objective(
                x, fx, grad_fx, hess_fx,
                /*compute_grad=*/false, /*compute_hess=*/false);
        }

        virtual void compute_objective(
            const Eigen::VectorXd& x,
            double& fx,
            Eigen::VectorXd& grad_fx) final
        {
            Eigen::SparseMatrix<double> hess_fx;
            return compute_objective(
                x, fx, grad_fx, hess_fx,
                /*compute_grad=*/true, /*compute_hess=*/false);
        }

        /// @returns the number of variables
        virtual int num_vars() const = 0;

        /// @returns \f$x_0\f$: the starting point for the optimization.
        virtual const Eigen::VectorXd& starting_point() = 0;

        /// @returns A vector of booleans indicating if a DoF is fixed.
        virtual const Eigen::VectorXb& is_dof_fixed() = 0;

        /// Determine if there is a collision between two configurations
        virtual bool has_collisions(
            const Eigen::VectorXd& xi, const Eigen::VectorXd& xj) const = 0;

        /// Compute the minimum distance among geometry
        virtual double compute_min_distance(const Eigen::VectorXd& x) const = 0;

        virtual bool is_barrier_problem() const { return false; }
        virtual bool is_constrained_problem() const { return false; }
    };

    /// Helper Functions for checking finite differences
    /// ------------------------------------------------
    Eigen::VectorXd eval_grad_objective_approx(
        OptimizationProblem& problem, const Eigen::VectorXd& x);
    Eigen::MatrixXd eval_hess_objective_approx(
        OptimizationProblem& problem, const Eigen::VectorXd& x);

} // namespace opt
} // namespace ccd
