#pragma once

#include <memory>

#include <solvers/optimization_solver.hpp>

namespace ccd {
namespace opt {

    class BarrierProblem : public virtual OptimizationProblem {
    public:
        virtual ~BarrierProblem() = default;

        /// Compute the objective function f(x)
        virtual double compute_objective(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& grad,
            Eigen::SparseMatrix<double>& hess,
            bool compute_grad = true,
            bool compute_hess = true) override;

        // Include thes lines to avoid issues with overriding inherited
        // functions with the same name.
        // (http://www.cplusplus.com/forum/beginner/24978/)
        using OptimizationProblem::compute_objective;

        /// Compute E(x) in f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k))
        virtual double compute_energy_term(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& grad,
            Eigen::SparseMatrix<double>& hess,
            bool compute_grad = true,
            bool compute_hess = true) = 0;

        /// Compute ∑_{k ∈ C} b(d(x_k)) in f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k))
        /// @returns number of active barriers
        virtual double compute_barrier_term(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& grad,
            Eigen::SparseMatrix<double>& hess,
            int& num_constraints,
            bool compute_grad = true,
            bool compute_hess = true) = 0;

        // --------------------------------------------------------------------
        // Convience functions
        // --------------------------------------------------------------------

        virtual double compute_energy_term(const Eigen::VectorXd& x) final
        {
            Eigen::VectorXd grad;
            Eigen::SparseMatrix<double> hess;
            return compute_energy_term(
                x, grad, hess,
                /*compute_grad=*/false, /*compute_hess=*/false);
        }

        virtual double compute_energy_term(
            const Eigen::VectorXd& x, Eigen::VectorXd& grad) final
        {
            Eigen::SparseMatrix<double> hess;
            return compute_energy_term(
                x, grad, hess,
                /*compute_grad=*/true, /*compute_hess=*/false);
        }

        virtual double compute_barrier_term(
            const Eigen::VectorXd& x, int& num_constraints) final
        {
            Eigen::VectorXd grad;
            Eigen::SparseMatrix<double> hess;
            return compute_barrier_term(
                x, grad, hess, num_constraints,
                /*compute_grad=*/false, /*compute_hess=*/false);
        }

        virtual double compute_barrier_term(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& grad,
            int& num_constraints) final
        {
            Eigen::SparseMatrix<double> hess;
            return compute_barrier_term(
                x, grad, hess, num_constraints,
                /*compute_grad=*/true, /*compute_hess=*/false);
        }

        virtual bool is_barrier_problem() const override { return true; }

        virtual double get_barrier_homotopy() const = 0;
        virtual void set_barrier_homotopy(const double eps) = 0;

        virtual double get_barrier_stiffness() const = 0;
        virtual void set_barrier_stiffness(const double kappa) = 0;
    };

    /// Helper Functions for checking finite differences
    /// ------------------------------------------------
    Eigen::VectorXd
    eval_grad_energy_approx(BarrierProblem& problem, const Eigen::VectorXd& x);
    Eigen::MatrixXd
    eval_hess_energy_approx(BarrierProblem& problem, const Eigen::VectorXd& x);

    // WARNING: You cannot create a eval_*_barrier_approx because barrier term
    // eval re-computes the constraint set which can lead to poor performance
    // and discontiuties when using finite differences.

} // namespace opt
} // namespace ccd
