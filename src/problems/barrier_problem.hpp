#pragma once

#include <memory>

#include <solvers/optimization_solver.hpp>

namespace ccd {
namespace opt {

    class BarrierProblem : public OptimizationProblem {
    public:
        virtual ~BarrierProblem() = default;

        /// Compute the objective function f(x)
        virtual void compute_objective(
            const Eigen::VectorXd& x,
            double& fx,
            Eigen::VectorXd& grad_fx,
            Eigen::SparseMatrix<double>& hess_fx,
            bool compute_grad = true,
            bool compute_hess = true) override;

        // Include thes lines to avoid issues with overriding inherited
        // functions with the same name.
        // (http://www.cplusplus.com/forum/beginner/24978/)
        using OptimizationProblem::compute_objective;

        /// Compute E(x) in f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k))
        virtual void compute_energy_term(
            const Eigen::VectorXd& x,
            double& Ex,
            Eigen::VectorXd& grad_Ex,
            Eigen::SparseMatrix<double>& hess_Ex,
            bool compute_grad = true,
            bool compute_hess = true) = 0;

        /// Compute ∑_{k ∈ C} b(d(x_k)) in f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k))
        /// @returns number of active barriers
        virtual int compute_barrier_term(
            const Eigen::VectorXd& x,
            double& Bx,
            Eigen::VectorXd& grad_Bx,
            Eigen::SparseMatrix<double>& hess_Bx,
            bool compute_grad = true,
            bool compute_hess = true) = 0;

        // --------------------------------------------------------------------
        // Convience functions
        // --------------------------------------------------------------------

        virtual void
        compute_energy_term(const Eigen::VectorXd& x, double& Ex) final
        {
            Eigen::VectorXd grad_Ex;
            Eigen::SparseMatrix<double> hess_Ex;
            return compute_energy_term(
                x, Ex, grad_Ex, hess_Ex,
                /*compute_grad=*/false, /*compute_hess=*/false);
        }

        virtual void compute_energy_term(
            const Eigen::VectorXd& x,
            double& Ex,
            Eigen::VectorXd& grad_Ex) final
        {
            Eigen::SparseMatrix<double> hess_Ex;
            return compute_energy_term(
                x, Ex, grad_Ex, hess_Ex,
                /*compute_grad=*/true, /*compute_hess=*/false);
        }

        virtual int
        compute_barrier_term(const Eigen::VectorXd& x, double& Bx) final
        {
            Eigen::VectorXd grad_Bx;
            Eigen::SparseMatrix<double> hess_Bx;
            return compute_barrier_term(
                x, Bx, grad_Bx, hess_Bx,
                /*compute_grad=*/false, /*compute_hess=*/false);
        }

        virtual int compute_barrier_term(
            const Eigen::VectorXd& x,
            double& Bx,
            Eigen::VectorXd& grad_Bx) final
        {
            Eigen::SparseMatrix<double> hess_Bx;
            return compute_barrier_term(
                x, Bx, grad_Bx, hess_Bx,
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
    Eigen::VectorXd
    eval_grad_barrier_approx(BarrierProblem& problem, const Eigen::VectorXd& x);
    Eigen::MatrixXd
    eval_hess_barrier_approx(BarrierProblem& problem, const Eigen::VectorXd& x);

} // namespace opt
} // namespace ccd
