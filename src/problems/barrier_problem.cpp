// Solve the optimization general_problem using Newton's Method with barriers
// for the constraints.
#include "barrier_problem.hpp"

#include <finitediff.hpp>

namespace ccd {
namespace opt {

    // Compute the objective function f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k))
    void BarrierProblem::compute_objective(
        const Eigen::VectorXd& x,
        double& fx,
        Eigen::VectorXd& grad_fx,
        Eigen::SparseMatrix<double>& hess_fx,
        bool compute_grad,
        bool compute_hess)
    {
        double Ex, Bx;
        Eigen::VectorXd grad_Ex, grad_Bx;
        Eigen::SparseMatrix<double> hess_Ex, hess_Bx;

        // Compute each term
        compute_energy_term(
            x, Ex, grad_Ex, hess_Ex, compute_grad, compute_hess);
        compute_barrier_term(
            x, Bx, grad_Bx, hess_Bx, compute_grad, compute_hess);

        double kappa = get_barrier_stiffness();

        fx = Ex + kappa * Bx;
        if (compute_grad) {
            grad_fx = grad_Ex + kappa * grad_Bx;
        }
        if (compute_hess) {
            hess_fx = hess_Ex + kappa * hess_Bx;
        }
    }

    Eigen::VectorXd
    eval_grad_energy_approx(BarrierProblem& problem, const Eigen::VectorXd& x)
    {
        auto f = [&](const Eigen::VectorXd& xk) -> double {
            double fx;
            problem.compute_energy_term(xk, fx);
            return fx;
        };
        Eigen::VectorXd grad(x.size());
        fd::finite_gradient(x, f, grad);
        return grad;
    }

    Eigen::MatrixXd
    eval_hess_energy_approx(BarrierProblem& problem, const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd hess;
        auto func_grad_f = [&](const Eigen::VectorXd& xk) -> Eigen::VectorXd {
            double fx;
            Eigen::VectorXd grad_fx;
            problem.compute_energy_term(xk, fx, grad_fx);
            return grad_fx;
        };

        fd::finite_jacobian(x, func_grad_f, hess, fd::AccuracyOrder::SECOND);
        return hess;
    }

} // namespace opt
} // namespace ccd
