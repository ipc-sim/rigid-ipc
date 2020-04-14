#include "optimization_problem.hpp"

#include <finitediff.hpp>
#include <logger.hpp>

namespace ccd {
namespace opt {

    Eigen::VectorXd eval_grad_objective_approx(
        OptimizationProblem& problem, const Eigen::VectorXd& x)
    {
        auto f = [&](const Eigen::VectorXd& xk) -> double {
            double fx;
            problem.compute_objective(xk, fx);
            return fx;
        };
        Eigen::VectorXd grad(x.size());
        fd::finite_gradient(x, f, grad);
        return grad;
    }

    Eigen::MatrixXd eval_hess_objective_approx(
        OptimizationProblem& problem, const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd hess;
        auto func_grad_f = [&](const Eigen::VectorXd& xk) -> Eigen::VectorXd {
            double fx;
            Eigen::VectorXd grad_fx;
            problem.compute_objective(xk, fx, grad_fx);
            return grad_fx;
        };

        fd::finite_jacobian(x, func_grad_f, hess, fd::AccuracyOrder::SECOND);
        return hess;
    }

} // namespace opt
} // namespace ccd
