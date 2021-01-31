#include "optimization_problem.hpp"

#include <finitediff.hpp>
#include <logger.hpp>

namespace ipc::rigid {

Eigen::VectorXd eval_grad_objective_approx(
    OptimizationProblem& problem, const Eigen::VectorXd& x)
{
    auto f = [&](const Eigen::VectorXd& xk) -> double {
        return problem.compute_objective(xk);
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
        Eigen::VectorXd grad;
        problem.compute_objective(xk, grad);
        return grad;
    };

    fd::finite_jacobian(x, func_grad_f, hess, fd::AccuracyOrder::SECOND);
    return hess;
}

} // namespace ipc::rigid
