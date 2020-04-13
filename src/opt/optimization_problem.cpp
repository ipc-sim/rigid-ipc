#include "optimization_problem.hpp"

#include <finitediff.hpp>
#include <logger.hpp>

namespace ccd {
namespace opt {

    Eigen::VectorXd eval_grad_f_approx(
        IUnconstraintedProblem& problem, const Eigen::VectorXd& x)
    {
        auto f = [&](const Eigen::VectorXd& xk) -> double {
            return problem.eval_f(xk);
        };
        Eigen::VectorXd grad(x.size());
        fd::finite_gradient(x, f, grad);
        return grad;
    }

    Eigen::MatrixXd eval_hess_f_approx(
        IUnconstraintedProblem& problem, const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd hess;
        auto func_grad_f = [&](const Eigen::VectorXd& xk) -> Eigen::VectorXd {
            return problem.eval_grad_f(xk);
        };

        fd::finite_jacobian(x, func_grad_f, hess, fd::AccuracyOrder::SECOND);
        return hess;
    }


} // namespace opt
} // namespace ccd
