#include "optimization_problem.hpp"

#include <autodiff/finitediff.hpp>
#include <logger.hpp>

namespace ccd {
namespace opt {

    Eigen::VectorXd eval_grad_f_approx(
        IUnconstraintedProblem& problem, const Eigen::VectorXd& x)
    {
        auto f = [&](const Eigen::VectorXd& xk) -> double {
            return problem.eval_f(xk);
        };
        Eigen::VectorXd grad;
        ccd::finite_gradient(x, f, grad);
        return grad;
    }

    Eigen::MatrixXd eval_hess_f_approx(
        IUnconstraintedProblem& problem, const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd hess;
        auto func_grad_f = [&](const Eigen::VectorXd& xk) -> Eigen::VectorXd {
            return problem.eval_grad_f(xk);
        };

        ccd::finite_jacobian(x, func_grad_f, hess, AccuracyOrder::SECOND);
        return hess;
    }

    Eigen::MatrixXd eval_jac_g_approx(
        IConstraintedProblem& problem, const Eigen::VectorXd& x)
    {

        Eigen::MatrixXd jac;
        auto func_g = [&](const Eigen::VectorXd& xk) -> Eigen::VectorXd {
            return problem.eval_g_set(xk, CstrSetFlag::KEEP_CSTR_SET);
        };
        ccd::finite_jacobian(x, func_g, jac, AccuracyOrder::SECOND);
        return jac;
    }

} // namespace opt
} // namespace ccd
