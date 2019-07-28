#include "optimization_problem.hpp"

#include <autodiff/finitediff.hpp>
#include <logger.hpp>

namespace ccd {
namespace opt {
    OptimizationProblem::OptimizationProblem(const std::string& name)
        : name_(name)
    {
    }

    OptimizationProblem::~OptimizationProblem() {}

    ////////////////////////////////////////////////////////////////////////////
    // Objective function and its derivatives.

    Eigen::VectorXd OptimizationProblem::eval_grad_f_approx(
        const Eigen::VectorXd& x)
    {
        Eigen::VectorXd grad;
        ccd::finite_gradient(x, func_f(), grad);
        return grad;
    }

    Eigen::MatrixXd OptimizationProblem::eval_hessian_f_approx(
        const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd hess;
        auto func_grad_f = [&](const Eigen::VectorXd& xk) -> Eigen::VectorXd {
            return eval_grad_f(xk);
        };

        ccd::finite_jacobian(x, func_grad_f, hess, AccuracyOrder::SECOND);
        return hess;
    }

    callback_f OptimizationProblem::func_f()
    {
        return [&](const Eigen::VectorXd& x) -> double { return eval_f(x); };
    }

    bool OptimizationProblem::compare_grad_f_approx(
        const Eigen::VectorXd& x, const Eigen::VectorXd& grad)
    {
        Eigen::VectorXd grad_approx = eval_grad_f_approx(x);
        return compare_gradient(grad, grad_approx);
    }

    bool OptimizationProblem::compare_hessian_f_approx(
        const Eigen::VectorXd& x, const Eigen::SparseMatrix<double>& hessian)
    {
        Eigen::MatrixXd hess_approx = eval_hessian_f_approx(x);
        return compare_jacobian(hessian.toDense(), hess_approx);
    }
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // Constraint function and its derivatives.

    Eigen::MatrixXd OptimizationProblem::eval_jac_g_approx(
        const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd jac;
        auto func_g = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            return eval_g(x);
        };
        ccd::finite_jacobian(x, func_g, jac);
        return jac;
    }

    bool OptimizationProblem::compare_jac_g_approx(
        const Eigen::VectorXd& x, const Eigen::MatrixXd& jac)
    {
        Eigen::MatrixXd jac_approx = eval_jac_g_approx(x);
        return compare_jacobian(jac, jac_approx);
    }
    bool OptimizationProblem::compare_jac_g_approx(
        const Eigen::VectorXd& x, const Eigen::MatrixXd& jac, double& diff_norm)
    {
        Eigen::MatrixXd jac_approx = eval_jac_g_approx(x);
        diff_norm = (jac_approx - jac).norm();
        bool same = compare_jacobian(jac, jac_approx);
        return same;
    }


} // namespace opt
} // namespace ccd
