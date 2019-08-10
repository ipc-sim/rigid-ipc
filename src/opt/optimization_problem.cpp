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
            return eval_g(x, /*update_constraints=*/false);
        };
        ccd::finite_jacobian(x, func_g, jac, AccuracyOrder::SECOND, 1e-8);
        return jac;
    }

    std::vector<Eigen::MatrixXd> OptimizationProblem::eval_hessian_g_approx(
        const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd finite_hess_i;
        int num_constraints = eval_jac_g(x).rows();
        std::vector<Eigen::MatrixXd> hessians;

        for (int i = 0; i < num_constraints; ++i) {
            auto f = [&](const Eigen::VectorXd& xk) -> Eigen::VectorXd {
                Eigen::MatrixXd actual_jac
                    = eval_jac_g(xk, /*update_constraints=*/false);
                return actual_jac.row(int(i));
            };
            Eigen::MatrixXd finite_hess_i;
            ccd::finite_jacobian(x, f, finite_hess_i, AccuracyOrder::SECOND, 1e-8);
            hessians.push_back(finite_hess_i);
        }
        return hessians;
    }

    bool OptimizationProblem::compare_jac_g_approx(
        const Eigen::VectorXd& x, const Eigen::MatrixXd& jac)
    {
        Eigen::MatrixXd jac_approx = eval_jac_g_approx(x);
        auto r = compare_jacobian(jac, jac_approx);
        if (r){
            return true;
        }
        else {
            return false;
        }
    }

    bool OptimizationProblem::compare_jac_g_approx(
        const Eigen::VectorXd& x, const Eigen::MatrixXd& jac, double& diff_norm)
    {
        Eigen::MatrixXd jac_approx = eval_jac_g_approx(x);
        diff_norm = (jac_approx - jac).norm();
        bool same = compare_jacobian(jac, jac_approx);
        return same;
    }
    bool OptimizationProblem::compare_hessian_g_approx(const Eigen::VectorXd& x,
        const std::vector<Eigen::SparseMatrix<double>>& hessians)
    {
        std::vector<Eigen::MatrixXd> hessians_approx = eval_hessian_g_approx(x);
        bool same = hessians_approx.size() == hessians.size();
        if (!same) {
            return false;
        }
        for (size_t i = 0; i < hessians_approx.size(); i++) {
            Eigen::MatrixXd hessian_i = hessians[i].toDense();
            same = compare_jacobian(hessian_i, hessians_approx[i]);
            if (!same) {
                return false;
            }
        }
        return true;
    }

} // namespace opt
} // namespace ccd
