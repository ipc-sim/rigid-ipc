#include "optimization_problem.hpp"

#include <autodiff/finitediff.hpp>

namespace ccd {
namespace opt {

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
    ////////////////////////////////////////////////////////////////////////////

    bool OptimizationProblem::validate_problem()
    {
        bool valid = true;
        valid &= num_vars == x0.rows();
        valid &= num_vars == x_lower.rows() || x_lower.size() == 0;
        valid &= num_vars == x_upper.rows() || x_upper.size() == 0;
        valid &= num_constraints == g_lower.rows() || g_lower.size() == 0;
        valid &= num_constraints == g_upper.rows() || g_upper.size() == 0;
        valid &= num_vars == is_dof_fixed.rows() || is_dof_fixed.size() == 0;

        if (!valid) {
            return false;
        }

        if (x_lower.size() == 0) {
            x_lower.resize(num_vars);
            x_lower.setConstant(NO_LOWER_BOUND); // no-lower-bound
        }
        if (x_upper.size() == 0) {
            x_upper.resize(num_vars);
            x_upper.setConstant(NO_UPPER_BOUND); // no-upper-bound
        }
        if (g_lower.size() == 0) {
            g_lower.resize(num_constraints);
            g_lower.setConstant(NO_LOWER_BOUND); // no-lower-bound
        }
        if (g_upper.size() == 0) {
            g_upper.resize(num_constraints);
            g_upper.setConstant(NO_UPPER_BOUND); // no-upper-bound
        }
        if (is_dof_fixed.size() == 0) {
            g_upper.resize(num_vars);
            g_upper.setZero(); // no-fixed-dof
        }
        return true;
    }

    // Check if all constraints are satisfied at a location.
    bool OptimizationProblem::are_constraints_satisfied(
        const Eigen::VectorXd& x, const double tol)
    {
        // TODO: Move this function to the constraint
        // return this->constraints.is_satisfied();
        Eigen::ArrayXd gx = eval_g(x).array();
        // TODO: Fix the lower and upper bounds of g for this to work again
        // return (this->g_lower.array() - 10 * tol <= gx).all()
        //     && (gx <= this->g_upper.array() + 10 * tol).all()
        //     && (this->x_lower.array() - 10 * tol <= x.array()).all()
        //     && (x.array() <= this->x_upper.array() + 10 * tol).all();
        // TODO: This does not work for the barrier constraint
        return (-10 * tol <= gx).all()
            && (this->x_lower.array() - 10 * tol <= x.array()).all()
            && (x.array() <= this->x_upper.array() + 10 * tol).all();
    }

} // namespace opt
} // namespace ccd
