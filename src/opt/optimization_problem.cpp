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

    Eigen::SparseMatrix<double> OptimizationProblem::eval_hessian_f_sparse(
        const Eigen::VectorXd& x)
    {
        return eval_hessian_f(x).sparseView();
    }

    Eigen::MatrixXd OptimizationProblem::eval_hessian_f_approx(
        const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd hess;
        ccd::finite_jacobian(x, func_grad_f(), hess);
        return hess;
    }

    // Evaluate the objective and its derivatives.
    void OptimizationProblem::eval_f_and_fdiff(
        const Eigen::VectorXd& x, double& value, Eigen::VectorXd& grad)
    {
        value = eval_f(x);
        grad = eval_grad_f(x);
    }

    // Evaluate the objective and its derivatives.
    void OptimizationProblem::eval_f_and_fdiff(const Eigen::VectorXd& x,
        double& value, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian)
    {
        value = eval_f(x);
        grad = eval_grad_f(x);
        hessian = eval_hessian_f(x);
    }

    void OptimizationProblem::eval_f_and_fdiff(const Eigen::VectorXd& x,
        double& value, Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hessian)
    {
        value = eval_f(x);
        grad = eval_grad_f(x);
        hessian = eval_hessian_f_sparse(x);
    }

    callback_f OptimizationProblem::func_f()
    {
        return [&](const Eigen::VectorXd& x) -> double { return eval_f(x); };
    }

    callback_grad_f OptimizationProblem::func_grad_f()
    {
        return [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            return eval_grad_f(x);
        };
    }

    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // Constraint function and its derivatives.

    void OptimizationProblem::eval_g(const Eigen::VectorXd& x,
        Eigen::VectorXd& gx, Eigen::SparseMatrix<double>& gx_jacobian,
        Eigen::VectorXi& gx_active)
    {
        gx = eval_g(x);
        gx_jacobian = eval_jac_g(x).sparseView();
        gx_active = Eigen::VectorXi::LinSpaced(gx.rows(), 0, int(gx.rows()));
    }

    void OptimizationProblem::eval_jac_g(
        const Eigen::VectorXd& x, Eigen::SparseMatrix<double>& jac_gx)
    {
        jac_gx = eval_jac_g(x).sparseView();
    }

    void OptimizationProblem::eval_g_and_gdiff(const Eigen::VectorXd& x,
        Eigen::VectorXd& gx, Eigen::MatrixXd& gx_jacobian,
        std::vector<Eigen::SparseMatrix<double>>& gx_hessian)
    {
        gx = eval_g(x);
        gx_jacobian = eval_jac_g(x);
        gx_hessian = eval_hessian_g(x);
    }

    callback_g OptimizationProblem::func_g()
    {
        return [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            return eval_g(x);
        };
    }

    ////////////////////////////////////////////////////////////////////////////

    bool OptimizationProblem::eval_intermediate_callback(
        const Eigen::VectorXd& /*x*/)
    {
        return true;
    }

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

    // Enable the line search mode. This functionality is not up to the child
    // class.
    void OptimizationProblem::enable_line_search_mode(
        const Eigen::VectorXd&) {};
    // Disable the line search mode.
    void OptimizationProblem::disable_line_search_mode() {};

} // namespace opt
} // namespace ccd
