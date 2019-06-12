#include "ad_hoc_problem.hpp"

#include <ccd/not_implemented_error.hpp>

namespace ccd {
namespace opt {

    AdHocProblem::AdHocProblem()
    {
        this->f = [](const Eigen::VectorXd&) -> double {
            throw NotImplementedError("Objective function not implemented!");
        };
        this->grad_f = [](const Eigen::VectorXd&) -> Eigen::VectorXd {
            throw NotImplementedError(
                "Gradient of the objective function not implemented!");
        };
        this->hessian_f = [](const Eigen::VectorXd&) -> Eigen::MatrixXd {
            throw NotImplementedError(
                "Hessian of the objective function not implemented!");
        };
        this->g = [](const Eigen::VectorXd&) -> Eigen::VectorXd {
            throw NotImplementedError("Constraint function not implemented!");
        };
        this->jac_g = [](const Eigen::VectorXd&) -> Eigen::MatrixXd {
            throw NotImplementedError(
                "Jacobian of the constraint function not implemented!");
        };
        this->hessian_g = [](const Eigen::VectorXd&)
            -> std::vector<Eigen::SparseMatrix<double>> {
            throw NotImplementedError("Second derivative of the constraint "
                                      "function not implemented!");
        };
    }

    AdHocProblem::AdHocProblem(int num_vars, int num_constraints)
        : AdHocProblem()
    {
        this->num_vars = num_vars;
        this->num_constraints = num_constraints;

        this->x0.resize(num_vars);
        this->x0.setConstant(0.0);

        this->x_lower.resize(this->num_vars);
        this->x_lower.setConstant(NO_LOWER_BOUND); // no-lower-bound

        this->x_upper.resize(this->num_vars);
        this->x_upper.setConstant(NO_UPPER_BOUND); // no-upper-bound

        this->g_lower.resize(this->num_constraints);
        this->g_lower.setConstant(NO_LOWER_BOUND); // no-lower-bound

        this->g_upper.resize(this->num_constraints);
        this->g_upper.setConstant(NO_UPPER_BOUND); // no-upper-bound

        this->fixed_dof.resize(this->num_vars, 1);
        this->fixed_dof.setConstant(false); // no-upper-bound
    }

    double AdHocProblem::eval_f(const Eigen::VectorXd& x) { return f(x); }

    Eigen::VectorXd AdHocProblem::eval_grad_f(const Eigen::VectorXd& x)
    {
        return grad_f(x);
    }

    Eigen::MatrixXd AdHocProblem::eval_hessian_f(const Eigen::VectorXd& x)
    {
        return hessian_f(x);
    }

    Eigen::VectorXd AdHocProblem::eval_g(const Eigen::VectorXd& x)
    {
        return g(x);
    };

    Eigen::MatrixXd AdHocProblem::eval_jac_g(const Eigen::VectorXd& x)
    {
        return jac_g(x);
    };

    std::vector<Eigen::SparseMatrix<double>> AdHocProblem::eval_hessian_g(
        const Eigen::VectorXd& x)
    {
        return hessian_g(x);
    };

} // namespace opt
} // namespace ccd
