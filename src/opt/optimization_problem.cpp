#include <ccd/not_implemented_error.hpp>
#include <opt/optimization_problem.hpp>

namespace ccd {
namespace opt {

    OptimizationProblem::~OptimizationProblem() {}
    callback_f OptimizationProblem::func_f()
    {
        callback_f f = [&](const Eigen::VectorXd& x) { return eval_f(x); };
        return f;
    }

    callback_grad_f OptimizationProblem::func_grad_f()
    {
        callback_grad_f f = [&](const Eigen::VectorXd& x) { return eval_grad_f(x); };
        return f;
    }

    callback_grad_f OptimizationProblem::func_g()
    {
        callback_grad_f f = [&](const Eigen::VectorXd& x) { return eval_g(x); };
        return f;
    }

    bool OptimizationProblem::validate_problem()
    {
        bool valid = true;
        valid &= num_vars == x0.rows();
        valid &= num_vars == x_lower.rows() || x_lower.size() == 0;
        valid &= num_vars == x_upper.rows() || x_upper.size() == 0;
        valid &= num_constraints == g_lower.rows() || g_lower.size() == 0;
        valid &= num_constraints == g_upper.rows() || g_upper.size() == 0;

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
        return true;
    }

    // Check if all constraints are satisfied at a location.
    bool OptimizationProblem::are_constraints_satisfied(
        const Eigen::VectorXd& x, const double tol)
    {
        Eigen::ArrayXd gx = eval_g(x).array();
        // return (this->g_lower.array() - 10 * tol <= gx).all()
        //     && (gx <= this->g_upper.array() + 10 * tol).all()
        //     && (this->x_lower.array() - 10 * tol <= x.array()).all()
        //     && (x.array() <= this->x_upper.array() + 10 * tol).all();
        return (-10 * tol <= gx).all()
            && (this->x_lower.array() - 10 * tol <= x.array()).all()
            && (x.array() <= this->x_upper.array() + 10 * tol).all();
    }

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
        this->hessian_g
            = [](const Eigen::VectorXd&) -> std::vector<Eigen::MatrixXd> {
            throw NotImplementedError("Second derivative of the constraint "
                                      "function not implemented!");
        };
    }
    double AdHocProblem::eval_f(const Eigen::VectorXd& x)
    {
        return f(x);
    }
    Eigen::VectorXd AdHocProblem::eval_grad_f(const Eigen::VectorXd& x)
    {
        return grad_f(x);
    }
    Eigen::MatrixXd AdHocProblem::eval_hessian_f(
        const Eigen::VectorXd& x)
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

    std::vector<Eigen::MatrixXd> AdHocProblem::eval_hessian_g(
        const Eigen::VectorXd& x)
    {
        return hessian_g(x);
    };

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

        this->fixed_dof.resize(this->num_vars);
        this->fixed_dof.setConstant(false); // no-upper-bound
    }

} // namespace opt
} // namespace ccd
