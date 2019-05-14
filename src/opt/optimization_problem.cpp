#include <ccd/not_implemented_error.hpp>
#include <opt/optimization_problem.hpp>

namespace ccd {
namespace opt {

    OptimizationProblem::OptimizationProblem()
        : num_vars()
        , num_constraints()
        , x0()
        , x_lower()
        , x_upper()
        , g_lower()
        , g_upper()
        , fixed_dof()
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

    // Resize fields accordingly
    OptimizationProblem::OptimizationProblem(int num_vars, int num_constraints)
        : OptimizationProblem()
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

    OptimizationProblem::OptimizationProblem(const Eigen::VectorXd& x0,
        const callback_f f, const callback_grad_f grad_f,
        const Eigen::VectorXd& x_lower, const Eigen::VectorXd& x_upper,
        const int num_constraints, const callback_g& g,
        const callback_jac_g& jac_g, const Eigen::VectorXd& g_lower,
        const Eigen::VectorXd& g_upper)
        : OptimizationProblem(int(x0.rows()), num_constraints)
    {
        this->x0 = x0;
        this->x_lower = x_lower;
        this->x_upper = x_upper;
        this->g_lower = g_lower;
        this->g_upper = g_upper;
        this->f = f;
        this->grad_f = grad_f;
        this->g = g;
        this->jac_g = jac_g;
    }

    OptimizationProblem::OptimizationProblem(const Eigen::VectorXd& x0,
        const callback_f f, const callback_grad_f grad_f,
        const callback_hessian_f& hessian_f, const Eigen::VectorXd& x_lower,
        const Eigen::VectorXd& x_upper, const int num_constraints,
        const callback_g& g, const callback_jac_g& jac_g,
        const callback_hessian_g& hessian_g, const Eigen::VectorXd& g_lower,
        const Eigen::VectorXd& g_upper)
        : OptimizationProblem(x0, f, grad_f, x_lower, x_upper, num_constraints,
            g, jac_g, g_lower, g_upper)
    {
        this->hessian_f = hessian_f;
        this->hessian_g = hessian_g;
    }

    bool OptimizationProblem::validate_problem()
    {
        bool valid = true;
        valid &= num_vars == x0.rows();
        valid &= num_vars == x_lower.rows() || x_lower.size() == 0;
        valid &= num_vars == x_upper.rows() || x_upper.size() == 0;
        valid &= num_constraints == g_lower.rows() || g_lower.size() == 0;
        valid &= num_constraints == g_upper.rows() || g_upper.size() == 0;
        valid &= (g == nullptr) == (jac_g == nullptr);
        valid &= f != nullptr;
        valid &= grad_f != nullptr;

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
        const Eigen::VectorXd& x, const double tol) const
    {
        Eigen::ArrayXd gx = this->g(x).array();
        // return (this->g_lower.array() - 10 * tol <= gx).all()
        //     && (gx <= this->g_upper.array() + 10 * tol).all()
        //     && (this->x_lower.array() - 10 * tol <= x.array()).all()
        //     && (x.array() <= this->x_upper.array() + 10 * tol).all();
        return (-10 * tol <= gx).all()
            && (this->x_lower.array() - 10 * tol <= x.array()).all()
            && (x.array() <= this->x_upper.array() + 10 * tol).all();
    }

} // namespace opt
} // namespace ccd
