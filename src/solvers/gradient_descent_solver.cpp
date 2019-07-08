/// Functions for optimizing functions.
#include "gradient_descent_solver.hpp"

#include <igl/slice.h>
#include <igl/slice_into.h>

#include <solvers/line_search.hpp>

#include <logger.hpp>

namespace ccd {
namespace opt {

    GradientDescentSolver::GradientDescentSolver()
        : absolute_tolerance(1e-5)
        , min_step_length(1e-12)
    {
    }

    GradientDescentSolver::~GradientDescentSolver() {}

    OptimizationResults GradientDescentSolver::solve(
        OptimizationProblem& problem)
    {
        Eigen::VectorXd x = problem.x0;
        Eigen::VectorXd gradient, gradient_free;
        Eigen::VectorXd delta_x = Eigen::VectorXd::Zero(problem.num_vars);
        double step_length = 1.0;

        // Simple function to make sure the optimization does not violate
        // constraints.
        auto constraint = [&problem](const Eigen::VectorXd& x) {
            return !std::isinf(problem.eval_f(x));
        };

        int iter = 0;
        std::string exit_reason = "exceeded the maximum allowable iterations";
        do {
            // Compute the gradient and hessian
            gradient = problem.eval_grad_f(x);
            // Remove rows of fixed dof
            igl::slice(gradient, free_dof, gradient_free);

            if (gradient_free.squaredNorm() <= absolute_tolerance) {
                exit_reason = "found a local optimum";
                break;
            }

            if (!compute_free_direction(gradient_free, delta_x)) {
                exit_reason = "direction computation failed";
                break;
            }

            // Perform a line search along Δx to stay in the feasible realm
            bool found_step_length;
            // step_length = std::min(1.0, 2 * step_length);
            step_length = 1.0;
            // problem.enable_line_search_mode(x + delta_x);
            found_step_length
                = constrained_line_search(x, delta_x, problem.func_f(),
                    gradient, constraint, step_length, min_step_length);
            // problem.disable_line_search_mode();

            if (!found_step_length) {
                spdlog::warn("solver=newtons_method iter={:d} failure=\"line "
                             "search using gradient descent\" failsafe=none",
                    iter + 1);
                exit_reason = "line-search failed";
                break;
            }

            x += step_length * delta_x; // Update x
            assert(constraint(x));

            problem.eval_intermediate_callback(x);
        } while (++iter <= max_iterations);

        spdlog::trace("solver=newtons_method total_iter={:d} "
                      "exit_reason=\"{:s}\" sqr_norm_grad={:g}",
            iter, exit_reason, gradient_free.squaredNorm());

        return OptimizationResults(x, problem.eval_f(x),
            gradient_free.squaredNorm() <= absolute_tolerance);
    }

    bool GradientDescentSolver::compute_free_direction(
        const Eigen::VectorXd& gradient_free, Eigen::VectorXd& delta_x)
    {
        Eigen::VectorXd delta_x_free;
        bool success = compute_direction(gradient_free, delta_x_free);
        if (success) {
            // Store the free dof back into delta_x
            igl::slice_into(delta_x_free, free_dof, delta_x);
        }
        return success;
    }

    bool GradientDescentSolver::compute_direction(
        const Eigen::VectorXd& gradient, Eigen::VectorXd& delta_x)
    {
        delta_x = -gradient;
        return true;
    }

} // namespace opt
} // namespace ccd
