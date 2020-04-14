/// Functions for optimizing functions.
#include "gradient_descent_solver.hpp"

#include <igl/slice.h>
#include <igl/slice_into.h>

#include <solvers/newton_solver.hpp> // init_free_dof()

#include <logger.hpp>

namespace ccd {
namespace opt {

    GradientDescentSolver::GradientDescentSolver()
        : absolute_tolerance(1e-5)
        , min_step_length(1e-12)
    {
    }

    void GradientDescentSolver::settings(const nlohmann::json& json)
    {
        max_iterations = json["max_iterations"].get<int>();
        absolute_tolerance = json["absolute_tolerance"].get<double>();
        min_step_length = json["min_step_length"].get<double>();
    }

    nlohmann::json GradientDescentSolver::settings() const
    {
        nlohmann::json json;
        json["max_iterations"] = max_iterations;
        json["absolute_tolerance"] = absolute_tolerance;
        json["min_step_length"] = min_step_length;
        return json;
    }

    void GradientDescentSolver::init_solve()
    {
        assert(problem_ptr != nullptr);
        free_dof = init_free_dof(problem_ptr->is_dof_fixed());
    }

    OptimizationResults GradientDescentSolver::solve(const Eigen::VectorXd& x0)
    {
        assert(problem_ptr != nullptr);

        Eigen::VectorXd x = x0;
        Eigen::VectorXd gradient, gradient_free;
        Eigen::VectorXd delta_x = Eigen::VectorXd::Zero(x.size());
        double step_length = 1.0;

        int iter = 0;
        std::string exit_reason = "exceeded the maximum allowable iterations";
        do {
            // Compute the gradient and hessian
            double fx;
            problem_ptr->compute_objective(x, fx, gradient);
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
            step_length = 1.0;
            found_step_length = line_search(x, delta_x, fx, step_length);

            if (!found_step_length) {
                spdlog::warn(
                    "solver={} iter={:d} failure=\"line "
                    "search using gradient descent\" failsafe=none",
                    name(), iter + 1);
                exit_reason = "line-search failed";
                break;
            }

            x += step_length * delta_x; // Update x

        } while (++iter <= max_iterations);

        spdlog::trace(
            "solver={} total_iter={:d} exit_reason=\"{:s}\" sqr_norm_grad={:g}",
            name(), iter, exit_reason, gradient_free.squaredNorm());

        double fx;
        problem_ptr->compute_objective(x, fx);
        return OptimizationResults(
            x, fx, gradient_free.squaredNorm() <= absolute_tolerance, true);
    }

    bool GradientDescentSolver::line_search(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const double fx,
        double& step_length)
    {
        bool success = false;

        int num_it = 0;
        while ((step_length * dir).norm() >= min_step_length) {
            Eigen::VectorXd xi = x + step_length * dir;

            bool no_collisions = !problem_ptr->has_collisions(x, xi);
            double fxi;
            problem_ptr->compute_objective(xi, fxi);

            num_it++;
            if (fxi < fx && no_collisions) {
                success = true;
                break; // while loop
            }

            step_length /= 2.0;
        }

        return success;
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
