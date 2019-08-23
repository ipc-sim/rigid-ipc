/// Functions for optimizing functions.
#include "gradient_descent_solver.hpp"

#include <igl/slice.h>
#include <igl/slice_into.h>

#include <solvers/line_search.hpp>

#include <logger.hpp>

namespace ccd {
namespace opt {

    GradientDescentSolver::GradientDescentSolver()
        : GradientDescentSolver("gradient_descent_solver")
    {
    }
    GradientDescentSolver::GradientDescentSolver(const std::string& name)
        : absolute_tolerance(1e-5)
        , min_step_length(1e-12)
        , name_(name)
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

    GradientDescentSolver::~GradientDescentSolver() {}

    OptimizationResults GradientDescentSolver::solve(IBarrierProblem& problem)
    {
        Eigen::VectorXd x = problem.starting_point();
        Eigen::VectorXd gradient, gradient_free;
        Eigen::VectorXd delta_x = Eigen::VectorXd::Zero(problem.num_vars());
        double step_length = 1.0;

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
            step_length = 1.0;
            found_step_length = line_search(
                problem, x, delta_x, problem.eval_f(x), step_length);

            if (!found_step_length) {
                spdlog::warn("solver=newtons_method iter={:d} failure=\"line "
                             "search using gradient descent\" failsafe=none",
                    iter + 1);
                exit_reason = "line-search failed";
                break;
            }

            x += step_length * delta_x; // Update x

        } while (++iter <= max_iterations);

        spdlog::trace("solver=newtons_method total_iter={:d} "
                      "exit_reason=\"{:s}\" sqr_norm_grad={:g}",
            iter, exit_reason, gradient_free.squaredNorm());

        return OptimizationResults(x, problem.eval_f(x),
            gradient_free.squaredNorm() <= absolute_tolerance);
    }

    bool GradientDescentSolver::line_search(IBarrierProblem& problem,
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const double fx,
        double& step_length)
    {

        double step_norm = (step_length * dir).norm();
        bool success = false;

        int num_it = 0;

        while (step_norm >= min_step_length) {

            Eigen::VectorXd xi = x + step_length * dir;
            double fxi = problem.eval_f(xi);

            bool min_rule = fxi < fx;
            bool cstr = !std::isinf(fxi);

            num_it += 1;
            if (min_rule && cstr) {
                success = true;
                break; // while loop
            }

            step_length /= 2.0;
            step_norm = (step_length * dir).norm();
        }

        return success;
    }

    void GradientDescentSolver::init_free_dof(Eigen::VectorXb is_dof_fixed)
    {
        free_dof = Eigen::VectorXi(is_dof_fixed.size() - is_dof_fixed.count());
        for (int i = 0, j = 0; i < is_dof_fixed.size(); i++) {
            if (!is_dof_fixed(i)) {
                free_dof(j++) = i;
            }
        }
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
