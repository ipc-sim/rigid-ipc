// Functions for optimizing functions.
// Includes Newton's method with and without constraints.
#include <cmath>
#include <iostream>

#include <Eigen/LU>
#include <Eigen/Sparse>

#include <opt/newtons_method.hpp>

namespace ccd {
namespace opt {

    // Performa Newton's Method to minimize the objective unconstrained.
    OptimizationResults newtons_method(const OptimizationProblem& problem,
        const SolverSettings& settings, const double mu)
    {
        Eigen::VectorXd x = problem.x0;

        int i = 0;
        Eigen::VectorXd g = problem.grad_f(x), delta_x;
        Eigen::MatrixXd H;
        double gamma;
        while (i++ <= settings.max_iter
            && g.squaredNorm() > settings.absolute_tolerance) {
            H = problem.hessian_f(x);
            H.diagonal().array() += mu;
            delta_x = H.lu().solve(-g);
            if (!line_search(x, delta_x, problem.f, gamma,
                    settings.line_search_tolerance)) {
                break; // Newton step unsuccessful
            }
            x += gamma * delta_x;  // Update x
            g = problem.grad_f(x); // Recompute the gradient
            // Save intermedtiate results
            settings.intermediate_cb(
                x, problem.f(x), Eigen::VectorXd::Zero(x.size()), 0, i);
        }

        return OptimizationResults(x, problem.f(x),
            i >= settings.max_iter
                || g.squaredNorm() <= settings.absolute_tolerance);
    }

    // Search along a search direction to find a scalar gamma in [0, 1] such
    // that f(x + gamma * dir) ≤ f(x).
    bool line_search(const Eigen::VectorXd& x, const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f, double& gamma,
        const double min_gamma)
    {
        return constrained_line_search(
            x, dir, f, [](const Eigen::VectorXd&) { return true; }, gamma,
            min_gamma);
    }

    // Search along a search direction to find a scalar gamma in [0, 1] such
    // that f(x + gamma * dir) ≤ f(x).
    bool constrained_line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const std::function<bool(const Eigen::VectorXd&)>& constraint,
        double& gamma, const double min_gamma)
    {
        gamma = 1.0;
        const double fx = f(x);
        while (gamma >= min_gamma) {
            if (f(x + gamma * dir) < fx && constraint(x + gamma * dir)) {
                return true;
            }
            gamma /= 2.0;
        }
        return false;
    }

} // namespace opt
} // namespace ccd
