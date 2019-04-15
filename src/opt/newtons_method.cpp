// Functions for optimizing functions.
// Includes Newton's method with and without constraints.
#include <cmath>
#include <iostream>

#include <Eigen/LU>
#include <Eigen/Sparse>

#include <opt/newtons_method.hpp>

namespace ccd {
namespace opt {

    bool line_search(const Eigen::VectorXd& x, const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f, const double fx,
        double& gamma)
    {
        gamma = 1.0;
        for (int i = 0; i <= 32; i++) {
            if (f(x + gamma * dir) < fx) {
                return true;
            }
            gamma /= 2.0;
        }
        return false;
    }

    // Search along a search direction to find a scalar gamma in [0, 1] such
    // that f(x + gamma * dir) ≤ f(x).
    double constrained_line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const std::function<bool(const Eigen::VectorXd&)>& constraint)
    {
        double gamma = 1.0;
        double fx = f(x);
        for (int i = 0; i <= 32; i++) {
            if (f(x + gamma * dir) < fx && constraint(x + gamma * dir))
                return gamma;
            gamma /= 2.0;
        }
        return 0.0;
    }

    // Perform a single step of Newton's Method to minimize the objective
    // unconstrained.
    bool newtons_method_step(
        const OptimizationProblem& problem, Eigen::VectorXd& x, const double mu)
    {
        Eigen::VectorXd g = problem.grad_f(x);
        Eigen::MatrixXd H = problem.hessian_f(x);
        H.diagonal().array() += mu;
        Eigen::VectorXd delta_x = H.lu().solve(-g);
        double gamma;
        if (!line_search(x, delta_x, problem.f, problem.f(x), gamma))
            return false;     // Newton step unsuccessful
        x += gamma * delta_x; // Store the return value for x_{n+1} in x
        return true;          // Newton step successful
    }

    // Performa Newton's Method to minimize the objective unconstrained.
    OptimizationResults newtons_method(const OptimizationProblem& problem,
        const SolverSettings& settings, const double mu)
    {
        Eigen::VectorXd x = problem.x0;

        int i = 0;
        while (i++ < settings.max_iter && newtons_method_step(problem, x, mu)
            && problem.grad_f(x).squaredNorm() > settings.absolute_tolerance)
            ;

        return OptimizationResults(x, problem.f(x), true);
    }

} // namespace opt
} // namespace ccd
