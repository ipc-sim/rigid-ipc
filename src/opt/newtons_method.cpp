// Functions for optimizing functions.
// Includes Newton's method with and without constraints.
#include <cmath>
#include <iostream>

#include <Eigen/LU>
#include <Eigen/Sparse>
#include <igl/slice.h>
#include <igl/slice_into.h>

#include <opt/newtons_method.hpp>

namespace ccd {
namespace opt {

    OptimizationResults newtons_method(OptimizationProblem& problem,
        const double absolute_tolerance, const double line_search_tolerance,
        const int max_iter, const double mu)
    {
        Eigen::VectorXd x = problem.x0;

        // Convert from the boolean vector to a vector of free dof indices
        Eigen::VectorXi free_dof(problem.num_vars - problem.fixed_dof.count());
        for (int i = 0, j = 0; i < problem.fixed_dof.size(); i++) {
            if (!problem.fixed_dof(i)) {
                free_dof(j++) = i;
            }
        }

        // Initalize the working variables
        Eigen::VectorXd g = problem.eval_grad_f(x),
                        delta_x = Eigen::VectorXd::Zero(problem.num_vars),
                        g_free; // subset of g for free degrees of freedom
        igl::slice(g, free_dof, 1, g_free); // Initialize g_free
        Eigen::MatrixXd H, H_free;
        double gamma;

        int iter = 0;
        while (iter <= max_iter && g_free.squaredNorm() > absolute_tolerance) {
            // Compute the full hessian
            H = problem.eval_hessian_f(x);
            H.diagonal().array() += mu;
            // Remove rows and columns of fixed dof
            igl::slice(H, free_dof, free_dof, H_free);

            // Store the free dof back into delta_x
            igl::slice_into(Eigen::VectorXd(H_free.lu().solve(-g_free)),
                free_dof, 1, delta_x);

            // Perform a line search along delta x to stay in the feasible realm
            if (!line_search(x, delta_x, problem.func_f(), gamma,
                    line_search_tolerance)) {
                break; // Newton step unsuccessful
            }

            x += gamma * delta_x; // Update x

            g = problem.eval_grad_f(x); // Recompute the gradient
            // Remove rows of fixed dof
            igl::slice(g, free_dof, 1, g_free);

            // Save intermedtiate results
            // settings.intermediate_cb(
            //     x, problem.f(x), Eigen::VectorXd::Zero(x.size()), 0, iter);

            iter++; // Increase iteration counter
        }

        // if (settings.verbosity > 0) {
        std::cout << "took " << iter << " iterations." << std::endl;
        // }

        return OptimizationResults(
            x, problem.eval_f(x), g.squaredNorm() <= absolute_tolerance);
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
