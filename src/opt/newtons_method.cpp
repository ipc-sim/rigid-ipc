// Functions for optimizing functions.
// Includes Newton's method with and without constraints.
#include <cmath>
#include <iostream>

#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <igl/slice.h>
#include <igl/slice_into.h>

#include <opt/newtons_method.hpp>

#include <profiler.hpp>

namespace ccd {
namespace opt {

    OptimizationResults newtons_method(OptimizationProblem& problem,
        const Eigen::VectorXi& free_dof, const double absolute_tolerance,
        const double line_search_tolerance, const int max_iter, bool verbose,
        const double mu)
    {
        // Initalize the working variables
        Eigen::VectorXd x = problem.x0;

        // subset of g for free degrees of freedom
        double value;
        Eigen::VectorXd gradient, gradient_free;
        Eigen::SparseMatrix<double> hessian, hessian_free;

        Eigen::VectorXd delta_x = Eigen::VectorXd::Zero(problem.num_vars);
        double gamma;

        // TODO: Can we use a better solver than LU?
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

        int iter = 0;
        problem.eval_intermediate_callback(x);

        do {
            // Compute the gradient and hessian
            problem.eval_f_and_fdiff(x, value, gradient, hessian);
            // Remove rows of fixed dof
            igl::slice(gradient, free_dof, 1, gradient_free);
            if (gradient_free.squaredNorm() <= absolute_tolerance) {
                break;
            }

            // Compute the full hessian and remove rows and columns of fixed dof
            igl::slice(hessian, free_dof, free_dof, hessian_free);
            // Add a small value to the diagonal to prevent singular matrices
            for (int i = 0; i < hessian_free.rows(); i++) {
                hessian_free.coeffRef(i, i) += mu;
            }

            PROFILE(
                solver.compute(hessian_free);
                if (solver.info() != Eigen::Success) {
                    std::cerr << "Sparse decomposition of the hessian failed!"
                              << std::endl;
                    break;
                }

                // Store the free dof back into delta_x
                igl::slice_into(Eigen::VectorXd(solver.solve(-gradient_free)),
                    free_dof, 1, delta_x);
                if (solver.info() != Eigen::Success) {
                    std::cerr << "Sparse solve for the update failed!"
                              << std::endl;
                    break;
                },
                ProfiledPoint::UPDATE_SOLVE);

            // Perform a line search along delta x to stay in the feasible realm
            if (!line_search(x, delta_x, problem.func_f(), gamma,
                    line_search_tolerance)) {
                break; // Newton step unsuccessful
            }

            x += gamma * delta_x; // Update x

            problem.eval_intermediate_callback(x);
        } while (++iter <= max_iter);

        if (verbose > 0) {
            std::cout << "took " << iter << " iterations." << std::endl;
        }

        return OptimizationResults(x, problem.eval_f(x),
            gradient_free.squaredNorm() <= absolute_tolerance);
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
