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
#ifdef PROFILE_FUNCTIONS
long number_of_update_solves = 0;
double time_spent_solving_for_update = 0;
#endif

namespace ccd {
namespace opt {

    OptimizationResults sparse_newtons_method(OptimizationProblem& problem,
        const Eigen::VectorXi& free_dof, const double absolute_tolerance,
        const double line_search_tolerance, const int max_iter, bool verbose,
        const double mu)
    {
        // Initalize the working variables
        Eigen::VectorXd x = problem.x0;
        Eigen::VectorXd delta_x = Eigen::VectorXd::Zero(problem.num_vars),
                        g_free; // subset of g for free degrees of freedom
        igl::slice(problem.eval_grad_f(x), free_dof, 1, g_free);
        Eigen::SparseMatrix<double> H_free;
        double gamma;

        int iter = 0;
        while (iter <= max_iter && g_free.squaredNorm() > absolute_tolerance) {
            // Compute the full hessian and remove rows and columns of fixed dof
            igl::slice(
                problem.eval_hessian_f_sparse(x), free_dof, free_dof, H_free);
            // Add a small value to the diagonal to prevent singular matrices
            for (int i = 0; i < H_free.rows(); i++) {
                H_free.coeffRef(i, i) += mu;
            }

#ifdef PROFILE_FUNCTIONS
            number_of_update_solves++;
            igl::Timer timer;
            timer.start();
#endif
            // TODO: Can we use a better solver than LU?
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(H_free);
            if (solver.info() != Eigen::Success) {
                std::cerr << "Sparse decomposition of the hessian failed!"
                          << std::endl;
                break;
            }
            // Store the free dof back into delta_x
            igl::slice_into(
                Eigen::VectorXd(solver.solve(-g_free)), free_dof, 1, delta_x);
            if (solver.info() != Eigen::Success) {
                std::cerr << "Sparse solve for the update failed!" << std::endl;
                break;
            }

#ifdef PROFILE_FUNCTIONS
            timer.stop();
            time_spent_solving_for_update += timer.getElapsedTime();
#endif

            // Perform a line search along delta x to stay in the feasible realm
            if (!line_search(x, delta_x, problem.func_f(), gamma,
                    line_search_tolerance)) {
                break; // Newton step unsuccessful
            }

            x += gamma * delta_x; // Update x

            // Recompute the gradient
            igl::slice(problem.eval_grad_f(x), free_dof, 1,
                g_free); // Remove rows of fixed dof


            problem.eval_intermediate_callback(x);

            iter++; // Increase iteration counter
        }

        if (verbose > 0) {
            std::cout << "took " << iter << " iterations." << std::endl;
        }

        return OptimizationResults(
            x, problem.eval_f(x), g_free.squaredNorm() <= absolute_tolerance);
    }

    OptimizationResults newtons_method(OptimizationProblem& problem,
        const Eigen::VectorXi& free_dof, const double absolute_tolerance,
        const double line_search_tolerance, const int max_iter, bool verbose,
        const double mu)
    {
        Eigen::VectorXd x = problem.x0;

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

#ifdef PROFILE_FUNCTIONS
            number_of_update_solves++;
            igl::Timer timer;
            timer.start();
#endif
            // Try to solve the hessian as a positive or negative semidefinite
            // matrix. I do not know if this is true. I do know that for
            // f(x) = 0.5 * ||U-U0||^2, hessian(f) = I.
            auto H_free_decomp = H_free.ldlt();
            if (H_free_decomp.info() == Eigen::NumericalIssue) {
                std::cerr
                    << "Hessian is not positive or negative semi-definite!"
                    << std::endl;
                // Store the free dof back into delta_x
                igl::slice_into(Eigen::VectorXd(H_free.lu().solve(-g_free)),
                    free_dof, 1, delta_x);
            } else {
                // Store the free dof back into delta_x
                igl::slice_into(Eigen::VectorXd(H_free_decomp.solve(-g_free)),
                    free_dof, 1, delta_x);
            }

#ifdef PROFILE_FUNCTIONS
            timer.stop();
            time_spent_solving_for_update += timer.getElapsedTime();
#endif

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

        if (verbose > 0) {
            std::cout << "took " << iter << " iterations." << std::endl;
        }

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
