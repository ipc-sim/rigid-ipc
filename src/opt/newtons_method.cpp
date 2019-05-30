// Functions for optimizing functions.
// Includes Newton's method with and without constraints.
#include <opt/newtons_method.hpp>

#include <cmath>
#include <iostream>

#include <Eigen/Sparse>
#include <igl/slice.h>
#include <igl/slice_into.h>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {
namespace opt {

    OptimizationResults newtons_method(OptimizationProblem& problem,
        const Eigen::VectorXi& free_dof, const double absolute_tolerance,
        const double line_search_tolerance, const int max_iter, const double mu)
    {
        // Initalize the working variables
        Eigen::VectorXd x = problem.x0;

        // subset of g for free degrees of freedom
        Eigen::VectorXd gradient, gradient_free;
        Eigen::SparseMatrix<double> hessian, hessian_free;

        Eigen::VectorXd delta_x = Eigen::VectorXd::Zero(problem.num_vars);
        double step_length = 1.0;

        // TODO: Can we use a better solver than LU?
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

        auto constraint = [&problem](const Eigen::VectorXd& x) {
            return problem.eval_f(x) < std::numeric_limits<double>::infinity();
        };

        int iter = 0;
        do {
            // Compute the gradient and hessian
            double _;
            problem.eval_f_and_fdiff(x, _, gradient, hessian);
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
            // step_length = std::min(1.0, 2 * step_length);
            step_length = 1.0;
            // problem.enable_line_search_mode(x + delta_x);
            bool found_step_length
                = constrained_line_search(x, delta_x, problem.func_f(),
                    gradient, constraint, step_length, line_search_tolerance);
            // problem.disable_line_search_mode();
            if (!found_step_length) {
                break; // Newton step unsuccessful
            }

            x += step_length * delta_x; // Update x
            assert(constraint(x));

            problem.eval_intermediate_callback(x);
        } while (++iter <= max_iter);

        std::ostringstream string_stream;
        if (gradient_free.squaredNorm() <= absolute_tolerance) {
            string_stream << "found a local optimum";
        } else if (iter > max_iter) {
            string_stream << "exceeded the maximum allowable iterations ("
                          << max_iter << ")";
        } else {
            string_stream << "line-search failed";
        }
        spdlog::debug(
            "method=newtons_method total_iter={:d} exit_code=\"{:s}\"", iter,
            string_stream.str());

        return OptimizationResults(x, problem.eval_f(x),
            gradient_free.squaredNorm() <= absolute_tolerance);
    }

    // Search along a search direction to find a scalar step_length in [0, 1]
    // such that f(x + step_length * dir) ≤ f(x).
    bool line_search(const Eigen::VectorXd& x, const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        double& step_length, const double min_step_length)
    {
        return line_search(x, dir, f, Eigen::VectorXd::Zero(dir.size()),
            step_length, min_step_length);
    }

    // Search along a search direction to find a scalar step_length in [0, 1]
    // such that f(x + step_length * dir) ≤ f(x).
    bool line_search(const Eigen::VectorXd& x, const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const Eigen::VectorXd& grad_fx, double& step_length,
        const double min_step_length)
    {
        return constrained_line_search(
            x, dir, f, grad_fx, [](const Eigen::VectorXd&) { return true; },
            step_length, min_step_length);
    }

    // Search along a search direction to find a scalar step_length in [0, 1]
    // such that f(x + step_length * dir) ≤ f(x).
    // TODO: Filter the dof that violate the constraints. These are the indices
    // i where ϕ([g(x)]_i) = ∞.
    bool constrained_line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const Eigen::VectorXd& grad_fx,
        const std::function<bool(const Eigen::VectorXd&)>& constraint,
        double& step_length, const double min_step_length)
    {
        const double fx = f(x); // Function value we want to beat

        // Wolfe conditions:
        // Armijo rule
        double wolfe_c1 = 1e-4;
        const double wolfe1 = wolfe_c1 * dir.transpose() * grad_fx;
        auto armijo_rule = [&]() {
            return f(x + step_length * dir) <= fx + step_length * wolfe1;
        };
        // Curvature condition
        // double wolfe_c2 = 0.9;
        // const double wolfe2 = -wolfe_c2 * dir.transpose() * grad_fx;
        auto curvature_consition = [&]() {
            // return -dir.transpose * grad_f(x + step_length * dir) <= wolfe2;
            return true;
        };

        while (step_length >= min_step_length) {
            if (armijo_rule() && curvature_consition()
                && constraint(x + step_length * dir)) {
                return true;
            }
            step_length /= 2.0;
        }
        return false;
    }

} // namespace opt
} // namespace ccd
