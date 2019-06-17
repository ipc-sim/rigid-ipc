// Functions for optimizing functions.
// Includes Newton's method with and without constraints.
#include "newton_solver.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/slice.h>
#include <igl/slice_into.h>

#include <solvers/line_search.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {
namespace opt {

    NewtonSolver::NewtonSolver()
        : free_dof()
        , absolute_tolerance(1e-5)
        , line_search_tolerance(1e-12)
        , max_iterations(3000)
    {
    }

    NewtonSolver::~NewtonSolver() {}

    // Initialize free_dof with indices of dof that are not fixed.
    void NewtonSolver::init_free_dof(Eigen::MatrixXb is_dof_fixed)
    {
        free_dof = Eigen::VectorXi(is_dof_fixed.size() - is_dof_fixed.count());
        for (int i = 0, j = 0; i < is_dof_fixed.size(); i++) {
            if (!is_dof_fixed(i)) {
                free_dof(j++) = i;
            }
        }
    }

    OptimizationResults NewtonSolver::solve(OptimizationProblem& problem)
    {
        // Initalize the working variables
        Eigen::VectorXd x = problem.x0;
        Eigen::VectorXd gradient, gradient_free;
        Eigen::SparseMatrix<double> hessian, hessian_free;
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
            double _;
            problem.eval_f_and_fdiff(x, _, gradient, hessian);

            // Remove rows of fixed dof
            igl::slice(gradient, free_dof, gradient_free);

            if (gradient_free.squaredNorm() <= absolute_tolerance) {
                exit_reason = "found a local optimum";
                break;
            }

            // Remove rows and columns of fixed dof of the hessian
            igl::slice(hessian, free_dof, free_dof, hessian_free);

            if (!compute_free_direction(
                    gradient_free, hessian_free, delta_x, true)) {
                exit_reason = "newton direction solve failed";
                break;
            }

            // Perform a line search along Δx to stay in the feasible realm
            bool found_step_length;
            // step_length = std::min(1.0, 2 * step_length);
            step_length = 1.0;
            // problem.enable_line_search_mode(x + delta_x);
            found_step_length
                = constrained_line_search(x, delta_x, problem.func_f(),
                    gradient, constraint, step_length, line_search_tolerance);
            // problem.disable_line_search_mode();

            if (!found_step_length) {
                // Revert to gradient descent if the newton direction
                // fails
                spdlog::warn("method=newtons_method iter={:d} failure=\"line "
                             "search using newton direction\" "
                             "failsafe=\"reverting to gradient descent\"",
                    iter + 1);
                igl::slice_into(-gradient_free, free_dof, 1, delta_x);

                // step_length = std::min(1.0, 2 * step_length);
                step_length = 1.0;
                // problem.enable_line_search_mode(x + delta_x);
                found_step_length = constrained_line_search(x, delta_x,
                    problem.func_f(), gradient, constraint, step_length,
                    line_search_tolerance);
                // problem.disable_line_search_mode();

                if (!found_step_length) {
                    spdlog::warn(
                        "method=newtons_method iter={:d} failure=\"line search "
                        "using gradient descent\" failsafe=none",
                        iter + 1);
                    exit_reason = "line-search failed";
                    break;
                }
            }

            x += step_length * delta_x; // Update x
            assert(constraint(x));

            problem.eval_intermediate_callback(x);
        } while (++iter <= max_iterations);

        spdlog::trace("method=newtons_method total_iter={:d} "
                      "exit_reason=\"{:s}\" sqr_norm_grad={:g}",
            iter, exit_reason, gradient_free.squaredNorm());

        return OptimizationResults(x, problem.eval_f(x),
            gradient_free.squaredNorm() <= absolute_tolerance);
    }

    // Solve for the newton direction for the limited degrees of freedom.
    // The fixed dof of x will have a delta_x of zero.
    bool NewtonSolver::compute_free_direction(
        const Eigen::VectorXd& gradient_free,
        const Eigen::SparseMatrix<double>& hessian_free,
        Eigen::VectorXd& delta_x, bool make_psd)
    {
        Eigen::VectorXd delta_x_free;
        if (!compute_direction(
                gradient_free, hessian_free, delta_x_free, make_psd)) {
            return false;
        }
        // Store the free dof back into delta_x
        igl::slice_into(delta_x_free, free_dof, delta_x);
        return true;
    }

    // Solve for the Newton direction (Δx = -H^{-1}∇f).
    // Return true if the solve was successful.
    bool NewtonSolver::compute_direction(const Eigen::VectorXd& gradient,
        const Eigen::SparseMatrix<double>& hessian, Eigen::VectorXd& delta_x,
        bool make_psd)
    {
        bool solve_success = false;

        // Profile the performance of the solve.
        PROFILE(
            // clang-format off
            // TODO: Can we use a better solver than LU?
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

            solver.compute(hessian);
            if (solver.info() != Eigen::Success) {
                spdlog::warn(
                    "method=newtons_method failure=\"sparse decomposition "
                    "of the hessian failed\" failsafe=none");
            } else {
                delta_x = solver.solve(-gradient);
                if (solver.info() != Eigen::Success) {
                    spdlog::warn("method=newtons_method failure=\"sparse solve "
                                 "for newton direction failed\" failsafe=none");
                } else {
                    solve_success = true;
                }
            },
            // clang-format on
            ProfiledPoint::UPDATE_SOLVE);

        if (solve_success && make_psd && delta_x.transpose() * gradient >= 0) {
            // If delta_x is not a descent direction then we want to modify the
            // hessian to be diagonally dominant with positive elements on the
            // diagonal (positive definite). We do this by adding μI to the
            // hessian. This can result in doing a step of gradient descent.
            Eigen::SparseMatrix<double> psd_hessian = hessian;
            double mu = make_matrix_positive_definite(psd_hessian);
            spdlog::debug("method=newtons_method failure=\"newton direction "
                          "not descent direction\" failsafe=\"H += μI\" μ={:g}",
                mu);
            solve_success
                = compute_direction(gradient, psd_hessian, delta_x, false);
            if (delta_x.transpose() * gradient >= 0) {
                spdlog::warn(
                    "method=newtons_method failure=\"newton direction not "
                    "descent direction\" failsafe=none dir_dot_grad={:g}",
                    delta_x.transpose() * gradient);
            }
        }

        return solve_success; // Was the solve successful
    }

    // Make the matrix positive definite (x^T A x > 0).
    double make_matrix_positive_definite(Eigen::SparseMatrix<double>& A)
    {
        // Conservative way of making A PSD by making it diagonally dominant
        // with all positive diagonal entries
        Eigen::SparseMatrix<double> I(A.rows(), A.rows());
        I.setIdentity();
        // Entries along the diagonal of A
        Eigen::VectorXd diag = Eigen::VectorXd::Zero(A.rows());
        // Sum of columns per row not including the diagonal entry
        // (∑_{i≠j}|a_{ij}|)
        Eigen::VectorXd sum_row = Eigen::VectorXd::Zero(A.rows());

        // Loop over elements adding them to the appropriate vector above
        for (int k = 0; k < A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it;
                 ++it) {
                if (it.row() == it.col()) { // Diagonal element
                    diag(it.row()) = it.value();
                } else { // Non-diagonal element
                    sum_row(it.row()) += abs(it.value());
                }
            }
        }
        // Take max to ensure all diagonal elements are dominant
        double mu = std::max((sum_row - diag).maxCoeff(), 0.0);
        A += mu * I;
        return mu;
    }

} // namespace opt
} // namespace ccd
