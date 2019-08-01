// Functions for optimizing functions.
#include "newton_solver.hpp"

#include <igl/slice.h>
#include <igl/slice_into.h>

#include <solvers/line_search.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {
namespace opt {

    NewtonSolver::NewtonSolver()
        : NewtonSolver("newton_solver")
    {
    }
    NewtonSolver::NewtonSolver(const std::string& name)
        : OptimizationSolver(name)
        , iteration_number(0)
        , absolute_tolerance(1e-5)
        , min_step_length(1e-12)
    {
    }

    NewtonSolver::~NewtonSolver() {}

    void NewtonSolver::settings(const nlohmann::json& json)
    {
        OptimizationSolver::settings(json);
        absolute_tolerance = json["absolute_tolerance"].get<double>();
        min_step_length = json["min_step_length"].get<double>();
    }

    nlohmann::json NewtonSolver::settings() const
    {
        nlohmann::json json = OptimizationSolver::settings();
        json["absolute_tolerance"] = absolute_tolerance;
        json["min_step_length"] = min_step_length;
        return json;
    }

    const static char* NEWTON_STEP_LOG
        = "solver=newton_solver iter={:d} action=begin_step";
    const static char* NEWTON_DIRECTION_LOG
        = "solver=newton_solver iter={:d} action=newton_line_search "
          "status=failure message=\"reverting to gradient descent\"";
    const static char* NEWTON_GRADIENT_LOG
        = "solver=newton_solver iter={:d} action=gradient_line_search "
          "status=failure message=\"no failsafe - break!\"";
    const static char* NEWTON_END_LOG = "solver=newton_solver action=END";
    const static char* NEWTON_BEGIND_LOG = "solver=newton_solver action=BEGIN";

    OptimizationResults NewtonSolver::solve(OptimizationProblem& problem)
    {
        NAMED_PROFILE_POINT(
            "newton_solver", SOLVER_STEP);
        NAMED_PROFILE_POINT(
            "newton_solver__eval_f", EVAL_F);
        NAMED_PROFILE_POINT(
            "newton_solver__slice", SLICE);
        NAMED_PROFILE_POINT(
            "newton_solver__compute_newton_direction", COMPUTE_DIRECTION);
        NAMED_PROFILE_POINT(
            "newton_solver__newton_line_search", NEWTON_LINE_SEARCH);
        NAMED_PROFILE_POINT(
            "newton_solver__gradient_line_search", GRADIENT_LINE_SEARCH);
        NAMED_PROFILE_POINT(
            "newton_solver__assert_constraint", ASSERT_CONSTRAINT);

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

        double eps = -1;
        if (problem.has_barrier_constraint()) {
            eps = problem.get_barrier_epsilon();
        }
        spdlog::trace(NEWTON_BEGIND_LOG);
        iteration_number = 0;
        std::string exit_reason = "exceeded the maximum allowable iterations";
        do {
            PROFILE_START(SOLVER_STEP);

            PROFILE_START(EVAL_F)
            // Compute the gradient and hessian
            double fx;
            problem.eval_f_and_fdiff(x, fx, gradient, hessian);
            PROFILE_END(EVAL_F)

            PROFILE_START(SLICE)
            // Remove rows of fixed dof
            igl::slice(gradient, free_dof, gradient_free);
            PROFILE_END(SLICE)

            spdlog::trace("{} eps={:e} x={} fx={:e} gradient_free={}",
                fmt::format(NEWTON_STEP_LOG, iteration_number + 1), eps,
                log::fmt_eigen(x), fx, log::fmt_eigen(gradient_free));

            if (gradient_free.squaredNorm() <= absolute_tolerance) {
                exit_reason = "found a local optimum";
                PROFILE_END(SOLVER_STEP);
                break;
            }

            PROFILE_START(SLICE)
            // Remove rows and columns of fixed dof of the hessian
            igl::slice(hessian, free_dof, free_dof, hessian_free);
            PROFILE_END(SLICE)

            PROFILE_START(COMPUTE_DIRECTION)
            bool found_direction = compute_free_direction(
                gradient_free, hessian_free, delta_x, true);
            PROFILE_END(COMPUTE_DIRECTION)

            if (!found_direction) {
                exit_reason = "newton direction solve failed";
                PROFILE_END(SOLVER_STEP);
                break;
            }

            step_length = std::min(1.0, 2.0 * step_length);
            // to restore when line-search fails
            double step_length_aux = step_length;

            PROFILE_START(NEWTON_LINE_SEARCH)
            // Perform a line search along Δx, and stay in the feasible realm
            std::string message = "";
            bool found_step_length
                = constrained_line_search(x, delta_x, problem.func_f(),
                    gradient, constraint, step_length, min_step_length);
            PROFILE_END(NEWTON_LINE_SEARCH)

            // Revert to gradient descent if the newton direction fails
            if (!found_step_length) {

                spdlog::warn(NEWTON_DIRECTION_LOG, iteration_number + 1);

                // replace delta_x by gradient direction
                PROFILE_START(SLICE)
                igl::slice_into(-gradient_free, free_dof, 1, delta_x);
                step_length = step_length_aux;
                PROFILE_END(SLICE)

                PROFILE_START(GRADIENT_LINE_SEARCH)
                found_step_length
                    = constrained_line_search(x, delta_x, problem.func_f(),
                        gradient, constraint, step_length, min_step_length);
                PROFILE_END(GRADIENT_LINE_SEARCH)

                if (!found_step_length) {
                    spdlog::warn(NEWTON_GRADIENT_LOG, iteration_number + 1);
                    exit_reason = "line-search failed";
                    PROFILE_END(SOLVER_STEP);
                    break;
                }
            }

            x += step_length * delta_x;
            PROFILE_START(ASSERT_CONSTRAINT)
            assert(constraint(x));
            PROFILE_END(ASSERT_CONSTRAINT)

            problem.eval_intermediate_callback(x);
            PROFILE_END(SOLVER_STEP);

        } while (++iteration_number <= max_iterations);

        spdlog::trace(
            "{} total_iter={:d} gradient_free_eq_norm={:e} message='{}'",
            NEWTON_END_LOG, iteration_number, gradient_free.squaredNorm(),
            exit_reason);

        return OptimizationResults(x, problem.eval_f(x),
            gradient_free.squaredNorm() <= absolute_tolerance);
    }

    // Solve for the newton direction for the limited degrees of freedom.
    // The fixed dof of x will have a delta_x of zero.
    bool NewtonSolver::compute_free_direction(
        const Eigen::VectorXd& gradient_free,
        const Eigen::SparseMatrix<double>& hessian_free,
        Eigen::VectorXd& delta_x,
        bool make_psd)
    {
        Eigen::VectorXd delta_x_free;
        bool success = compute_direction(
            gradient_free, hessian_free, delta_x_free, make_psd);
        if (success) {
            // Store the free dof back into delta_x
            igl::slice_into(delta_x_free, free_dof, delta_x);
        }
        return success;
    }

    // Solve for the Newton direction (Δx = -H^{-1}∇f).
    // Return true if the solve was successful.
    bool NewtonSolver::compute_direction(const Eigen::VectorXd& gradient,
        const Eigen::SparseMatrix<double>& hessian,
        Eigen::VectorXd& delta_x,
        bool make_psd)
    {
        bool solve_success = false;

        // TODO: Can we use a better solver than LU?
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(hessian);
        if (solver.info() != Eigen::Success) {
            spdlog::warn(
                "solver=newton iter={:d} failure=\"sparse decomposition of the hessian failed\" failsafe=none",
                iteration_number);
        } else {
            delta_x = solver.solve(-gradient);
            if (solver.info() != Eigen::Success) {
                spdlog::warn(
                    "solver=newton iter={:d} failure=\"sparse solve for newton direction failed\" failsafe=none",
                    iteration_number);
            } else {
                solve_success = true;
            }
        }

        if (solve_success && make_psd && delta_x.transpose() * gradient >= 0) {
            // If delta_x is not a descent direction then we want to modify the
            // hessian to be diagonally dominant with positive elements on the
            // diagonal (positive definite). We do this by adding μI to the
            // hessian. This can result in doing a step of gradient descent.
            Eigen::SparseMatrix<double> psd_hessian = hessian;
            double mu = make_matrix_positive_definite(psd_hessian);
            spdlog::debug(
                "solver=newton iter={:d} failure=\"newton direction not descent direction\" failsafe=\"H += μI\" μ={:g}",
                iteration_number, mu);
            solve_success
                = compute_direction(gradient, psd_hessian, delta_x, false);
            if (delta_x.transpose() * gradient >= 0) {
                spdlog::warn(
                    "solver=newton iter={:d} failure=\"newton direction not descent direction\" failsafe=none dir_dot_grad={:g}",
                    iteration_number, delta_x.transpose() * gradient);
            }
        }

        return solve_success;
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
