// Functions for optimizing functions.
#include "newton_solver.hpp"

#include <igl/slice.h>
#include <igl/slice_into.h>

#include <constants.hpp>
#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {
namespace opt {

    NewtonSolver::NewtonSolver()
        : max_iterations(1000)
        , iteration_number(0)
    {
    }

    void NewtonSolver::settings(const nlohmann::json& json)
    {
        max_iterations = json["max_iterations"].get<int>();
    }

    nlohmann::json NewtonSolver::settings() const
    {
        nlohmann::json json;
        json["max_iterations"] = max_iterations;
        return json;
    }

    void NewtonSolver::init_solve(const Eigen::VectorXd& x0)
    {
        assert(problem_ptr != nullptr);
        reset_stats();
        free_dof = init_free_dof(problem_ptr->is_dof_fixed());
    }

    Eigen::VectorXi init_free_dof(const Eigen::VectorXb& is_dof_fixed)
    {
        Eigen::VectorXi free_dof(is_dof_fixed.size() - is_dof_fixed.count());
        for (int i = 0, j = 0; i < is_dof_fixed.size(); i++) {
            if (!is_dof_fixed(i)) {
                free_dof(j++) = i;
            }
        }
        return free_dof;
    }

    std::string NewtonSolver::stats()
    {
        return fmt::format(
            "total_newton_steps={:d} total_ls_steps={:d} count_fx={:d} "
            "count_grad={:d} count_hess={:d} count_ccd={:d}",
            newton_iterations, ls_iterations, num_fx, num_grad_fx,
            num_hessian_fx, num_collision_check);
    }

    void NewtonSolver::reset_stats()
    {
        num_fx = 0;
        num_grad_fx = 0;
        num_hessian_fx = 0;
        num_collision_check = 0;
        ls_iterations = 0;
        newton_iterations = 0;
    }

    OptimizationResults NewtonSolver::solve(const Eigen::VectorXd& x0)
    {
        assert(problem_ptr != nullptr);
        // Initalize the working variables
        Eigen::VectorXd x = x0;
        Eigen::VectorXd gradient, gradient_free;
        Eigen::SparseMatrix<double> hessian, hessian_free;

        double step_length = 1.0;

        spdlog::debug("solver={} action=BEGIN", name());

        std::string exit_reason = "exceeded the maximum allowable iterations";

        Eigen::VectorXd direction(problem_ptr->num_vars());
        Eigen::VectorXd grad_direction(problem_ptr->num_vars());
        Eigen::VectorXd direction_free(free_dof.size());

        bool success = false;

        for (iteration_number = 0; iteration_number < max_iterations;
             iteration_number++) {

            double fx = problem_ptr->compute_objective(x, gradient, hessian);

            newton_iterations++;
            num_fx++;
            num_grad_fx++;
            num_hessian_fx++;

            // Remove rows and cols of fixed dof
            igl::slice(gradient, free_dof, gradient_free);
            igl::slice(hessian, free_dof, free_dof, hessian_free);
            compute_direction(
                gradient_free, hessian_free, direction_free, /*make_psd=*/true);

            ///////////////////////////////////////////////////////////////////
            // Line search over newton direction
            // get grad direction for lineseach
            grad_direction.setZero();
            igl::slice_into(gradient_free, free_dof, grad_direction);
            direction.setZero();
            igl::slice_into(direction_free, free_dof, direction);

            // check for newton termination
            if (converged(gradient_free, direction_free)) {
                exit_reason = "found a local optimum with newton dir";
                success = true;
                break;
            }

            step_length = 1.0;
            bool found_newton_step =
                line_search(x, direction, fx, grad_direction, step_length);
            ///////////////////////////////////////////////////////////////////

            ///////////////////////////////////////////////////////////////////
            // When newton direction fails, revert to gradient descent
            if (!found_newton_step) {
                spdlog::warn(
                    "solver={} iter={:d} failure=\"newton line-search\" "
                    "failsafe=\"gradient descent\"",
                    name(), iteration_number + 1);

                ///////////////////////////////////////////////////////////////
                // Line search over -gradient direction
                // replace delta_x by gradient direction
                direction.setZero();
                direction_free = -gradient_free;
                igl::slice_into(direction_free, free_dof, 1, direction);

                // check for newton termination again
                if (converged(gradient_free, direction_free)) {
                    exit_reason = "found a local optimum with -grad dir";
                    success = true;
                    break;
                }
                step_length = 1;
                bool found_gradient_step =
                    line_search(x, direction, fx, grad_direction, step_length);
                ///////////////////////////////////////////////////////////////

                // When gradient direction fails, exit
                if (!found_gradient_step) {
                    spdlog::error(
                        "solver={} iter={:d} failure=\"gradient line-search\" "
                        "failsafe=\"none\"",
                        name(), iteration_number + 1);
                    exit_reason = "line-search failed";
                    break;
                }
            }
            ///////////////////////////////////////////////////////////////////

            spdlog::debug(
                "solver={} iter={:d} step_length={:g}", name(),
                iteration_number, step_length);

            Eigen::VectorXd xk = x + step_length * direction;
            assert(!problem_ptr->has_collisions(x, xk));

            post_step_update(x, xk);

            x = xk;
        } // end for loop

        spdlog::info(
            "solver={} action=END total_iter={:d} exit_reason=\"{}\"", name(),
            iteration_number, exit_reason);

        return OptimizationResults(
            x, problem_ptr->compute_objective(x), success, true);
    }

    bool NewtonSolver::line_search(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const double fx,
        const Eigen::VectorXd& grad_fx,
        double& alpha)
    {
        NAMED_PROFILE_POINT("line_search", LINE_SEARCH);
        PROFILE_START(LINE_SEARCH);

        bool success = false;
        int num_it = 0;
        double lower_bound = line_search_lower_bound() / -grad_fx.dot(dir);
        assert(std::isfinite(lower_bound));

        while (alpha > lower_bound) {
            num_it++;        // Count the number of iterations
            ls_iterations++; // Count the gloabal number of iterations

            // Compute the next variable
            Eigen::VectorXd xi = x + alpha * dir;

            // Check for collisions between newton updates
            num_collision_check++; // Count the number of collision checks
            if (!problem_ptr->has_collisions(x, xi)) {
                num_fx++; // Count the number of objective computations
                if (problem_ptr->compute_objective(xi) < fx) {
                    success = true;
                    break; // while loop
                }
            }

            // Try again with a smaller alpha
            alpha /= 2.0;
        }

        PROFILE_MESSAGE(
            LINE_SEARCH,
            fmt::format(
                "success,{},it,{},dir,{:10e}", success, num_it, dir.norm()))
        PROFILE_END(LINE_SEARCH);

        if (!success) {
            spdlog::warn(
                "solver={} iter={:d} failure=\"line-search α ≤ {:g}\"", name(),
                iteration_number, lower_bound);
        }

        return success;
    }

    bool NewtonSolver::compute_direction(
        const Eigen::VectorXd& gradient,
        const Eigen::SparseMatrix<double>& hessian,
        Eigen::VectorXd& direction,
        bool make_psd)
    {
        // Solve for the Newton direction (Δx = -H⁻⅟∇f).
        // Return true if the solve was successful.
        bool solve_success = false;

        // TODO: Can we use a better solver than LU?
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(hessian);
        if (solver.info() == Eigen::Success) {
            direction = solver.solve(-gradient);
            if (solver.info() == Eigen::Success) {
                solve_success = true;
            } else {
                spdlog::warn(
                    "solver={} iter={:d} failure=\"sparse solve for newton "
                    "direction\" failsafe=\"gradient descent\"",
                    name(), iteration_number);
            }
        } else {
            spdlog::warn(
                "solver={} iter={:d} failure=\"sparse decomposition of the "
                "hessian\" failsafe=\"gradient descent\"",
                name(), iteration_number);
        }

        if (!solve_success) {
            direction = -gradient;
        }

        if (solve_success && make_psd && direction.dot(gradient) > 0) {
            // If delta_x is not a descent direction then we want to modify the
            // hessian to be diagonally dominant with positive elements on the
            // diagonal (positive definite). We do this by adding μI to the
            // hessian. This can result in doing a step of gradient descent.
            Eigen::SparseMatrix<double> psd_hessian = hessian;
            double mu = make_matrix_positive_definite(psd_hessian);
            spdlog::debug(
                "solver={} iter={:d} failure=\"newton direction not descent "
                "direction\" failsafe=\"H += μI\" μ={:g}",
                name(), iteration_number, mu);
            solve_success =
                compute_direction(gradient, psd_hessian, direction, false);
            double dir_dot_grad = direction.dot(gradient);
            if (dir_dot_grad > 0) {
                spdlog::warn(
                    "solver={} iter={:d} failure=\"adjusted newton "
                    "direction not descent direction\" failsafe=\"gradient "
                    "descent\" dir_dot_grad={:g}",
                    name(), iteration_number, dir_dot_grad);
                direction = -gradient;
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
