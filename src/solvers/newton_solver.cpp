// Functions for optimizing functions.
#include "newton_solver.hpp"

#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/writeOBJ.h>

#include <constants.hpp>
#include <geometry/distance.hpp>
#include <logger.hpp>
#include <physics/simulation_problem.hpp>
#include <profiler.hpp>

namespace ccd {
namespace opt {

    NewtonSolver::NewtonSolver()
        : max_iterations(1000)
        , iteration_number(0)
        , convergence_criteria(ConvergenceCriteria::ENERGY)
    {
    }

    void NewtonSolver::settings(const nlohmann::json& json)
    {
        max_iterations = json["max_iterations"].get<int>();
        convergence_criteria =
            json["convergence_criteria"].get<ConvergenceCriteria>();
        reset_stats();
    }

    nlohmann::json NewtonSolver::settings() const
    {
        nlohmann::json json;
        json["max_iterations"] = max_iterations;
        json["convergence_criteria"] = convergence_criteria;
        return json;
    }

    void NewtonSolver::init_solve(const Eigen::VectorXd& x0)
    {
        assert(problem_ptr != nullptr);
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
            "total_newton_steps={:d} total_ls_steps={:d} "
            "num_newton_ls_fails={:d} num_grad_ls_fails={:d} count_fx={:d} "
            "count_grad={:d} count_hess={:d} count_ccd={:d}",
            newton_iterations, ls_iterations, num_newton_ls_fails,
            num_grad_ls_fails, num_fx, num_grad_fx, num_hessian_fx,
            num_collision_check);
    }

    void NewtonSolver::reset_stats()
    {
        num_fx = 0;
        num_grad_fx = 0;
        num_hessian_fx = 0;
        num_collision_check = 0;
        ls_iterations = 0;
        newton_iterations = 0;
        num_newton_ls_fails = 0;
        num_grad_ls_fails = 0;
    }

    bool NewtonSolver::converged() const
    {
        switch (convergence_criteria) {
        case ConvergenceCriteria::VELOCITY: {
            Eigen::MatrixXd V_prev = problem_ptr->world_vertices(x);
            Eigen::MatrixXd V = problem_ptr->world_vertices(x + direction);
            // igl::writeOBJ(
            //     fmt::format("x{:04d}.obj", iteration_number), V,
            //     dynamic_cast<physics::SimulationProblem*>(problem_ptr)
            //         ->faces());
            double step_max_speed = (V - V_prev).lpNorm<Eigen::Infinity>()
                / problem_ptr->timestep();

            double tol = Constants::NEWTON_VELOCITY_CONVERGENCE_TOL
                * problem_ptr->world_bbox_diagonal();

            spdlog::info(
                "solve={} iter={:d} step_max_speed={:g} tol={:g}", //
                name(), iteration_number, step_max_speed, tol);
            return step_max_speed <= tol;
        }
        case ConvergenceCriteria::ENERGY: {
            spdlog::debug(
                "solve={} iter={:d} step_energy={:g} tol={:g}", //
                name(), iteration_number,
                abs(gradient_free.dot(direction_free)),
                Constants::NEWTON_ENERGY_CONVERGENCE_TOL);
            return abs(gradient_free.dot(direction_free))
                <= Constants::NEWTON_ENERGY_CONVERGENCE_TOL;
        }
        }
        throw NotImplementedError("Invalid convergence criteria option!");
    }

    OptimizationResults NewtonSolver::solve(const Eigen::VectorXd& x0)
    {
        assert(problem_ptr != nullptr);
        // Initialize the working variables
        x_prev = x0;
        x = x0;

        // igl::writeOBJ(
        //     fmt::format("x{:04d}.obj", 0), problem_ptr->world_vertices(x),
        //     dynamic_cast<physics::SimulationProblem*>(problem_ptr)->faces());

        double step_length = 1.0;

        spdlog::debug("solver={} action=BEGIN", name());

        std::string exit_reason = "exceeded the maximum allowable iterations";

        direction.setZero(problem_ptr->num_vars());
        grad_direction.setZero(problem_ptr->num_vars());

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
            igl::slice_into(gradient_free, free_dof, grad_direction);
            igl::slice_into(direction_free, free_dof, direction);

            // check for newton termination
            if (iteration_number > 0 && converged()) {
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
                // If I forced it to take a step when it should have converged
                if (iteration_number > 0) {
                    num_newton_ls_fails++;
                    spdlog::warn(
                        "solver={} iter={:d} failure=\"newton line-search\" "
                        "failsafe=\"gradient descent\"",
                        name(), iteration_number + 1);
                }

                ///////////////////////////////////////////////////////////////
                // Line search over -gradient direction
                // replace delta_x by gradient direction
                direction_free = -gradient_free;
                direction = -grad_direction;

                // check for newton termination again
                if (iteration_number > 0 && converged()) {
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
                    // If I forced it to take a step when it should have
                    // converged
                    if (iteration_number == 0 && converged()) {
                        // Do not consider this a failure
                        spdlog::error(
                            "solver={} failure=\"converged without taking a "
                            "step\"",
                            name());
                        exit_reason = "found a local optimum with -grad dir";
                        break;
                    }
                    num_grad_ls_fails++;
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

            x_prev = x;
            x += step_length * direction;
            assert(!problem_ptr->has_collisions(x_prev, x));

            post_step_update();
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
        double& step_length)
    {
        NAMED_PROFILE_POINT("line_search", LINE_SEARCH);
        PROFILE_START(LINE_SEARCH);

        bool success = false;
        int num_it = 0;
        double lower_bound = line_search_lower_bound() / -grad_fx.dot(dir);
        // double lower_bound = line_search_lower_bound();

        while (std::isfinite(lower_bound) && step_length > lower_bound) {
            num_it++;        // Count the number of iterations
            ls_iterations++; // Count the gloabal number of iterations

            // Compute the next variable
            Eigen::VectorXd xi = x + step_length * dir;

            // Check for collisions between newton updates
            num_collision_check++; // Count the number of collision checks
            if (!problem_ptr->has_collisions(x, xi)) {
                num_fx++; // Count the number of objective computations
                if (problem_ptr->compute_objective(xi) < fx) {
                    success = true;
                    break; // while loop
                }
            }

            // Try again with a smaller step_length
            step_length /= 2.0;
        }

        PROFILE_MESSAGE(
            LINE_SEARCH,
            fmt::format(
                "success,{},it,{},dir,{:10e}", success, num_it, dir.norm()))
        PROFILE_END(LINE_SEARCH);

        if (!success && iteration_number > 0) {
            if (!std::isfinite(lower_bound)) {
                spdlog::warn(
                    "solver={} iter={:d} failure=\"line-search ∇f(x)⋅dir=0\"",
                    name(), iteration_number + 1, lower_bound);
            } else {
                spdlog::warn(
                    "solver={} iter={:d} failure=\"line-search α ≤ {:g} / "
                    "{:g}\"",
                    name(), iteration_number + 1, line_search_lower_bound(),
                    -grad_fx.dot(dir));
            }
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

        // TODO: Can we use a better solver than Eigen's LDLT?
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(hessian);
        if (solver.info() == Eigen::Success) {
            direction = solver.solve(-gradient);
            if (solver.info() == Eigen::Success) {
                solve_success = true;
            } else {
                spdlog::warn(
                    "solver={} iter={:d} failure=\"sparse solve for newton "
                    "direction\" failsafe=\"gradient descent\"",
                    name(), iteration_number + 1);
            }
        } else {
            spdlog::warn(
                "solver={} iter={:d} failure=\"sparse decomposition of the "
                "hessian\" failsafe=\"gradient descent\"",
                name(), iteration_number + 1);
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
                name(), iteration_number + 1, mu);
            solve_success =
                compute_direction(gradient, psd_hessian, direction, false);
            double dir_dot_grad = direction.dot(gradient);
            if (dir_dot_grad > 0) {
                spdlog::warn(
                    "solver={} iter={:d} failure=\"adjusted newton "
                    "direction not descent direction\" failsafe=\"gradient "
                    "descent\" dir_dot_grad={:g}",
                    name(), iteration_number + 1, dir_dot_grad);
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
