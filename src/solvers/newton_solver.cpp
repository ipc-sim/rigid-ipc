// Functions for optimizing functions.
#include "newton_solver.hpp"

#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/writeOBJ.h>

#include <constants.hpp>
#include <logger.hpp>
#include <profiler.hpp>
#include <solvers/line_search.hpp> // sample_search_direction

// #define USE_GRADIENT_DESCENT

namespace ccd {
namespace opt {

    NewtonSolver::NewtonSolver()
        : max_iterations(1000)
        , iteration_number(0)
        , convergence_criteria(ConvergenceCriteria::ENERGY)
        , energy_conv_tol(Constants::DEFAULT_NEWTON_ENERGY_CONVERGENCE_TOL)
        , velocity_conv_tol(Constants::DEFAULT_NEWTON_VELOCITY_CONVERGENCE_TOL)
        , m_line_search_lower_bound(Constants::DEFAULT_LINE_SEARCH_LOWER_BOUND)
    {
        linear_solver = polysolve::LinearSolver::create("", "");
    }

    void NewtonSolver::settings(const nlohmann::json& json)
    {
        max_iterations = json["max_iterations"];
        convergence_criteria = json["convergence_criteria"];
        energy_conv_tol = json["energy_conv_tol"];
        velocity_conv_tol = json["velocity_conv_tol"];
        is_velocity_conv_tol_abs = json["is_velocity_conv_tol_abs"];
        m_line_search_lower_bound = json["line_search_lower_bound"];

        linear_solver_settings = json["linear_solver"];
        try {
            linear_solver = polysolve::LinearSolver::create(
                linear_solver_settings["name"], "");
        } catch (const std::runtime_error& err) {
            spdlog::error(
                "{}! Using Eigen::SimplicialLDLT instead.", err.what());
            linear_solver_settings["name"] = "Eigen::SimplicialLDLT";
            linear_solver = polysolve::LinearSolver::create(
                linear_solver_settings["name"], "");
        }
        linear_solver->setParameters(linear_solver_settings);

        reset_stats();
    }

    nlohmann::json NewtonSolver::settings() const
    {
        nlohmann::json settings;
        settings["max_iterations"] = max_iterations;
        settings["convergence_criteria"] = convergence_criteria;
        settings["linear_solver"] = linear_solver_settings;
        settings["energy_conv_tol"] = energy_conv_tol;
        settings["velocity_conv_tol"] = velocity_conv_tol;
        settings["is_velocity_conv_tol_abs"] = is_velocity_conv_tol_abs;
        return settings;
    }

    void NewtonSolver::init_solve(const Eigen::VectorXd& x0)
    {
        assert(problem_ptr != nullptr);
    }

    nlohmann::json NewtonSolver::stats() const
    {
        return { { "total_newton_steps", newton_iterations },
                 { "total_ls_steps", ls_iterations },
                 { "num_newton_ls_fails", num_newton_ls_fails },
                 { "num_grad_ls_fails", num_grad_ls_fails },
                 { "count_fx", num_fx },
                 { "count_grad", num_grad_fx },
                 { "count_hess", num_hessian_fx },
                 { "count_ccd", num_collision_check },
                 { "total_regularizations", regularization_iterations } };
    }

    std::string NewtonSolver::stats_string() const
    {
        return fmt::format(
            "total_newton_steps={:d} total_ls_steps={:d} "
            "num_newton_ls_fails={:d} num_grad_ls_fails={:d} count_fx={:d} "
            "count_grad={:d} count_hess={:d} count_ccd={:d} "
            "total_regularizations={:d}",
            newton_iterations, ls_iterations, num_newton_ls_fails,
            num_grad_ls_fails, num_fx, num_grad_fx, num_hessian_fx,
            num_collision_check, regularization_iterations);
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
        regularization_iterations = 0;
    }

    bool NewtonSolver::energy_converged() const
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

            // TODO: Renable this with a better check for static objects
            double tol = velocity_conv_tol;
            if (!is_velocity_conv_tol_abs) {
                tol *= problem_ptr->world_bbox_diagonal();
            }

            spdlog::info(
                "solver={} iter={:d} step_max_speed={:g} tol={:g}", //
                name(), iteration_number, step_max_speed, tol);

            double newton_step_energy =
                sqrt(abs(gradient_free.dot(direction_free)));
            double gradient_step_energy =
                sqrt(abs(gradient_free.dot(gradient_free)));
            double gradient_mass_step_energy = sqrt(
                abs(grad_direction.transpose() * problem_ptr->mass_matrix()
                    * grad_direction));

            // GREP_NCONV,1,0.130871,1.02885e-05,1.15864e-06,2.46467e-06,0.00434816,0.00399746,9.82746e-06
            /*
            std::cout << fmt::format(
                             "{},{:d},{:g},{:g},{:g},{:g},{:g},{:g},{:g}",
                             step_max_speed <= tol ? "GREP_CONV" : "GREP_NCONV",
                             iteration_number, step_max_speed,
                             newton_step_energy, gradient_step_energy,
                             gradient_mass_step_energy, direction_free.norm(),
                             direction_free.lpNorm<Eigen::Infinity>(),
                             sqrt(
                                 abs(direction.transpose()
                                     * problem_ptr->mass_matrix() * direction)))
                      << std::endl;
            */
            return step_max_speed <= tol;
        }
        case ConvergenceCriteria::ENERGY: {
            double step_energy = abs(gradient_free.dot(direction_free));
            // double step_energy = abs(gradient_free.dot(gradient_free));
            double tol = energy_conv_tol;
            spdlog::info(
                "solver={} iter={:d} step_energy={:g} tol={:g}", //
                name(), iteration_number, step_energy, tol);
            return step_energy <= tol;
        }
        }
        throw NotImplementedError("Invalid convergence criteria option!");
    }

    bool NewtonSolver::converged() const
    {
        return energy_converged()
            && problem_ptr->are_equality_constraints_satisfied(x);
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
        double regulariztion_coeff = 0;

        spdlog::debug("solver={} action=BEGIN", name());

        std::string exit_reason = "exceeded the maximum allowable iterations";

        direction.setZero(problem_ptr->num_vars());
        grad_direction.setZero(problem_ptr->num_vars());

        bool success = false;

        for (iteration_number = 0; iteration_number < max_iterations;
             iteration_number++) {
            double fx = problem_ptr->compute_objective(x, gradient, hessian);

            num_fx++;
            num_grad_fx++;
            num_hessian_fx++;

            // Remove rows and cols of fixed DoF
            Eigen::VectorXi free_dof = problem_ptr->free_dof();
            igl::slice(gradient, free_dof, gradient_free);
            igl::slice(hessian, free_dof, free_dof, hessian_free);

#ifdef USE_GRADIENT_DESCENT
            direction_free = -gradient_free;
#else
            bool solve_success = compute_regularized_direction(
                fx, gradient_free, hessian_free, direction_free,
                regulariztion_coeff);
            if (!solve_success) {
                exit_reason = "regularization failed";
                break;
            }
#endif

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

            step_length = 1;
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
                        name(), iteration_number);
                }

                ///////////////////////////////////////////////////////////////
                // Line search over -gradient direction
                // replace delta_x by gradient direction
                direction_free = -gradient_free;
                direction = -grad_direction;

                // WARNING: Disable convergence using negative gradient
                // check for newton termination again
                // if (iteration_number > 0 && converged()) {
                //     exit_reason = "found a local optimum with -grad dir";
                //     success = true;
                //     break;
                // }
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
                        spdlog::warn(
                            "solver={} failure=\"converged without taking a "
                            "step\"",
                            name());
                        exit_reason = "found a local optimum with -grad dir";
                        success = true;
                        break;
                    }
                    num_grad_ls_fails++;
                    spdlog::error(
                        "solver={} iter={:d} failure=\"gradient line-search\" "
                        "failsafe=\"none\"",
                        name(), iteration_number);
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
            newton_iterations++; // Only count complete steps

            post_step_update();
        } // end for loop

        spdlog::info(
            "solver={} action=END total_iter={:d} exit_reason=\"{}\"", name(),
            iteration_number, exit_reason);

        return OptimizationResults(
            x, problem_ptr->compute_objective(x), success, true,
            iteration_number);
    }

    bool NewtonSolver::line_search(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const double fx,
        const Eigen::VectorXd& grad_fx,
        double& step_length)
    {
        NAMED_PROFILE_POINT("NewtonSolver::line_search", LINE_SEARCH);
        PROFILE_START(LINE_SEARCH);

        bool success = false;
        int num_it = 0;
        // double lower_bound = line_search_lower_bound() / -grad_fx.dot(dir);
        double lower_bound =
            line_search_lower_bound() / sqrt(abs(grad_fx.dot(dir)));
        lower_bound = std::min(lower_bound, 1e-1);
        // double lower_bound = line_search_lower_bound() / dir.squaredNorm();
        // double lower_bound = line_search_lower_bound();

        // Filter the step length so that x to x + α * Δx is collision free for
        // α ≤ step_length.
        num_collision_check++; // Count the number of collision checks
        double max_step_size =
            std::min(problem_ptr->compute_earliest_toi(x, x + dir), 1.0);
        step_length = std::min(step_length, max_step_size);
        // while (problem_ptr->has_collisions(x, x + step_length * dir)) {
        //     spdlog::error(
        //         "max_step_size={:g} has collisions halving the step!",
        //         step_length);
        //     step_length *= 0.8;
        // }
        if (step_length < lower_bound) {
            spdlog::error(
                "solver={} iter={:d} failure=\"initial step_length (α={:g}) is"
                " less than lower bound ({:g})\"",
                name(), iteration_number, step_length, lower_bound);
        }

        double fxi = std::numeric_limits<double>::infinity();
        while (std::isfinite(lower_bound) && step_length >= lower_bound) {
            num_it++;        // Count the number of iterations
            ls_iterations++; // Count the gloabal number of iterations

            // Compute the next variable
            Eigen::VectorXd xi = x + step_length * dir;

            // NOTE: We do not need to check for collisions because we filtered
            // the step length.
            // Check for collisions between newton updates
            // if (!problem_ptr->has_collisions(x, xi)) {
            fxi = problem_ptr->compute_objective(xi);
            num_fx++; // Count the number of objective computations
            if (fxi < fx) {
                success = true;
                break; // while loop
            }
            // }

            // Try again with a smaller step_length
            step_length /= 2.0;
        }

        PROFILE_MESSAGE(
            LINE_SEARCH, "success,it,dir",
            fmt::format("{},{:d},{:10e}", success, num_it, dir.norm()));
        PROFILE_END(LINE_SEARCH);

        if (!success && iteration_number > 0 && max_step_size >= lower_bound) {
            if (!std::isfinite(lower_bound)) {
                spdlog::warn(
                    "solver={} iter={:d} failure=\"line-search ∇f(x)⋅dir=0\"",
                    name(), iteration_number);
            } else {
                spdlog::warn(
                    "solver={} iter={:d} failure=\"line-search α ≤ {:g} / {:g} "
                    "= {:g}; f(x + αΔx)-f(x)={:g}; α_max={:g}\"",
                    name(), iteration_number, line_search_lower_bound(),
                    sqrt(abs(grad_fx.dot(dir))), lower_bound, fxi - fx,
                    max_step_size);
                sample_search_direction(
                    x, dir,
                    [&](const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
                        double fx = problem_ptr->compute_objective(x, grad);
                        Eigen::VectorXi free_dof = problem_ptr->free_dof();
                        Eigen::VectorXd grad_free;
                        igl::slice(grad, free_dof, grad_free);
                        grad.setZero();
                        igl::slice_into(grad_free, free_dof, grad);
                        return fx;
                    },
                    max_step_size);
            }
        }

        spdlog::debug(
            "solver={} iter={:d} line_search_success={} step_length={:g}",
            name(), iteration_number, success, step_length);

        return success;
    }

    bool NewtonSolver::compute_regularized_direction(
        double& fx,
        Eigen::VectorXd& gradient,
        Eigen::SparseMatrix<double>& hessian,
        Eigen::VectorXd& direction,
        double& coeff)
    {
        bool success = false;
        while (!success) {
            double regularized_fx = fx;
            Eigen::VectorXd regularized_gradient = gradient;
            Eigen::SparseMatrix<double> regularized_hessian = hessian;

            if (coeff > 0) {
                // regularized_fx += coeff / 2 * (x - x_prev).squaredNorm();
                // regularized_gradient += coeff * (x - x_prev);
                Eigen::SparseMatrix<double> I(hessian.rows(), hessian.cols());
                I.setIdentity();
                regularized_hessian += coeff * I;
                regularization_iterations++;
            }

            success = compute_direction(
                regularized_gradient, regularized_hessian, direction,
                /*make_psd=*/false);
            success = success && regularized_gradient.dot(direction) <= 0;

            // Update coefficient adaptivly when the solve fails
            if (success) {
                coeff /= 2;
                if (coeff < 1e-8) {
                    coeff = 0;
                }
                fx = regularized_fx;
                gradient = regularized_gradient;
                hessian = regularized_hessian;
            } else {
                coeff = std::max(2 * coeff, 1e-8);
                if (!std::isfinite(coeff)) {
                    spdlog::error(
                        "solver={} iter={:d} failure=\"regularization failed "
                        "(coeff={:g})\" failsafe=\"none\"",
                        name(), iteration_number, coeff);
                    return false;
                }
                spdlog::warn(
                    "solver={} iter={:d} failure=\"solve failed (∇f⋅Δx={:g}, "
                    "||H||_∞={:g}); increasing regularization coeff={:g}\"",
                    name(), iteration_number, gradient_free.dot(direction_free),
                    Eigen::MatrixXd(hessian).lpNorm<Eigen::Infinity>(), coeff);
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
        PROFILE_POINT("NewtonSolver::compute_direction:linear_solve");
        PROFILE_START();

        // Check if the hessian is positive semi-definite.
        // Eigen::LLT<Eigen::MatrixXd> LLT_H((Eigen::MatrixXd(hessian)));
        // if (LLT_H.info() == Eigen::NumericalIssue) {
        //     spdlog::warn(
        //         "solver={} iter={:d} failure=\"possibly non semi-positive "
        //         "definite Hessian\"",
        //         name(), iteration_number);
        // }

        // Solve for the Newton direction (Δx = -H⁻¹∇f).
        // Return true if the solve was successful.
        bool solve_success = false;

        // if (hessian.rows() <= 1200) { // <= 200 bodies
        //     Eigen::MatrixXd dense_hessian(hessian);
        //     direction = dense_hessian.ldlt().solve(-gradient);
        //     solve_success = true;
        // } else {
        linear_solver->analyzePattern(hessian, hessian.rows());
        linear_solver->factorize(hessian);
        nlohmann::json info;
        linear_solver->getInfo(info);
        // TODO: This check only works for direct Eigen solvers
        if (!info.contains("solver_info") || info["solver_info"] == "Success") {
            // TODO: Do we have a better initial guess for iterative
            // solvers?
            direction = Eigen::VectorXd::Zero(gradient.size());
            linear_solver->solve(-gradient, direction);
            linear_solver->getInfo(info);
            if (!info.contains("solver_info")
                || info["solver_info"] == "Success") {
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
        // }

        PROFILE_END();

        // Check solve residual
        if (solve_success) {
            double solve_residual = (hessian * direction + gradient).norm();
            if (solve_residual > 1e-8) {
                spdlog::warn(
                    "solver={} iter={:d} "
                    "failure=\"linear solve residual ({:g}) > 1e-10\" "
                    "failsafe=\"gradient descent\"",
                    name(), iteration_number, solve_residual);
                solve_success = false;
            }
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
            spdlog::warn(
                "solver={} iter={:d} failure=\"newton direction not descent "
                "direction\" failsafe=\"H += μI\" μ={:g}",
                name(), iteration_number, mu);
            solve_success =
                compute_direction(gradient, psd_hessian, direction, false);
            double dir_dot_grad = direction.dot(gradient);
            if (dir_dot_grad > 0) {
                spdlog::error(
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

    void NewtonSolver::post_step_update()
    {
        if (energy_converged()
            && !problem_ptr->are_equality_constraints_satisfied(x)) {
            spdlog::info(
                "solver={} iter={:d} msg=\"updated augmented Lagrangian\"",
                name(), iteration_number);
            problem_ptr->update_augmented_lagrangian(x);
        }
    }

} // namespace opt
} // namespace ccd
