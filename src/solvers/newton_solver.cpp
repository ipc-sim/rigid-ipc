// Functions for optimizing functions.
#include "newton_solver.hpp"

#include <igl/slice.h>
#include <igl/slice_into.h>

#include <constants.hpp>
#include <logger.hpp>
#include <profiler.hpp>

#if defined(DEBUG_LINESEARCH) || defined(DEBUG_COLLISIONS)
#include <io/serialize_json.hpp>
#include <nlohmann/json.hpp>

#endif

namespace ccd {
namespace opt {

    NewtonSolver::NewtonSolver()
        : NewtonSolver("newton_solver")
    {
    }
    NewtonSolver::NewtonSolver(const std::string& name)
        : max_iterations(1000)
        , iteration_number(0)
        , name_(name)

    {
    }

    NewtonSolver::~NewtonSolver() {}

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

    OptimizationResults NewtonSolver::solve(IBarrierProblem& problem)
    {
        static int global_it = 0;
        // Initalize the working variables
        Eigen::VectorXd x = problem.starting_point();
        Eigen::VectorXd gradient, gradient_free;
        Eigen::SparseMatrix<double> hessian, hessian_free;

        double step_length = 1.0;

        spdlog::debug("solver=newton action=BEGIN");

        double tolerance = std::max(e_b_, c_ * m_ / t_);
        std::string exit_reason = "exceeded the maximum allowable iterations";

        Eigen::VectorXd direction(problem.num_vars());
        Eigen::VectorXd grad_direction(problem.num_vars());
        Eigen::VectorXd direction_free(free_dof.size());

        bool success = false;

        for (iteration_number = 0; iteration_number < max_iterations;
             iteration_number++) {

            double fx;
            problem.eval_f_and_fdiff(x, fx, gradient, hessian);

            debug_newton_iterations += 1;
            debug_num_fx += 1;
            debug_num_grad_fx += 1;
            debug_num_hessian_fx += 1;

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
            if (abs(direction_free.dot(gradient_free)) <= tolerance) {
                exit_reason = "found a local optimum with newton dir";
                success = true;
                break;
            }

            // TODO: Test not resetting this to 1.0
            // step_length = std::min(2 * step_length, 1.0);
            step_length = 1.0;
            bool found_newton_step = line_search(
                problem, x, direction, fx, grad_direction, step_length,
                /*log_failure=*/true);
            ///////////////////////////////////////////////////////////////////

            // When newton direction fails, revert to gradient descent
            if (!found_newton_step) {
                spdlog::warn(
                    "solver=newton iter={:d} failure=\"newton line-search\" "
                    "failsafe=\"gradient descent\"",
                    iteration_number + 1);

                ///////////////////////////////////////////////////////////////
                // Line search over -gradient direction
                // replace delta_x by gradient direction
                direction.setZero();
                direction_free = -gradient_free;
                igl::slice_into(direction_free, free_dof, 1, direction);

                // check for newton termination again
                if (abs(direction_free.dot(gradient_free)) <= tolerance) {
                    exit_reason = "found a local optimum with -grad dir";
                    success = true;
                    break;
                }
                step_length = 1;
                bool found_gradient_step = line_search(
                    problem, x, direction, fx, grad_direction, step_length,
                    /*log_failure=*/true);
                ///////////////////////////////////////////////////////////////

                // When gradient direction fails, exit
                if (!found_gradient_step) {
                    spdlog::error(
                        "solver=newton iter={:d} failure=\"gradient "
                        "line-search\" failsafe=\"none\"",
                        iteration_number + 1);
                    exit_reason = "line-search failed";
                    break;
                }
            }

            spdlog::debug(
                "solver=newton iter={:d} step_length={:g}", iteration_number,
                step_length);

            Eigen::VectorXd xk = x + step_length * direction;
            assert(!problem.has_collisions(x, xk));
            x = xk;
        } // end for loop

        spdlog::debug(
            "solver=newton action=END total_iter={:d} "
            "exit_reason=\"{}\"",
            iteration_number, exit_reason);

        return OptimizationResults(x, problem.eval_f(x), success, true);
    }

    void NewtonSolver::init_free_dof(Eigen::VectorXb is_dof_fixed)
    {
        free_dof = Eigen::VectorXi(is_dof_fixed.size() - is_dof_fixed.count());
        for (int i = 0, j = 0; i < is_dof_fixed.size(); i++) {
            if (!is_dof_fixed(i)) {
                free_dof(j++) = i;
            }
        }
    }

    std::string NewtonSolver::debug_stats()
    {
        return fmt::format(
            "total_newton_steps,{},total_ls_steps,{},count_fx,{}"
            ",count_grad,{},count_hess,{},count_ccd,{}",
            debug_newton_iterations, debug_ls_iterations, debug_num_fx,
            debug_num_grad_fx, debug_num_hessian_fx, debug_num_collision_check);
    }

    void NewtonSolver::debug_reset_stats()
    {
        debug_num_fx = 0;
        debug_num_grad_fx = 0;
        debug_num_hessian_fx = 0;
        debug_num_collision_check = 0;
        debug_ls_iterations = 0;
        debug_newton_iterations = 0;
    }

    bool NewtonSolver::line_search(
        IBarrierProblem& problem,
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const double fx,
        const Eigen::VectorXd& grad_fx,
        double& alpha,
        bool log_failure)
    {

        static int global_it = 0;
        PROFILE_POINT("line_search");
        PROFILE_START();

        bool success = false;

        std::stringstream debug;
        debug << "global_it,it,step_length,collisions,f(x),fx0,f'*step,>lb, "
                 "barrier_eps\n";

        int num_it = 0;
        global_it += 1;

        double lower_bound = std::min(1E-12, c_ * e_b_ / 10.0);
        const double eps = problem.get_barrier_epsilon();

        assert(isfinite(-grad_fx.dot(dir)));

        while (-grad_fx.dot(dir) * alpha > lower_bound) {
            debug_ls_iterations += 1;
            Eigen::VectorXd xi = x + alpha * dir;

            bool no_collisions = !problem.has_collisions(x, xi);
            double fxi = problem.eval_f(xi);

            debug_num_collision_check += 1;
            debug_num_fx += 1;

            debug << fmt::format(
                "{},{},{:.18e},{},{:.18e},{:.18e},{:.18e},{:.18e}, {:.8e}\n",
                global_it, num_it, alpha,
                (no_collisions ? "no-collisions" : "has-collision"), fxi, fx,
                -grad_fx.dot(dir) * alpha, lower_bound, eps);

            num_it += 1;
            if (fxi < fx && no_collisions) {
                success = true;
                break; // while loop
            }

            alpha /= 2.0;
        }

        if (log_failure && !success) {
            std::cout
                << fmt::format(
                       "linesearch fail num_it={} `-grad_fx.dot(dir)`={:.18e} "
                       "`-grad_fx.dot(dir) * alpha `={:.18e} "
                       "lower_bound={:.8e} alpha={:.8e} fx={:.8e} t={:.8e}",
                       num_it, -grad_fx.dot(dir), -grad_fx.dot(dir) * alpha,
                       lower_bound, alpha, fx, t_)
                << std::endl;
            std::string fout = fmt::format(
                "{}/linesearch_{}.csv", DATA_OUTPUT_DIR, ccd::logger::now());
            std::ofstream myfile;
            myfile.open(fout);
            myfile << debug.str();
            myfile.close();
            spdlog::debug("saved failure to `{}`", fout);
        }

        // clang-format off
        PROFILE_MESSAGE(,
            fmt::format(
                "success,{},it,{},dir,{:10e}", success, num_it, dir.norm()))
        // clang-format on
        PROFILE_END();
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
                    "solver=newton iter={:d} failure=\"sparse solve for newton "
                    "direction\" failsafe=\"gradient descent\"",
                    iteration_number);
            }
        } else {
            spdlog::warn(
                "solver=newton iter={:d} failure=\"sparse decomposition of the "
                "hessian\" failsafe=\"gradient descent\"",
                iteration_number);
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
                "solver=newton iter={:d} failure=\"newton direction "
                "not descent direction\" failsafe=\"H += μI\" μ={:g}",
                iteration_number, mu);
            solve_success =
                compute_direction(gradient, psd_hessian, direction, false);
            double dir_dot_grad = direction.dot(gradient);
            if (dir_dot_grad > 0) {
                spdlog::warn(
                    "solver=newton iter={:d} failure=\"adjusted newton "
                    "direction not descent direction\" failsafe=\"gradient "
                    "descent\" dir_dot_grad={:g}",
                    iteration_number, dir_dot_grad);
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
