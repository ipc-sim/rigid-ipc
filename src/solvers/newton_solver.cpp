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
        : NewtonSolver("newton_solver")
    {
    }
    NewtonSolver::NewtonSolver(const std::string& name)
        : absolute_tolerance(Constants::NEWTON_ABSOLUTE_TOLERANCE)
        , min_step_length(1e-12)
        , max_iterations(1000)
        , iteration_number(0)
        , name_(name)

    {
    }

    NewtonSolver::~NewtonSolver() {}

    void NewtonSolver::settings(const nlohmann::json& json)
    {
        max_iterations = json["max_iterations"].get<int>();
        absolute_tolerance = json["absolute_tolerance"].get<double>();
        min_step_length = json["min_step_length"].get<double>();
    }

    nlohmann::json NewtonSolver::settings() const
    {
        nlohmann::json json;
        json["max_iterations"] = max_iterations;
        json["absolute_tolerance"] = absolute_tolerance;
        json["min_step_length"] = min_step_length;
        return json;
    }

    const static char* NEWTON_DIRECTION_LOG
        = "solver=newton_solver iter={:d} action=newton_line_search "
          "status=failure message=\"reverting to gradient descent\"";
    const static char* NEWTON_GRADIENT_LOG
        = "solver=newton_solver iter={:d} action=gradient_line_search "
          "status=failure message=\"no failsafe - break!\"";
    const static char* NEWTON_END_LOG = "solver=newton_solver action=END";
    const static char* NEWTON_BEGIND_LOG = "solver=newton_solver action=BEGIN";

    OptimizationResults NewtonSolver::solve(IBarrierProblem& problem)
    {
        NAMED_PROFILE_POINT("newton_solver", SOLVER_STEP);
        NAMED_PROFILE_POINT("newton_solver__eval_f_and_fdiff", EVAL_F);
        NAMED_PROFILE_POINT(
            "newton_solver__compute_newton_direction", COMPUTE_DIRECTION);
        NAMED_PROFILE_POINT(
            "newton_solver__newton_line_search", NEWTON_LINE_SEARCH);
        NAMED_PROFILE_POINT(
            "newton_solver__gradient_line_search", GRADIENT_LINE_SEARCH);
        NAMED_PROFILE_POINT(
            "newton_solver__assert_constraint", ASSERT_CONSTRAINT);

        // Initalize the working variables
        Eigen::VectorXd x = problem.starting_point();
        Eigen::VectorXd gradient, gradient_free;
        Eigen::SparseMatrix<double> hessian, hessian_free;

        double step_length = 1.0;

        spdlog::trace(NEWTON_BEGIND_LOG);
        iteration_number = 0;
        std::string exit_reason = "exceeded the maximum allowable iterations";
        do {
            PROFILE_START(SOLVER_STEP);

            PROFILE_START(EVAL_F) // Compute the gradient and hessian
            double fx;
            problem.eval_f_and_fdiff(x, fx, gradient, hessian);
            PROFILE_END(EVAL_F)

            // Remove rows of fixed dof
            igl::slice(gradient, free_dof, gradient_free);

            if (gradient_free.squaredNorm() <= absolute_tolerance) {
                exit_reason = "found a local optimum";
                PROFILE_END(SOLVER_STEP);
                break;
            }

            // Remove rows and columns of fixed dof of the hessian
            igl::slice(hessian, free_dof, free_dof, hessian_free);

            PROFILE_START(COMPUTE_DIRECTION)

            Eigen::VectorXd direction
                = Eigen::VectorXd::Zero(problem.num_vars());
            Eigen::VectorXd direction_free(free_dof.size());
            bool found_direction = compute_direction(
                gradient_free, hessian_free, direction_free, /*make_psd=*/true);
            igl::slice_into(direction_free, free_dof, direction);

            PROFILE_END(COMPUTE_DIRECTION)

            if (!found_direction) {
                exit_reason = "newton direction solve failed";
                PROFILE_END(SOLVER_STEP);
                break;
            }

            step_length = 1; // std::min(1.0, 2.0 * step_length);
            //            double step_length_aux
            //                = step_length; // to restore when line-search
            //                fails

            PROFILE_START(NEWTON_LINE_SEARCH)
            bool found_step_length
                = line_search(problem, x, direction, fx, step_length);
            PROFILE_END(NEWTON_LINE_SEARCH)

            // Revert to gradient descent if the newton direction fails
            if (!found_step_length) {
                step_length = 1; // step_length_aux;

                spdlog::warn(NEWTON_DIRECTION_LOG, iteration_number + 1);

                // replace delta_x by gradient direction
                igl::slice_into(-gradient_free, free_dof, 1, direction);

                PROFILE_START(GRADIENT_LINE_SEARCH)
                bool found_step_length2
                    = line_search(problem, x, direction, fx, step_length, true);
                PROFILE_END(GRADIENT_LINE_SEARCH)

                if (!found_step_length2) {
                    spdlog::error(NEWTON_GRADIENT_LOG, iteration_number + 1);
                    exit_reason = "line-search failed";
                    PROFILE_END(SOLVER_STEP);

                    // Remove this line one we are confident this doesn't happen
                    // Keep this line commented until we know the affect of
                    // failing
                    // throw std::logic_error("gradient line-search failed");
                    break;
                }
            }

            x += step_length * direction;
            PROFILE_START(ASSERT_CONSTRAINT)
            assert(!std::isinf(problem.eval_f(x)));
            PROFILE_END(ASSERT_CONSTRAINT)

            PROFILE_END(SOLVER_STEP);

        } while (++iteration_number <= max_iterations);

        spdlog::trace(
            "{} total_iter={:d} gradient_free_eq_norm={:e} message='{}'",
            NEWTON_END_LOG, iteration_number, gradient_free.squaredNorm(),
            exit_reason);

        return OptimizationResults(x, problem.eval_f(x),
            gradient_free.squaredNorm() <= absolute_tolerance);
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

    bool NewtonSolver::line_search(IBarrierProblem& problem,
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const double fx,
        double& step_length,
        bool log_failure)
    {

        static int global_it = 0;
        PROFILE_POINT("line_search");
        PROFILE_START();

        double step_norm = (step_length * dir).norm();
        bool success = false;

        std::stringstream debug;
        debug << "global_it,it,step_length,step_norm,f(x),fx0\n";

        int num_it = 0;

        global_it += 1;
        // Multiprecision fx = problem.eval_mp_f(x);
        while (step_norm > Constants::LINESEARCH_MIN_STEP_NORM) {

            Eigen::VectorXd xi = x + step_length * dir;
            // Multiprecision fxi = problem.eval_mp_f(xi);
            double fxi = problem.eval_f(xi);

            bool min_rule = fxi < fx;
            bool cstr = !std::isinf(fxi);
            // bool cstr = !std::isinf(fxi.to_double());

            // debug <<
            // fmt::format("{},{},{:.18e},{:.18e},{:.18e},{:.18e}\n",
            // global_it, num_it,
            //     step_length, step_norm, fxi.to_double(),
            //     fx.to_double());
            debug << fmt::format("{},{},{:.18e},{:.18e},{:.18e},{:.18e}\n",
                global_it, num_it, step_length, step_norm, fxi, fx);

            num_it += 1;
            if (min_rule && cstr) {
                success = true;
                assert(!std::isinf(problem.eval_f(xi)));
                break; // while loop
            }

            step_length /= 2.0;
            step_norm = (step_length * dir).norm();
        }

        if (log_failure && !success) {
            std::string fout = fmt::format(
                "{}/linesearch_{}.csv", DATA_OUTPUT_DIR, ccd::log::now());
            std::ofstream myfile;
            myfile.open(fout);
            myfile << debug.str();
            myfile.close();
            spdlog::debug("saved failure to `{}`", fout);
        }

        PROFILE_MESSAGE(,
            fmt::format(
                "success,{},it,{},dir,{:10e}", success, num_it, dir.norm()))
        PROFILE_END();
        return success;
    }

    bool NewtonSolver::compute_direction(const Eigen::VectorXd& gradient,
        const Eigen::SparseMatrix<double>& hessian,
        Eigen::VectorXd& direction,
        bool make_psd)
    {
        // Solve for the Newton direction (Δx = -H^{-1}∇f).
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
                    "solver=newton iter={:d} failure=\"sparse solve for newton direction failed\" failsafe=none",
                    iteration_number);
            }
        } else {
            spdlog::warn(
                "solver=newton iter={:d} failure=\"sparse decomposition of the hessian failed\" failsafe=none",
                iteration_number);
        }

        if (solve_success && make_psd
            && direction.transpose() * gradient >= 0) {
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
                = compute_direction(gradient, psd_hessian, direction, false);
            if (direction.transpose() * gradient >= 0) {
                spdlog::warn(
                    "solver=newton iter={:d} failure=\"newton direction not descent direction\" failsafe=none dir_dot_grad={:g}",
                    iteration_number, direction.transpose() * gradient);
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
