// Solve the optimization general_problem using Newton's Method with barriers
// for the constraints.

#include "barrier_solver.hpp"

#include <opt/barrier.hpp>
#include <utils/eigen_ext.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {
namespace opt {

    BarrierSolver::BarrierSolver()
        : min_barrier_epsilon(1e-5)
        , inner_solver_type(BarrierInnerSolver::NEWTON)
        , num_outer_iterations_(0)
    {
        newton_inner_solver = NewtonSolver();
        bfgs_inner_solver = BFGSSolver();
        gradient_descent_inner_solver = GradientDescentSolver();
    }

    OptimizationSolver& BarrierSolver::get_inner_solver()
    {
        switch (inner_solver_type) {
        case BarrierInnerSolver::NEWTON:
            return newton_inner_solver;
        case BarrierInnerSolver::BFGS:
            return bfgs_inner_solver;
        default:
            return gradient_descent_inner_solver;
        }
    }

    void BarrierSolver::clear()
    {
        num_outer_iterations_ = 0;
        general_problem_ptr = nullptr;
        barrier_problem_ptr.reset();
    }
    void BarrierSolver::init(OptimizationProblem& original_problem)
    {
        assert(original_problem.has_barrier_constraint());
        general_problem_ptr = &original_problem;
        barrier_problem_ptr
            = std::make_unique<BarrierProblem>(original_problem);

        OptimizationSolver& inner_solver = get_inner_solver();

        // Convert from the boolean vector to a vector of free dof indices
        inner_solver.init_free_dof(barrier_problem_ptr->is_dof_fixed());

        // Calculate the maximum number of iteration allowable
        inner_solver.max_iterations = int(max_iterations
            / ceil(-log2(min_barrier_epsilon) + log2(barrier_epsilon())));

        barrier_problem_ptr->eval_intermediate_callback(
            barrier_problem_ptr->x0);
        num_outer_iterations_ = 0;
    }

    OptimizationResults BarrierSolver::step_solve()
    {
        assert(general_problem_ptr != nullptr);
        assert(barrier_problem_ptr != nullptr);

        OptimizationResults results;
        OptimizationSolver& inner_solver = get_inner_solver();

        // Log the epsilon and the newton method will log the number of
        // iterations.
        spdlog::trace("solver=barrier ϵ={:g}",
            general_problem_ptr->get_barrier_epsilon());

        // Optimize for a fixed epsilon
        results = inner_solver.solve(*barrier_problem_ptr);
        // Save the original problems objective
        results.minf = general_problem_ptr->eval_f(results.x);

        // Steepen the barrier
        double eps = barrier_epsilon();
        general_problem_ptr->set_barrier_epsilon(eps / 2);

        // Start next iteration from the ending optimal position
        barrier_problem_ptr->x0 = results.x;

        // TODO: This should check if the barrier constraints are satisfied.
        results.success = results.minf >= 0
            && barrier_problem_ptr->eval_f(results.x)
                < std::numeric_limits<double>::infinity();

        results.finished = barrier_epsilon() <= min_barrier_epsilon;

        num_outer_iterations_ += 1;
        return results;
    }

    OptimizationResults BarrierSolver::solve(
        OptimizationProblem& original_problem)
    {
        init(original_problem);

        OptimizationResults results;
        do {
            results = step_solve();
        } while (barrier_epsilon() > min_barrier_epsilon);

        return results;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// BARRIER PROBLEM
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BarrierProblem::BarrierProblem(OptimizationProblem& general_problem)
        : OptimizationProblem("BarrierProblem")
        , general_problem(&general_problem)
    {
        num_vars = general_problem.num_vars;
        num_constraints = 0;
        x0 = general_problem.x0;
    }

    const Eigen::VectorXb& BarrierProblem::is_dof_fixed()
    {
        return general_problem->is_dof_fixed();
    }

    double BarrierProblem::eval_f(const Eigen::VectorXd& x)
    {
        double f_uk = general_problem->eval_f(x);
        double gx = general_problem->eval_g(x).sum();
        return f_uk + gx;
    }

    Eigen::VectorXd BarrierProblem::eval_grad_f(const Eigen::VectorXd& x)
    {
        Eigen::VectorXd f_uk_gradient = general_problem->eval_grad_f(x);
        Eigen::MatrixXd dgx = general_problem->eval_jac_g(x);

        f_uk_gradient += dgx.colwise().sum().transpose();

        return f_uk_gradient;
    }

    Eigen::SparseMatrix<double> BarrierProblem::eval_hessian_f(
        const Eigen::VectorXd& x)
    {
        Eigen::SparseMatrix<double> f_uk_hessian
            = general_problem->eval_hessian_f(x);
        std::vector<Eigen::SparseMatrix<double>> ddgx
            = general_problem->eval_hessian_g(x);

        PROFILE(
            // clang-format off
            for (const auto& ddgx_i : ddgx) {
                f_uk_hessian += ddgx_i;
            },
            // clang-format on
            ProfiledPoint::SUMMING_HESSIAN)

        return f_uk_hessian;
    }

    void BarrierProblem::eval_f_and_fdiff(const Eigen::VectorXd& x,
        double& f_uk,
        Eigen::VectorXd& f_uk_gradient,
        Eigen::SparseMatrix<double>& f_uk_hessian)
    {
        general_problem->eval_f_and_fdiff(x, f_uk, f_uk_gradient, f_uk_hessian);

        Eigen::VectorXd gx;
        Eigen::MatrixXd dgx;
        std::vector<Eigen::SparseMatrix<double>> ddgx;
        general_problem->eval_g_and_gdiff(x, gx, dgx, ddgx);

        f_uk += gx.sum();
        f_uk_gradient += dgx.colwise().sum().transpose();

        PROFILE(
            // clang-format off
            for (const auto& ddgx_i : ddgx) {
                f_uk_hessian += ddgx_i;
            },
            // clang-format on
            ProfiledPoint::SUMMING_HESSIAN)
    }

    void BarrierProblem::enable_line_search_mode(const Eigen::VectorXd& max_x)
    {
        general_problem->enable_line_search_mode(max_x);
    }

    void BarrierProblem::disable_line_search_mode()
    {
        general_problem->disable_line_search_mode();
    }

    bool BarrierProblem::eval_intermediate_callback(const Eigen::VectorXd& x)
    {
        return general_problem->eval_intermediate_callback(x);
    }

} // namespace opt
} // namespace ccd
