// Solve the optimization general_problem using Newton's Method with barriers
// for the constraints.

#include "barrier_solver.hpp"

#include <opt/barrier.hpp>
#include <utils/eigen_ext.hpp>

#include <solvers/solver_factory.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {
namespace opt {

    BarrierSolver::BarrierSolver()
        : BarrierSolver("barrier_solver")
    {
    }

    BarrierSolver::BarrierSolver(const std::string& name)
        : min_barrier_epsilon(1e-5)
        , max_iterations(1000)
        , num_outer_iterations_(0)
        , name_(name)

    {
    }
    void BarrierSolver::set_problem(IBarrierGeneralProblem& problem){
        general_problem_ptr = &problem;
    }

    void BarrierSolver::settings(const nlohmann::json& json)
    {
        max_iterations = json["max_iterations"].get<int>();
        min_barrier_epsilon = json["min_barrier_epsilon"].get<double>();
        inner_solver_ptr = SolverFactory::factory().get_barrier_inner_solver(
            json["inner_solver"].get<std::string>());
    }

    nlohmann::json BarrierSolver::settings() const
    {
        nlohmann::json json;
        json["max_iterations"] = max_iterations;
        json["min_barrier_epsilon"] = min_barrier_epsilon;
        json["inner_solver"] = get_inner_solver().name();
        return json;
    }

    void BarrierSolver::init_solve()
    {
        num_outer_iterations_ = 0;
        assert(general_problem_ptr != nullptr);
        barrier_problem_ptr.reset();

        barrier_problem_ptr
            = std::make_unique<BarrierProblem>(*general_problem_ptr);

        IBarrierOptimizationSolver& inner_solver = get_inner_solver();

        // Convert from the boolean vector to a vector of free dof indices
        inner_solver.init_free_dof(barrier_problem_ptr->is_dof_fixed());

        num_outer_iterations_ = 0;
    }

    OptimizationResults BarrierSolver::step_solve()
    {
        assert(general_problem_ptr != nullptr);
        assert(barrier_problem_ptr != nullptr);

        OptimizationResults results;
        IBarrierOptimizationSolver& inner_solver = get_inner_solver();

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

    OptimizationResults BarrierSolver::solve()
    {
        init_solve();

        OptimizationResults results;
        do {
            auto msg = fmt::format(
                "class=BarrierSolver function=solve.step epsilon={}",
                barrier_epsilon());

            results = step_solve();

        } while (barrier_epsilon() > min_barrier_epsilon);

        return results;
    }

    Eigen::VectorXd BarrierSolver::get_grad_kkt() const
    {
        return barrier_problem_ptr->eval_grad_f(barrier_problem_ptr->x0);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// BARRIER PROBLEM
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BarrierProblem::BarrierProblem(IBarrierGeneralProblem& general_problem)
        : general_problem(&general_problem)
    {
        num_vars_ = general_problem.num_vars();

        x0 = general_problem.starting_point();
    }

    const Eigen::VectorXb& BarrierProblem::is_dof_fixed()
    {
        return general_problem->is_dof_fixed();
    }

    double BarrierProblem::eval_f(const Eigen::VectorXd& x)
    {
        return eval_f_(x, /*update_cstr_set=*/true);
    }
    double BarrierProblem::eval_f_(
        const Eigen::VectorXd& x, const bool update_constraint_set)
    {
        NAMED_PROFILE_POINT("barrier_problem__eval_f", EVAL_F)
        NAMED_PROFILE_POINT("barrier_problem__eval_g", EVAL_G)

        PROFILE_START(EVAL_F)
        double f_uk = general_problem->eval_f(x);
        PROFILE_END(EVAL_F)

        PROFILE_START(EVAL_G)
        auto gx_ = general_problem->eval_g_(x, update_constraint_set);
        double gx = gx_.sum();
        PROFILE_MESSAGE(EVAL_G,
            fmt::format("epsilon,{:10e},gx_sum,{:10e}",
                general_problem->get_barrier_epsilon(), gx))

        PROFILE_END(EVAL_G)

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
        Eigen::SparseMatrix<double> f_uk_hessian;
        f_uk_hessian = general_problem->eval_hessian_f(x);
        std::vector<Eigen::SparseMatrix<double>> ddgx
            = general_problem->eval_hessian_g(x);

        for (const auto& ddgx_i : ddgx) {
            f_uk_hessian += ddgx_i;
        }
        return f_uk_hessian;
    }

    void BarrierProblem::eval_f_and_fdiff(const Eigen::VectorXd& x,
        double& f_uk,
        Eigen::VectorXd& f_uk_gradient,
        Eigen::SparseMatrix<double>& f_uk_hessian)
    {
        NAMED_PROFILE_POINT("barrier_problem__eval_f_and_fdiff", EVAL_F)
        NAMED_PROFILE_POINT("barrier_problem__eval_g_and_gdiff", EVAL_G)

        PROFILE_START(EVAL_F)
        general_problem->eval_f_and_fdiff(x, f_uk, f_uk_gradient, f_uk_hessian);
        PROFILE_END(EVAL_F)

        Eigen::VectorXd gx;
        Eigen::MatrixXd dgx;
        std::vector<Eigen::SparseMatrix<double>> ddgx;

        PROFILE_START(EVAL_G)
        general_problem->eval_g_and_gdiff(x, gx, dgx, ddgx);
        PROFILE_END(EVAL_G)

        f_uk += gx.sum();
        f_uk_gradient += dgx.colwise().sum().transpose();

        for (const auto& ddgx_i : ddgx) {
            f_uk_hessian += ddgx_i;
        }
    }

    void BarrierProblem::eval_f_and_fdiff(
        const Eigen::VectorXd& x, double& f_uk, Eigen::VectorXd& f_uk_gradient)
    {
        NAMED_PROFILE_POINT("barrier_problem__eval_f_and_fdiff", EVAL_F)
        NAMED_PROFILE_POINT("barrier_problem__eval_g_and_gdiff", EVAL_G)

        PROFILE_START(EVAL_F)
        general_problem->eval_f_and_fdiff(x, f_uk, f_uk_gradient);
        PROFILE_END(EVAL_F)

        Eigen::VectorXd gx;
        Eigen::MatrixXd dgx;
        std::vector<Eigen::SparseMatrix<double>> ddgx;

        PROFILE_START(EVAL_G)
        general_problem->eval_g_and_gdiff(x, gx, dgx, ddgx);
        PROFILE_END(EVAL_G)

        f_uk += gx.sum();
    }

} // namespace opt
} // namespace ccd
