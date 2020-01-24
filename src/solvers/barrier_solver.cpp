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
        : tinit(1.0)
        , t(1.0)
        , m(1.0)
        , c(0.01)
        , e_b(1e-5)
        , t_inc(2)
        , num_outer_iterations_(0)

        , name_(name)

    {
    }
    void BarrierSolver::set_problem(IBarrierGeneralProblem& problem)
    {
        general_problem_ptr = &problem;
    }

    void BarrierSolver::settings(const nlohmann::json& json)
    {
        // termination criteria variables
        e_b = json["e_b"].get<double>();
        m = json["m"].get<double>();
        tinit = json["t_init"].get<double>();
        c = json["c"].get<double>();
        t_inc = json["t_inc"].get<double>();

        inner_solver_ptr = SolverFactory::factory().get_barrier_inner_solver(
            json["inner_solver"].get<std::string>());
    }
    void BarrierSolver::inner_solver_settings(const nlohmann::json& json)
    {
        inner_solver_ptr->settings(json);
    }

    nlohmann::json BarrierSolver::settings() const
    {
        nlohmann::json json;
        json["e_b"] = e_b;
        json["t_inc"] = t_inc;
        json["t_init"] = tinit;
        json["m"] = m;
        json["c"] = c;
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

        // Find a good initial value for `t`
        int num_active_barriers;
        Eigen::VectorXd grad_B = barrier_problem_ptr->eval_grad_B(
            barrier_problem_ptr->x0, num_active_barriers);
        if (false /*num_active_barriers > 0*/) {
            //            Eigen::VectorXd grad_E
            //                =
            //                barrier_problem_ptr->eval_grad_E(barrier_problem_ptr->x0);
            //            t = tinit = -grad_B.dot(grad_E) / grad_E.dot(grad_E);
        } else {
            t = tinit;
        }

        spdlog::debug(
            "solver_init eps_barrier={} m={} t={} e_b={} c={} t_inc={}",
            barrier_epsilon(), m, t, e_b, c, t_inc);
        barrier_problem_ptr->t = t;

        debug.str("");
        debug_max_constraints = 0;
    }

    OptimizationResults BarrierSolver::step_solve()
    {
        assert(general_problem_ptr != nullptr);
        assert(barrier_problem_ptr != nullptr);

        OptimizationResults results;
        IBarrierOptimizationSolver& inner_solver = get_inner_solver();

        spdlog::debug(
            "\tsolve_step BEGIN it={} eps_barrier={} m={} t={} e_b={} c={}",
            num_outer_iterations_, barrier_epsilon(), m, t, e_b, c);

        // propagate variables
        barrier_problem_ptr->t = t;
        inner_solver_ptr->e_b(e_b);
        inner_solver_ptr->c(c);
        inner_solver_ptr->t(t);
        inner_solver_ptr->m(m);

        int num_constraints = 0;
        barrier_problem_ptr->eval_grad_B(
            barrier_problem_ptr->x0, num_constraints);
        debug_max_constraints
            = std::max(debug_max_constraints, num_constraints);

        results = inner_solver.solve(*barrier_problem_ptr);

#ifdef DEBUG_LINESEARCH
        Eigen::VectorXd xdiff = barrier_problem_ptr->x0 - results.x;
        Eigen::MatrixXd vdiff
            = barrier_problem_ptr->debug_vertices(barrier_problem_ptr->x0)
            - barrier_problem_ptr->debug_vertices(results.x);

        double min_dist = barrier_problem_ptr->debug_min_distance(results.x);
        double min_dist_diff;
        if (min_dist > -1) {
            min_dist_diff = barrier_problem_ptr->debug_min_distance(
                                barrier_problem_ptr->x0)
                - min_dist;
        } else {
            min_dist_diff = -1;
        }
        std::string min_dist_str
            = min_dist > -1 ? fmt::format("{:.10e}", min_dist) : "NA";
        std::string min_dist_diff_str
            = min_dist_diff > -1 ? fmt::format("{:.10e}", min_dist_diff) : "NA";

        double Ex = general_problem_ptr->eval_f(results.x);
        double Bx = general_problem_ptr->eval_g(results.x).sum();
        debug << fmt::format("{},{},{:.10e},{:.10e},{},{},{:.10e}\n",
            num_outer_iterations_, t, xdiff.norm(), vdiff.norm(), min_dist_str,
            min_dist_diff_str, Ex, Bx);
#endif

        results.minf = general_problem_ptr->eval_f(results.x);
        // Start next iteration from the ending optimal position
        barrier_problem_ptr->x0 = results.x;
        t *= t_inc;

        num_outer_iterations_ += 1;
        spdlog::debug("\tsolve_step END it={} epsilon={} m / t ={} e_b={}",
            num_outer_iterations_, barrier_epsilon(), m / t, e_b);
        return results;
    }

    OptimizationResults BarrierSolver::solve()
    {
        init_solve();
        inner_solver_ptr->debug_reset_stats();

        OptimizationResults results;
        do {
            results = step_solve();
        } while (m / t > e_b); // barrier_epsilon() > min_barrier_epsilon

        int num_constraints = 0;
        barrier_problem_ptr->eval_grad_B(
            barrier_problem_ptr->x0, num_constraints);
        debug_max_constraints
            = std::max(debug_max_constraints, num_constraints);

        // make one last iteration with exactly eb
        t = m / e_b;
        const double t_used = t;
        results = step_solve();

        double min_dist = barrier_problem_ptr->debug_min_distance(results.x);
        std::cout << fmt::format("GREP_ME,t,{},min_dist,{},num_constraints,{}",
                         t_used, min_dist, debug_max_constraints)
                  << std::endl;

#ifdef DEBUG_LINESEARCH
        std::cout
            << fmt::format(
                   "tinit={} tend={} m={}  e_b={} c={} eps_barrier={} t_inc={}",
                   tinit, t / t_inc, m, e_b, c, barrier_epsilon(), t_inc)
            << std::endl;
        std::cout << "outer_it, t, var_diff_norm, vertex_diff_norm, "
                     "min_distance, min_distance_diff, E(x), B(x)"
                  << std::endl;
        std::cout << debug.str() << std::flush;
#endif
        spdlog::info(
            "BARRIER_STATS,c,{},tinit,{},tinc,{},num_outer_iterations,{},{}", c,
            tinit, t_inc, num_outer_iterations_,
            inner_solver_ptr->debug_stats());

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

    double BarrierProblem::get_termination_threshold() const
    {
        return inner_solver_threshold;
    }

    const Eigen::VectorXb& BarrierProblem::is_dof_fixed()
    {
        return general_problem->is_dof_fixed();
    }

    double BarrierProblem::eval_f(const Eigen::VectorXd& x)
    {
        NAMED_PROFILE_POINT("barrier_problem__eval_f", EVAL_F)
        NAMED_PROFILE_POINT("barrier_problem__eval_g", EVAL_G)

        PROFILE_START(EVAL_F)
        double f_uk = general_problem->eval_f(x);
        PROFILE_END(EVAL_F)

        PROFILE_START(EVAL_G)
        auto gx_ = general_problem->eval_g(x);
        double gx = gx_.sum();
        PROFILE_MESSAGE(EVAL_G,
            fmt::format("epsilon,{:10e},gx_sum,{:10e}",
                general_problem->get_barrier_epsilon(), gx))

        PROFILE_END(EVAL_G)

        //        return t * f_uk + gx;
        return f_uk + gx / t;
    }

    bool BarrierProblem::has_collisions(
        const Eigen::VectorXd& sigma_i, const Eigen::VectorXd& sigma_j) const
    {
        return general_problem->has_collisions(sigma_i, sigma_j);
    }

    Eigen::VectorXd BarrierProblem::eval_grad_f(const Eigen::VectorXd& x)
    {

        Eigen::VectorXd f_uk_gradient = general_problem->eval_grad_f(x);
        Eigen::MatrixXd dgx = general_problem->eval_jac_g(x);

        Eigen::MatrixXd g_uk_gradient;
        g_uk_gradient = dgx.colwise().sum().transpose();
        // return t * f_uk_gradient + g_uk_gradient;
        return f_uk_gradient + g_uk_gradient / t;
    }

    Eigen::SparseMatrix<double> BarrierProblem::eval_hessian_f(
        const Eigen::VectorXd& x)
    {
        Eigen::SparseMatrix<double> f_uk_hessian;
        f_uk_hessian = general_problem->eval_hessian_f(x);
        std::vector<Eigen::SparseMatrix<double>> ddgx
            = general_problem->eval_hessian_g(x);

        Eigen::SparseMatrix<double> g_uk_hessian;
        g_uk_hessian.resize(f_uk_hessian.rows(), f_uk_hessian.cols());
        g_uk_hessian.setZero();
        for (const auto& ddgx_i : ddgx) {
            g_uk_hessian += ddgx_i;
        }
        // return t * f_uk_hessian + g_uk_hessian;
        return f_uk_hessian + g_uk_hessian / t;
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

        //        f_uk = t * f_uk + gx.sum();
        //        f_uk_gradient = t * f_uk_gradient +
        //        dgx.colwise().sum().transpose();

        //        f_uk_hessian = t * f_uk_hessian;
        //        for (const auto& ddgx_i : ddgx) {
        //            f_uk_hessian += ddgx_i;
        //        }

        f_uk = f_uk + gx.sum() / t;
        f_uk_gradient = f_uk_gradient + dgx.colwise().sum().transpose() / t;

        for (const auto& ddgx_i : ddgx) {
            f_uk_hessian += ddgx_i / t;
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

        //        f_uk = t * f_uk + gx.sum();
        //        f_uk_gradient = t * f_uk_gradient +
        //        dgx.colwise().sum().transpose();
        f_uk = f_uk + gx.sum() / t;
        f_uk_gradient = f_uk_gradient + dgx.colwise().sum().transpose() / t;
    }

    Eigen::VectorXd BarrierProblem::eval_grad_E(const Eigen::VectorXd& xk)
    {
        return general_problem->eval_grad_f(xk);
    }

    Eigen::VectorXd BarrierProblem::eval_grad_B(
        const Eigen::VectorXd& xk, int& num_active_b)
    {
        Eigen::MatrixXd dgx = general_problem->eval_jac_g(xk);
        num_active_b = dgx.rows();
        return dgx.colwise().sum().transpose();
    }

    Multiprecision BarrierProblem::eval_mp_f(const Eigen::VectorXd& /*x*/)
    {

        //        double f_uk = general_problem->eval_f(x);
        //        auto gx_ = general_problem->eval_mp_g(x);
        //        Multiprecision gx = gx_.sum();

        //        return t * f_uk + gx;
        return 0;
    }

} // namespace opt
} // namespace ccd
