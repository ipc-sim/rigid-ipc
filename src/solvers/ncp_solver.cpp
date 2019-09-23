#include "ncp_solver.hpp"

#include <iostream>

#include <igl/slice.h>
#include <igl/slice_into.h>

#include <solvers/lcp_solver.hpp>
#include <solvers/line_search.hpp>

#include <constants.hpp>
#include <logger.hpp>

namespace ccd {
namespace opt {

    // !important: this needs to be define in the enum namespace
    NLOHMANN_JSON_SERIALIZE_ENUM(NcpUpdate,
        { { NcpUpdate::LINEARIZED, "linearized" },
            { NcpUpdate::G_GRADIENT, "g_gradient" } })

    NLOHMANN_JSON_SERIALIZE_ENUM(LCPSolver,
        { { LCPSolver::LCP_NEWTON, "lcp_newton" },
            { LCPSolver::LCP_MOSEK, "lcp_mosek" },
            { LCPSolver::LCP_GAUSS_SEIDEL, "lcp_gauss_seidel" } })

    NCPSolver::~NCPSolver() {}
    NCPSolver::NCPSolver()
        : NCPSolver("ncp_solver")
    {
    }

    NCPSolver::NCPSolver(const std::string& name)
        : do_line_search(true)
        , solve_for_active_cstr(true)
        , convergence_tolerance(1e-6)
        , update_type(NcpUpdate::LINEARIZED)
        , lcp_solver(LCPSolver::LCP_GAUSS_SEIDEL)
        , max_iterations(1000)
        , num_outer_iterations_(0)
        , name_(name)

    {
        Asolver
            = std::make_shared<Eigen::SparseLU<Eigen::SparseMatrix<double>>>();
    }

    void NCPSolver::settings(const nlohmann::json& json)
    {
        max_iterations = json["max_iterations"].get<int>();
        do_line_search = json["do_line_search"].get<bool>();
        solve_for_active_cstr = json["solve_for_active_cstr"].get<bool>();
        convergence_tolerance = json["convergence_tolerance"].get<double>();
        update_type = json["update_type"].get<NcpUpdate>();
        lcp_solver = json["lcp_solver"].get<LCPSolver>();
    }

    nlohmann::json NCPSolver::settings() const
    {
        nlohmann::json json;
        json["max_iterations"] = max_iterations;
        json["do_line_search"] = do_line_search;
        json["solve_for_active_cstr"] = solve_for_active_cstr;
        json["convergence_tolerance"] = convergence_tolerance;
        json["update_type"] = update_type;
        json["lcp_solver"] = lcp_solver;
        return json;
    }

    bool NCPSolver::solve_ncp(const Eigen::SparseMatrix<double>& f_A,
        const Eigen::VectorXd& f_b,
        INCPProblem& opt_problem,
        Eigen::VectorXd& x_opt,
        Eigen::VectorXd& lambda_opt)
    {
        A = f_A;
        b = f_b;

        problem_ptr_ = &opt_problem;
        compute_initial_solution();

        OptimizationResults result;
        for (int i = 0; i < max_iterations; ++i) {
            result = step_solve();
            if (result.finished) {
                break;
            }
        }

        x_opt = xi;
        lambda_opt = lambda_i;
        return result.success;
    }

    void NCPSolver::set_problem(INCPProblem& problem)
    {
        problem_ptr_ = &problem;
    }

    OptimizationResults NCPSolver::solve()
    {
        auto results = solve(/*use_grad=*/true);
        bool valid = std::isfinite(results.x.maxCoeff());
        if (results.finished && valid) {
            return results;
        }
        spdlog::info("fallback to use normals");
        return solve(/*use_grad=*/false);
    }

    OptimizationResults NCPSolver::solve(const bool use_grad)
    {
        m_use_gradient = use_grad;
        init_solve();
        debug.str("");

        OptimizationResults result;
        for (int i = 0; i < max_iterations; ++i) {
            result = step_solve();
            if (result.finished) {
                break;
            }
        }

        if (result.finished) {
            spdlog::debug("solver=ncp_solver action=solve status=success");
        } else {
            spdlog::error(
                "solver=ncp_solver action=solve status=failed message='max_iterations reached' it={}",
                num_outer_iterations_);
        }

        return result;
    }

    void NCPSolver::init_solve()
    {
        assert(problem_ptr_ != nullptr);
        num_outer_iterations_ = 0;
        compute_linear_system(*problem_ptr_);
        compute_initial_solution();
    }

    void NCPSolver::compute_initial_solution()
    {
        Asolver->compute(A);
        if (Asolver->info() != Eigen::Success) {
            std::cerr << "LU failed - ncp_solver.cpp" << std::endl;
            assert(0);
        }

        // 1. Solve assumming constraints are not violated
        xi = Ainv(b);
        if (m_use_gradient) {
            problem_ptr_->eval_g(xi, g_xi, jac_g_xi);
        } else {
            problem_ptr_->eval_g_normal(xi, g_xi, jac_g_xi);
        }

        zero_out_fixed_dof(problem_ptr_->is_dof_fixed(), jac_g_xi);

        lambda_i.resize(g_xi.rows());
        lambda_i.setZero();
    }

    void NCPSolver::compute_linear_system(INCPProblem& opt_problem)
    {
        Eigen::VectorXd x0 = Eigen::VectorXd::Zero(opt_problem.num_vars());
        A = opt_problem.eval_hessian_f(x0);
        b = -opt_problem.eval_grad_f(x0);
    }

    OptimizationResults NCPSolver::step_solve()
    {

        OptimizationResults result;
        num_outer_iterations_ += 1;

        if ((g_xi.array() >= 0.0).all()) {
            result.finished = true;
            result.success = true;
            result.x = xi;
            result.minf = problem_ptr_->eval_f(xi);
            spdlog::trace("solver=ncp_solver it={} step=collisions_solved}",
                num_outer_iterations_);
            return result;
        }

        Eigen::VectorXd delta_i;
        solve_lcp(xi, g_xi, jac_g_xi, delta_i, lambda_i);

        double alpha = 1.0, alpha_prev = 1.0;
        double step_norm = (alpha * delta_i).norm();
        double total_vol = g_xi.sum();

        bool line_search_success = false;

        if (do_line_search) {
            while (step_norm >= convergence_tolerance) {
                Eigen::VectorXd x_next = xi + delta_i * alpha;
                Eigen::VectorXd g_next = problem_ptr_->eval_g(x_next);

                bool in_unfeasible = (g_next.array() < 0.0).any();
                double total_vol_next = g_next.sum();
                spdlog::trace(
                    "solver=ncp_solver it={} action=line_search gamma={} step_norm={} unfeasible={} vol={:.10e} vol_0={:10e}",
                    num_outer_iterations_, alpha, step_norm, in_unfeasible,
                    total_vol_next, total_vol);

                if (in_unfeasible && total_vol_next > total_vol) {
                    line_search_success = true;
                    step_norm = (alpha * delta_i).norm();
                    if (step_norm < convergence_tolerance) {
                        alpha = alpha_prev;
                    }
                    break;
                }
                alpha_prev = alpha;
                alpha = alpha / 2.0;
                step_norm = (alpha * delta_i).norm();
            } /// end of line search

            if (!line_search_success) {
                alpha = alpha_prev;
            }
        }

        while ((delta_i * alpha).norm() > 1e-10) {
            Eigen::VectorXd x_next = xi + delta_i * alpha;
            Eigen::VectorXd g_next = problem_ptr_->eval_g(x_next);
            if (g_next.sum() > g_xi.sum()) {
                break;
            }
            alpha = alpha / 2.0;
        }

        xi = xi + delta_i * alpha;

        if (num_outer_iterations_ == Constants::NCP_FALLBACK_ITERATIONS) {
            spdlog::warn("starting to use normal instead of gradients");
        }
        if (!m_use_gradient
            || num_outer_iterations_ > Constants::NCP_FALLBACK_ITERATIONS) {
            problem_ptr_->eval_g_normal(xi, g_xi, jac_g_xi);
        } else {
            problem_ptr_->eval_g(xi, g_xi, jac_g_xi);
        }

        zero_out_fixed_dof(problem_ptr_->is_dof_fixed(), jac_g_xi);
        result.finished = false;
        result.success = false;
        result.x = xi;
        result.minf = problem_ptr_->eval_f(xi);
        return result;
    }

    void NCPSolver::solve_lcp(const Eigen::VectorXd& xi,
        const Eigen::VectorXd& g_xi,
        const Eigen::MatrixXd& jac_g_xi,
        Eigen::VectorXd& delta_x,
        Eigen::VectorXd& lambda_i) const
    {
        // 2.1 Linearize the problem and solve for primal variables (xᵢ₊₁) and
        // dual variables (λᵢ) Linearization:
        //      A xᵢ₊₁ = b + ∇g(xᵢ)ᵀ λᵢ
        //      0 ≤ λᵢ ⟂ g(xᵢ) + ∇g(xᵢ) Δx ≥ 0
        // Update:
        //      xᵢ₊₁ = xᵢ + Δx
        // Δx:
        //      g_gradient update) Δx = A⁻¹ [∇g(xᵢ)]ᵀ λᵢ
        //      linearized update) Δx = A⁻¹ [∇g(xᵢ)]ᵀ λᵢ + A⁻¹b - xᵢ
        //
        // We want to take our problem to the form
        //      s = q + N (Mλᵢ + p)
        //      0 ≤ λᵢ ⟂ s ≥ 0
        // where
        //      q = g(xᵢ)
        //      N = ∇g(xᵢ)
        //      M = A⁻¹ [∇g(xᵢ)]ᵀ
        //      p = Δx - A⁻¹ [∇g(xᵢ)]ᵀ λᵢ
        uint dof = uint(xi.rows());

        // p
        Eigen::VectorXd p(dof);
        if (update_type == NcpUpdate::G_GRADIENT) {
            p.setZero();
        } else { // update_type == NcpUpdate::LINEARIZED
            p = Ainv(b) - xi;
        }

        lambda_i.setZero();

        uint num_constraints = uint(g_xi.rows());
        Eigen::MatrixXd M(dof, num_constraints);
        // Compute M = A⁻¹ [∇g(xᵢ)]ᵀ as a series of linear system solves
        for (uint ci = 0; ci < num_constraints; ++ci) {
            assert(M.cols() > ci);
            assert(uint(jac_g_xi.rows()) > ci);
            M.col(ci) = Ainv(jac_g_xi.row(int(ci)));
        }

        // lcp_solve(q, N, M, p)
        lcp_solve(g_xi, jac_g_xi, M, p, lcp_solver, lambda_i);
        // Δx = A⁻¹ [∇g(xᵢ)]ᵀ λᵢ + (A⁻¹b - xᵢ) * (update_type == LINEARIZED)
        delta_x = M * lambda_i + p;
    }

    Eigen::VectorXd NCPSolver::Ainv(const Eigen::VectorXd& x) const
    {
        auto y = Asolver->solve(x);
        if (Asolver->info() != Eigen::Success) {
            spdlog::error("Linear solve failed - ncp_solver.cpp");
            std::cerr << "Linear solve failed - ncp_solver.cpp" << std::endl;
            assert(0);
        }
        return y;
    }

    Eigen::VectorXd NCPSolver::get_grad_kkt() const
    {
        return A * xi - (b + jac_g_xi.transpose() * lambda_i);
    }

    void zero_out_fixed_dof(
        const Eigen::VectorXb& is_fixed, Eigen::MatrixXd& jac)
    {
        assert(jac.cols() == is_fixed.rows());
        for (int i = 0; i < is_fixed.rows(); ++i) {
            if (is_fixed(i)) {
                jac.col(i).setZero();
            }
        }
    }

} // namespace opt
} // namespace ccd
