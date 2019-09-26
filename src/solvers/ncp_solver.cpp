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
    NLOHMANN_JSON_SERIALIZE_ENUM(NCPUpdate,
        { { NCPUpdate::LINEARIZED, "linearized" },
            { NCPUpdate::G_GRADIENT, "g_gradient" } })

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
        , update_type(NCPUpdate::LINEARIZED)
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
        update_type = json["update_type"].get<NCPUpdate>();
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

    void NCPSolver::set_problem(INCPProblem& problem)
    {
        problem_ptr_ = &problem;
    }

    // Solve the internal problem
    OptimizationResults NCPSolver::solve()
    {
        auto results = solve(/*use_grad=*/true);
        bool valid = std::isfinite(results.x.maxCoeff());
        if (!results.finished || !valid) {
            spdlog::info("solver=ncp_solver action=solve "
                         "message='fallback to use normals'");
            results = solve(/*use_grad=*/false);
        }
        return results;
    }

    // Solve the internal problem using either gradient or normal
    OptimizationResults NCPSolver::solve(const bool use_grad)
    {
        m_use_gradient = use_grad;
        init_solve();

        OptimizationResults results;
        for (int i = 0; i < max_iterations && !results.finished; ++i) {
            results = step_solve();
        }

        if (results.finished) {
            spdlog::debug("solver=ncp_solver action=solve status=success");
        } else {
            spdlog::error("solver=ncp_solver action=solve status=failed "
                          "message='max_iterations reached' it={}",
                num_outer_iterations_);
        }

        return results;
    }

    void NCPSolver::init_solve()
    {
        assert(problem_ptr_ != nullptr);
        num_outer_iterations_ = 0;
        compute_linear_system(*problem_ptr_);
        compute_initial_solution();
    }

    void NCPSolver::compute_linear_system(INCPProblem& opt_problem)
    {
        Eigen::VectorXd x0 = Eigen::VectorXd::Zero(opt_problem.num_vars());
        A = opt_problem.eval_hessian_f(x0);
        b = -opt_problem.eval_grad_f(x0);
    }

    void NCPSolver::compute_initial_solution()
    {
        Asolver->compute(A);
        if (Asolver->info() != Eigen::Success) {
            spdlog::error(
                "solver=ncp_solver failure='LU decomposition of A failed'");
            assert(0);
        }

        // Solve assumming constraints are not violated
        xi = Ainv(b);
        if (m_use_gradient) {
            problem_ptr_->eval_g(xi, g_xi, jac_g_xi);
        } else {
            problem_ptr_->eval_g_normal(xi, g_xi, jac_g_xi);
        }

        zero_out_fixed_dof(problem_ptr_->is_dof_fixed(), jac_g_xi);

        // Initialize slack variables to zero to satisfy complementarity
        lambda_i.resize(g_xi.rows());
        lambda_i.setZero();
    }

    OptimizationResults NCPSolver::step_solve()
    {
        num_outer_iterations_ += 1;

        Eigen::VectorXd delta_xi = solve_lcp();

        double alpha = 1.0, alpha_prev = 1.0;
        double step_norm = (alpha * delta_xi).norm();
        double total_vol = g_xi.sum();

        bool line_search_success = false;

        if (do_line_search) {
            while (step_norm >= convergence_tolerance) {
                Eigen::VectorXd x_next = xi + alpha * delta_xi;
                Eigen::VectorXd g_next = problem_ptr_->eval_g(x_next);

                bool in_infeasible = (g_next.array() < 0.0).any();
                double total_vol_next = g_next.sum();
                spdlog::trace(
                    "solver=ncp_solver it={} action=line_search gamma={} "
                    "step_norm={} unfeasible={} vol={:.10e} vol_0={:10e}",
                    num_outer_iterations_, alpha, step_norm, in_infeasible,
                    total_vol_next, total_vol);

                if (in_infeasible && total_vol_next > total_vol) {
                    line_search_success = true;
                    step_norm = (alpha * delta_xi).norm();
                    if (step_norm < convergence_tolerance) {
                        alpha = alpha_prev;
                    }
                    break;
                }
                alpha_prev = alpha;
                alpha = alpha / 2.0;
                step_norm = (alpha * delta_xi).norm();
            } // end of line search

            if (!line_search_success) {
                alpha = alpha_prev;
            }
        }

        while ((alpha * delta_xi).norm() > 1e-10) {
            Eigen::VectorXd x_next = xi + alpha * delta_xi;
            Eigen::VectorXd g_next = problem_ptr_->eval_g(x_next);
            if (g_next.sum() > g_xi.sum()) {
                break;
            }
            alpha = alpha / 2.0;
        }

        xi = xi + alpha * delta_xi;

        if (num_outer_iterations_ == Constants::NCP_FALLBACK_ITERATIONS) {
            spdlog::warn(
                "solver=ncp_solver it={:d} "
                "message='starting to use normal instead of gradients'",
                num_outer_iterations());
        }
        if (!m_use_gradient
            || num_outer_iterations() > Constants::NCP_FALLBACK_ITERATIONS) {
            problem_ptr_->eval_g_normal(xi, g_xi, jac_g_xi);
        } else {
            problem_ptr_->eval_g(xi, g_xi, jac_g_xi);
        }

        zero_out_fixed_dof(problem_ptr_->is_dof_fixed(), jac_g_xi);
        if(g_xi.size() != 0){
            spdlog::debug("solve=ncp_solver it={:d} abs_tol={:g} min(g(x))={:g}",
                num_outer_iterations(), Constants::NCP_ABS_TOL, g_xi.minCoeff());
        }
        bool success = (g_xi.array() >= -Constants::NCP_ABS_TOL).all();
        return OptimizationResults(
            xi, problem_ptr_->eval_f(xi), success, success);
    }

    Eigen::VectorXd NCPSolver::solve_lcp()
    {
        // Linearize the problem and solve for primal variables (xᵢ₊₁) and dual
        // variables (λᵢ)
        //
        // Linearization:
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
        switch (update_type) {
        case NCPUpdate::G_GRADIENT:
            p.setZero();
            break;
        case NCPUpdate::LINEARIZED:
            p = Ainv(b) - xi;
            break;
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
        return /*delta_x = */ M * lambda_i + p;
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
