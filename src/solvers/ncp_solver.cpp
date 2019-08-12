#include "ncp_solver.hpp"

#include <iostream>

#include <igl/slice.h>
#include <igl/slice_into.h>

#include <solvers/lcp_solver.hpp>
#include <solvers/line_search.hpp>

#include <logger.hpp>

namespace ccd {
namespace opt {

    // !important: this needs to be define in the enum namespace
    NLOHMANN_JSON_SERIALIZE_ENUM(NcpUpdate,
        { { NcpUpdate::LINEARIZED, "linearized" },
            { NcpUpdate::G_GRADIENT, "g_gradients" } })

    NLOHMANN_JSON_SERIALIZE_ENUM(LCPSolver,
        { { LCPSolver::LCP_MOSEK, "lcp_mosek" },
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
        IVolumeProblem& opt_problem,
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

    void NCPSolver::set_problem(IVolumeProblem& problem)
    {
        problem_ptr_ = &problem;
    }

    OptimizationResults NCPSolver::solve()
    {
        init_solve();

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
            spdlog::warn(
                "solver=ncp_solver action=solve status=failed message='max_iterations reached'");
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
        problem_ptr_->eval_g(xi, g_xi, jac_g_xi, g_active);
        zero_out_fixed_dof(problem_ptr_->is_dof_fixed(), jac_g_xi);

        lambda_i.resize(g_xi.rows());
        lambda_i.setZero();
    }

    void NCPSolver::compute_linear_system(IVolumeProblem& opt_problem)
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
            spdlog::trace(
                "solver=ncp_solver it={} step=collisions_solved xi={}",
                num_outer_iterations_, log::fmt_eigen(xi));
            return result;
        }

        Eigen::VectorXd delta_i;
        spdlog::trace("solver=ncp_solver it={} action=solve_lcp status=BEGIN",
            num_outer_iterations_);

        solve_lcp(xi, g_xi, jac_g_xi, g_active, delta_i, lambda_i);

        spdlog::trace(
            "solver=ncp_solver it={} action=solve_lcp status=END delta_x={} jac_g_xi={} lambda_i={}",
            num_outer_iterations_, log::fmt_eigen(delta_i),
            log::fmt_eigen(jac_g_xi), log::fmt_eigen(lambda_i));

        double alpha = 1.0, alpha_prev = 1.0;
        double step_norm = (alpha * delta_i).norm();
        double total_vol = g_xi.sum();

        bool line_search_success = false;
        spdlog::trace("solver=ncp_solver it={} action=line_search status=BEGIN",
            num_outer_iterations_);

        /// write file for Denis v.v
        //        int N = 1000;
        //        std::string filename = DATA_OUTPUT_DIR
        //            + fmt::format("/solve_it_{}.csv", num_outer_iterations_);
        //        spdlog::debug("writting into {}", filename);
        //        std::ofstream myfile;
        //        myfile.open(filename);
        //        myfile << fmt::format("delta_i,
        //        {}\n",log::fmt_eigen(delta_i)); myfile << "s,sum V, |Ax -(b
        //        +J^T lambda)|, g_xi, J\n"; for (int it = 0; it < N; it++) {
        //            double s = (it + 1) / double(N);
        //            Eigen::VectorXd x_next = xi + s * delta_i;
        //            problem_ptr_->eval_g(x_next, g_xi, jac_g_xi, g_active);

        //            double residual
        //                = (A * x_next - (b + jac_g_xi.transpose() *
        //                lambda_i)).norm();
        //            Eigen::VectorXd g_xi_active;
        //            igl::slice(g_xi, g_active, g_xi_active);

        //            uint dof = uint(xi.rows());
        //            Eigen::SparseMatrix<double> jac_g_xi_active;
        //            igl::slice(jac_g_xi, g_active,
        //                Eigen::VectorXi::LinSpaced(dof, 0, int(dof)),
        //                jac_g_xi_active);
        //            myfile << fmt::format("{:10e},{:10e}, {:10e}, {}, {}\n",
        //            s,
        //                g_xi.sum(), residual, log::fmt_eigen(g_xi, 10),
        //                log::fmt_eigen(jac_g_xi, 10));
        //        }
        //        myfile.close();

        /// line search

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
                        spdlog::trace(
                            "solver=ncp_solver it={} action=line_search status=SUCCESS converged=YES",
                            num_outer_iterations_);
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

        xi = xi + delta_i * alpha;
        problem_ptr_->eval_g(xi, g_xi, jac_g_xi, g_active);
        zero_out_fixed_dof(problem_ptr_->is_dof_fixed(), jac_g_xi);
        result.finished = false;
        result.success = false;
        result.x = xi;
        result.minf = problem_ptr_->eval_f(xi);
        return result;
    }

    void NCPSolver::solve_lcp(const Eigen::VectorXd& xi,
        const Eigen::VectorXd& g_xi,
        const Eigen::SparseMatrix<double>& jac_g_xi,
        const Eigen::VectorXi& g_active,
        Eigen::VectorXd& delta_x,
        Eigen::VectorXd& lambda_i) const
    {
        // 2.1 Linearize the problem and solve for \alpha
        // Linearization:
        //      A x_{i+1} = b + jac_g(x_i)^T \lambda
        //      0 <= \alpha_i \perp g(x_i) + jac_g(x_i) delta_x >=0
        // Update:
        //      x_{i+1} = x_{i} + delta_x
        // Delta x:
        //      a) delta_x = A^{-1} jac_x(x_i)^T \alpha_i
        //      b) delta_x = A^{-1} jac_x(x_i)^T \alpha_i + A^{-1}b - x_i
        //
        // We want to take our problem to the form
        //      s = q + N [M alpha + p]
        //      0 <= alpha \perp s >=0
        // clang-format off
        //      s = ...q...+ ....N......[...........M........ \alpha   + ......p......]
        //      s = g(x_i) + jac_g(x_i) [ A^{-1} jac_x(x_i)^T \alpha_i + A^{-1}b - x_i]
        // clang-format on
        uint dof = uint(xi.rows());

        Eigen::VectorXd p(dof);
        if (update_type == NcpUpdate::G_GRADIENT) {
            p.setZero();
        } else {
            p = Ainv(b) - xi;
        }

        lambda_i.setZero();

        uint num_constraints = uint(g_xi.rows());
        if (solve_for_active_cstr) {
            uint num_active_constraints = uint(g_active.rows());
            spdlog::trace(
                "solver=ncp_solver it={} action=solve_ncp active_cstr={}/{}",
                num_outer_iterations_, num_active_constraints, num_constraints);

            // create lcp problem for ACTIVE constraints only
            // ----------------------------------------------
            Eigen::VectorXd g_xi_active;
            igl::slice(g_xi, g_active, g_xi_active);
            assert(g_xi_active.rows() == num_active_constraints);

            Eigen::SparseMatrix<double> jac_g_xi_active;
            igl::slice(jac_g_xi, g_active,
                Eigen::VectorXi::LinSpaced(dof, 0, int(dof)), jac_g_xi_active);
            assert(jac_g_xi_active.rows() == int(num_active_constraints));
            assert(jac_g_xi_active.cols() == int(dof));

            Eigen::MatrixXd M_active(dof, num_active_constraints);
            for (uint ci = 0; ci < num_active_constraints; ++ci) {
                M_active.col(ci) = Ainv(jac_g_xi_active.row(int(ci)));
            }

            Eigen::VectorXd alpha_i_active(num_active_constraints);
            lcp_solve(g_xi_active, jac_g_xi_active, M_active, p, lcp_solver,
                alpha_i_active);

            delta_x = M_active * alpha_i_active + p;

            lambda_i.setZero();
            igl::slice_into(alpha_i_active, g_active, lambda_i);

        } else {
            Eigen::MatrixXd M(dof, num_constraints);
            for (uint ci = 0; ci < num_constraints; ++ci) {
                assert(M.cols() > ci);
                assert(uint(jac_g_xi.rows()) > ci);
                M.col(ci) = Ainv(jac_g_xi.row(int(ci)));
            }
            lcp_solve(g_xi, jac_g_xi, M, p, lcp_solver, lambda_i);
            delta_x = M * lambda_i + p;
        }
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

    //    void NCPSolver::eval_f(
    //        const Eigen::MatrixXd& /*points*/, Eigen::VectorXd& /*fx*/)
    //    {
    //        // TODO!!!!
    //    }

    void zero_out_fixed_dof(
        const Eigen::VectorXb& is_fixed, Eigen::SparseMatrix<double>& jac)
    {
        assert(jac.cols() == is_fixed.rows());
        int r, c;
        for (int k = 0; k < jac.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(jac, k); it;
                 ++it) {
                r = it.row(); // row index
                c = it.col(); // col index
                if (is_fixed(c)) {
                    jac.coeffRef(r, c) = 0.0;
                }
            }
        }
    }

} // namespace opt
} // namespace ccd
