#include "ncp_solver.hpp"

#include <iostream>

#include <igl/slice.h>
#include <igl/slice_into.h>

#include <opt/lcp_solver.hpp>
#include <opt/newtons_method.hpp> //> line search

#include <logger.hpp>

namespace ccd {
namespace opt {

    NCPSolver::~NCPSolver() {}
    NCPSolver::NCPSolver()
        : keep_in_unfeasible(true)
        , check_convergence(true)
        , check_volume_increase(true)
        , solve_for_active_cstr(true)
        , convergence_tolerance(1e-6)
        , update_type(NcpUpdate::LINEARIZED)
        , lcp_solver(LCPSolver::LCP_GAUSS_SEIDEL)
        , max_iterations(100)
    {
        Asolver
            = std::make_shared<Eigen::SparseLU<Eigen::SparseMatrix<double>>>();
    }

    bool NCPSolver::solve_ncp(const Eigen::SparseMatrix<double>& f_A,
        const Eigen::VectorXd& f_b, OptimizationProblem& opt_problem,
        Eigen::VectorXd& x_opt, Eigen::VectorXd& alpha_opt)
    {
        A = f_A;
        b = f_b;
        problem = &opt_problem;

        bool success = compute();
        x_opt = xi;
        alpha_opt = alpha_i;
        return success;
    }

    OptimizationResults NCPSolver::solve(OptimizationProblem& opt_problem)
    {
        compute_linear_system(opt_problem);
        problem = &opt_problem;

        OptimizationResults result;
        result.success = compute();
        result.x = xi;
        result.minf = opt_problem.eval_f(xi);
        return result;
    }

    void NCPSolver::compute_linear_system(OptimizationProblem& opt_problem)
    {
        Eigen::VectorXd x0 = Eigen::VectorXd::Zero(opt_problem.num_vars);
        A = opt_problem.eval_hessian_f_sparse(x0);
        b = -opt_problem.eval_grad_f(x0);
    }

    bool NCPSolver::compute()
    {
        initialize();
        problem->eval_intermediate_callback(xi);

        // 2. solve constraints with successive linearizations
        for (int i = 0; i < max_iterations; ++i) {
            bool collisions_solved = (g_xi.array() >= 0).all();
            if (collisions_solved) {
                break;
            }

            Eigen::VectorXd delta_x;
            solve_lcp(delta_x);

            // 2.2 Do line-search on delta to ensure volume is increased (-->0)
            double gamma = 1.0;
            auto eval_g = [&](const Eigen::VectorXd& y) {
                return -problem->eval_g(y).sum();
            };

            // TODO: This is sending a zero for the gradient which is not
            // correct
            if (check_volume_increase
                && !line_search(xi, delta_x, eval_g, gamma)) {
            }

            if (keep_in_unfeasible) {
                move_to_unfeasible_domain(delta_x, gamma);
            }

            update_candidate(delta_x, gamma);
        }
        return (g_xi.array() >= 0).all();
    }

    void NCPSolver::initialize()
    {
        // We solve the NCP problem
        //      Ax = b + \nabla g(x)^T \alpha
        //      0 <= alpha \perp g(x) >=0

        Asolver->compute(A);
        if (Asolver->info() != Eigen::Success) {
            std::cerr << "LU failed - ncp_solver.cpp" << std::endl;
            assert(0);
        }

        // Solve assumming constraints are not violated
        xi = Ainv(b);
        problem->eval_g(xi, g_xi, jac_g_xi, g_active);
        problem->eval_intermediate_callback(xi);

        alpha_i.resize(problem->num_constraints);
        alpha_i.setZero();
    }

    void NCPSolver::update_candidate(Eigen::VectorXd& delta_x, double& gamma)
    {
        xi = xi + delta_x * gamma;
        problem->eval_g(xi, g_xi, jac_g_xi, g_active);
        problem->eval_intermediate_callback(xi);
    }

    Eigen::VectorXd NCPSolver::Ainv(const Eigen::VectorXd& x)
    {
        auto y = Asolver->solve(x);
        if (Asolver->info() != Eigen::Success) {
            spdlog::error("Linear solve failed - ncp_solver.cpp");
            std::cerr << "Linear solve failed - ncp_solver.cpp" << std::endl;
            assert(0);
        }
        return y;
    }

    void NCPSolver::solve_lcp(Eigen::VectorXd& delta_x)
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

        uint num_constraints = uint(g_xi.rows());
        if (solve_for_active_cstr) {
            uint num_active_constraints = uint(g_active.rows());
            spdlog::trace("solver=ncp_solver active_cstr={}/{}",
                num_active_constraints, num_constraints);

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

            alpha_i.setZero();
            igl::slice_into(alpha_i_active, g_active, alpha_i);

        } else {
            Eigen::MatrixXd M(dof, num_constraints);
            for (uint ci = 0; ci < num_constraints; ++ci) {
                assert(M.cols() > ci);
                assert(uint(jac_g_xi.rows()) > ci);
                M.col(ci) = Ainv(jac_g_xi.row(int(ci)));
            }
            lcp_solve(g_xi, jac_g_xi, M, p, lcp_solver, alpha_i);
            delta_x = M * alpha_i + p;
        }
    }

    void NCPSolver::move_to_unfeasible_domain(
        Eigen::VectorXd& delta_x, double& gamma)
    {
        Eigen::VectorXd x_next;
        double gamma_prev = gamma;
        for (int j = 0; j < 32; j++) {
            x_next = xi + delta_x * gamma;
            g_xi = problem->eval_g(x_next);

            // if any is negative
            if (!(g_xi.array() >= 0).all()) {
                break;
            }
            gamma_prev = gamma; // feasible gamma
            gamma /= 2.0;
        }

        // check if point satisfy equality constraint, then
        // find the point that should have also satisfied g(x)>0
        if (check_convergence) {
            problem->eval_jac_g(x_next, jac_g_xi);
            Eigen::VectorXd eq
                = A * x_next - (b + jac_g_xi.transpose() * alpha_i);
            double eq_norm = eq.squaredNorm();
            if (eq_norm < convergence_tolerance) {
                gamma = gamma_prev;
            }
            spdlog::trace("solver=ncp_solver step=check_convergence "
                          "||eq||^2={} passed={}",
                eq_norm, eq_norm < convergence_tolerance);
        }
    }

} // namespace opt
} // namespace ccd
