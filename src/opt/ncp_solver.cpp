#include "ncp_solver.hpp"

#include <iostream>
#include <opt/lcp_solver.hpp>
#include <opt/newtons_method.hpp> //> line search

#include <logger.hpp>
namespace ccd {
namespace opt {

    NCPSolver::~NCPSolver() {}
    NCPSolver::NCPSolver()
        : keep_in_unfeasible(true)
        , check_convergence(false)
        , convergence_tolerance(1e-6)
        , update_type(NcpUpdate::LINEARIZED)
        , lcp_solver(LCPSolver::LCP_GAUSS_SEIDEL)
        , max_iterations(100)
    {
        Asolver
            = std::make_shared<Eigen::SparseLU<Eigen::SparseMatrix<double>>>();
    }

    NCPSolver::NCPSolver(const bool keep_in_unfeasible,
        const bool check_convergence, const double convergence_tolerance,
        const NcpUpdate update_type, const LCPSolver lcp_solver,
        const int max_iter)
        : keep_in_unfeasible(keep_in_unfeasible)
        , check_convergence(check_convergence)
        , convergence_tolerance(convergence_tolerance)
        , update_type(update_type)
        , lcp_solver(lcp_solver)
        , max_iterations(max_iter)
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
            auto eval_g = [&](const Eigen::VectorXd& y) {
                return -problem->eval_g(y).sum();
            };
            double gamma = 1.0;
            // TODO: This is sending a zero for the gradient which is not
            // correct
            if (!line_search(xi, delta_x, eval_g,
                    Eigen::VectorXd::Zero(delta_x.rows()), gamma)) {
                // ?????
            }

            if (keep_in_unfeasible) {
                move_to_unfeasible_domain(delta_x, gamma);
            }

            update_candidate(delta_x, gamma);
            problem->eval_intermediate_callback(xi);
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
        g_xi = problem->eval_g(xi);
        jac_g_xi = problem->eval_jac_g(xi);

        alpha_i.resize(uint(g_xi.rows()));
        alpha_i.setZero();

        //        callback(xi, alpha_i, 1.0);
    }

    Eigen::VectorXd NCPSolver::Ainv(const Eigen::VectorXd& x)
    {
        auto y = Asolver->solve(x);
        if (Asolver->info() != Eigen::Success) {
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
        uint num_constraints = uint(g_xi.rows());

        Eigen::VectorXd p(dof);
        if (update_type == NcpUpdate::G_GRADIENT) {
            p.setZero();
        } else {
            p = Ainv(b) - xi;
        }

        Eigen::MatrixXd M(dof, num_constraints);
        assert(num_constraints == jac_g_xi.rows());
        assert(dof == jac_g_xi.cols());

        // computes A^{-1} jac_x(x_i)^T for each row of jac_x
        for (uint ci = 0; ci < num_constraints; ++ci) {
            M.col(ci) = Ainv(jac_g_xi.row(ci));
        }

        lcp_solve(/*q=*/g_xi, /*N=*/jac_g_xi, M, p, lcp_solver, alpha_i);

        delta_x = M * alpha_i + p;
    }

    void NCPSolver::move_to_unfeasible_domain(
        Eigen::VectorXd& delta_x, double& gamma)
    {
        Eigen::VectorXd x_next;
        double gamma_prev = gamma;
        for (int j = 0; j < 32; j++) {
            x_next = xi + delta_x * gamma;
            g_xi = problem->eval_g(x_next);
            //            if (j == 0) {
            //                std::cout << fmt::format(" gxi_0={}\n", g_xi.sum()
            //                > 0);
            //            }
            if (!(g_xi.array() >= 0).all()) {
                break;
            }
            gamma_prev = gamma; // feasible gamma
            gamma /= 2.0;
        }

        // check if point satisfy equality constraint, then
        // find the point that should have also satisfied g(x)>0
        jac_g_xi = problem->eval_jac_g(x_next);
        // TODO: we need to be able to test this when recomputing the
        // collision set
        //                eq = A * x_next - (b + jac_g_xi.transpose() *
        //                alpha_i); std::cout << fmt::format("
        //                gamma_n={:.3e} eq={:.3e} gxi_n={}\n",
        //                    gamma, eq.squaredNorm(), g_xi.sum() > 0);

        //                if (eq.squaredNorm() < convergence_tolerance)
        //                {
        //                    gamma = gamma_prev;
        //                }
    }

    void NCPSolver::update_candidate(Eigen::VectorXd& delta_x, double& gamma)
    {
        xi = xi + delta_x * gamma;
        g_xi = problem->eval_g(xi);
        jac_g_xi = problem->eval_jac_g(xi);

        //        callback(xi, alpha_i, gamma);
    }

    bool solve_ncp(const Eigen::SparseMatrix<double>& A,
        const Eigen::VectorXd& b, OptimizationProblem& problem,
        const int max_iter, const callback_intermediate_ncp& callback,
        const NcpUpdate update_type, const LCPSolver lcp_solver,
        Eigen::VectorXd& xi, Eigen::VectorXd& alpha_i,
        const bool keep_in_unfeasible, const bool check_convergence,
        const double convergence_tolerance)
    {

        auto solver
            = std::make_unique<NCPSolver>(keep_in_unfeasible, check_convergence,
                convergence_tolerance, update_type, lcp_solver, max_iter);
        bool sucess = solver->solve_ncp(A, b, problem, xi, alpha_i);

        return sucess;
    }

} // namespace opt
} // namespace ccd
