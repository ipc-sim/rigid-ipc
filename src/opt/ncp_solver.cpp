#include "ncp_solver.hpp"

#include <iostream>
#include <opt/lcp_solver.hpp>
#include <opt/newtons_method.hpp> //> line search

#include <logger.hpp>
namespace ccd {
namespace opt {

    bool solve_ncp(const Eigen::SparseMatrix<double>& A,
        const Eigen::VectorXd& b,
        const std::function<Eigen::VectorXd(const Eigen::VectorXd& x)>& g,
        const std::function<Eigen::MatrixXd(const Eigen::VectorXd& x)>& jac_g,
        const int max_iter, const callback_intermediate_ncp& callback,
        const NcpUpdate update_type, const LCPSolver lcp_solver,
        Eigen::VectorXd& xi, Eigen::VectorXd& alpha_i,
        const bool check_convergence, const bool check_convergence_unfeasible,
        const double convergence_tolerance)
    {
        Eigen::SparseLU<Eigen::SparseMatrix<double>> Asolver;
        // We solve the NCP problem
        //      Ax = b + \nabla g(x)^T \alpha
        //      0 <= alpha \perp g(x) >=0

        Asolver.compute(A);
        if (Asolver.info() != Eigen::Success) {
            std::cerr << "LU failed - ncp_solver.cpp" << std::endl;
            assert(0);
        }

        // wrapper for the solve to check for success
        auto Ainv = [&Asolver](const Eigen::VectorXd x) {
            auto y = Asolver.solve(x);
            if (Asolver.info() != Eigen::Success) {
                std::cerr << "Linear solve failed - ncp_solver.cpp"
                          << std::endl;
                assert(0);
            }
            return y;
        };

        // 1. solve assumming constraints are not violated
        xi = Ainv(b);
        Eigen::VectorXd g_xi = g(xi);
        Eigen::MatrixXd jac_g_xi = jac_g(xi);

        uint num_constraints = uint(g_xi.rows());
        uint dof = uint(xi.rows());
        alpha_i.resize(num_constraints);
        alpha_i.setZero();

        callback(xi, alpha_i, 1.0);

        // 2. solve constraints with successive linearizations
        for (int i = 0; i < max_iter; ++i) {
            // Step 2 ends when all constraints are satisfied
            // and [optional] when the equality condition converges

            // TODO: we need to be able to test this when recomputing the
            // collision set
            //            Eigen::VectorXd eq = A * xi - (b +
            //            jac_g_xi.transpose() * alpha_i);
            if ((g_xi.array() >= 0).all()) {
                //              && (eq.squaredNorm() < convergence_tolerance ||
                //              !check_convergence)
                break;
            }

            // number of constraints can change when we update the collision set
            num_constraints = uint(g_xi.rows());
            alpha_i.resize(num_constraints);
            alpha_i.setZero();

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

            Eigen::VectorXd delta_x = M * alpha_i + p;

            // 2.2 Do line-search on delta to ensure volume is increased (-->0)
            auto eval_g = [&](const Eigen::VectorXd& y) { return -g(y).sum(); };
            double gamma = 1.0;
            if (!line_search(xi, delta_x, eval_g, gamma)) {
                // ?????
            }

            // 2.3 [optinal] Enforce staying in the unfeasible domain until
            // we converge to a point that solve the equality constraint.

            // TODO: we would like to do this check with a FIXED set of
            // collisions
            if (check_convergence_unfeasible) {
                Eigen::VectorXd x_next;
                double gamma_prev = gamma;
                std::cout << fmt::format("[{}] gamma_0={:.3e}", i, gamma);
                for (int j = 0; j < 32; j++) {
                    x_next = xi + delta_x * gamma;
                    g_xi = g(x_next);
                    if (j == 0) {
                        std::cout << fmt::format(" gxi_0={}\n", g_xi.sum() > 0);
                    }
                    if (!(g_xi.array() >= 0).all()) {
                        break;
                    }
                    gamma_prev = gamma; // feasible gamma
                    gamma /= 2.0;
                }

                // check if point satisfy equality constraint, then
                // find the point that should have also satisfied g(x)>0
                jac_g_xi = jac_g(x_next);
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

            // 2.4 Update candidate:
            xi = xi + delta_x * gamma;
            g_xi = g(xi);
            jac_g_xi = jac_g(xi);

            callback(xi, alpha_i, gamma);
        }

        return (g_xi.array() >= 0).all();
    }

} // namespace opt
} // namespace ccd
