#include "lcp_solver.hpp"

#if BUILD_WITH_MOSEK
#include <igl/mosek/mosek_quadprog.h>
#include <opt/OptimizationProblem.hpp> // no-upper-bound
#endif
#include <iostream>

#include <logger.hpp>

namespace ccd {
namespace opt {

    bool lcp_solve(const Eigen::VectorXd& q, const Eigen::MatrixXd& N,
        const Eigen::MatrixXd& M, const Eigen::VectorXd& p,
        const LCPSolver solver, Eigen::VectorXd& x)
    {
        bool success = false;
        switch (solver) {
        case LCP_GAUSS_SEIDEL:
            success = lcp_gauss_seidel(q, N, M, p, x);
            break;
        case LCP_MOSEK:
            //  s = q + N*(M * x + p) --> N* M x + (N*p + q)
            success = lcp_mosek(N * M, N * p + q, x);
            break;
        }
        if (success) {
            Eigen::VectorXd s = q + N * (M * x + p);
            double err = x.transpose() * s;
            spdlog::debug("lcp_solve={} x={} s={} x^Ts={}", LCPSolverNames[solver], ccd::log::fmt_eigen(x),
                ccd::log::fmt_eigen(s), err);
        }
        return success;
    }

    double get_fischer_burmeister(Eigen::VectorXd& x, Eigen::VectorXd& s)
    {
        double fb = 0;
        assert(x.size() == s.size());

        for (int i = 0; i < x.size(); ++i) {
            double ent = x[i] + s[i] - sqrt(x[i] * x[i] + s[i] * s[i]);
            fb += ent * ent;
        }
        return sqrt(fb);
    }

    bool lcp_gauss_seidel(const Eigen::VectorXd& q, const Eigen::MatrixXd& N,
        const Eigen::MatrixXd& M, const Eigen::VectorXd& p, Eigen::VectorXd& x)
    {
        // LCP Problem:
        //      s = q + N*(M * x + p)
        //      0 <= x \perp s >=0

        const long num_constraints = q.rows();
        const long dof = p.rows();

        assert(N.cols() == dof);
        assert(N.rows() == num_constraints);
        assert(M.cols() == num_constraints);
        assert(M.rows() == dof);
        assert(p.rows() == dof);
        assert(q.rows() == num_constraints);

        // ----------------------
        // Gauss-Seidel steps
        // ----------------------

        const uint num_gs_steps = 1000;
        Eigen::VectorXd dk(dof);
        x.resize(num_constraints);
        x.setZero();

        double FB = -1.0;
        for (uint jj = 0; jj < num_gs_steps; ++jj) {
            // update each alpha
            for (uint ci = 0; ci < num_constraints; ++ci) {
                dk = p;
                for (uint cj = 0; cj < num_constraints; cj++) {
                    if (cj == ci) {
                        continue;
                    }
                    dk += M.col(cj) * x(cj);
                }
                double ndk = q[ci] + N.row(ci).dot(dk);
                // WARN: I'm not sure this is supposed to be >=0
                if (ndk >= 0) { // no contact impulse.
                    x(ci) = 0.0;
                    continue;
                }
                // TODO: the expression  N.col(ci).dot(M.col(ci)); is
                // constant in the GS iteration
                x(ci) = -ndk / N.row(ci).dot(M.col(ci));
            }
            Eigen::VectorXd s = N * (M * x + p) + q;
            FB = get_fischer_burmeister(x, s);

            if (FB < 1e-4) {
                break;
            }
        }
        if (FB >= 1e-4) {
            std::cerr << "gs_convergence=" << FB
                      << " num_constraints=" << num_constraints << std::endl;
            return false;
        }
        return true;
    }

    bool lcp_mosek(
        const Eigen::MatrixXd& M, const Eigen::VectorXd& q, Eigen::VectorXd& x)
    {
        // LCP Problem:
        //      s = (M * x + q)
        //      0 <= x \perp s >=0
        // Equivalent EQ porblem:
        //      Min  x^T M x + x^T q
        //      s.t  (Mx + q) >= 0 ----> -q <= Mx <= inf
        //                  x >= 0 ----> 0 <= x <= inf

        //
        // Minimize: 1/2 * x^T * Q * x + c^T * x + cf
        //
        // Subject to: lc ≤ Ax ≤ uc
        //             lx ≤ x ≤ ux
        const int x_dof = int(x.rows());
        const int q_dof = int(q.rows());
        assert(M.cols() == x_dof);
        assert(M.rows() == q_dof);

        Eigen::SparseMatrix<double> A = M.sparseView();
        const Eigen::VectorXd& lc = -q;
        Eigen::VectorXd uc(q_dof);
        uc.setConstant(NO_UPPER_BOUND);

        Eigen::VectorXd lx(x_dof), ux(x_dof);
        lx.setZero();
        ux.setConstant(NO_UPPER_BOUND);

        igl::mosek::MosekData mosek_data;
        bool success = igl::mosek::mosek_quadprog(
            A, q, /*cf=*/0.0, A, lc, uc, lx, ux, mosek_data, x);

        return success;
    }

} // namespace opt
} // namespace ccd
