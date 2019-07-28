#include "lcp_solver.hpp"

#if BUILD_WITH_MOSEK
#include <igl/mosek/mosek_quadprog.h>
#include <opt/optimization_problem.hpp> // no-upper-bound
#endif
#include <iostream>

#include <logger.hpp>

namespace ccd {
namespace opt {

    bool lcp_solve(const Eigen::VectorXd& gxi,
        const Eigen::SparseMatrix<double>& jac_gxi,
        const Eigen::MatrixXd& tilde_jac_gxi,
        const Eigen::VectorXd& tilde_b,
        const LCPSolver solver,
        Eigen::VectorXd& alpha)
    {
        bool success = false;
        switch (solver) {
        case LCP_GAUSS_SEIDEL:
            success
                = lcp_gauss_seidel(gxi, jac_gxi, tilde_jac_gxi, tilde_b, alpha);
            break;
#if BUILD_WITH_MOSEK
        case LCP_MOSEK:
            //  s = jac_gxi * (tilde_jac_gxi * alpha + tilde_b) + gx
            //  s = jac_gxi * tilde_jac_gxi * alpha + (jac_gxi* tilde_b + gxi)
            success = lcp_mosek(
                jac_gxi * tilde_jac_gxi, jac_gxi * tilde_b + gxi, alpha);
            break;
#endif
        }

        Eigen::VectorXd s = jac_gxi * tilde_jac_gxi * alpha + (jac_gxi* tilde_b + gxi);
        double err = alpha.transpose() * s;
        spdlog::trace("solver=lcp_solver lcp_solver={} x^Ts={}",
            LCPSolverNames[solver], err);

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

    bool lcp_gauss_seidel(const Eigen::VectorXd& gxi,
        const Eigen::SparseMatrix<double>& jac_gxi,
        const Eigen::MatrixXd& tilde_jac_gxi,
        const Eigen::VectorXd& tilde_b,
        Eigen::VectorXd& alpha)
    {
        // LCP Problem:
        //      s = q + N*(M * x + p)
        //      0 <= x \perp s >=0

        const long num_constraints = gxi.rows();
        const long dof = tilde_b.rows();

        assert(jac_gxi.cols() == dof);
        assert(jac_gxi.rows() == num_constraints);
        assert(tilde_jac_gxi.cols() == num_constraints);
        assert(tilde_jac_gxi.rows() == dof);
        assert(tilde_b.rows() == dof);
        assert(gxi.rows() == num_constraints);

        // ----------------------
        // Gauss-Seidel steps
        // ----------------------

        const uint num_gs_steps = 10000;
        Eigen::VectorXd dk(dof);
        alpha.resize(num_constraints);
        alpha.setZero();

        double FB = -1.0;
        for (uint jj = 0; jj < num_gs_steps; ++jj) {
            // update each alpha
            for (uint ci = 0; ci < num_constraints; ++ci) {
                dk = tilde_b;
                for (uint cj = 0; cj < num_constraints; ++cj) {
                    if (cj == ci) {
                        continue;
                    }
                    dk += tilde_jac_gxi.col(cj) * alpha(cj);
                }
                double ndk = gxi[ci] + jac_gxi.row(int(ci)).dot(dk);
                // WARN: I'm not sure this is supposed to be >=0
                if (ndk >= 0) { // no contact impulse.
                    alpha(ci) = 0.0;
                    continue;
                }
                // TODO: the expression  N.col(ci).dot(Ntilde.col(ci)); is
                // constant in the GS iteration
                alpha(ci)
                    = -ndk / jac_gxi.row(int(ci)).dot(tilde_jac_gxi.col(ci));
            }
            Eigen::VectorXd s
                = jac_gxi * (tilde_jac_gxi * alpha + tilde_b) + gxi;
            FB = get_fischer_burmeister(alpha, s);

            if (FB < 1e-10) {
                break;
            }
        }
        spdlog::trace("solver=GaussSeidel gs_convergence={} num_constraints={}",
            FB, num_constraints);
        if (FB >= 1e-10) {
            spdlog::warn(
                "solver=GaussSeidel gs_convergence={} num_constraints={}", FB,
                num_constraints);
            return false;
        }
        return true;
    }

#if BUILD_WITH_MOSEK
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
#endif

} // namespace opt
} // namespace ccd
