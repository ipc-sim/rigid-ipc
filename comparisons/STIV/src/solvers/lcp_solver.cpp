#include "lcp_solver.hpp"

#if BUILD_WITH_MOSEK
#include <igl/mosek/mosek_quadprog.h>
#include <opt/optimization_problem.hpp> // no-upper-bound
#endif
#include <iostream>

#include <constants.hpp>
#include <logger.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace opt {

    static const char* LCPSolverNames[]
        = { "LCP_GAUSS_SEIDEL", "LCP_MOSEK", "LCP_NEWTON" };

    bool lcp_solve(const Eigen::VectorXd& gxi,
        const Eigen::MatrixXd& jac_gxi,
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
        case LCP_MOSEK:
#if BUILD_WITH_MOSEK
            //  s = jac_gxi * (tilde_jac_gxi * alpha + tilde_b) + gx
            //  s = jac_gxi * tilde_jac_gxi * alpha + jac_gxi * tilde_b + gxi
            success = lcp_mosek(
                jac_gxi * tilde_jac_gxi, jac_gxi * tilde_b + gxi, alpha);
            break;
#else
            throw NotImplementedError("Mosek is not loaded");
#endif
        case LCP_NEWTON:
            success = lcp_newton(
                jac_gxi * tilde_jac_gxi, jac_gxi * tilde_b + gxi, alpha);
            break;
        }

        Eigen::VectorXd s
            = jac_gxi * tilde_jac_gxi * alpha + (jac_gxi * tilde_b + gxi);
        double err = alpha.transpose() * s;
        spdlog::debug("solver=lcp_solver lcp_solver={} xᵀs={:g}",
            LCPSolverNames[solver], err);

        return success;
    }

    double get_fischer_burmeister(Eigen::VectorXd& x, Eigen::VectorXd& s)
    {
        assert(x.size() == s.size());
        Eigen::VectorXd e
            = x + s - (x.array().pow(2) + s.array().pow(2)).sqrt().matrix();
        return sqrt(e.sum());
    }

    bool lcp_gauss_seidel(const Eigen::VectorXd& gxi, // q
        const Eigen::MatrixXd& jac_gxi,               // N
        const Eigen::MatrixXd& tilde_jac_gxi,         // M
        const Eigen::VectorXd& tilde_b,               // p
        Eigen::VectorXd& alpha)
    {
        // LCP Problem:
        //      s = q + N (Mx + p)
        //      0 ≤ x ⟂ s ≥ 0
        // where
        //      q = g(xᵢ)
        //      N = ∇g(xᵢ)
        //      M = Ñ
        //      q = b̃

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

        Eigen::VectorXd dk(dof);
        alpha.resize(num_constraints);
        alpha.setZero();

        double FB = -1.0;
        for (uint jj = 0; jj < Constants::GUASS_SEIDEL_MAX_ITER; ++jj) {
            // update each alpha
            for (uint ci = 0; ci < uint(num_constraints); ++ci) {
                dk = tilde_b;
                for (uint cj = 0; cj < uint(num_constraints); ++cj) {
                    if (cj == ci) {
                        continue;
                    }
                    dk += tilde_jac_gxi.col(int(cj)) * alpha(int(cj));
                }
                double ndk = gxi[int(ci)] + jac_gxi.row(int(ci)).dot(dk);

                // WARN: I'm not sure this is supposed to be >=0
                if (ndk >= 0) { // no contact impulse.
                    alpha(int(ci)) = 0.0;
                    continue;
                }
                // TODO: the expression  N.col(ci).dot(Ntilde.col(ci)); is
                // constant in the GS iteration
                double dividend
                    = jac_gxi.row(int(ci)).dot(tilde_jac_gxi.col(int(ci)));
                double thr = std::abs(dividend)
                    / std::max(std::abs(dividend), std::abs(ndk));
                if (thr > 1e-12) {
                    alpha(int(ci)) = -ndk / dividend;
                } else {
                    alpha(int(ci)) = 0.0;
                }
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
        //      s = Mx + q
        //      0 ≤ x ⟂ s ≥ 0
        // Equivalent QP:
        //      min   xᵀ M x + xᵀq
        //      s.t.  (Mx + q) ≥ 0 → -q ≤ Mx
        //                   x ≥ 0 →  0 ≤ x
        //
        // Standard QP formulation:
        //      min   1/2 * x^T * Q * x + c^T * x + cf
        //      s.t.  lc ≤ Ax ≤ uc
        //            lx ≤  x ≤ ux

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
        if (success) {
            spdlog::debug("solver=mosek_lcp action=solve status=success");
        } else {
            spdlog::error("solver=mosek_lcp action=solve status=failed");
        }
        return success;
    }
#endif

} // namespace opt
} // namespace ccd
