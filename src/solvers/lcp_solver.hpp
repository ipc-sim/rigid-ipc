#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <opt/optimization_problem.hpp>

namespace ccd {
namespace opt {

    /**
     * We handle the LCP(A, b), which seeks vectors s and x which satisfy the
     * following constraints
     * \f{align*}
     *     s = A x + b \\
     *     0 \leq x \perp s \geq 0
     * \f}
     */
    enum LCPSolver {
        /// Solve the LCP using the Guass-Seidel method
        LCP_GAUSS_SEIDEL,
        /// Solve the LCP as a QP using Mosek
        LCP_MOSEK,
        /// Sove the LCP using minimum-map Newton's method
        LCP_NEWTON
    };

    /**
     * @brief Solve the LCP using a given solver.
     */
    bool lcp_solve(const Eigen::VectorXd& gxi,
        const Eigen::MatrixXd& jac_gxi,
        const Eigen::MatrixXd& tilde_jac_gxi,
        const Eigen::VectorXd& tilde_b,
        const LCPSolver solver,
        Eigen::VectorXd& alpha);

    /**
     * @brief Solve the LCP using Gauss-Seidel
     *
     * To use our gauss-seidel implementation we require to express
     * the problem as
     *      s = q + N (Mx + p)
     *
     * @param gxi \f$q\f$
     * @param jac_gxi \f$N\f$
     * @param tilde_jac_gxi \f$M\f$
     * @param tilde_b \f$p\f$
     * @param alpha \f$x\f$
     */
    bool lcp_gauss_seidel(const Eigen::VectorXd& gxi,
        const Eigen::MatrixXd& jac_gxi,
        const Eigen::MatrixXd& tilde_jac_gxi,
        const Eigen::VectorXd& tilde_b,
        Eigen::VectorXd& alpha);

#if BUILD_WITH_MOSEK
    bool lcp_mosek(
        const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x);
#endif

    bool lcp_newton(const Eigen::MatrixXd& A_dense,
        const Eigen::VectorXd& b,
        Eigen::VectorXd& x);

} // namespace opt
} // namespace ccd
