#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ccd {
namespace opt {

    /**
     * We handle the LCP problem LCP(A, b), which seeks vectors
     * s and x which satisfy the following constraints
     * \begin{align*}
     *     s = A x + b \\\\
     *     0 \leq x \perp s \geq 0
     * \end{align*}
     */
    enum LCPSolver {
        /// Solve the LCP using the Guass-Seidel method
        LCP_GAUSS_SEIDEL,
        /// Solve the LCP as a QP using Mosek
        LCP_MOSEK
    };
    static const char* LCPSolverNames[] = { "LCP_GAUSS_SEIDEL", "LCP_MOSEK" };

    bool lcp_solve(const Eigen::VectorXd& gxi,
        const Eigen::SparseMatrix<double>& jac_gxi,
        const Eigen::MatrixXd& tilde_jac_gxi, const Eigen::VectorXd& ci,
        const LCPSolver solver, Eigen::VectorXd& alpha);

    /**
     * solves the lcp problem using gauss-seidel
     *
     * To use our gauss-seidel implementation we require to express
     * the problem as
     *      s = q + N*(M * x + p)
     * */
    bool lcp_gauss_seidel(const Eigen::VectorXd& gxi,
        const Eigen::SparseMatrix<double>& jac_gxi,
        const Eigen::MatrixXd& tilde_jac_gxi, const Eigen::VectorXd& ci,
        Eigen::VectorXd& alpha);

    bool lcp_mosek(
        const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x);

} // namespace opt
} // namespace ccd
