#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ccd {
namespace opt {

    /**
     *  We handle the LCP problem LCP(A, b), which seeks vectors
     *  s and x which satisfy the following constraints
     *      s = A * x + b
     *      0 <= x \perp s >=0
     */
    enum LCPSolver { LCP_GAUSS_SEIDEL, LCP_MOSEK };
    static const char* LCPSolverNames[] = { "LCP_GAUSS_SEIDEL", "LCP_MOSEK" };

    bool lcp_solve(const Eigen::VectorXd& q, const Eigen::MatrixXd& N,
        const Eigen::MatrixXd& M, const Eigen::VectorXd& p,
        const LCPSolver solver, Eigen::VectorXd& x);

    /**
     * solves the lcp problem using gauss-seidel
     *
     * To use our gauss-seidel implementation we require to express
     * the problem as
     *      s = q + N*(M * x + p)
     * */
    bool lcp_gauss_seidel(const Eigen::VectorXd& q, const Eigen::MatrixXd& N,
        const Eigen::MatrixXd& M, const Eigen::VectorXd& p, Eigen::VectorXd& x);

    bool lcp_mosek(
        const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x);

} // namespace opt
} // namespace ccd
