#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ccd {
namespace opt {

    /**
     *  We handle the LCP problem LCP(M,q), which seeks vectors
     *  s and x which satisfy the following constraints
     *      s = N (M^T x + p) + q
     *      x >=0, s >=0, x^T s = 0
     */
    enum LCPSolver { LCP_GAUSS_SEIDEL, LCP_MOSEK };

    void lcp_solve(const Eigen::VectorXd& q, const Eigen::MatrixXd& N,
        const Eigen::VectorXd& p, const Eigen::MatrixXd& M,
        const LCPSolver solver, Eigen::VectorXd& x);

    ///> solves the lcp problem using gauss-seidel
    void lcp_gauss_seidel(const Eigen::VectorXd& q, const Eigen::MatrixXd& N,
        const Eigen::VectorXd& p, const Eigen::MatrixXd& M, Eigen::VectorXd& x);

} // namespace opt
} // namespace ccd
