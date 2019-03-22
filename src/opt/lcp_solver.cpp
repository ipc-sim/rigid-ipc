#include "lcp_solver.hpp"

#include <iostream>

namespace ccd {
namespace opt {
    /**
     *  We handle the LCP problem LCP(M,q), which seeks vectors
     *  s and x which satisfy the following constraints
     *      s = N (M x + p) + q
     *      x >=0, s >=0, x^T s = 0
     */
    void lcp_solve(const Eigen::VectorXd& q, const Eigen::MatrixXd& N,
        const Eigen::VectorXd& p, const Eigen::MatrixXd& M,
        const LCPSolver solver, Eigen::VectorXd& x)
    {
        switch (solver) {
        case LCP_GAUSS_SEIDEL:
            return lcp_gauss_seidel(q, N, p, M, x);
        case LCP_MOSEK:
            return;
        }
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

    void lcp_gauss_seidel(const Eigen::VectorXd& q, const Eigen::MatrixXd& N,
        const Eigen::VectorXd& p, const Eigen::MatrixXd& M, Eigen::VectorXd& x)
    {

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
        }

    }

} // namespace opt
} // namespace ccd
