#include "lcp_solver.hpp"

#include <Eigen/Core>

#include <logger.hpp>
#include <utils/eigen_ext.hpp>

namespace ccd {
namespace opt {
    static const double eps = 2.2204e-16;

    inline Eigen::VectorXd fischer(
        const Eigen::VectorXd& x, const Eigen::VectorXd& s)
    {
        return (x.array().pow(2) + s.array().pow(2)).sqrt().matrix() - x - s;
    }

    inline Eigen::VectorXd phi_lambda(
        const Eigen::VectorXd& a, const Eigen::VectorXd& b, double lambda)
    {
        return (lambda * fischer(a, b))
            - ((1 - lambda) * a.array().max(0) * b.array().max(0)).matrix();
    }

    inline double psi_lambda(
        const Eigen::VectorXd& a, const Eigen::VectorXd& b, double l)
    {
        Eigen::VectorXd phi_l = phi_lambda(a, b, l);
        return 0.5 * phi_l.transpose() * phi_l;
    }

    bool lcp_newton(const Eigen::MatrixXd& A_dense,
        const Eigen::VectorXd& b,
        Eigen::VectorXd& x)
    {
        const Eigen::SparseMatrix<double> A = A_dense.sparseView();
        size_t num_vars = b.rows();
        Eigen::VectorXd x0 = Eigen::VectorXd::Zero(b.rows());
        int max_iter = 3000;  // num_vars / 2;
        double tol_rel = 0.0; // 1e-4;
        double tol_abs = 0.0; // 10 * eps;

        double lambda = 1;

        //--- Here comes a bunch of magic constants ---------------------------
        double h = 1e-7;      // Fixed constant used to evaluate the directional
                              // detivative
        double alpha = 0.5;   // Step reduction parameter for projected Armijo
                              // backtracking line search
        double beta = 0.001;  // Sufficent decrease parameter for projected
                              // Armijo backtracking line search
        double gamma = 1e-28; // Perturbation values used to fix near singular
                              // points in derivative
        double rho = eps; // Descent direction test parameter used to test if
                          // the Newton direction does a good enough job.

        double old_err;
        double err = std::numeric_limits<double>::infinity();
        x = x0;
        for (int i = 0; i < max_iter; i++) {
            Eigen::VectorXd s = A * x + b; // Slack variables

            //--- Test all stopping criteria used
            //---------------------------------

            Eigen::VectorXd phi = phi_lambda(x, s, lambda);
            old_err = err;
            err = 0.5 * phi.transpose() * phi;

            if (abs(err - old_err) / abs(old_err) < tol_rel) {
                return true;
            }
            if (err < tol_abs) {
                return true;
            }

            Eigen::VectorXd delta_x = Eigen::VectorXd::Zero(num_vars);

            assert(x.rows() == phi.rows());
            Eigen::VectorXb S(x.rows()); // Bitmask for singular indices
            std::vector<size_t> I;       // Bitmask for non-singular indices
            for (size_t j = 0; j < x.rows(); j++) {
                bool is_singular = abs(phi(j)) < gamma && abs(x(j)) < gamma;
                S(j) = is_singular;
                if (!is_singular) {
                    I.push_back(j);
                }
            }

            // Perturbation: works on full system
            Eigen::VectorXd px = x;
            Eigen::VectorXd dir
                = x.unaryExpr([](double xi) { return (0 < xi) - (xi < 0); });
            // dir(dir==0) = 1;
            for (size_t j = 0; j < dir.rows(); j++) {
                if (dir(j) == 0) {
                    dir(j) = 1;
                }
            }
            // px(S==1) = gamma * dir(S==1);
            for (size_t j = 0; j < S.rows(); j++) {
                if (S(j)) {
                    px(j) = gamma * dir(j);
                }
            }

            Eigen::VectorXd p
                = (px.array() / (px.array().pow(2) + s.array().pow(2)).sqrt())
                - 1;
            Eigen::VectorXd q
                = (s.array() / (px.array().pow(2) + s.array().pow(2)).sqrt())
                - 1;
            Eigen::SparseMatrix<double> J = Eigen::SparseDiagonal<double>(p)
                    * Eigen::SparseDiagonal<double>(
                          Eigen::VectorXd::Ones(num_vars))
                + Eigen::SparseDiagonal<double>(q) * A;
            // Billups solution actually pick the gradient direction as the
            // descent direction if the Newton system cannot be solved. Look in
            // his thesis at page 85.
            //
            // //if rcond(J) < eps
            // //  dx = -J'*phi;          // Too bad we pick gradient direction
            // //end
            //
            // Instead one may use the pseudo inverse if badly conditioned
            //
            // if (rcond(J) < eps * 2) {
            //     dx = pinv(J) * (-phi); // badly conditioned
            // } else {
            //     dx = J \ (-phi); // well conditioned
            // }
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(J);
            if (solver.info() == Eigen::Success) {
                delta_x = solver.solve(-phi);
                if (solver.info() != Eigen::Success) {
                    spdlog::warn(
                        "solver=newton iter={:d} failure=\"sparse solve for newton direction failed\" failsafe=none",
                        i);
                    delta_x = -phi;
                }
            } else {
                spdlog::warn(
                    "solver=newton iter={:d} failure=\"sparse decomposition of the hessian failed\" failsafe=none",
                    i);
                delta_x = -phi; // Revert to gradient descent
            }
            //
            // Or one can hope for the best and try GMRES
            //
            //[dx ~]  = gmres( J, (-phi), restart);
            //

            // Test if the search direction is smaller than numerical precision.
            // That is if it is too close to zero.
            if (delta_x.array().abs().maxCoeff() < eps) {
                // Rather than just giving up we may just use the gradient
                // direction instead. However, I am lazy here! dx = nabla_phi'
                return false;
            }

            Eigen::VectorXd nabla_phi = phi.transpose() * J;
            if (nabla_phi.norm() < tol_abs) {
                return false;
            }

            // Test if our search direction is a 'sufficient' descent direction
            if ((nabla_phi.transpose() * delta_x)(0)
                > (-rho * (delta_x.transpose() * delta_x))(0)) {
                // Rather than just giving up we may
                // just use the gradient direction
                // instead. However, I am lazy here!
                //  dx = nabla_phi'
                return false;
            }

            //--- Armijo backtracking combined with a projected line-search
            //-------
            double tau = 1.0; // Current step length
            double f_0 = err;
            double grad_f = beta * (nabla_phi.transpose() * delta_x)(0);
            Eigen::VectorXd x_k = x;

            while (true) {

                x_k = (x + delta_x * tau).array().max(0);
                Eigen::VectorXd y_k = A * x_k + b;
                Eigen::VectorXd phi_k = phi_lambda(y_k, x_k, lambda);
                double f_k = 0.5 * (phi_k.transpose() * phi_k)(0);

                // Perform Armijo codition to see if we got a sufficient
                // decrease
                if (f_k <= f_0 + tau * grad_f) {
                    break;
                }

                // Test if time-step became too small
                if (tau * tau < gamma) {
                    return false;
                }

                tau = alpha * tau;
            }

            // Update iterate with result from Armijo backtracking
            x = x_k;
        }

        return false;
    }
} // namespace opt
} // namespace ccd
