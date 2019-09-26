#include "lcp_solver.hpp"

#include <Eigen/Core>

#include <logger.hpp>
#include <utils/eigen_ext.hpp>

namespace ccd {
namespace opt {

    /// @brief Compute the Fischer error measure.
    inline Eigen::VectorXd fischer(
        const Eigen::VectorXd& x, const Eigen::VectorXd& s)
    {
        return (x.array().pow(2) + s.array().pow(2)).sqrt().matrix() - x - s;
    }

    /// @brief Compute the sign of the input.
    inline double sign(double x) { return (0 < x) - (x < 0); }

    // Solve the LCP using Fischer-Newton:
    //      s = Ax + b
    //      0 ≤ x ⟂ s ≥ 0
    bool lcp_newton(const Eigen::MatrixXd& A_dense,
        const Eigen::VectorXd& b,
        Eigen::VectorXd& x)
    {
        const Eigen::SparseMatrix<double> A = A_dense.sparseView();
        int num_vars = b.rows();
        Eigen::VectorXd x0 = Eigen::VectorXd::Zero(num_vars);

        // Function constants
        const double EPS = 2.2204e-16; // Machine epsilon
        const int MAX_ITER = 3000;
        const double TOL_REL = 1e-5;
        const double TOL_ABS = 0e-10;
        // Perturbation values used to fix near singular points in derivative
        const double GAMMA = 1e-28;

        double old_err;
        double err = std::numeric_limits<double>::infinity();
        x = x0;
        for (int i = 0; i < MAX_ITER; i++) {
            Eigen::VectorXd s = A * x + b; // Slack variables

            // Test all stopping criteria used
            Eigen::VectorXd phi = fischer(x, s);
            old_err = err;
            err = phi.squaredNorm() / 2;
            // Test relative error
            if (abs(err - old_err) / abs(old_err) < TOL_REL) {
                return true;
            }
            // Test absolute error
            if (err < TOL_ABS) {
                return true;
            }

            // Bitmask for singular indices
            assert(x.rows() == phi.rows());
            auto S = phi.array().abs() < GAMMA && x.array().abs() < GAMMA;

            // Perturbation: works on full system
            Eigen::VectorXd px = x;
            Eigen::VectorXd dir = x.unaryExpr(&sign);
            // dir(dir==0) = 1;
            dir += (dir.array() == 0).matrix().cast<double>();

            // px(S==1) = GAMMA * dir(S==1);
            for (size_t j = 0; j < S.rows(); j++) {
                if (S(j)) {
                    px(j) = GAMMA * dir(j);
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

            Eigen::VectorXd nabla_phi = J.transpose() * phi;
            // Test if we have dropped into a local minimia if so we are stuck
            if (nabla_phi.norm() < TOL_ABS) {
                spdlog::warn("solver=ficher_newton_lcp iter={:d} "
                             "failure=\"||∇ϕ|| < TOL_ABS\" ||∇ϕ||={:g} "
                             "failsafe=none",
                    i, nabla_phi.norm());
                return false;
            }

            // Solve Δx = -J⁻¹ϕ
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(J);
            Eigen::VectorXd delta_x;
            if (solver.info() == Eigen::Success) {
                delta_x = solver.solve(-phi);
                if (solver.info() != Eigen::Success) {
                    spdlog::warn(
                        "solver=ficher_newton_lcp iter={:d} "
                        "failure=\"sparse solve for newton direction failed\" "
                        "failsafe=\"revert to gradient descent\"",
                        i);
                    delta_x = -nabla_phi; // Revert to gradient descent
                }
            } else {
                spdlog::warn(
                    "solver=ficher_newton_lcp iter={:d} "
                    "failure=\"sparse decomposition of the hessian failed\" "
                    "failsafe=\"revert to gradient descent\"",
                    i);
                delta_x = -nabla_phi; // Revert to gradient descent
            }

            // Test if the search direction is smaller than numerical precision.
            if (delta_x.array().abs().maxCoeff() < EPS) {
                spdlog::warn(
                    "solver=ficher_newton_lcp iter={:d} "
                    "failure=\"search direction too small\" max(Δx)={:g} "
                    "failsafe=none",
                    i, delta_x.array().abs().maxCoeff());
                return false;
            }

            // Test if our search direction is a 'sufficient' descent direction
            double descent_magnitude = (nabla_phi.transpose() * delta_x)(0);
            if (descent_magnitude > -EPS * descent_magnitude) {
                spdlog::warn(
                    "solver=ficher_newton_lcp iter={:d} "
                    "failure=\"search direction not descent direction\" "
                    "(∇ϕ)ᵀΔx={:g} failsafe=\"revert to gradient descent\"",
                    i, descent_magnitude);
                delta_x = -nabla_phi; // Revert to gradient descent
            }

            // Armijo backtracking combined with a projected line-search
            double alpha = 1.0; // Current step length
            double f_0 = err;
            // Sufficent decrease parameter for Armijo backtracking line search
            double beta = 0.001;
            double grad_f = beta * descent_magnitude;
            Eigen::VectorXd x_k;

            while (true) {
                x_k = (x + alpha * delta_x).array().max(0); // project it
                Eigen::VectorXd s_k = A * x_k + b;
                double f_k = fischer(x_k, s_k).squaredNorm() / 2;

                // Perform Armijo codition to see if we got a sufficient
                // decrease
                if (f_k <= f_0 + alpha * grad_f) {
                    break;
                }

                // Test if step length became too small
                if (alpha * alpha < GAMMA) {
                    spdlog::warn(
                        "solver=ficher_newton_lcp iter={:d} "
                        "failure=\"step length too small\" step_length={:g}",
                        i, alpha);
                    return false;
                }

                alpha /= 2;
            }

            // Update iterate with result from Armijo backtracking
            x = x_k;
        }

        spdlog::warn("solver=ficher_newton_lcp "
                     "failure=\"too many iterations\" MAX_ITER={:d}",
            MAX_ITER);
        return false;
    }
} // namespace opt
} // namespace ccd
