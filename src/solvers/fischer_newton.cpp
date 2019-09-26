#include "lcp_solver.hpp"

#include <Eigen/Core>

#include <constants.hpp>
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

    // inline Eigen::SparseMatrix<double> grad_fischer(
    //     const Eigen::VectorXd& x, const Eigen::VectorXd& s)
    // {
    // }

    /// @brief Compute the sign of the input.
    inline double sign(double x) { return (0 < x) - (x < 0); }

    // Solve the LCP:
    //      s = Ax + b
    //      0 ≤ x ⟂ s ≥ 0
    // by minimizing the Fisher error:
    //      f(x) = (x² + s²) - x - s
    // using Newton's method.
    bool lcp_newton(const Eigen::MatrixXd& A_dense,
        const Eigen::VectorXd& b,
        Eigen::VectorXd& x)
    {
        const Eigen::SparseMatrix<double> A = A_dense.sparseView();
        int num_vars = b.rows();
        assert(A.rows() == num_vars);
        Eigen::VectorXd x0 = Eigen::VectorXd::Zero(num_vars);

        // Function constants
        const double EPS = std::numeric_limits<double>::epsilon();

        double old_err;
        double err = std::numeric_limits<double>::infinity();
        x = x0;
        for (int i = 0; i < Constants::FISCHER_MAX_ITER; i++) {
            Eigen::VectorXd s = A * x + b; // Slack variables

            // Test all stopping criteria used
            Eigen::VectorXd phi = fischer(x, s);
            old_err = err;
            err = phi.squaredNorm() / 2;
            // Test relative error
            if (abs(err - old_err) / abs(old_err)
                < Constants::FISCHER_REL_TOL) {
                spdlog::debug("solver=fischer_newton_lcp iter={} "
                              "status=success rel_tol={} rel_err={}",
                    i, Constants::FISCHER_REL_TOL,
                    abs(err - old_err) / abs(old_err));
                return true;
            }
            // Test absolute error
            if (err < Constants::FISCHER_ABS_TOL) {
                spdlog::debug("solver=fischer_newton_lcp iter={} "
                              "status=success abs_tol={} abs_err={}",
                    i, Constants::FISCHER_ABS_TOL, err);
                return true;
            }

            // Bitmask for singular indices
            auto S = phi.array().abs() < Constants::FISCHER_SINGULAR_TOL
                && x.array().abs() < Constants::FISCHER_SINGULAR_TOL;

            // Perturbation: works on full system
            Eigen::VectorXd px = x;
            Eigen::VectorXd dir = x.unaryExpr(&sign);
            // dir(dir==0) = 1;
            dir += (dir.array() == 0).matrix().cast<double>();

            // px(S==1) = Constants::FISCHER_SINGULAR_TOL * dir(S==1);
            for (size_t j = 0; j < S.rows(); j++) {
                if (S(j)) {
                    px(j) = Constants::FISCHER_SINGULAR_TOL * dir(j);
                }
            }

            Eigen::ArrayXd denom
                = (px.array().pow(2) + s.array().pow(2)).sqrt();
            Eigen::VectorXd p = (px.array() / denom) - 1;
            Eigen::VectorXd q = (s.array() / denom) - 1;
            Eigen::SparseMatrix<double> J = Eigen::SparseDiagonal<double>(p)
                + Eigen::SparseDiagonal<double>(q) * A;

            Eigen::VectorXd nabla_phi = J.transpose() * phi;
            // Test if we have dropped into a local minimia if so we are stuck
            if (nabla_phi.norm() < Constants::FISCHER_ABS_TOL) {
                spdlog::warn("solver=ficher_newton_lcp iter={:d} "
                             "failure='||∇ϕ|| < {:g}' ||∇ϕ||={:g} "
                             "failsafe=none",
                    i, Constants::FISCHER_ABS_TOL, nabla_phi.norm());
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
                        "failure='sparse solve for newton direction failed' "
                        "failsafe='revert to gradient descent'",
                        i);
                    delta_x = -nabla_phi; // Revert to gradient descent
                }
            } else {
                spdlog::warn(
                    "solver=ficher_newton_lcp iter={:d} "
                    "failure='sparse decomposition of the hessian failed' "
                    "failsafe='revert to gradient descent'",
                    i);
                delta_x = -nabla_phi; // Revert to gradient descent
            }

            // Test if the search direction is smaller than numerical precision.
            if (delta_x.array().abs().maxCoeff() < EPS) {
                spdlog::warn(
                    "solver=ficher_newton_lcp iter={:d} "
                    "failure='search direction too small' max(|Δx|)={:g} "
                    "failsafe=none",
                    i, delta_x.array().abs().maxCoeff());
                return false;
            }

            // Test if our search direction is a 'sufficient' descent direction
            double descent_magnitude = (nabla_phi.transpose() * delta_x)(0);
            if (descent_magnitude > -EPS * descent_magnitude) {
                spdlog::warn(
                    "solver=ficher_newton_lcp iter={:d} "
                    "failure='search direction not descent direction' "
                    "(∇ϕ)ᵀΔx={:g} failsafe='revert to gradient descent'",
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
                if (alpha * alpha < Constants::FISCHER_SINGULAR_TOL) {
                    spdlog::warn(
                        "solver=ficher_newton_lcp iter={:d} "
                        "failure='step length too small' step_length={:g}",
                        i, alpha);
                    return false;
                }

                alpha /= 2;
            }

            // Update iterate with result from Armijo backtracking
            x = x_k;
        }

        spdlog::warn("solver=ficher_newton_lcp "
                     "failure='too many iterations' MAX_ITER={:d}",
            Constants::FISCHER_MAX_ITER);
        return false;
    }
} // namespace opt
} // namespace ccd
