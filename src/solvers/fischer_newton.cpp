#include "lcp_solver.hpp"

#include <Eigen/Core>

#include <constants.hpp>
#include <logger.hpp>
#include <utils/eigen_ext.hpp>

namespace ccd {
namespace opt {

    inline Eigen::VectorXd fischer(
        const Eigen::VectorXd& x, const Eigen::VectorXd& s);
    inline double fischer_error(const Eigen::VectorXd& phi);
    Eigen::SparseMatrix<double> jac_fischer(const Eigen::VectorXd& x,
        const Eigen::VectorXd& s,
        const Eigen::SparseMatrix<double>& A);
    inline Eigen::VectorXd grad_fischer_error(
        const Eigen::VectorXd& phi, const Eigen::SparseMatrix<double>& J);
    Eigen::VectorXd perturb(
        const Eigen::VectorXd& x, const Eigen::VectorXd& phi);
    bool is_converged(double old_err, double err, int iter);
    Eigen::VectorXd compute_search_direction(const Eigen::VectorXd& phi,
        const Eigen::SparseMatrix<double>& J,
        const Eigen::VectorXd& grad_psi,
        const int iter);
    // Armijo backtracking combined with a projected line-search
    bool projected_line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const double fx,
        const Eigen::VectorXd& grad_fx,
        double& step_length,
        int iter,
        const double armijo_rule_coeff = 0.001);

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
        // Use a sparse matrix
        const Eigen::SparseMatrix<double> A = A_dense.sparseView();

        int num_vars = b.rows();
        assert(A.rows() == num_vars);

        // x₀ = 0
        x = Eigen::VectorXd::Zero(num_vars);

        double old_err, err = std::numeric_limits<double>::infinity();
        for (int i = 0; i < Constants::FISCHER_MAX_ITER; i++) {
            Eigen::VectorXd s = A * x + b; // slack variable

            // Test all stopping criteria used
            Eigen::VectorXd phi = fischer(x, s);
            old_err = err;
            err = fischer_error(phi); // ψ

            if (is_converged(old_err, err, i)) {
                return true;
            }

            Eigen::VectorXd px = perturb(x, phi);
            Eigen::SparseMatrix<double> J = jac_fischer(px, s, A);
            Eigen::VectorXd grad_psi = grad_fischer_error(phi, J);

            // Test if we have dropped into a local minimia if so we are stuck
            if (grad_psi.norm() < Constants::FISCHER_ABS_TOL) {
                spdlog::warn("solver=ficher_newton_lcp iter={:d} "
                             "failure='||∇ψ|| < {:g}' ||∇ψ||={:g} "
                             "failsafe=none",
                    i, Constants::FISCHER_ABS_TOL, grad_psi.norm());
                return false;
            }

            Eigen::VectorXd delta_x
                = compute_search_direction(phi, J, grad_psi, i);

            // Update iterate with result from Armijo backtracking
            double step_length = 1.0;
            bool success = projected_line_search(x, delta_x,
                [&A, &b](const Eigen::VectorXd& x) {
                    Eigen::VectorXd s = A * x + b;
                    Eigen::VectorXd phi = fischer(x, s);
                    return fischer_error(phi);
                },
                /*fx=*/err, grad_psi, step_length, i);
            if (!success) {
                delta_x = -grad_psi;
                success = projected_line_search(x, delta_x,
                    [&A, &b](const Eigen::VectorXd& x) {
                        Eigen::VectorXd s = A * x + b;
                        Eigen::VectorXd phi = fischer(x, s);
                        return fischer_error(phi);
                    },
                    /*fx=*/err, grad_psi, step_length, i);
                if (!success) {
                    return false;
                }
            }

            x += step_length * delta_x;
        }

        spdlog::warn("solver=ficher_newton_lcp "
                     "failure='too many iterations' MAX_ITER={:d}",
            Constants::FISCHER_MAX_ITER);
        return false;
    }

    /// @brief Compute the Fischer error measure.
    inline Eigen::VectorXd fischer(
        const Eigen::VectorXd& x, const Eigen::VectorXd& s)
    {
        //             _________
        // ϕ(x, s) = ⎷(x² + s²) - x - s
        //
        return (x.array().pow(2) + s.array().pow(2)).sqrt().matrix() - x - s;
    }

    inline double fischer_error(const Eigen::VectorXd& phi)
    {
        // ψ(x, s) = ½ϕ(x, s)ᵀϕ(x, s)
        return phi.squaredNorm() / 2;
    }

    Eigen::SparseMatrix<double> jac_fischer(const Eigen::VectorXd& x,
        const Eigen::VectorXd& s,
        const Eigen::SparseMatrix<double>& A)
    {
        //                ___________             _x̲ᵢ̲_+̲_s̲ᵢ̲a̲ᵢ̲_
        // ∇ϕᵢ(x, s) = ∇⎷(xᵢ² + sᵢ²) - xᵢ - sᵢ =  xᵢ² + s²    - aᵢ - 1
        Eigen::ArrayXd denom = (x.array().pow(2) + s.array().pow(2)).sqrt();
        Eigen::VectorXd p = (x.array() / denom) - 1;
        Eigen::VectorXd q = (s.array() / denom) - 1;
        return Eigen::SparseDiagonal<double>(p)
            + Eigen::SparseDiagonal<double>(q) * A;
    }

    inline Eigen::VectorXd grad_fischer_error(
        const Eigen::VectorXd& phi, const Eigen::SparseMatrix<double>& J)
    {
        // ∇ψ(x, s) = ∇ϕ(x, s)ᵀϕ(x, s) = Jᵀϕ
        return J.transpose() * phi;
    }

    /// @brief Compute the sign of the input.
    inline double sign(double x) { return (0 < x) - (x < 0); }

    Eigen::VectorXd perturb(
        const Eigen::VectorXd& x, const Eigen::VectorXd& phi)
    {
        // Bitmask for singular indices
        Eigen::ArrayXb is_singular
            = phi.array().abs() < Constants::FISCHER_SINGULAR_TOL
            && x.array().abs() < Constants::FISCHER_SINGULAR_TOL;

        // Perturbation: works on full system
        Eigen::VectorXd px = x; // perturbed x
        Eigen::VectorXd x_sign = x.unaryExpr(&sign);
        // x_sign(x_sign==0) = 1;
        x_sign += (x_sign.array() == 0).matrix().cast<double>();
        // px(S==1) = Constants::FISCHER_SINGULAR_TOL * dir(S==1);
        for (size_t i = 0; i < x.rows(); i++) {
            if (is_singular(i)) {
                px(i) = Constants::FISCHER_SINGULAR_TOL * x_sign(i);
            }
        }

        return px;
    }

    bool is_converged(double old_err, double err, int iter)
    {
        // Test relative error
        if (abs(err - old_err) / abs(old_err) < Constants::FISCHER_REL_TOL) {
            spdlog::debug("solver=fischer_newton_lcp iter={} "
                          "status=success rel_tol={} rel_err={}",
                iter, Constants::FISCHER_REL_TOL,
                abs(err - old_err) / abs(old_err));
            return true;
        }
        // Test absolute error
        if (err < Constants::FISCHER_ABS_TOL) {
            spdlog::debug("solver=fischer_newton_lcp iter={} "
                          "status=success abs_tol={} abs_err={}",
                iter, Constants::FISCHER_ABS_TOL, err);
            return true;
        }
        return false;
    }

    Eigen::VectorXd compute_search_direction(const Eigen::VectorXd& phi,
        const Eigen::SparseMatrix<double>& J,
        const Eigen::VectorXd& grad_psi,
        const int iter)
    {
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
                    iter);
                delta_x = -grad_psi; // Revert to gradient descent
            }
        } else {
            spdlog::warn("solver=ficher_newton_lcp iter={:d} "
                         "failure='sparse decomposition of the hessian failed' "
                         "failsafe='revert to gradient descent'",
                iter);
            delta_x = -grad_psi; // Revert to gradient descent
        }

        return delta_x;
    }

    // Armijo backtracking combined with a projected line-search
    bool projected_line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const double fx,
        const Eigen::VectorXd& grad_fx,
        double& step_length,
        int iter,
        const double armijo_rule_coeff)
    {
        const double EPS = std::numeric_limits<double>::epsilon();

        // Test if the search direction is smaller than numerical precision.
        if (dir.array().abs().maxCoeff() < EPS) {
            spdlog::warn("solver=ficher_newton_lcp iter={:d} "
                         "failure='search direction too small' max(|Δx|)={:g} "
                         "failsafe=none",
                iter, dir.array().abs().maxCoeff());
            return false;
        }

        // Test if our search direction is a 'sufficient' descent direction
        double descent_magnitude = (grad_fx.transpose() * dir)(0);
        if (descent_magnitude > -EPS * descent_magnitude) {
            spdlog::warn("solver=ficher_newton_lcp iter={:d} "
                         "failure='search direction not descent direction' "
                         "(∇ϕ)ᵀΔx={:g} failsafe='revert to gradient descent'",
                iter, descent_magnitude);
            return false;
        }

        step_length = 1.0; // Current step length

        // Sufficent decrease parameter for Armijo backtracking line search
        double armijo_term = armijo_rule_coeff * descent_magnitude;

        Eigen::VectorXd x_k;
        do {
            x_k = (x + step_length * dir).array().max(0); // project it

            // Perform Armijo codition to see if we got a sufficient decrease
            if (f(x_k) <= fx + step_length * armijo_term) {
                return true;
            }

            step_length /= 2;
        } while ((step_length * dir).norm() >= Constants::FISCHER_SINGULAR_TOL);

        spdlog::warn("solver=ficher_newton_lcp iter={:d} "
                     "failure='step length too small' step_length={:g}",
            iter, step_length);
        return false;
    }

} // namespace opt
} // namespace ccd
