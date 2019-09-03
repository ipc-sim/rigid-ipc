/**
 * Functions for optimizing functions.
 * Includes Newton's method.
 */

#pragma once

#include <Eigen/Core>

#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>
#include <solvers/optimization_solver.hpp>

namespace ccd {

/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    class NewtonSolver : public virtual IBarrierOptimizationSolver {
    public:
        NewtonSolver();
        NewtonSolver(const std::string& name);
        virtual ~NewtonSolver() override;

        // From IOptimizationSolver
        // -------------------------------------------------
        /**
         * @brief Perform Newton's Method to minimize the objective, \f$f(x)\f$,
         * of the problem unconstrained.
         *
         * @param[in] problem  The optimization problem to minimize
         *                     unconstrained.
         *
         * @return The results of the optimization including the minimizer,
         * minimum, and if the optimization was successful.
         */
        virtual OptimizationResults solve(IBarrierProblem& problem) override;

        // From IBarrierOptimizationSolver
        // -------------------------------------------------
        const std::string& name() const override { return name_; }
        void settings(const nlohmann::json& /*json*/) override;
        nlohmann::json settings() const override;
        void init_free_dof(Eigen::VectorXb is_dof_fixed) override;

        // -------------------------------------------------
        /**
         * @brief Solve for the Newton direction
         *        (\f$\Delta x = -H^{-1} \nabla f \f$).
         *
         * @param[in]  gradient  Gradient of the objective function.
         * @param[in]  hessian   Hessian of the objective function.
         * @param[out] delta_x   Output newton direction.
         * @param[in]  make_psd  If delta_x is not a descent direction, then
         * make the hessian positive semi-definite.
         *
         * @return Returns true if the solve was successful.
         */
        bool compute_direction(const Eigen::VectorXd& gradient,
            const Eigen::SparseMatrix<double>& hessian,
            Eigen::VectorXd& delta_x,
            bool make_psd = false);

        void c(const double value) override { c_ = value; }
        void e_b(const double value) override { e_b_ = value; }
        void t(const double value) override { t_ = value; }
        void m(const double value) override { m_ = value; }

        int max_iterations;

        double c_;
        double e_b_;
        double t_;
        double m_;

    protected:
        bool line_search(IBarrierProblem& problem,
            const Eigen::VectorXd& x,
            const Eigen::VectorXd& dir,
            const double fx,
            const Eigen::VectorXd& grad_fx,
            double& step_length,
            bool log_failure = false);

        Eigen::VectorXi free_dof; ///< @breif Indices of the free degrees.
        int iteration_number;     ///< @brief The current iteration number.

        std::string name_;
    };

    /**
     * @brief Make the matrix positive definite (\f$x^T A x > 0\$).
     *
     * @param A The matrix to make positive definite.
     *
     * @return The scale of the update to the diagonal.
     */
    double make_matrix_positive_definite(Eigen::SparseMatrix<double>& A);

} // namespace opt
} // namespace ccd
