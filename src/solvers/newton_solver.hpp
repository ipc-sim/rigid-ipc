#pragma once

#include <Eigen/Core>

#include <constants.hpp>
#include <solvers/optimization_solver.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    class NewtonSolver : public virtual OptimizationSolver {
    public:
        NewtonSolver();
        virtual ~NewtonSolver() = default;

        /// Initialize the state of the solver using the settings saved in JSON
        virtual void settings(const nlohmann::json& params) override;
        /// Export the state of the solver using the settings saved in JSON
        virtual nlohmann::json settings() const override;

        /// An identifier for the solver class
        static std::string solver_name() { return "newton_solver"; }
        /// An identifier for this solver
        virtual std::string name() const override
        {
            return NewtonSolver::solver_name();
        }

        /// Initialize the solver with a problem to solve
        virtual void set_problem(OptimizationProblem& problem) override
        {
            this->problem_ptr = &problem;
        }

        /// Initialize the solver state for a new solve
        virtual void init_solve(const Eigen::VectorXd& x0) override;

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
        virtual OptimizationResults solve(const Eigen::VectorXd& x0) override;

        /// Perform a single step of solving the optimization problem
        virtual OptimizationResults step_solve() override
        {
            throw NotImplementedError(
                "Taking a single newton step is not implemented!");
        };

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
        virtual bool compute_direction(
            const Eigen::VectorXd& gradient,
            const Eigen::SparseMatrix<double>& hessian,
            Eigen::VectorXd& delta_x,
            bool make_psd = false);

        virtual std::string stats() override;

        int max_iterations;

    protected:
        virtual bool
        converged(const Eigen::VectorXd& grad, const Eigen::VectorXd& dir) const
        {
            return abs(dir.dot(grad)) <= Constants::NEWTON_ABSOLUTE_TOLERANCE;
        }

        virtual void
        post_step_update(const Eigen::VectorXd& xi, const Eigen::VectorXd& xj)
        {
            // Maybe a child might want to do something here
        }

        virtual bool line_search(
            const Eigen::VectorXd& x,
            const Eigen::VectorXd& dir,
            const double fx,
            const Eigen::VectorXd& grad_fx,
            double& step_length);

        virtual double line_search_lower_bound() const
        {
            return Constants::LINE_SEARCH_LOWER_BOUND;
        }

        /// @brief Pointer to the problem to solve.
        OptimizationProblem* problem_ptr;

        Eigen::VectorXi free_dof; ///< @brief Indices of the free degrees.
        int iteration_number;     ///< @brief The current iteration number.

    private:
        void reset_stats();

        int num_fx = 0;
        int num_grad_fx = 0;
        int num_hessian_fx = 0;
        int num_collision_check = 0;
        int ls_iterations = 0;
        int newton_iterations = 0;
    };

    /**
     * @brief Make the matrix positive definite (\f$x^T A x > 0\$).
     *
     * @param A The matrix to make positive definite.
     *
     * @return The scale of the update to the diagonal.
     */
    double make_matrix_positive_definite(Eigen::SparseMatrix<double>& A);

    /// Construct indicies of free DoF
    Eigen::VectorXi init_free_dof(const Eigen::VectorXb& is_dof_fixed);

} // namespace opt
} // namespace ccd
