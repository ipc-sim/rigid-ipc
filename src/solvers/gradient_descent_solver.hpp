#pragma once

#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>
#include <solvers/barrier_solver.hpp>

namespace ccd {

/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    class GradientDescentSolver : public virtual BarrierInnerSolver {
    public:
        GradientDescentSolver();
        virtual ~GradientDescentSolver() = default;

        /// Initialize the state of the solver using the settings saved in JSON
        void settings(const nlohmann::json& params) override;
        /// Export the state of the solver using the settings saved in JSON
        nlohmann::json settings() const override;

        static std::string solver_name() { return "gradient_descent_solver"; }
        /// An identifier for this solver
        virtual std::string name() const override
        {
            return GradientDescentSolver::solver_name();
        }

        /// Initialize the solver with a problem to solve
        virtual void set_problem(OptimizationProblem& problem) override
        {
            this->problem_ptr = &problem;
        }

        /// Initialize the solver state for a new solve
        virtual void init_solve() override;

        /**
         * @brief Perform gradient descent to minimize the objective,
         * \f$f(x)\f$, of the problem unconstrained.
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

        double absolute_tolerance; ///< @brief Convergence tolerance.
        double min_step_length;    ///< @brief Minimum step length.
        int max_iterations;

    protected:
        /**
         * @brief Compute the direction for the limited degrees of freedom.
         *
         * The fixed dof of x will have a delta_x of zero.
         *
         * @param[in]  gradient_free  Gradient of the free dof of the objective
         *                            function.
         * @param[out] delta_x        Output full dof newton direction.
         *
         * @return Returns true if the computation was successful.
         */
        bool compute_free_direction(
            const Eigen::VectorXd& gradient_free, Eigen::VectorXd& delta_x);

        /**
         * @brief Solve for the direction (\f$\Delta x = -\nabla f \f$).
         *
         * @param[in]  gradient  Gradient of the objective function.
         * @param[out] delta_x   Output newton direction.
         *
         * @return Returns true if the computation was successful.
         */
        static bool compute_direction(
            const Eigen::VectorXd& gradient, Eigen::VectorXd& delta_x);

        bool line_search(
            const Eigen::VectorXd& x,
            const Eigen::VectorXd& dir,
            const double fx,
            double& step_length);

        OptimizationProblem* problem_ptr;
        Eigen::VectorXi free_dof; ///< @breif Indices of the free degrees.
    };

} // namespace opt
} // namespace ccd
