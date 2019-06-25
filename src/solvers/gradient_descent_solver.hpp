/// Functions for optimizing functions.

#pragma once

#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>
#include <solvers/newton_solver.hpp>

namespace ccd {

/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    class GradientDescentSolver : public NewtonSolver {
    public:
        GradientDescentSolver();
        virtual ~GradientDescentSolver() override;

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
        virtual OptimizationResults solve(
            OptimizationProblem& problem) override;

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
    };

} // namespace opt
} // namespace ccd
