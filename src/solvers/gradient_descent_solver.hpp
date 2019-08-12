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

    class GradientDescentSolver : public IBarrierOptimizationSolver {
    public:
        GradientDescentSolver();
        GradientDescentSolver(const std::string& name);

        virtual ~GradientDescentSolver() override;

        // From IBarrierOptimizationSolver
        // --------------------------------
        const std::string& name() const override { return name_; }
        void settings(const nlohmann::json& /*json*/) override;
        nlohmann::json settings() const override;
        void init_free_dof(Eigen::VectorXb is_dof_fixed) override;

        // From IOptimizationSolver
        // --------------------------------
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

        Eigen::VectorXi free_dof; ///< @breif Indices of the free degrees.
        std::string name_;
    };

} // namespace opt
} // namespace ccd
