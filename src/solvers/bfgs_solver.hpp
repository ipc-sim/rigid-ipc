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

    class BFGSSolver : public OptimizationSolver {
    public:
        BFGSSolver();
        BFGSSolver(const std::string& name);
        virtual ~BFGSSolver() override;

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
        virtual OptimizationResults solve(
            OptimizationProblem& problem) override;

        void settings(const nlohmann::json& /*json*/) override;
        nlohmann::json settings() const override;

        double absolute_tolerance; ///< @brief Convergence tolerance.
        double min_step_length;    ///< @brief Minimum step length.
    };

} // namespace opt
} // namespace ccd
