#pragma once

#include <Eigen/Core>

#include <opt/OptimizationProblem.hpp>
#include <opt/OptimizationResults.hpp>
#include <opt/SolverSettings.hpp>

namespace ccd {
/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    enum QPSolver {
        OSQP, ///< @brief Use OSQP to solve the qudratic program.
        MOSEK ///< @brief Use Mosek to solve the qudratic program.
    };

    /**
     * @brief Optimize the displacments using linearized constraints
     *
     * @param[in,out] problem Optimization problem to solve
     * @return The results of the optimization
     */
    OptimizationResults solve_problem_with_linearized_constraints(
        const OptimizationProblem& problem, const SolverSettings& settings);

} // namespace opt

} // namespace ccd
