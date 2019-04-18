#pragma once

#include <Eigen/Core>

#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>
#include <opt/solver_settings.hpp>

namespace ccd {
/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

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
