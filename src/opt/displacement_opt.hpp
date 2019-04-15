/**
 * Methods for optimizing the displacments with a non-linear interference volume
 * constraint.
 */

#pragma once

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>
#include <opt/solver.hpp>

namespace ccd {

/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    /**
     * @brief Runs the optimization problem for the given method and initial
     * value
     * @param[in,out] problem : Optimization problem (possibly modified for
     *                          validation)
     * @param[in] U0          : Initial displacements
     * @param[in] settings    : Solver settings including tolerances and method
     * @return Optimization Result of the optimization
     * */
    OptimizationResults displacement_optimization(OptimizationProblem& problem,
        const Eigen::MatrixX2d& U0, SolverSettings& settings);

    /**
     * @brief Solves the KKT conditions of the Optimization Problem
     *      (U - Uk) = \nabla g(U)
     *      s.t V(U) >= 0
     */
    OptimizationResults solve_ncp_displacement_optimization(
        OptimizationProblem& problem, SolverSettings& settings);

} // namespace opt
} // namespace ccd
