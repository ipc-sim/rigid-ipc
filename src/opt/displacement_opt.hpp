/**
 * Methods for optimizing the displacments with a non-linear interference volume
 * constraint.
 */

#pragma once

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>
#include <opt/solver.hpp>
#include <opt/ncp_solver.hpp>

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

    class NCPDisplacementOptimization {
    public:
        NCPDisplacementOptimization();
        OptimizationResults solve(OptimizationProblem& problem);

        // ---------------
        // Settings
        // ---------------
        int max_iterations;
        NcpUpdate update_method;
        LCPSolver lcp_solver;
        bool keep_in_unfeasible;
        bool check_convergence;
        double convegence_tolerance;

    };
} // namespace opt
} // namespace ccd
