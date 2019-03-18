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
     * @brief Creates an Optimization problem for displacements C(U) ≤ 0.
     *      Min ||U - Uk||^2
     *      s.t V(U) >= 0
     * @param[in] V                : Vertices
     * @param[in] U                : Displacments
     * @param[in] E                : Edges
     * @param[in] volume_epsilon   : Epsilon value for volume computation.
     * @param[in] detection_method : Method of detecting collisions.
     * @param[out] problem         : The optimization problem.
     */
    void setup_displacement_optimization_problem(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const DetectionMethod ccd_detection_method,
        OptimizationProblem& problem);

    /**
     * @brief Run the entire CCD pipeline.
     *
     * @param[in] V                : Vertices
     * @param[in] U                : Displacments
     * @param[in] E                : Edges
     * @param[in] detection_method : Method of detecting collisions.
     * @param[out] ee_impacts      : Vector of all edge-edge impacts.
     * @param[out] edge_impact_map : Earliest impact per edge.
     */
    void detect_collisions(const Eigen::MatrixX2d& V, const Eigen::MatrixX2d& U,
        const Eigen::MatrixX2i& E, const ccd::DetectionMethod detection_method,
        EdgeEdgeImpacts& ee_impacts, Eigen::VectorXi& edge_impact_map);

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
        const Eigen::MatrixX2d& U0, std::vector<Eigen::MatrixX2d>& u_history,
        std::vector<double>& f_history, std::vector<double>& g_history,
        SolverSettings& settings);

    /**
     * @brief Save JSON file of optimization objectives per iteration.
     *
     * @param[in] method      : Optimization method used.
     * @param[in] objectives  : Vector of objective values per iteration.
     * @param[in] constraints : Vector of sum of constraint values per
     *                          iteration.
     */
    void export_intermediate(const OptimizationMethod method,
        const std::vector<double>& objectives,
        const std::vector<double>& constraints);

} // namespace opt
} // namespace ccd
