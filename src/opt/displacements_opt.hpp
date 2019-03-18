/**
 * Methods for optimizing the displacments with a non-linear interference volume
 * constraint.
 */

#pragma once

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>
#include <opt/minimize.hpp>

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
     * @param[in] V     : Vertices
     * @param[in] U     : Displacments
     * @param[in] E     : Edges
     * @param[in] volume_epsilon    : Epsilon value for volume computation.
     * @param[in] detection_method  : Method of detecting collisions.
     *
     * @param[out] problem          : The optimization problem.
     */
    void setup_displacement_optimization(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const DetectionMethod ccd_detection_method,
        OptimizationProblem& problem);

    /**
     * @brief Runs the optimization problem for the given method and initial
     * value
     * @param[in] opt_method:   Optimization method
     * @param[in] U0:           Initial displacements
     * @param[in,out] problem:  Optimization problem (possibly modified for
     * validation)
     * @param[out] u_history:   Optimization steps for displacements
     * @param[out] f_history:   Optimization steps for functional
     * @param[out] g_history:   Optimization steps for constraints
     *
     * @return Optimization Result of the optimization
     * */
    OptimizationResult displacements_optimization(
        const OptimizationMethod& opt_method, const Eigen::MatrixX2d& U0,
        OptimizationProblem& problem, std::vector<Eigen::MatrixX2d>& u_history,
        std::vector<double>& f_history, std::vector<double>& g_history);

    /**
     * @brief Optimize the displacments with the volume constraint C(U) ≤ 0.
     *
     * @param[in] V Vertices
     * @param[in] U Displacments
     * @param[in] E Edges
     * @param[in] volume_epsilon Epsilon value for volume computation.
     * @param[in] detection_method Method of detecting collisions.
     * @param[in] opt_method Method of non-linear constrained optimization.
     * @param[in] max_iter Maximum number of iteration to take.
     * @param[in,out] Uopt Initial value for the optimization and storage for
     * the optimal displacments.
     * @return The value of the objective for the minimum found.
     */
    double displacements_optimization(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const DetectionMethod ccd_detection_method,
        const OptimizationMethod opt_method, const unsigned max_iter,
        Eigen::MatrixX2d& Uopt);

    /**
     * @brief Run the entire CCD pipeline.
     *
     * @param[in] V Vertices
     * @param[in] U Displacments
     * @param[in] E Edges
     * @param[in] detection_method Method of detecting collisions.
     * @param[out] ee_impacts Vector of all edge-edge impacts.
     * @param[out] edge_impact_map Earliest impact per edge.
     */
    void detect_collisions(const Eigen::MatrixX2d& V, const Eigen::MatrixX2d& U,
        const Eigen::MatrixX2i& E, const ccd::DetectionMethod detection_method,
        EdgeEdgeImpacts& ee_impacts, Eigen::VectorXi& edge_impact_map);

    void export_intermediate(const OptimizationMethod method,
        const std::vector<double>& objectives,
        const std::vector<double>& constraints);
} // namespace opt

} // namespace ccd
