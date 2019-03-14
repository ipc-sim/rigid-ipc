/**
 * Methods for optimizing the displacments with a non-linear interference volume
 * constraint.
 */

#pragma once

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>

namespace ccd {

/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    /// @brief Methods available for non-linear constrained optimization.
    enum OptimizationMethod {
        /// @brief Method of Moving Asymptotes (NLopt)
        MMA,
        /// @brief Sequential Least-Squares Quadratic Programming (NLopt)
        SLSQP,
        /// @brief Interior-Point Method (Ipopt)
        IP,
        /// @breif Linearize the constraints and solve the QP (OSQP/Mosek)
        LINEARIZED_CONSTRAINTS,
        /// @breif Nonlinear Complementarity Problem
        NCP
    };

    /**
     * @brief Optimize the displacments with the volume constraint C(U) ≤ 0.
     *
     * @param[in] V Vertices
     * @param[in] U Displacments
     * @param[in] E Edges
     * @param[in] volume_epsilon Epsilon value for volume computation.
     * @param[in] ccd_detection_method Method of detecting collisions.
     * @param[in] opt_method Method of non-linear constrained optimization.
     * @param[in] max_iter Maximum number of iteration to take.
     * @param[in,out] Uopt Initial value for the optimization and storage for
     * the optimal displacments.
     * @return The value of the objective for the minimum found.
     */
    bool solve_problem(const Eigen::MatrixX2d& V, const Eigen::MatrixX2d& U,
        const Eigen::MatrixX2i& E, const double volume_epsilon,
        const DetectionMethod ccd_detection_method,
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

    /**
     * @brief Save JSON file of optimization objectives per iteration.
     *
     * @param[in] method Optimization method used.
     * @param[in] objectives Vector of objective values per iteration.
     * @param[in] objectives Vector of sum of constraint values per iteration.
     */
    void export_intermediate(const OptimizationMethod method,
        const std::vector<double>& objectives,
        const std::vector<double>& constraints);
} // namespace opt

} // namespace ccd
