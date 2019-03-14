#pragma once

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>

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
     * @brief Optimize the displacments using a linearized constraint.
     *
     * @param[in] V Vertices
     * @param[in] U Displacments
     * @param[in] E Edges
     * @param[in] VOLUME_EPSILON Epsilon value for volume computation.
     * @param[in] DETECTION_METHOD Method of detecting collisions.
     * @param[in] opt_method Method of non-linear constrained optimization.
     * @param[in] max_iter Maximum number of iteration to take.
     * @param[in,out] Uopt Initial value for the optimization and storage for
     * the optimal displacments.
     * @return The value of the objective for the minimum found.
     */
    bool solve_problem_with_linearized_constraints(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double VOLUME_EPSILON, const DetectionMethod CCD_DETECTION_METHOD,
        const int max_iter, Eigen::MatrixX2d& Uopt);

} // namespace opt

} // namespace ccd
