#pragma once

#include <Eigen/Core>
#include <nlopt.hpp>

#include <ccd/collision_detection.hpp>

namespace ccd {
/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    /**
     * @brief Optimize the displacments using NLopt.
     *
     * @param[in] V Matrix of vertex positions.
     * @param[in] U Matrix of initial displacments.
     * @param[in] E Matrix of edge pairs as vertex indices.
     * @param[in] volume_epsilon Epsilon for volume computations.
     * @param[in] ccd_detection_method Method of detecting collisions.
     * @param[in] opt_method NLopt algorithm to use for optimization.
     * @param[in] max_iter Maximum number of optimization iterations.
     * @param[in,out] Uopt Inital value for optimization and optimal
     *                     displacments to return at termination.
     */
    void displacements_optimization_nlopt(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const DetectionMethod ccd_detection_method,
        const nlopt::algorithm opt_method, const unsigned max_iter,
        Eigen::MatrixX2d& Uopt);

}

}
