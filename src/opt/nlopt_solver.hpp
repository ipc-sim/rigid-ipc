/// Solve a optimization problem with NLopt.
#pragma once

#include <Eigen/Core>
#include <nlopt.hpp>

#include <ccd/collision_detection.hpp>

namespace ccd {
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
     * @return The value of the objective for the minimum found.
     */
    bool solve_problem_with_nlopt(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const DetectionMethod ccd_detection_method,
        const nlopt::algorithm opt_method, const unsigned max_iter,
        Eigen::MatrixX2d& Uopt);

    /**
     * @brief Output the reason for NLopt's termination.
     *
     * @param[in] opt Optimizer object.
     * @param[in] reuslt Optimization result to parse.
     */
    void print_nlopt_termination_reason(
        const nlopt::opt& opt, const nlopt::result result);

} // namespace opt
} // namespace ccd
