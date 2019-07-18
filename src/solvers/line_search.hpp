/**
 * Functions for optimizing functions.
 * Includes line search to find a step length to reduce a function.
 */

#pragma once

#include <Eigen/Core>

namespace ccd {
namespace opt {

    /**
     * @brief Search along a search direction to find a scalar \f$\step_length
     * \in [0, 1]\f$ such that \f$f(x + \step_length \vec{dir}) \leq f(x)\f$.
     *
     * @param[in] x                Starting point for the line search.
     * @param[in] dir              Direction to search along.
     * @param[in] f                Function of x to minimize.
     * @param[out] step_length     Scalar coefficent of the direction to step.
     * @param[in] min_step_length  Minimum value of step_length before the line
     *                             search fails.
     *
     * @return True if the line search was successful, false otherwise.
     */
    bool line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        double& step_length,
        const double min_step_length = 1e-10);

    /**
     * @brief Search along a search direction to find a scalar \f$\step_length
     * \in [0, 1]\f$ such that \f$f(x + \step_length \vec{dir}) \leq f(x)\f$.
     *
     * @param[in] x                Starting point for the line search.
     * @param[in] dir              Direction to search along.
     * @param[in] f                Function of x to minimize.
     * @param[in] grad_fx          The precomputed value of ∇f(x).
     * @param[out] step_length     Scalar coefficent of the direction to step.
     * @param[in] min_step_length  Minimum value of step_length before the line
     *                             search fails.
     *
     * @return True if the line search was successful, false otherwise.
     */
    bool line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const Eigen::VectorXd& grad_fx,
        double& step_length,
        const double min_step_length = 1e-10,
        const double armijo_rule_coeff = 0);

    /**
     * @brief Search along a search direction to find a scalar \f$\step_length
     * \in [0, 1]\f$ such that \f$f(x + \step_length \vec{dir}) \leq f(x)\f$.
     *
     * @param[in] x                 Starting point for the line search.
     * @param[in] dir               Direction to search along.
     * @param[in] f                 Function of x to minimize.
     * @param[in] grad_fx          The precomputed value of ∇f(x).
     * @param[in] constraint        Constraint on x such that constraint(x) must
     *                              be true.
     * @param[out] step_length      Scalar coefficent of the direction to step.
     * @param[in] min_step_length   Minimum value of step_length before the line
     * search fails.
     *
     * @return True if the line search was successful, false otherwise.
     */
    bool constrained_line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const Eigen::VectorXd& grad_fx,
        const std::function<bool(const Eigen::VectorXd&)>& constraint,
        double& step_length,
        const double min_step_length = 1e-10,
        const double armijo_rule_coeff = 0);

    /**
     * @brief Log values along a search direction.
     *
     * @param[in] x                 Starting point for the line search.
     * @param[in] dir               Direction to search along.
     * @param[in] f                 Function of x to sample.
     * @param[in] grad_f            Gradient of f to sample.
     */
    void sample_search_direction(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f);

} // namespace opt
} // namespace ccd
