/**
 * Functions for optimizing functions.
 * Includes line search to find a step length to reduce a function.
 */

#pragma once

#include <Eigen/Core>

namespace ccd {
namespace opt {

    /**
     * @brief Search along a search direction to find a scalar \f$\alpha
     * \in [0, 1]\f$ such that \f$f(x + \alpha \Delta x) \leq f(x)\f$.
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
    bool line_search(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        double& step_length,
        const double min_step_length = 1e-10);

    /**
     * @brief Search along a search direction to find a scalar \f$\alpha
     * \in [0, 1]\f$ such that \f$f(x + \alpha \Delta x) < f(x)\f$.
     *
     * @param[in] x                Starting point for the line search.
     * @param[in] dir              Direction to search along.
     * @param[in] f                Function of x to minimize.
     * @param[in] grad_fx          The precomputed value of \f$\nabla f(x)\f$.
     * @param[out] step_length     Scalar coefficent of the direction to step.
     * @param[in] min_step_length  Minimum value of step_length before the line
     *                             search fails.
     *
     * @return True if the line search was successful, false otherwise.
     */
    bool line_search(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const Eigen::VectorXd& grad_fx,
        double& step_length,
        const double min_step_length = 1e-10,
        const double armijo_rule_coeff = 0);

    /**
     * @brief Search along a search direction to find a scalar \f$\alpha
     * \in [0, 1]\f$ such that \f$f(x + \alpha \Delta x) < f(x)\f$.
     *
     * @param[in] x                Starting point for the line search.
     * @param[in] dir              Direction to search along \f$(\Delta x)\f$.
     * @param[in] f                Function of x to minimize.
     * @param[in] grad_fx          The precomputed value of \f$\nabla f(x)\f$.
     * @param[in] constraint       Constraint on x such that constraint(x) must
     *                             be true.
     * @param[out] step_length     Scalar coefficent of the direction to step
     *                             \f$(\alpha)\f$.
     * @param[in] min_step_length  Minimum value of step_length before the line
     *                             search fails.
     *
     * @return True if the line search was successful, false otherwise.
     */
    bool constrained_line_search(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const Eigen::VectorXd& grad_fx,
        const std::function<bool(const Eigen::VectorXd&)>& constraint,
        double& step_length,
        const double min_step_length = 1e-10,
        const double armijo_rule_coeff = 0);

} // namespace opt
} // namespace ccd
