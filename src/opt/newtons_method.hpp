/**
 * Functions for optimizing functions.
 * Includes Newton's method with and without constraints.
 */

#pragma once

#include <Eigen/Core>

#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>

namespace ccd {

/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    /**
     * @brief Perform Newton's Method to minimize the objective, \f$f(x)\f$, of
     * the problem unconstrained.
     *
     * @param[in] problem  The optimization problem to minimize unconstrained.
     * @param[in] mu       A small value to add to the hessian diagonal to
     *                     prevent it from being singular.
     *
     * @return The results of the optimization including the minimizer, minimum,
     *         and if the optimization was successful.
     */

    OptimizationResults newtons_method(OptimizationProblem& problem,
        const double absolute_tolerance, const double line_search_tolerance,
        const int max_iter, const double mu = 1e-5);

    /**
     * @brief Search along a search direction to find a scalar \f$\gamma \in [0,
     * 1]\f$ such that \f$f(x + \gamma \vec{dir}) \leq f(x)\f$.
     *
     * @param[in] x Starting point for the line search
     * @param[in] dir Direction to search along
     * @param[in] f Function of x to minimize
     * @param[in] fx The precomputed value of f(x)
     * @param[out] gamma Scalar coefficent of the direction to step
     * @param[in] min_gamma Minimum value of gamma before the line search fails
     * @return True if the line search was successful, false otherwise.
     */
    bool line_search(const Eigen::VectorXd& x, const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f, double& gamma,
        const double min_gamma = 1e-10);

    /**
     * @brief Search along a search direction to find a scalar \f$\gamma \in [0,
     * 1]\f$ such that \f$f(x + \gamma \vec{dir}) \leq f(x)\f$.
     *
     * @param[in] x Starting point for the line search
     * @param[in] dir Direction to search along
     * @param[in] f Function of x to minimize
     * @param[in] constraint Constrain on x such that constraint(x) must be
     *                       true
     * @param[out] gamma Scalar coefficent of the direction to step
     * @param[in] min_gamma Minimum value of gamma before the line search fails
     * @return True if the line search was successful, false otherwise.
     */
    bool constrained_line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const std::function<bool(const Eigen::VectorXd&)>& constraint,
        double& gamma, const double min_gamma = 1e-10);

} // namespace opt
} // namespace ccd
