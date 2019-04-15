/**
 * Functions for optimizing functions.
 * Includes Newton's method with and without constraints.
 */

#pragma once

#include <Eigen/Core>

#include <opt/OptimizationProblem.hpp>
#include <opt/OptimizationResults.hpp>
#include <opt/SolverSettings.hpp>

namespace ccd {

/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    bool line_search(const Eigen::VectorXd& x, const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f, const double fx,
        double& gamma);

    /**
     * @brief Search along a search direction to find a scalar \f$\gamma \in [0,
     * 1]\f$ such that \f$f(x + \gamma \vec{dir}) \leq f(x)\f$.
     *
     * @param[in] x Starting point for the line search.
     * @param[in] dir Direction to search along.
     * @param[in] f Function of x to minimize.
     * @param[in] constraint Constrain on x such that constraint(x) must be
     * true.
     * @return A scalar \f$\gamma\f$ that scales the search direction in an
     *      optimal way.
     */
    double constrained_line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const std::function<bool(const Eigen::VectorXd&)>& constraint);

    /**
     * @brief Perform a single step of Newton's Method to minimize the
     *        objective, \f$f(x)\f$, of the problem unconstrained.
     *
     * @param[in] problem The optimization problem to minimize unconstrained.
     * @param[in] mu      A small value to add to the hessian diagonal to
     *                    prevent it from being singular.
     * @param[in, out] x  The starting value for optimization and the place to
     *                    store the resulting updated value.
     *
     * @return A boolean for if the step was successful.
     */
    bool newtons_method_step(const OptimizationProblem& problem,
        Eigen::VectorXd& x, const double mu = 1e-5);

    /**
     * @brief Perform Newton's Method to minimize the objective, \f$f(x)\f$, of
     * the problem unconstrained.
     *
     * @param[in] problem  The optimization problem to minimize unconstrained.
     * @param[in] settings The settings of the optimization (tolerances and
     *                     maximum number of iterations).
     * @param[in] mu       A small value to add to the hessian diagonal to
     *                     prevent it from being singular.
     *
     * @return The results of the optimization including the minimizer, minimum,
     *         and if the optimization was successful.
     */
    OptimizationResults newtons_method(const OptimizationProblem& problem,
        const SolverSettings& settings, const double mu = 1e-5);

} // namespace opt
} // namespace ccd
