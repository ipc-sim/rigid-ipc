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
        const Eigen::VectorXi& free_dof, const double absolute_tolerance,
        const double line_search_tolerance, const int max_iter);

    // Solve for the Newton direction (Δx = -H^{-1}∇f).
    // Return true if the solve was successful.
    bool solve_for_newton_direction(const Eigen::VectorXd& gradient,
        const Eigen::SparseMatrix<double>& hessian, Eigen::VectorXd& delta_x,
        bool make_psd = false);

    // Solve for the newton direction for the limited degrees of freedom.
    // The fixed dof of x will have a delta_x of zero.
    bool solve_for_free_newton_direction(const Eigen::VectorXd& gradient_free,
        const Eigen::SparseMatrix<double>& hessian_free,
        const Eigen::VectorXi& free_dof, Eigen::VectorXd& delta_x,
        bool make_psd = false);

    // Make the matrix positive definite (x^T A x > 0).
    double make_matrix_positive_definite(Eigen::SparseMatrix<double>& A);

} // namespace opt
} // namespace ccd
