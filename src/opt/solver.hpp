/**
 *  We handle optimization problems of the form
 *      MIN     f(x)      x ∈ Rⁿ
 *
 *   s.t.       g_L ≤ g(x) ≤ g_U
 *              x_L ≤  x   ≤ x_U
 *
 */
#pragma once

#include <opt/OptimizationProblem.hpp>
#include <opt/OptimizationResults.hpp>
#include <opt/SolverSettings.hpp>

namespace ccd {
namespace opt {

    /**
     * @brief Solve an optimization problem
     * We handle optimization problems of the form
     * \begin{align}
     *      &\min f(x) & x \in \mathbb{R}^n \\\\
     *      &\text{subject to} & g_L \leq g(x) \leq g_U \\\\
     *      &                   & x_L \leq x \leq x_U \\\\
     * \end{align}
     */
    OptimizationResults solve_problem(
        OptimizationProblem& problem, const SolverSettings& settings);

} // namespace opt
} // namespace ccd
