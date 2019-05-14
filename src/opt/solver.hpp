/**
 *  We handle optimization problems of the form
 *      MIN     f(x)      x ∈ Rⁿ
 *
 *   s.t.       g_L ≤ g(x) ≤ g_U
 *              x_L ≤  x   ≤ x_U
 *
 */
#pragma once

#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>
#include <opt/solver_settings.hpp>

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
        OptimizationProblem& problem, SolverSettings& settings);

    class OptimizationSolver {
    public:
        OptimizationSolver();
        virtual ~OptimizationSolver();
        virtual OptimizationResults solve(OptimizationProblem& problem) = 0;
    };

} // namespace opt
} // namespace ccd
