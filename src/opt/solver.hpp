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

    class OptimizationSolver {
    public:
        OptimizationSolver();
        virtual ~OptimizationSolver();
        virtual OptimizationResults solve(OptimizationProblem& problem) = 0;
    };

} // namespace opt
} // namespace ccd
