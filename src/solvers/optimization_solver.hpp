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

namespace ccd {
namespace opt {

    class OptimizationSolver {
    protected:
        Eigen::VectorXi free_dof; ///< @breif Indices of the free degrees.

    public:
        OptimizationSolver();
        OptimizationSolver(const int max_iterations);
        virtual ~OptimizationSolver();

        virtual void init_free_dof(Eigen::MatrixXb is_dof_fixed);

        virtual OptimizationResults solve(OptimizationProblem& problem) = 0;

        int max_iterations; ///< @brief Maximum number of iteration.
    };

} // namespace opt
} // namespace ccd
