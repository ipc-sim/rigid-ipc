/**
 * @brief Solve the optimization problem using Newton's Method with barriers for
 * the constraints.
 */

#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>
#include <opt/solver_settings.hpp>

namespace ccd {
namespace opt {

    /**
     * @brief Perform Newton's Method with a barrier  to minimize a function
     * f(x) + barrier(g(x)).
     *
     * @param[in] problem The optimization problem to solve.
     * @param[in] settings The solver settings to use.
     * @return The results of the optimization including x* and minf.
     */
    OptimizationResults solve_problem_with_barrier_newton(
        OptimizationProblem& problem, SolverSettings& settings);

    /**
     * @brief Apply the constraints of the general problem as a barrier added to
     * the objective of the original problem.
     *
     * @param[in] general_problem Problem to convert to a problem with barrier
     *                            on the objective
     * @param[in] epsilon The barrier epsilon controlling the steepness. This
     *                    value is stored by reference in the barrier functions.
     * @paramp[out] barrier_problem Problem to store the modified objective of
     *                              the general problem.
     */
    void setup_barrier_problem(OptimizationProblem& general_problem,
        double& epsilon, OptimizationProblem& barrier_problem);

} // namespace opt
} // namespace ccd
