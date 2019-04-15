/**
 * @brief Solve the optimization problem using Newton's Method with barriers for
 * the constraints.
 */

#include <opt/OptimizationProblem.hpp>
#include <opt/OptimizationResults.hpp>
#include <opt/SolverSettings.hpp>

namespace ccd {
namespace opt {

    void setup_barrier_problem(const OptimizationProblem& general_problem,
        double& epsilon, OptimizationProblem& barrier_problem);

    // Performa Newton's Method to minimize a function f(x).
    OptimizationResults solve_problem_with_barrier_newton(
        const OptimizationProblem& problem, const SolverSettings& settings);

} // namespace opt
} // namespace ccd
