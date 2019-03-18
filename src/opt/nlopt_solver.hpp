/// Solve a optimization problem with NLopt.
#pragma once

#include <nlopt.hpp>

#include <opt/OptimizationProblem.hpp>
#include <opt/OptimizationResults.hpp>
#include <opt/SolverSettings.hpp>

namespace ccd {
namespace opt {

    /**
     * @brief Optimize the displacments using NLopt
     *
     * @param[in] algorithm NLopt algorithm to use for optimization
     * @param[in,out] problem Optimization problem to solve
     * @return The results of the optimization
     */
    OptimizationResults solve_problem_with_nlopt(
        OptimizationProblem& problem, const SolverSettings& settings);

    /**
     * @brief Computes NLopt's objective for an OptimizationProblem
     *
     * @param[in] x Optimization variable for which to compute \f$f(x)\f$ and
     *              \f$\nabla f(x)\f$.
     * @param[out] grad Location to store \f$\nabla f(x)\f$.
     * @param[in] data A pointer to the OptimizationProblem.
     * @return Returns the value \f$f(x)\f$.
     */
    double nlopt_objective(
        const std::vector<double>& x, std::vector<double>& grad, void* data);
    /**
     * @brief Computes NLopt's constraints for an OptimizationProblem
     *
     * These are NLopt inequality constraints where constraints(x) ≤ 0.
     *
     * @param[in] m Number of constraints
     * @param[out] results Storage for the constraint values
     * @param[in] n Number of varaibles
     * @param[in] x The optimization variables
     * @param[out] grad Storage for the constraint Jacobian values
     * @param[in] data A pointer to the OptimizationProblem.
     */
    void nlopt_inequality_constraints(unsigned m, double* results, unsigned n,
        const double* x, double* grad, void* data);

    /**
     * @brief Output the reason for NLopt's termination
     *
     * @param[in] opt Optimizer object
     * @param[in] reuslt Optimization result to parse
     */
    void print_nlopt_termination_reason(
        const nlopt::opt& opt, const nlopt::result result);

} // namespace opt
} // namespace ccd
