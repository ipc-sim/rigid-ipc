/// Solve a optimization problem with NLopt.
#pragma once

#include <nlopt.hpp>

#include <array>

#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>

#include <solvers/optimization_solver.hpp>

namespace ccd {
namespace opt {

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

    enum class NLOptAlgorithm {
        /// \f$\Delta x = A^{-1} jac_x(x_i)^T \alpha_i\f$
        G_GRADIENT,
        /// \f$\Delta x = A^{-1} jac_x(x_i)^T \alpha_i + A^{-1}b - x_i\f$
        LINEARIZED
    };
    static const std::array<nlopt::algorithm, 2> NLOptAlgorithm
        = { { nlopt::LD_MMA, nlopt::LD_SLSQP } };

    static const char* NLOptAlgorithmNames[] = { "MMA", "SLSQP" };

    class NLOptSolver : public OptimizationSolver {
    public:
        NLOptSolver();
        ~NLOptSolver() override;
        OptimizationResults solve(OptimizationProblem& problem) override;

        // Settings
        // -----------
        nlopt::algorithm algorithm;
        double absolute_tolerance;
        double relative_tolerance;
        int max_iterations;
        double max_time;
        bool verbose;
    };

} // namespace opt
} // namespace ccd
