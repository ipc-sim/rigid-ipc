/**
 *  We handle optimization problems of the form
 *      MIN     f(x)      x ∈ Rⁿ
 *
 *   s.t.       g_L ≤ g(x) ≤ g_U
 *              x_L ≤  x   ≤ x_U
 *
 */
#pragma once

#include <nlohmann/json.hpp>

#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>

namespace ccd {
namespace opt {

    /// Basic Interface class for non-spetiallized solvers
    class IBasicOptimizationSolver {
    public:
        virtual ~IBasicOptimizationSolver() = default;
        virtual OptimizationResults solve(IConstraintedProblem& problem) = 0;
    };

    // Interface class for optimization solvers used by BarrierSolver
    class IBarrierOptimizationSolver {
    public:
        virtual ~IBarrierOptimizationSolver() = default;
        virtual OptimizationResults solve(IBarrierProblem& problem) = 0;
        virtual const std::string& name() const = 0;
        virtual void settings(const nlohmann::json& json) = 0;
        virtual nlohmann::json settings() const = 0;
        virtual void init_free_dof(Eigen::VectorXb is_dof_fixed) = 0;

        virtual void c(const double) {}
        virtual void e_b(const double) {}
        virtual void t(const double) {}
        virtual void m(const double) {}

#ifdef DEBUG_LINESEARCH
        virtual void debug_stats(){}
#endif
    };

    // Interface class for optimizatino solvers used by State
    class IStateOptimizationSolver {
    public:
        virtual ~IStateOptimizationSolver() = default;

        virtual OptimizationResults solve() = 0;
        virtual void init_solve() = 0;
        virtual OptimizationResults step_solve() = 0;

        virtual void settings(const nlohmann::json& json) = 0;
        virtual nlohmann::json settings() const = 0;
        virtual const std::string& name() const = 0;

        virtual bool has_inner_solver() = 0;
        virtual const IBarrierOptimizationSolver& inner_solver() = 0;

        // UI debugging
        virtual Eigen::VectorXd get_grad_kkt() const = 0;
        virtual int num_outer_iterations() const = 0;
    };

} // namespace opt
} // namespace ccd
