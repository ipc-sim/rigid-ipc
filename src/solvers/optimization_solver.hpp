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

#include <logger.hpp>
#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>
#include <utils/not_implemented_error.hpp>

namespace ipc::rigid {

class OptimizationSolver {
public:
    virtual ~OptimizationSolver() {};

    /// Initialize the state of the solver using the settings saved in JSON
    virtual void settings(const nlohmann::json& params) = 0;
    /// Export the state of the solver using the settings saved in JSON
    virtual nlohmann::json settings() const = 0;

    /// An identifier for this solver
    virtual std::string name() const = 0;

    /// Initialize the solver with a problem to solve
    virtual void set_problem(OptimizationProblem& problem) = 0;

    /// Initialize the solver state for a new solve
    virtual void init_solve(const Eigen::VectorXd& x0) = 0;
    /// Solve the saved optimization problem to completion
    virtual OptimizationResults solve(const Eigen::VectorXd& x0) = 0;
    /// Perform a single step of solving the optimization problem
    virtual OptimizationResults step_solve() = 0;

    virtual std::string stats_string() const { return ""; }
    virtual nlohmann::json stats() const { return nlohmann::json(); }

    virtual bool has_inner_solver() const { return false; }
    virtual OptimizationSolver& inner_solver()
    {
        throw NotImplementedError(
            fmt::format("{} does not have and inner_solver", name()));
    }
};

} // namespace ipc::rigid
