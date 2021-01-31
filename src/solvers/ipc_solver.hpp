#pragma once

#include <problems/barrier_problem.hpp>
#include <solvers/newton_solver.hpp>
#include <utils/not_implemented_error.hpp>

namespace ipc::rigid {

class IPCSolver : public virtual NewtonSolver {
public:
    IPCSolver();
    virtual ~IPCSolver() = default;

    /// Initialize the state of the solver using the settings saved in JSON
    virtual void settings(const nlohmann::json& params) override;
    /// Export the state of the solver using the settings saved in JSON
    virtual nlohmann::json settings() const override;

    /// An identifier for the solver class
    static std::string solver_name() { return "ipc_solver"; }
    /// An identifier for this solver
    virtual std::string name() const override
    {
        return IPCSolver::solver_name();
    }

    /// Initialize the solver state for a new solve
    virtual void init_solve(const Eigen::VectorXd& x0) override;
    /// Solve the saved optimization problem to completion
    virtual OptimizationResults solve(const Eigen::VectorXd& x0) override;

    virtual std::string stats_string() const override;
    virtual nlohmann::json stats() const override;

protected:
    /// Adaptivly update the barrier stiffness at the end of each step
    void post_step_update() override;

    BarrierProblem* barrier_problem_ptr()
    {
        if (problem_ptr == nullptr) {
            return nullptr;
        }
        assert(problem_ptr->is_barrier_problem());
        return dynamic_cast<BarrierProblem*>(problem_ptr);
    }

    ///////////////////////////////////////////////////////////////////////
    // User simulation parameters

    /// @brief Minimum barrier stiffness scale.
    double min_barrier_stiffness_scale;

    /// @brief Activation distance of adaptive barrier stiffness.
    double dhat_epsilon;

    ///////////////////////////////////////////////////////////////////////
    // Computed values

    /// @brief The max value for adaptive barrier stiffness.
    double max_barrier_stiffness;

    /// @brief The minimum distance of the previous iteration.
    double prev_min_distance;

private:
    int num_kappa_updates = 0;
};

} // namespace ipc::rigid
