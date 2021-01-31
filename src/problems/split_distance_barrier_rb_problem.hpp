#pragma once

#include <problems/distance_barrier_rb_problem.hpp>

namespace ipc::rigid {

/// @brief A diststance barrier rigid body problem using a split time
/// stepping method.
class SplitDistanceBarrierRBProblem : public DistanceBarrierRBProblem {
public:
    SplitDistanceBarrierRBProblem();
    virtual ~SplitDistanceBarrierRBProblem() = default;

    virtual bool settings(const nlohmann::json& params) override;
    nlohmann::json settings() const override;

    /// The name of the class as a string
    static std::string problem_name()
    {
        return "split_distance_barrier_rb_problem";
    }

    /// The name of the class as a string
    virtual std::string name() const override
    {
        return SplitDistanceBarrierRBProblem::problem_name();
    }

    ////////////////////////////////////////////////////////////
    // Rigid Body Problem

    /// @brief Solve for and take a time step.
    /// @param[out] had_collision True if the step had collisions.
    /// @param[out] has_intersections True if the step resulted in
    ///     intersections.
    /// @param[in] solve_collision True if collisions should be resolved.
    void simulation_step(
        bool& had_collisions,
        bool& has_intersections,
        bool solve_collisions = true) override;

    ////////////////////////////////////////////////////////////
    // Barrier Problem

    /// @brief Compute \f$E(x)\f$ in
    /// \f$f(x) = E(x) + \kappa \sum_{k \in C} b(d(x_k))\f$
    double compute_energy_term(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad = true,
        bool compute_hess = true) override;

    // Include thes lines to avoid issues with overriding inherited
    // functions with the same name.
    // (http://www.cplusplus.com/forum/beginner/24978/)
    using BarrierProblem::compute_barrier_term;
    using BarrierProblem::compute_energy_term;

protected:
    /// Update the stored poses and inital value for the solver.
    void update_dof() override;

    /// Update problem using current status of bodies.
    void update_constraints() override;

    /// Take a step from the current DOF to the provided ones.
    bool take_step(const Eigen::VectorXd& x) override;

    /// Apply restitution to solve for the current velocities.
    void solve_velocities();

    /// Unconstrained time-stepping method
    std::shared_ptr<TimeStepper> m_time_stepper;

    /// Rigid body poses at end of time-step
    PosesD poses_t1;

    /// @brief Original impacts used for velocity resitution
    /// @todo Replace this with the std::vector version
    Impacts original_impacts;
};

} // namespace ipc::rigid
