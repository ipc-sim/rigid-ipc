#pragma once

#include <nlohmann/json.hpp>

#include <opt/collision_constraint.hpp>
#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>

#include <solvers/optimization_solver.hpp>

namespace ccd {
namespace physics {

    /// @brief interface class for simulation specific methods
    class ISimulationProblem {
    public:
        virtual ~ISimulationProblem() = default;
        enum CollisionCheck { EXACT = 0, CONSERVATIVE };

        virtual opt::CollisionConstraint& constraint() = 0;
        virtual const opt::CollisionConstraint& constraint() const = 0;
        virtual opt::IStateOptimizationSolver& solver() = 0;

        virtual std::string name() = 0;
        virtual void settings(const nlohmann::json& params) = 0;
        virtual nlohmann::json settings() const = 0;
        virtual nlohmann::json state() const = 0;
        virtual void state(const nlohmann::json& s) = 0;

        /// @brief  does a single simulation step. Returns true if there is
        /// a collision
        virtual bool simulation_step(const double time_step) = 0;

        /// @brief moves status to given positions
        virtual bool
        take_step(const Eigen::VectorXd& positions, const double time_step) = 0;

        /// @brief Check for intersections at the end of the time-step.
        virtual bool has_intersections() const = 0;

        /// \brief update optimization problem using current status.
        virtual void update_constraint() = 0;
        virtual opt::OptimizationResults solve_constraints() = 0;
        virtual void init_solve() = 0;
        virtual opt::OptimizationResults step_solve() = 0;

        /// -----------------------------------------------------------------
        virtual Eigen::MatrixXd vertices() const = 0;
        virtual const Eigen::MatrixXi& edges() const = 0;
        virtual const Eigen::MatrixXi& faces() const = 0;
        virtual Eigen::MatrixXd velocities() const = 0;
        virtual Eigen::VectorXi group_ids() const = 0;

        virtual const Eigen::MatrixXb& particle_dof_fixed() const = 0;

        virtual bool is_rb_problem() const { return false; };
    };

} // namespace physics
} // namespace ccd
