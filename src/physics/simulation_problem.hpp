#pragma once

#include <nlohmann/json.hpp>

#include <opt/collision_constraint.hpp>
#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>

#include <solvers/optimization_solver.hpp>

namespace ccd {
namespace physics {

    /// @brief interface class for simulation specific methods
    /// Simulation problems are not optimization problems because they do not
    /// necessarily have an objective function to optimize.
    class SimulationProblem {
    public:
        virtual ~SimulationProblem() = default;
        enum CollisionCheck { EXACT = 0, CONSERVATIVE };

        virtual std::string name() const = 0;

        virtual opt::CollisionConstraint& constraint() = 0;
        virtual const opt::CollisionConstraint& constraint() const = 0;
        virtual opt::OptimizationSolver& solver() = 0;

        virtual nlohmann::json settings() const = 0;
        virtual void settings(const nlohmann::json& params) = 0;

        virtual nlohmann::json state() const = 0;
        virtual void state(const nlohmann::json& s) = 0;

        virtual double timestep() const = 0;
        virtual void timestep(double timestep) = 0;

        /// @brief Compute the step but does not take it
        /// @returns true if there is a collision
        virtual bool simulation_step() = 0;

        /// @brief moves status to given positions
        virtual bool take_step(const Eigen::VectorXd& x) = 0;

        /// @brief Check for intersections at the end of the time-step.
        virtual bool has_intersections() const = 0;

        /// @brief update optimization problem using current status
        virtual void update_constraint() = 0;
        virtual opt::OptimizationResults solve_constraints() = 0;
        virtual void init_solve() = 0;
        virtual opt::OptimizationResults step_solve() = 0;

        /// -----------------------------------------------------------------
        virtual Eigen::MatrixXd vertices() const = 0;
        virtual const Eigen::MatrixXi& edges() const = 0;
        virtual const Eigen::MatrixXi& faces() const = 0;
        virtual Eigen::MatrixXd velocities() const = 0;
        virtual const Eigen::VectorXi& group_ids() const = 0;

        virtual const Eigen::MatrixXb& particle_dof_fixed() const = 0;

        virtual int dim() const = 0;

        virtual bool is_rb_problem() const { return false; };
    };

} // namespace physics
} // namespace ccd
