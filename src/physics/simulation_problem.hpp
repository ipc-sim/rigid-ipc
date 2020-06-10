#pragma once

#include <nlohmann/json.hpp>

#include <opt/collision_constraint.hpp>
#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>

#include <solvers/optimization_solver.hpp>

namespace ccd {
namespace physics {

    /**
     * @brief Interface class for simulation specific methods.
     *
     * Simulation problems are not optimization problems because they do not
     * necessarily have an objective function to optimize.
     */
    class SimulationProblem {
    public:
        virtual ~SimulationProblem() = default;
        enum CollisionCheck { EXACT = 0, CONSERVATIVE };

        virtual std::string name() const = 0;

        virtual opt::CollisionConstraint& constraint() = 0;
        virtual const opt::CollisionConstraint& constraint() const = 0;
        virtual opt::OptimizationSolver& solver() = 0;

        /// Get the settings of the simulation
        virtual nlohmann::json settings() const = 0;
        /// Set the settings of the simulation
        virtual void settings(const nlohmann::json& params) = 0;

        /// Get the state of the simulation
        virtual nlohmann::json state() const = 0;
        /// Set the state of the simulation
        virtual void state(const nlohmann::json& s) = 0;

        virtual double timestep() const = 0;        ///< Get the timestep size
        virtual void timestep(double timestep) = 0; ///< Set the timestep size

        /**
         * @brief Takes a step in the simulation
         * @param had_collision True if the step had collisions.
         * @param has_intersections True if the resulting bodies are
         *                          intersecting.
         * @param solve_collision True if collisions should be solved.
         */
        virtual void simulation_step(
            bool& had_collision,
            bool& has_intersections,
            bool solve_collision = true) = 0;

        /// @brief Check for intersections at the end of the time-step.
        virtual bool has_intersections() const = 0;

        virtual Eigen::MatrixXd vertices() const = 0;
        virtual const Eigen::MatrixXi& edges() const = 0;
        virtual const Eigen::MatrixXi& faces() const = 0;
        virtual Eigen::MatrixXd velocities() const = 0;
        virtual const Eigen::VectorXi& group_ids() const = 0;

        virtual const Eigen::MatrixXb& particle_dof_fixed() const = 0;

        virtual int dim() const = 0; ///< Spatial dimension (e.g. 3 for 3D)

        virtual bool is_rb_problem() const { return false; };
    };

} // namespace physics
} // namespace ccd
