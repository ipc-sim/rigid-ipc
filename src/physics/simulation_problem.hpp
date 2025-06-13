#pragma once

#include <nlohmann/json.hpp>

#include <opt/collision_constraint.hpp>
#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>

#include <solvers/optimization_solver.hpp>

namespace ipc::rigid {

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

    virtual CollisionConstraint& constraint() = 0;
    virtual const CollisionConstraint& constraint() const = 0;
    virtual OptimizationSolver& solver() = 0;

    /// Get the settings of the simulation
    virtual nlohmann::json settings() const = 0;
    /// Set the settings of the simulation
    virtual bool settings(const nlohmann::json& params) = 0;

    /// Get the state of the simulation
    virtual nlohmann::json state() const = 0;
    /// Set the state of the simulation
    virtual void state(const nlohmann::json& s) = 0;

    virtual double timestep() const = 0;        ///< Get the timestep size
    virtual void timestep(double timestep) = 0; ///< Set the timestep size

    /// @brief Takes a step in the simulation
    /// @param[out] had_collisions True if the step had collisions.
    /// @param[out] has_intersections True if the resulting bodies are
    ///                               intersecting.
    /// @param[in] solve_collisions True if collisions should be solved.
    virtual void simulation_step(
        bool& had_collisions,
        bool& has_intersections,
        bool solve_collisions = true) = 0;

    virtual int dim() const = 0; ///< Spatial dimension (e.g. 3 for 3D)
    virtual long num_vertices() const = 0;
    virtual long num_edges() const = 0;
    virtual long num_faces() const = 0;
    virtual size_t num_bodies() const = 0;

    virtual Eigen::MatrixXd vertices() const = 0;
    virtual const std::vector<int>& codim_vertices_to_vertices() const = 0;
    virtual const Eigen::MatrixXi& edges() const = 0;
    virtual const std::vector<int>& codim_edges_to_edges() const = 0;
    virtual const Eigen::MatrixXi& faces() const = 0;
    virtual Eigen::MatrixXd vertices(size_t i) const { return vertices(); }
    virtual const Eigen::MatrixXi& edges(size_t i) const { return edges(); }
    virtual const std::vector<int>& codim_edges_to_edges(size_t i) const
    {
        return codim_edges_to_edges();
    }
    virtual const Eigen::MatrixXi& faces(size_t i) const { return faces(); }
    virtual Eigen::MatrixXd velocities() const = 0;
    virtual const Eigen::VectorXi& group_ids() const = 0;

    virtual const MatrixXb& vertex_dof_fixed() const = 0;

    virtual bool is_rb_problem() const { return false; };

    virtual int num_contacts() const = 0;

    /// Compute the minimum distance among geometry
    virtual double compute_min_distance() const = 0;

    OptimizationResults opt_result;
};

} // namespace ipc::rigid
