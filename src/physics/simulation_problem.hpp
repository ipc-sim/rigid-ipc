#pragma once

#include <nlohmann/json.hpp>

#include <opt/collision_constraint.hpp>
#include <opt/optimization_problem.hpp>

namespace ccd {
namespace physics {

    /// @brief interface class for simulation specific methods
    class SimulationProblem : public opt::OptimizationProblem {
    public:
        enum CollisionCheck { EXACT = 0, CONSERVATIVE };

        SimulationProblem(const std::string& name);

        virtual void init(const nlohmann::json& params) = 0;

        virtual const opt::CollisionConstraint& constraint() = 0;

        /// @brief  does a single simulation step. Returns true if there is
        /// a collision
        virtual bool simulation_step(const double time_step) = 0;

        /// @bried return candidate position of vertices, assuming no collisions
        virtual Eigen::MatrixXd vertices_next(const double time_step) = 0;

        /// @brief moves status to given positions
        virtual bool take_step(
            const Eigen::VectorXd& positions, const double time_step)
            = 0;

        /// \brief update optimization problem using current status.
        virtual void update_constraint() = 0;

        /// \brief position of vertices at the END of the step
        virtual Eigen::MatrixXd vertices() = 0;

        /// \brief position of vertices at the BEGINING of the step
        virtual Eigen::MatrixXd vertices_prev() = 0;

        /// \brief velocity of vertices at the END of the step
        /// as_delta = true, will return vertices / time_step;
        virtual Eigen::MatrixXd velocities(
            const bool as_delta, const double time_step)
            = 0;

        /// \brief collision forced applied to fix END position of vertices
        /// as_delta = true, will return Fc M^-1 / (time_step^2);
        virtual Eigen::MatrixXd collision_force(
            const bool as_delta, const double time_step)
            = 0;

        virtual const Eigen::MatrixXi& edges() = 0;

        virtual const Eigen::MatrixXb& particle_dof_fixed() = 0;
        virtual const Eigen::VectorXd& gravity() = 0;

        virtual void create_sample_points(
                const Eigen::MatrixXd& xy_points, Eigen::MatrixXd& sample_points) {
            throw NotImplementedError("create_sample_points not implemented");
        }
    };

} // namespace physics
} // namespace ccd
