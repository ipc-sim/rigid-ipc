#pragma once

#include <nlohmann/json.hpp>

#include <opt/optimization_problem.hpp>
#include <opt/collision_constraint.hpp>

namespace ccd {
namespace physics {

    /// @brief interface class for simulation specific methods
    class SimulationProblem : public opt::OptimizationProblem {
    public:
        SimulationProblem(const std::string& name);

        virtual void init(const nlohmann::json& params)
            = 0;

        virtual const opt::CollisionConstraint& constraint()=0;

        /// @brief  does a single simulation step. Returns true if there is
        /// a collision
        virtual bool simulation_step(const double time_step) = 0;

        /// @brief moves status to given positions
        virtual bool take_step(
            const Eigen::VectorXd& positions, const double time_step)
            = 0;

        /// \brief update optimization problem using current status.
        virtual void update_constraint() = 0;

        virtual Eigen::MatrixXd vertices() = 0;
        virtual const Eigen::MatrixXi& edges() = 0;

        virtual const Eigen::MatrixXb& particle_dof_fixed()=0;

    };

} // namespace physics
} // namespace ccd
