#include "SimState.hpp"

#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>

#include <io/read_rb_scene.hpp>
#include <logger.hpp>
#include <opt/barrier_constraint.hpp>

namespace ccd {

SimState::SimState()
    : m_timestep_size(0.1)
{

}

void SimState::load_scene(const std::string& filename)
{
    using nlohmann::json;
    std::ifstream input(filename);
    json scene = json::parse(input);

    if (io::is_rb_scene(scene)) {
        std::vector<physics::RigidBody> rbs;
        io::read_rb_scene(scene, rbs);
        m_problem.init(rbs, m_ccd_constraint);
        m_num_simulation_steps = 0.0;
    }
}

void SimState::simulation_step()
{
    m_step_had_collision = m_problem.simulation_step(m_timestep_size);

    if (m_step_had_collision) {
        Eigen::VectorXd gamma_1 = solve_collision();
        m_step_has_collision = m_problem.take_step(gamma_1, m_timestep_size);

        if (m_step_has_collision) {
            spdlog::error("simulation_step it={} collisions=unsolved",
                m_num_simulation_steps);
        } else {
            spdlog::trace("simulation_step it={} collisions=resolved",
                m_num_simulation_steps);
        }
    }
    m_num_simulation_steps += 1;
}

Eigen::MatrixXd SimState::solve_collision()
{
    m_problem.update_constraint();
    m_problem.validate_problem();
    auto result = m_ccd_solver.solve(m_problem);

    return result.x;
}
} // namespace ccd
