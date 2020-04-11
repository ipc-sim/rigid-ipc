#pragma once

#include <nlohmann/json.hpp>

#include <memory> // shared_ptr

#include <physics/simulation_problem.hpp>
#include <solvers/optimization_solver.hpp>

namespace ccd {

class SimState {
public:
    SimState();

    bool load_scene(const std::string& filename);
    bool reload_scene();
    bool load_simulation(const nlohmann::json& args);
    bool init(const nlohmann::json& args);

    void simulation_step();
    bool solve_collision();
    void collision_resolution_step();

    void save_simulation(const std::string& filename);
    void save_simulation_step();

    void run_simulation(const std::string& fout);

    const nlohmann::json& get_config() { return args; }
    nlohmann::json get_active_config();

    // CCD
    // ----------------------------------------------
    std::shared_ptr<physics::ISimulationProblem> problem_ptr;
    double m_timestep_size;

    bool m_step_had_collision;     ///< last step had a collision
    bool m_step_has_collision;     ///< last step failed to solve collisions
    bool m_step_has_intersections; ///< last step resulted in intersections
    bool m_solve_collisions;    ///< solve collisions automatically on the step
    int m_num_simulation_steps; ///< counts simulation steps
    int m_max_simulation_steps;

    std::string scene_file;

    nlohmann::json args;

protected:
    bool m_dirty_constraints;
    std::vector<Eigen::MatrixXd> vertices_sequence;
    std::vector<nlohmann::json> state_sequence;
};

} // namespace ccd
