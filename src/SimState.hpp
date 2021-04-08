#pragma once

#include <igl/Timer.h>
#include <nlohmann/json.hpp>

#include <memory> // shared_ptr

#include <physics/simulation_problem.hpp>
#include <solvers/optimization_solver.hpp>

namespace ipc::rigid {

class SimState {
public:
    SimState();

    bool load_scene(const std::string& filename, const std::string& patch = "");
    bool reload_scene();
    bool load_simulation(const nlohmann::json& args);
    bool init(const nlohmann::json& args);

    void simulation_step();

    bool save_simulation(const std::string& filename);
    void save_simulation_step();

    bool save_obj_sequence(const std::string& dir_name);
    bool save_gltf(const std::string& filename);

    void run_simulation(const std::string& fout);

    const nlohmann::json& get_config() { return args; }
    nlohmann::json get_active_config();

    // ----------------------------------------------
    std::shared_ptr<SimulationProblem> problem_ptr;

    bool m_step_had_collision;     ///< last step had a collision
    bool m_step_has_collision;     ///< last step failed to solve collisions
    bool m_step_has_intersections; ///< last step resulted in intersections
    bool m_solve_collisions;    ///< solve collisions automatically on the step
    int m_num_simulation_steps; ///< counts simulation steps
    int m_max_simulation_steps; ///< maximum number of time-steps to take
    int m_checkpoint_frequency; ///< time-steps between checkpoints

    std::string scene_file;

    nlohmann::json args;

    std::vector<nlohmann::json> state_sequence;
    std::vector<double> step_timings;
    std::vector<int> solver_iterations;
    std::vector<int> num_contacts;
    std::vector<double> step_minimum_distances;

protected:
    igl::Timer step_timer;
    size_t initial_rss;

    bool m_dirty_constraints;
};

} // namespace ipc::rigid
