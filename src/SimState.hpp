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
    void init(const nlohmann::json& args);

    void simulation_step();
    bool solve_collision();
    void collision_resolution_step();

    void get_collision_gradient(Eigen::MatrixXd& fx);

    void save_simulation(const std::string& filename);

    const nlohmann::json& get_config() { return args; }
    nlohmann::json get_active_config();

    // CCD
    // ----------------------------------------------
    std::shared_ptr<physics::SimulationProblem> problem_ptr;
    std::shared_ptr<opt::OptimizationSolver> ccd_solver_ptr;
    double m_timestep_size;

    bool m_step_had_collision;  ///< last step had a collision
    bool m_step_has_collision;  ///< last step failed to solve collisions
    bool m_solve_collisions;    ///< solve collisions automatically on the step
    int m_num_simulation_steps; ///< counts simulation steps

    std::string scene_file;

    nlohmann::json args;

    // for visualization
    Eigen::MatrixXd grid_V;
    Eigen::MatrixXi grid_F;

protected:
    bool m_dirty_constraints;
    void init_background_grid(const nlohmann::json& args);
    std::vector<Eigen::MatrixXd> vertices_sequence;
};

} // namespace ccd
