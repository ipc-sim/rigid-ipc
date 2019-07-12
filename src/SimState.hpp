#pragma once

#include <nlohmann/json.hpp>

#include <memory> // shared_ptr

#include <solvers/optimization_solver.hpp>
#include <physics/simulation_problem.hpp>

namespace ccd {

class SimState {
public:
    SimState();

    void load_scene(const std::string& filename);
    void reload_scene();
    void init(const nlohmann::json& args);

    void simulation_step();
    bool solve_collision();
    void collision_resolution_step();

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
protected:
    bool m_dirty_constraints;
};

} // namespace ccd
