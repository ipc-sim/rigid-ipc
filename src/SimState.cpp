#include "SimState.hpp"

#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>

#include <io/read_rb_scene.hpp>
#include <logger.hpp>

#include <physics/problem_factory.hpp>
#include <solvers/solver_factory.hpp>

namespace ccd {

SimState::SimState()
    : m_timestep_size(0.1)
    , m_solve_collisions(true)
    , m_dirty_constraints(false)
{
}

void SimState::load_scene(const std::string& filename)
{
    using nlohmann::json;
    std::ifstream input(filename);
    if (input.good()) {
        scene_file = filename;
        json scene = json::parse(input);
        init(scene);
    }
}

void SimState::reload_scene() { return load_scene(scene_file); }

void SimState::init(const nlohmann::json& args_in)
{
    using namespace nlohmann;

    // clang-format off
    args = R"({
        "scene_type":"rigid_body_problem",
        "collision_solver":"barrier_solver",
        "rigid_body_problem":{
              "rigid_bodies": [],
              "constraint": "barrier_constraint",
              "update_constraint_set": true,
              "use_chain_functional": false,
              "gravity":[0.0,0.0,0.0],
              "collision_eps": 0.1
         },
        "particles_problem":{
            "vertices":[],
            "edges":[],
            "velocities":[],
            "x_fixed":[],
            "y_fixed":[],
            "use_mass_matrix": true,
            "constraint": "barrier_constraint",
            "update_constraint_set": true,
            "gravity":[0.0,0.0],
            "collision_eps": 0.1
         },
        "barrier_solver": {},
        "timestep_size": 0.1
    })"_json;
    // clang-format on

    args.merge_patch(args_in);
    m_timestep_size = args["timestep_size"].get<double>();
    problem_ptr
        = physics::ProblemFactory::factory().get_problem(args["scene_type"]);
    problem_ptr->init(args[args["scene_type"].get<std::string>()]);

    ccd_solver_ptr
        = opt::SolverFactory::factory().get_solver(args["collision_solver"]);
    ccd_solver_ptr->clear();

    m_num_simulation_steps = 0;
    m_dirty_constraints = true;
}

void SimState::simulation_step()
{
    m_step_has_collision = false;
    m_step_had_collision = problem_ptr->simulation_step(m_timestep_size);
    m_dirty_constraints = true;

    if (m_solve_collisions && m_step_had_collision) {
        solve_collision();
    }
    m_num_simulation_steps += 1;
}

bool SimState::solve_collision()
{
    if (m_dirty_constraints){
        problem_ptr->update_constraint();
        m_dirty_constraints = false;
    }
    auto result = ccd_solver_ptr->solve(*problem_ptr);
    m_step_has_collision = problem_ptr->take_step(result.x, m_timestep_size);

    if (m_step_has_collision) {
        spdlog::error("simulation_step it={} collisions=unsolved",
            m_num_simulation_steps);
    } else {
        spdlog::trace("simulation_step it={} collisions=resolved",
            m_num_simulation_steps);
    }
    return m_step_has_collision;
}
void SimState::collision_resolution_step(){
    if (m_dirty_constraints){
        problem_ptr->update_constraint();
        ccd_solver_ptr->init(*problem_ptr);
        m_dirty_constraints = false;

    }
    auto result = ccd_solver_ptr->step_solve();
    // TODO: use results.finished
    problem_ptr->take_step(result.x, m_timestep_size);
}
} // namespace ccd
