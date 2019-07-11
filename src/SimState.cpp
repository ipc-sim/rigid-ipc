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
{
}

void SimState::load_scene(const std::string& filename)
{
    using nlohmann::json;
    std::ifstream input(filename);
    if (input.good()){
        scene_file = filename;
        json scene = json::parse(input);
        init(scene);
    }
}

void SimState::reload_scene(){
    return load_scene(scene_file);
}

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
              "use_chain_functional": false
         },
        "particles_problem":{
            "vertices":[],
            "edges":[],
            "velocities":[],
            "x_fixed":[],
            "y_fixed":[],
            "use_mass_matrix": true,
            "constraint": "barrier_constraint",
            "update_constraint_set": true
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

    m_num_simulation_steps = 0;
}

void SimState::simulation_step()
{
    m_step_has_collision = false;
    m_step_had_collision = problem_ptr->simulation_step(m_timestep_size);

    if (m_step_had_collision) {
        Eigen::VectorXd gamma_1 = solve_collision();
        m_step_has_collision = problem_ptr->take_step(gamma_1, m_timestep_size);

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
    problem_ptr->update_constraint();
    auto result = ccd_solver_ptr->solve(*problem_ptr);

    return result.x;
}
} // namespace ccd
