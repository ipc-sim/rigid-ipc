#include "SimState.hpp"

#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>

#include <io/read_rb_scene.hpp>
#include <io/serialize_json.hpp>
#include <logger.hpp>

#include <physics/problem_factory.hpp>
#include <solvers/solver_factory.hpp>

#include <utils/regular_2d_grid.hpp>

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
        "timestep_size": 0.1,
        "viewport_bbox": {"min":[0,0,0],"max":[0,0,0]}
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


    // scale background grid
    Eigen::VectorXd vmin, vmax;
    io::from_json(args["viewport_bbox"]["min"], vmin);
    io::from_json(args["viewport_bbox"]["max"], vmax);
    if ((vmax - vmin).squaredNorm() == 0.0){
        // default size to fix data
        vmin = problem_ptr->vertices().colwise().minCoeff();
        vmax = problem_ptr->vertices().colwise().maxCoeff();
    }

    regular_2d_grid(50, /*triangle=*/false, grid_V, grid_F);
    // center at 0.0
    grid_V.col(0).array() -= 0.5;
    grid_V.col(1).array() -= 0.5;

    double scale = (vmax - vmin).norm();
    grid_V *= scale;
    grid_V.rowwise() += problem_ptr->vertices().colwise().mean();


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

void SimState::get_collision_functional_isolines(Eigen::VectorXd& fx){
    if (m_dirty_constraints){
        problem_ptr->update_constraint();
        ccd_solver_ptr->init(*problem_ptr);
        m_dirty_constraints = false;
    }

    Eigen::MatrixXd Xk;
    problem_ptr->create_sample_points(grid_V, Xk);
    ccd_solver_ptr->eval_f(Xk, fx);

}
} // namespace ccd
