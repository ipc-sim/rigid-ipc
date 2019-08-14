#include "SimState.hpp"

#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>

#include <io/read_rb_scene.hpp>
#include <io/serialize_json.hpp>
#include <logger.hpp>

#include <opt/constraint_factory.hpp>
#include <problems/problem_factory.hpp>
#include <solvers/solver_factory.hpp>

#include <profiler.hpp>
#include <utils/regular_2d_grid.hpp>

namespace ccd {

SimState::SimState()
    : m_timestep_size(0.1)
    , m_solve_collisions(true)
    , m_max_simulation_steps(-1)
    , m_dirty_constraints(false)

{
}

bool SimState::load_scene(const std::string& filename)
{
    using nlohmann::json;
    std::ifstream input(filename);
    PROFILER_CLEAR()

    if (input.good()) {
        scene_file = filename;
        json scene = json::parse(input, nullptr, false);
        if (scene.is_discarded()) {
            spdlog::error("Invalid Json file");
            return false;
        }
        if (scene.find("args") != scene.end()) {
            return load_simulation(scene);
        } else {
            return init(scene);
        }
    }
    return false;
}

bool SimState::reload_scene() { return load_scene(scene_file); }

bool SimState::load_simulation(const nlohmann::json& args)
{
    // load original setup
    bool success = init(args["args"]);
    if (!success) {
        return false;
    }

    // now reload simulation history
    vertices_sequence.clear();
    for (auto& jv : args["animation"]["vertices_sequence"]) {
        Eigen::MatrixXd v;
        io::from_json(jv, v);
        vertices_sequence.push_back(v);
    };
    state_sequence = args["animation"]["state_sequence"]
                         .get<std::vector<nlohmann::json>>();
    problem_ptr->state(state_sequence.back());
    return true;
}

bool SimState::init(const nlohmann::json& args_in)
{
    using namespace nlohmann;

    // clang-format off
    args = R"({
        "max_iterations":-1,
        "scene_type":"distance_barrier_rb_problem",
        "rigid_body_problem":{
              "rigid_bodies": [],
              "coefficient_restitution":0.0,
              "gravity":[0.0,0.0,0.0],
              "collision_eps": 0.0
         },
        "particles_problem":{
            "vertices":[],
            "edges":[],
            "velocities":[],
            "x_fixed":[],
            "y_fixed":[],
            "use_mass_matrix": true,
            "gravity":[0.0,0.0],
            "collision_eps": 0.0
         },
        "barrier_solver": {
            "inner_solver": "newton_solver",
            "min_barrier_epsilon":1e-5,
            "max_iterations": 0
        },
        "gradient_descent_solver":{
            "absolute_tolerance": 1e-5,
            "min_step_length": 1e-12,
            "max_iterations": 3000
        },
        "newton_solver":{
           "absolute_tolerance": 1e-5,
           "min_step_length": 1e-10,
           "max_iterations": 3000
        },
        "bfgs_solver":{
            "absolute_tolerance": 1e-5,
            "min_step_length": 1e-12,
            "max_iterations": 3000
        },
        "ncp_solver":{
            "max_iterations": 1000,
            "do_line_search": false,
            "solve_for_active_cstr": true,
            "convergence_tolerance":1e-6,
            "update_type": "linearized",
            "lcp_solver": "lcp_gauss_seidel"
        },
        "time_barrier_constraint":{
            "initial_epsilon":"min_toi",
            "custom_initial_epsilon":1.0,
            "detection_method": "hash_grid",
            "custom_hashgrid_cellsize":-1
        },
       "distance_barrier_constraint":{
           "custom_initial_epsilon":0.5,
           "detection_method": "hash_grid",
           "use_distance_hashgrid": true,
           "active_constraint_scale" : 1.5,
           "custom_hashgrid_cellsize":-1
       },
       "volume_constraint":{
           "detection_method": "hash_grid",
           "volume_epsilon": 1e-6,
           "custom_hashgrid_cellsize":-1
       },
        "timestep_size": 0.1,
        "viewport_bbox": {"min":[0,0],"max":[0,0]}
    })"_json;
    // clang-format on

    // check that incomming json doesn't have any unkown keys
    // to avoid stupid bugs
    auto patch = json::diff(args, args_in);
    bool valid = true;
    for (auto& op : patch) {
        if (op["op"].get<std::string>().compare("add") == 0) {
            auto new_path = json::json_pointer(op["path"].get<std::string>());
            if (args_in[new_path.parent_pointer()].is_array()) {
                valid = true;
            } else {
                valid = false;
                spdlog::error("Unknown key in json path={}",
                    op["path"].get<std::string>());
            }
        }
    }
    if (!valid) {
        return false;
    }

    args.merge_patch(args_in);
    m_max_simulation_steps = args["max_iterations"].get<int>();
    m_timestep_size = args["timestep_size"].get<double>();

    auto problem_name = args["scene_type"].get<std::string>();

    problem_ptr = ProblemFactory::factory().get_problem(problem_name);
    problem_ptr->settings(args);

    m_num_simulation_steps = 0;
    m_dirty_constraints = true;

    init_background_grid(args);

    vertices_sequence.clear();
    vertices_sequence.push_back(problem_ptr->vertices());

    state_sequence.clear();
    state_sequence.push_back(problem_ptr->state());

    return true;
}

nlohmann::json SimState::get_active_config()
{
    nlohmann::json active_args;
    active_args["timestep_size"] = m_timestep_size;
    active_args["scene_type"] = problem_ptr->name();

    active_args[problem_ptr->name()] = problem_ptr->settings();
    active_args[problem_ptr->constraint().name()]
        = problem_ptr->constraint().settings();
    active_args[problem_ptr->solver().name()]
        = problem_ptr->solver().settings();
    if (problem_ptr->solver().has_inner_solver()) {
        active_args[problem_ptr->solver().inner_solver().name()]
            = problem_ptr->solver().inner_solver().settings();
    }

    return active_args;
}

void SimState::run_simulation(const std::string& fout)
{

    PROFILE_MAIN_POINT("run_simulation")
    PROFILE_START()

    spdlog::info("starting simulation {}", scene_file);
    m_solve_collisions = true;
    for (int i = 0; i < m_max_simulation_steps; ++i) {
        simulation_step();
    }
    save_simulation(fout);
    spdlog::info("simulation results saved to {}", fout);

    PROFILE_END()
    LOG_PROFILER(scene_file);
}

void SimState::simulation_step()
{
    m_num_simulation_steps += 1;
    m_step_has_collision = false;
    m_step_had_collision = problem_ptr->simulation_step(m_timestep_size);
    m_dirty_constraints = true;

    if (m_step_had_collision) {
        spdlog::debug("sim_state action=simulation_step status=had_collision");
    }

    if (m_solve_collisions && m_step_had_collision) {
        solve_collision();
    }
    vertices_sequence.push_back(problem_ptr->vertices());
    state_sequence.push_back(problem_ptr->state());
}

bool SimState::solve_collision()
{
    PROFILE_POINT("sim_state_solve_collisions");

    if (m_dirty_constraints) {
        problem_ptr->update_constraint();
        m_dirty_constraints = false;
    }

    ccd::opt::OptimizationResults result;

    PROFILE_START();

    result = problem_ptr->solve_constraints();

    PROFILE_END();

    m_step_has_collision = problem_ptr->take_step(result.x, m_timestep_size);

    if (m_step_has_collision) {
        spdlog::error(
            "sim_state action=solve_collisions sim_it={} status=collisions_unsolved",
            m_num_simulation_steps);
    } else {
        spdlog::debug(
            "sim_state action=solve_collisions sim_it={} status=collisions_solved",
            m_num_simulation_steps);
    }
    return m_step_has_collision;
}

void SimState::collision_resolution_step()
{
    if (m_dirty_constraints) {
        problem_ptr->update_constraint();
        problem_ptr->init_solve();
        m_dirty_constraints = false;
    }

    ccd::opt::OptimizationResults result;
    result = problem_ptr->step_solve();

    // TODO: use results.finished
    m_step_has_collision = problem_ptr->take_step(result.x, m_timestep_size);
    spdlog::debug(
        "sim_state action=collision_resolution_step status=collisions_{}",
        m_step_has_collision ? "unsolved" : "solved");
}

void SimState::save_simulation(const std::string& filename)
{
    nlohmann::json results;
    results["args"] = args;
    results["active_args"] = get_active_config();
    std::vector<nlohmann::json> vs;
    for (auto& v : vertices_sequence) {
        vs.push_back(io::to_json(v));
    }
    results["animation"] = nlohmann::json();
    results["animation"]["vertices_sequence"] = vs;
    results["animation"]["state_sequence"] = state_sequence;
    results["animation"]["edges"] = io::to_json(problem_ptr->edges());

    std::ofstream o(filename);
    o << std::setw(4) << results << std::endl;
}

void SimState::get_collision_gradient(Eigen::MatrixXd& fgrad)
{
    if (m_dirty_constraints) {
        fgrad.resize(0, 0);
        return;
    }

    //    fgrad = ccd_solver_ptr->get_grad_kkt();
    problem_ptr->unflatten_dof(fgrad);
}

void SimState::init_background_grid(const nlohmann::json& args_)
{

    // scale background grid
    Eigen::VectorXd vmin, vmax;
    io::from_json(args_["viewport_bbox"]["min"], vmin);
    io::from_json(args_["viewport_bbox"]["max"], vmax);
    Eigen::VectorXd d = (vmax - vmin);
    if (d.squaredNorm() == 0.0) {
        // default size to fit data
        vmin = problem_ptr->vertices().colwise().minCoeff();
        vmax = problem_ptr->vertices().colwise().maxCoeff();
        // add some padding
        d = (vmax - vmin);
        vmin -= d / 2.0;
        vmax += d / 2.0;
        d = (vmax - vmin);
    }

    double ratio = d[0] / d[1];
    regular_2d_grid(int(std::ceil(50 * ratio)), 50, grid_V, grid_F);

    // center at 0.0
    grid_V.col(0).array() -= 0.5;
    grid_V.col(1).array() -= 0.5;

    grid_V.col(0) *= d[0];
    grid_V.col(1) *= d[1];
    grid_V.rowwise() += (vmax + vmin).transpose() / 2.0;
}
} // namespace ccd
