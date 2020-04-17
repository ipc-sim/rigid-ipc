#include "SimState.hpp"

#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <nlohmann/json.hpp>

#include <io/read_rb_scene.hpp>
#include <io/serialize_json.hpp>
#include <logger.hpp>

#include <problems/problem_factory.hpp>

#include <profiler.hpp>
#include <utils/regular_2d_grid.hpp>

namespace ccd {

SimState::SimState()
    : m_timestep_size(0.1)
    , m_step_had_collision(false)
    , m_step_has_collision(false)
    , m_step_has_intersections(false)
    , m_solve_collisions(true)
    , m_num_simulation_steps(0)
    , m_max_simulation_steps(-1)
    , m_dirty_constraints(false)
{
}

bool SimState::load_scene(const std::string& filename)
{
    PROFILER_CLEAR()

    std::string ext = boost::filesystem::extension(filename);
    boost::algorithm::to_lower(ext); // modifies ext

    nlohmann::json scene;
    if (ext == ".mjcf") {
        // TODO: Add converter from MCJF to JSON
        // scene = ...
        spdlog::error("MuJoCo file format not supported yet", ext);
        return false;
    } else if (ext == ".json") {
        std::ifstream input(filename);
        if (input.good()) {
            scene = nlohmann::json::parse(input, nullptr, false);
        } else {
            spdlog::error("Unable to open json file: {}", filename);
            return false;
        }
    } else {
        spdlog::error("Unknown scene file format: {}", ext);
        return false;
    }

    scene_file = filename;
    if (scene.is_discarded()) {
        spdlog::error("Invalid Json file");
        return false;
    }

    // Check if this is a saved simulation file
    if (scene.find("args") != scene.end()) {
        return load_simulation(scene);
    } else {
        return init(scene);
    }
    return false;
}

bool SimState::reload_scene() { return load_scene(scene_file); }

bool SimState::load_simulation(const nlohmann::json& input_args)
{
    // load original setup
    bool success = init(input_args["args"]);
    if (!success) {
        return false;
    }

    // now reload simulation history
    vertices_sequence.clear();
    for (auto& jv : input_args["animation"]["vertices_sequence"]) {
        Eigen::MatrixXd v;
        io::from_json(jv, v);
        vertices_sequence.push_back(v);
    };
    state_sequence = input_args["animation"]["state_sequence"]
                         .get<std::vector<nlohmann::json>>();
    m_num_simulation_steps = int(state_sequence.size());
    problem_ptr->state(state_sequence.back());
    return true;
}

bool SimState::init(const nlohmann::json& args_in)
{
    using namespace nlohmann;

    // Default args
    args = R"({
        "max_iterations": -1,
        "max_time": -1,
        "timestep_size": 0.1,
        "scene_type": "distance_barrier_rb_problem",
        "rigid_body_problem": {
            "rigid_bodies": [],
            "coefficient_restitution": 0.0,
            "gravity": [0.0, 0.0, 0.0],
            "collision_eps": 0.0,
            "time_stepper": "default"
        },
        "barrier_solver": {
            "inner_solver": "DEPRECATED",
            "e_b": 1e-6,
            "t_inc": 100,
            "t_init": 100,
            "m": 1,
            "c": 0.1
        },
        "newton_solver": {
            "max_iterations": 3000
        },
        "ipc_solver": {
            "max_iterations": 3000,
            "dhat_epsilon": 1e-9
        },
        "ncp_solver": {
            "max_iterations": 1000,
            "do_line_search": false,
            "solve_for_active_cstr": true,
            "convergence_tolerance": 1e-6,
            "update_type": "linearized",
            "lcp_solver": "lcp_gauss_seidel"
        },
        "distance_barrier_constraint": {
            "detection_method": "hash_grid",
            "trajectory_type": "screwing",
            "custom_initial_epsilon": 1e-2,
            "min_distance": 1e-10,
            "active_constraint_scale": 1.01,
            "barrier_type": "poly_log"
        },
        "volume_constraint": {
            "detection_method": "hash_grid",
            "trajectory_type": "screwing",
            "volume_epsilon": 1e-6,
            "custom_hashgrid_cellsize": -1,
            "time_epsilon": 1e-4
        },
        "viewport_bbox": {
            "min": [0, 0],
            "max": [0, 0]
        }
    })"_json;

    // check that incomming json doesn't have any unkown keys to avoid stupid
    // bugs
    auto patch = json::diff(args, args_in);
    for (auto& op : patch) {
        if (op["op"].get<std::string>().compare("add") == 0) {
            auto new_path = json::json_pointer(op["path"].get<std::string>());
            if (!args_in[new_path.parent_pointer()].is_array()) {
                spdlog::error(
                    "Unknown key in json path={}",
                    op["path"].get<std::string>());
            }
        }
    }

    args.merge_patch(args_in);
    m_max_simulation_steps = args["max_iterations"].get<int>();
    m_timestep_size = args["timestep_size"].get<double>();
    double max_time = args["max_time"].get<double>();
    if (max_time >= 0) {
        assert(m_max_simulation_steps == -1);
        m_max_simulation_steps = int(ceil(max_time / m_timestep_size));
    }

    auto problem_name = args["scene_type"].get<std::string>();
    problem_ptr = ProblemFactory::factory().get_problem(problem_name);
    problem_ptr->settings(args);

    m_num_simulation_steps = 0;
    m_dirty_constraints = true;

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
    active_args[problem_ptr->constraint().name()] =
        problem_ptr->constraint().settings();
    active_args[problem_ptr->solver().name()] =
        problem_ptr->solver().settings();
    if (problem_ptr->solver().has_inner_solver()) {
        active_args[problem_ptr->solver().inner_solver().name()] =
            problem_ptr->solver().inner_solver().settings();
    }

    return active_args;
}

void SimState::run_simulation(const std::string& fout)
{
    PROFILE_MAIN_POINT("run_simulation")
    PROFILE_START()

    spdlog::info("Starting simulation {}", scene_file);
    spdlog::info("Running {} iterations", m_max_simulation_steps);
    m_solve_collisions = true;
    for (int i = 0; i < m_max_simulation_steps; ++i) {
        simulation_step();
        save_simulation_step();
        spdlog::info(
            "Finished it={} sim_step={}", i + 1, m_num_simulation_steps);
    }
    save_simulation(fout);
    spdlog::info("Simulation results saved to {}", fout);

    PROFILE_END()
    LOG_PROFILER(scene_file);
}

void SimState::simulation_step()
{
    m_num_simulation_steps += 1;
    m_step_has_collision = false;
    m_step_has_intersections = false;
    // Take unconstrained step
    m_step_had_collision = problem_ptr->simulation_step(m_timestep_size);
    m_dirty_constraints = true;

    if (m_step_had_collision) {
        spdlog::debug("sim_state action=simulation_step status=had_collision");
    }

    if (m_solve_collisions && m_step_had_collision) {
        // Solve with collision constraints
        this->solve_collision();
    } else {
        m_step_has_intersections = problem_ptr->has_intersections();
    }
}

void SimState::save_simulation_step()
{
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

    // Check for intersections at the final poses.
    m_step_has_intersections =
        problem_ptr->take_step(result.x, m_timestep_size);
    // TODO: Check for collisions along the entire trajectory
    m_step_has_collision = false; // ...

    if (m_step_has_intersections) {
        spdlog::error(
            "sim_state action=solve_collisions sim_it={} "
            "status=has_intersections",
            m_num_simulation_steps);
    } else {
        spdlog::debug(
            "sim_state action=solve_collisions sim_it={} "
            "status=no_intersections",
            m_num_simulation_steps);
    }
    spdlog::debug(
        "sim_state action=solve_collisions sim_it={} "
        "status={}",
        m_num_simulation_steps,
        m_step_has_collision ? "linearized_collisions_unsolved"
                             : "collisions_solved");
    return m_step_has_intersections;
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
    // Check for intersections here
    m_step_has_intersections =
        problem_ptr->take_step(result.x, m_timestep_size);
    // TODO: Check for collisions along the entire time-step.
    m_step_has_collision = false; // ...
    spdlog::debug(
        "sim_state action=collision_resolution_step collisions={} "
        "intersections={}",
        m_step_has_collision ? "unsolved" : "solved",
        m_step_has_intersections ? "unsolved" : "solved");
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
    // TODO: Consider adding faces
    results["animation"]["group_id"] = io::to_json(problem_ptr->group_ids());

    std::ofstream o(filename);
    o << std::setw(4) << results << std::endl;
}

} // namespace ccd
