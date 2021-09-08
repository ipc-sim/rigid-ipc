#include "SimState.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <string>

#include <ghc/fs_std.hpp> // filesystem
#include <igl/Timer.h>
#include <igl/write_triangle_mesh.h>
#include <nlohmann/json.hpp>

#include <constants.hpp>
#include <io/read_rb_scene.hpp>
#include <io/serialize_json.hpp>
#include <io/write_gltf.hpp>
#include <io/write_obj.hpp>
#include <physics/rigid_body_problem.hpp>
#include <problems/problem_factory.hpp>
#include <utils/get_rss.hpp>
#include <utils/regular_2d_grid.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ipc::rigid {

SimState::SimState()
    : m_step_had_collision(false)
    , m_step_has_collision(false)
    , m_step_has_intersections(false)
    , m_solve_collisions(true)
    , m_num_simulation_steps(0)
    , m_max_simulation_steps(-1)
    , m_checkpoint_frequency(100)
    , m_dirty_constraints(false)
{
    initial_rss = getCurrentRSS();
}

void to_lower(std::string& s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return std::tolower(c);
    });
}

bool SimState::load_scene(const std::string& filename, const std::string& patch)
{
    PROFILER_CLEAR();
    initial_rss = getCurrentRSS();

    std::string ext = fs::path(filename).extension().string();
    to_lower(ext); // modifies ext

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
            if (patch.size()) {
                scene.merge_patch(nlohmann::json::parse(patch));
            }
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
    state_sequence = input_args["animation"]["state_sequence"]
                         .get<std::vector<nlohmann::json>>();
    if (input_args.find("stats") != input_args.end()) {
        const auto& stats = input_args["stats"];
        step_timings = stats["step_timings"].get<std::vector<double>>();
        solver_iterations = stats["solver_iterations"].get<std::vector<int>>();
        num_contacts = stats["num_contacts"].get<std::vector<int>>();
        step_minimum_distances =
            stats["step_minimum_distances"].get<std::vector<double>>();
    }
    m_num_simulation_steps = int(state_sequence.size()) - 1;
    problem_ptr->state(state_sequence.back());
    return true;
}

bool SimState::init(const nlohmann::json& args_in)
{
    using namespace nlohmann;

    // Default args (any null values will be initalized below)
    args = R"({
        "max_iterations": -1,
        "max_time": -1,
        "timestep": 0.01,
        "scene_type": "distance_barrier_rb_problem",
        "solver": "ipc_solver",
        "rigid_body_problem": {
            "rigid_bodies": [],
            "coefficient_restitution": 0.0,
            "coefficient_friction": 0.0,
            "gravity": [0.0, 0.0, 0.0],
            "collision_eps": 0.0,
            "time_stepper": "default",
            "do_intersection_check": false
        },
        "homotopy_solver": {
            "inner_solver": "DEPRECATED",
            "e_b": 1e-6,
            "t_inc": 100,
            "t_init": 100,
            "m": 1,
            "c": 0.1
        },
        "newton_solver": {
            "max_iterations": 3000,
            "convergence_criteria": "velocity",
            "energy_conv_tol": null,
            "velocity_conv_tol": null,
            "is_velocity_conv_tol_abs": false,
            "line_search_lower_bound": null,
            "linear_solver": {
                "name": "Eigen::SimplicialLDLT",
                "max_iter": 1000,
                "tolerance": 1e-10,
                "pre_max_iter": 1,
                "mtype": 2
            }
        },
        "ipc_solver": {
            "dhat_epsilon": 1e-9,
            "min_barrier_stiffness_scale": null
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
            "detection_method": "bvh",
            "trajectory_type": "piecewise_linear",
            "initial_barrier_activation_distance": 1e-3,
            "minimum_separation_distance": 0,
            "barrier_type": "ipc"
        },
        "friction_constraints": {
            "static_friction_speed_bound": 1e-3,
            "iterations": 1
        },
        "volume_constraint": {
            "detection_method": "hash_grid",
            "trajectory_type": "piecewise_linear",
            "volume_epsilon": 1e-6,
            "custom_hashgrid_cellsize": -1,
            "time_epsilon": 1e-4
        },
        "viewport_bbox": {
            "min": [0, 0],
            "max": [0, 0]
        }
    })"_json;
    // Fill in default values
    args["newton_solver"]["energy_conv_tol"] =
        Constants::DEFAULT_NEWTON_ENERGY_CONVERGENCE_TOL;
    args["newton_solver"]["velocity_conv_tol"] =
        Constants::DEFAULT_NEWTON_VELOCITY_CONVERGENCE_TOL;
    args["newton_solver"]["line_search_lower_bound"] =
        Constants::DEFAULT_LINE_SEARCH_LOWER_BOUND;
    args["ipc_solver"]["min_barrier_stiffness_scale"] =
        Constants::DEFAULT_MIN_BARRIER_STIFFNESS_SCALE;

    // Share the newton solver settings with IPC
    json newton_settings = args["newton_solver"];    // make a copy of newton
    newton_settings.merge_patch(args["ipc_solver"]); // apply ipc to newton
    args["ipc_solver"] = newton_settings; // set ipc to updated newton

    // check that incomming json doesn't have any unkown keys to avoid stupid
    // bugs
    auto patch = json::diff(args, args_in);
    for (auto& op : patch) {
        if (op["op"].get<std::string>().compare("add") == 0) {
            auto new_path = json::json_pointer(op["path"].get<std::string>());
            if (!args_in[new_path.parent_pointer()].is_array()
                && new_path.back().size() != 0
                && new_path.back().front() != '#') {
                spdlog::error(
                    "Unknown key in json path={}",
                    op["path"].get<std::string>());
            }
        }
    }

    args.merge_patch(args_in);

    // Share the newton solver settings with IPC
    newton_settings = args["newton_solver"];         // make a copy of newton
    newton_settings.merge_patch(args["ipc_solver"]); // apply ipc to newton
    args["ipc_solver"] = newton_settings; // set ipc to updated newton

    auto problem_name = args["scene_type"].get<std::string>();
    auto tmp_problem_ptr = ProblemFactory::factory().get_problem(problem_name);
    if (tmp_problem_ptr == nullptr) {
        spdlog::error("Invalid scene_type \"{}\"!", problem_name);
        return false;
    }
    problem_ptr = tmp_problem_ptr;

    bool success = problem_ptr->settings(args);
    if (!success) {
        return false;
    }

    m_max_simulation_steps = args["max_iterations"].get<int>();
    problem_ptr->timestep(args["timestep"].get<double>());
    double max_time = args["max_time"].get<double>();
    if (max_time >= 0) {
        assert(m_max_simulation_steps == -1); // is default value?
        m_max_simulation_steps = int(ceil(max_time / problem_ptr->timestep()));
    }

    m_num_simulation_steps = 0;
    m_dirty_constraints = true;

    state_sequence.clear();
    state_sequence.push_back(problem_ptr->state());
    step_timings.clear();
    solver_iterations.clear();
    num_contacts.clear();
    step_minimum_distances.clear();

    return true;
}

nlohmann::json SimState::get_active_config()
{
    nlohmann::json active_args;
    active_args["timestep"] = problem_ptr->timestep();
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

void print_progress_bar(
    int cur_iter,
    int max_iter,
    double time_passed,
    int bar_width = 70,
    bool clear_end = false)
{
    float progress = std::max(0.0f, std::min(1.0f, cur_iter / float(max_iter)));
    int pos = bar_width * progress;
    int len = fmt::format("{:d}", max_iter).length();
    fmt::print(
        "{3:2d}%|{0:█>{1}}{0: >{2}}| {4: >{6}d}/{5:d} ({7:.1f}s/{8:.1f}s)\r", //
        "", pos, std::max(bar_width - pos, 0), int(progress * 100.0), cur_iter,
        max_iter, len, time_passed,
        cur_iter > 0 ? (time_passed * max_iter) / cur_iter
                     : std::numeric_limits<double>::infinity());
    if (progress >= 1) {
        if (clear_end) {
            fmt::print("{0: {1}}", "", bar_width * 2);
        } else {
            std::cout << std::endl;
        }
    }
    std::cout.flush();
}

void SimState::run_simulation(const std::string& fout)
{
    PROFILE_MAIN_POINT("run_simulation");
    PROFILE_START();

    fs::path fout_path(fout);
    fs::create_directories(fout_path.parent_path());
    std::string chkpt_base =
        (fout_path.parent_path() / fout_path.stem()).string();

    if (m_max_simulation_steps <= 0) {
        m_max_simulation_steps = 1000;
    }

    spdlog::info("Starting simulation {}", scene_file);
    spdlog::info("Running {} iterations", m_max_simulation_steps);

    igl::Timer timer;
    timer.start();

    m_solve_collisions = true;
    print_progress_bar(0, m_max_simulation_steps, 0);
    for (int i = 0; i < m_max_simulation_steps; ++i) {
        simulation_step();
        save_simulation_step();
        spdlog::info(
            "Finished it={} sim_step={}", i + 1, m_num_simulation_steps);

        // Checkpoint the simulation every m_checkpoint_frequency time-steps
        if ((i + 1) % m_checkpoint_frequency == 0
            && (i + 1) < m_max_simulation_steps) {
            std::string chkpt_fout = fmt::format(
                "{}-chkpt{:05d}.json", chkpt_base, m_num_simulation_steps);
            save_simulation(chkpt_fout);
            spdlog::info("Simulation checkpoint saved to {}", chkpt_fout);
        }
        print_progress_bar(
            i + 1, m_max_simulation_steps, timer.getElapsedTime());
    }

    timer.stop();
    fmt::print(
        "Simulation finished (total_runtime={:g}s average_fps={:g})\n",
        timer.getElapsedTime(),
        m_max_simulation_steps / timer.getElapsedTime());

    save_simulation(fout);
    spdlog::info("Simulation results saved to {}", fout);
    fs::path gltf_filename(fout);
    gltf_filename.replace_extension(".glb");
    save_gltf(gltf_filename.string());
    spdlog::info("Animation saved to {}", gltf_filename.string());

    PROFILE_END();
    LOG_PROFILER(scene_file);
}

void SimState::simulation_step()
{
    m_num_simulation_steps += 1;
    m_step_has_collision = false;
    m_step_has_intersections = false;

    step_timer.start();
    problem_ptr->simulation_step(
        m_step_had_collision, m_step_has_intersections, m_solve_collisions);
    step_timer.stop();

    if (m_step_had_collision) {
        spdlog::debug("sim_state action=simulation_step status=had_collision");
    }

    if (m_step_has_intersections) {
        spdlog::error(
            "sim_state action=simulation_step sim_it={} "
            "status=has_intersections",
            m_num_simulation_steps);
    } else {
        spdlog::debug(
            "sim_state action=simulation_step sim_it={} "
            "status=no_intersections",
            m_num_simulation_steps);
    }
    // spdlog::debug(
    //     "sim_state action=solve_collisions sim_it={} "
    //     "status={}",
    //     m_num_simulation_steps,
    //     m_step_has_collision ? "linearized_collisions_unsolved"
    //                          : "collisions_solved");
}

void SimState::save_simulation_step()
{
    PROFILE_POINT("SimState::save_simulation_step");
    PROFILE_START();

    state_sequence.push_back(problem_ptr->state());
    step_timings.push_back(step_timer.getElapsedTime());
    solver_iterations.push_back(problem_ptr->opt_result.num_iterations);
    num_contacts.push_back(problem_ptr->num_contacts());
    step_minimum_distances.push_back(problem_ptr->compute_min_distance());

    PROFILE_END();
}

bool SimState::save_simulation(const std::string& filename)
{
    PROFILE_POINT("SimState::save_simulation");
    PROFILE_START();

    nlohmann::json results;
    results["args"] = args;
    results["animation"] = nlohmann::json();
    results["animation"]["state_sequence"] = state_sequence;

    nlohmann::json stats;
    stats["dim"] = problem_ptr->dim();
    stats["num_bodies"] = problem_ptr->num_bodies();
    stats["num_vertices"] = problem_ptr->num_vertices();
    stats["num_edges"] = problem_ptr->num_edges();
    stats["num_faces"] = problem_ptr->num_faces();
    stats["num_timesteps"] = m_max_simulation_steps;
    stats["memory"] = getPeakRSS() - initial_rss;
    stats["step_timings"] = step_timings;
    stats["solver_iterations"] = solver_iterations;
    stats["num_contacts"] = num_contacts;
    stats["step_minimum_distances"] = step_minimum_distances;
    stats["solve_stats"] = problem_ptr->solver().stats();
    results["stats"] = stats;

    std::ofstream file(filename);
    if (!file) {
        PROFILE_END();
        return false;
    }

    file << results.dump();

    PROFILE_END();
    return true;
}

bool SimState::save_obj_sequence(const std::string& dir_name)
{
    // Create the output directory if it does not exist
    fs::path dir_path(dir_name);
    fs::create_directories(dir_path);

    bool success = true;
    for (int i = 0; i < state_sequence.size(); i++) {
        const auto& state = state_sequence[i];
        problem_ptr->state(state);
        write_obj(
            (dir_path / fmt::format("{:05d}.obj", i)).string(), *problem_ptr,
            false);
    }

    problem_ptr->state(state_sequence.back());

    return success;
}

bool SimState::save_gltf(const std::string& filename)
{
    std::vector<PosesD> poses(state_sequence.size());

    bool success = true;
    for (int i = 0; i < state_sequence.size(); i++) {
        const auto& state = state_sequence[i];
        std::vector<nlohmann::json> jrbs = state["rigid_bodies"];
        for (int j = 0; j < jrbs.size(); j++) {
            VectorMax3d position;
            VectorMax3d rotation;
            from_json(jrbs[j]["position"], position);
            from_json(jrbs[j]["rotation"], rotation);
            poses[i].emplace_back(position, rotation);
        }
    }

    std::shared_ptr<RigidBodyProblem> rbp =
        std::dynamic_pointer_cast<RigidBodyProblem>(problem_ptr);

    return write_gltf(
        filename, rbp->m_assembler, poses, problem_ptr->timestep());
}

} // namespace ipc::rigid
