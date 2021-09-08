// Simulation executable

#include <CLI/CLI.hpp>

#include <tbb/global_control.h>
#include <tbb/task_scheduler_init.h>
#include <thread>

#include <ghc/fs_std.hpp> // filesystem

#include <SimState.hpp>
#ifdef RIGID_IPC_WITH_OPENGL
#include <viewer/UISimState.hpp>
#endif
#include <logger.hpp>
#include <profiler.hpp>

int main(int argc, char* argv[])
{
    using namespace ipc::rigid;
    set_logger_level(spdlog::level::info);

    CLI::App app("run simulation");

#ifdef RIGID_IPC_WITH_OPENGL
    bool with_viewer = true;
#else
    bool with_viewer = false;
#endif
    app.add_flag("--gui,!--ngui", with_viewer, "use the OpenGL GUI")
        ->default_val(with_viewer);

    std::string scene_path = "";
    app.add_option(
        "scene_path,-i,-s,--scene-path", scene_path,
        "JSON file with input scene");

    std::string output_dir = "";
    app.add_option(
        "output_dir,-o,--output-path", output_dir,
        "directory for results (ngui only)");

    std::string output_name = "sim.json";
    app.add_option(
           "-f,--output-name", output_name,
           "name for simulation file (ngui only)")
        ->default_val(output_name);

    int num_steps = -1;
    app.add_option(
        "--num-steps", num_steps, "number of time-steps (ngui only)");

    int checkpoint_freq = -1;
    app.add_option(
        "--chkpt,--checkpoint-frequency", checkpoint_freq,
        "number of time-steps between checkpoints (ngui only)");

    spdlog::level::level_enum loglevel = spdlog::level::info;
    app.add_option("--log,--loglevel", loglevel, "log level")
        ->default_val(loglevel)
        ->transform(CLI::CheckedTransformer(
            SPDLOG_LEVEL_NAMES_TO_LEVELS, CLI::ignore_case));

    int nthreads = tbb::task_scheduler_init::default_num_threads();
    app.add_option("--nthreads", nthreads, "maximum number of threads to use")
        ->default_val(nthreads);

    std::string patch = "";
    app.add_option("--patch", patch, "patch to input file (ngui only)")
        ->default_val(patch);

    CLI11_PARSE(app, argc, argv);

    set_logger_level(loglevel);

    if (nthreads <= 0) {
        nthreads = tbb::task_scheduler_init::default_num_threads();
    }

    if (nthreads > tbb::task_scheduler_init::default_num_threads()) {
        spdlog::warn(
            "Attempting to use more threads than available ({:d} > {:d})!",
            nthreads, tbb::task_scheduler_init::default_num_threads());
    }

    tbb::global_control thread_limiter(
        tbb::global_control::max_allowed_parallelism, nthreads);

    if (with_viewer) {
#ifdef RIGID_IPC_WITH_OPENGL
        UISimState ui;
        ui.launch(scene_path);
#else
        exit(app.exit(CLI::Error(
            "gui",
            "Unable to use GUI mode because OpenGL is disable in CMake!")));
#endif
    } else {
        if (scene_path.empty()) {
            exit(app.exit(CLI::Error(
                "scene_path", "Must provide a scene path in ngui mode!")));
        }

        if (output_dir.empty()) {
            exit(app.exit(CLI::Error(
                "output_dir",
                "Must provide a output directory in ngui mode!")));
        }

        // Create the output directory if it does not exist
        fs::create_directories(fs::path(output_dir));
        PROFILER_OUTDIR(output_dir);
        std::string fout = fmt::format("{}/{}", output_dir, output_name);

        SimState sim;

        bool success = sim.load_scene(scene_path, patch);
        if (!success) {
            return 1;
        }

        if (num_steps > 0) {
            sim.m_max_simulation_steps = num_steps;
        }

        if (checkpoint_freq > 0) {
            sim.m_checkpoint_frequency = checkpoint_freq;
        }

        sim.run_simulation(fout);
    }
}
